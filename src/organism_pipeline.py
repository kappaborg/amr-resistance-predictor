"""Organism-general pipeline — the same method, any pathogen (generalization proof).

Runs QC -> AMRFinderPlus -> MLST -> feature matrix -> multi-drug labels -> per-drug models
(rules baseline / logreg / calibrated XGBoost, VME<=target operating point) for a configured
organism, writing tag-based outputs. Reuses the K. pneumoniae pipeline's pure functions unchanged.

    python -m src.organism_pipeline --organism ecoli
    python -m src.organism_pipeline --organism ecoli --models-only   # skip annotate/mlst if cached
"""
from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from sklearn.frozen import FrozenEstimator
from sklearn.calibration import CalibratedClassifierCV
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import brier_score_loss
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from xgboost import XGBClassifier

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from src.data.qc_metadata import fetch_qc  # noqa: E402
from src.evaluation.metrics import full_report  # noqa: E402
from src.models.train_panel import pick_threshold  # noqa: E402

ENV_BIN = "/opt/homebrew/anaconda3/envs/amr-resistance-predictor/bin"
AMR_BIN, MLST_BIN = f"{ENV_BIN}/amrfinder", f"{ENV_BIN}/mlst"
import os  # noqa: E402
MLST_ENV = {**os.environ, "PATH": f"{ENV_BIN}:{os.environ.get('PATH', '')}"}
AMR_API = "https://www.bv-brc.org/api/genome_amr/"
LAB = 'eq(evidence,%22Laboratory%20Method%22)'
GENOMES = REPO / "data/raw/genomes"
AMR_CACHE = REPO / "data/interim/amrfinder"
MLST_CACHE = REPO / "data/interim/mlst"

# Drug-specific known-gene rules (organism-independent mechanisms; grlA_ = S. aureus parC homolog).
RULES = {
    "ciprofloxacin": ["gyrA_", "parC_", "grlA_", "gyrB_", "parE_", "qnr", "aac(6')-Ib-cr", "qepA"],
    "gentamicin": ["aac(3)", "ant(2'')-Ia", "armA", "rmt"],
    "trimethoprim_sulfamethoxazole": ["sul1", "sul2", "sul3", "dfrA", "dfrB", "dfrG"],
    "ampicillin": ["blaTEM", "blaSHV", "blaOXA-1", "blaCARB", "blaPSE", "blaCTX-M", "blaCMY", "blaDHA"],
    "ceftriaxone": ["blaCTX-M", "blaCMY", "blaDHA", "blaSHV-12", "blaSHV-2", "blaOXA-48", "blaKPC", "blaNDM"],
    # S. aureus panel
    "oxacillin": ["mecA", "mecC", "mecB"],                       # MRSA — PBP2a
    "cefoxitin": ["mecA", "mecC"],                               # S. aureus: mecA screen (not AmpC)
    "erythromycin": ["erm", "msr(A)", "mph(C)", "mef", "lsa"],   # macrolide (MLSb + efflux)
    "clindamycin": ["erm", "lnu", "vga", "sal(A)", "lsa"],       # lincosamide (MLSb)
    # A. baumannii panel (carbapenemases incl. OXA-23/24/58/51-like)
    "meropenem": ["blaKPC", "blaNDM", "blaVIM", "blaIMP", "blaOXA-23", "blaOXA-24", "blaOXA-58",
                  "blaOXA-51", "blaOXA-40", "blaOXA-48"],
    "imipenem": ["blaKPC", "blaNDM", "blaVIM", "blaIMP", "blaOXA-23", "blaOXA-24", "blaOXA-58",
                 "blaOXA-51", "blaOXA-40", "blaOXA-48"],
    "amikacin": ["aac(6')-I", "armA", "rmt", "aph(3')"],
    # --- Round-4 expansion drugs (known-gene baselines; ⚕ review before use) ---
    # ACQUIRED ESBL/carbapenemase only — blaPDC excluded: it is the INTRINSIC P. aeruginosa AmpC
    # (present in ~100% of genomes), so it cannot discriminate R from S; resistance is via its
    # derepression/overexpression, which a determinant-presence baseline cannot see (⚕).
    "ceftazidime": ["blaCTX-M", "blaSHV-12", "blaSHV-2", "blaPER", "blaVEB", "blaGES", "blaKPC",
                    "blaVIM", "blaNDM", "blaIMP", "blaOXA-10"],
    # aph(3') excluded — matches the intrinsic P. aeruginosa aph(3')-IIb (~99%), which does not confer
    # tobramycin resistance; tobramycin resistance = AAC(6')/ANT(2'')/ANT(4') + 16S methyltransferases.
    "tobramycin": ["aac(6')", "ant(2'')-Ia", "ant(4')", "armA", "rmt"],
    "chloramphenicol": ["cat", "catA", "catB", "cmlA", "cml", "floR"],
    "tetracycline": ["tet(A)", "tet(B)", "tet(C)", "tet(D)", "tet(G)", "tet(K)", "tet(L)", "tet(M)",
                     "tet(O)", "tet(W)", "tet(X)"],
    "vancomycin": ["vanA", "vanB", "vanC", "vanD", "vanM", "vanN"],
    "penicillin": ["pbp1a", "pbp2b", "pbp2x", "blaZ", "mecA"],   # S. pneumoniae PBP mosaics (hard); ⚕
}

# Organism-specific rule-gene OVERRIDES — the honest known-gene baseline differs by organism for some
# drugs, so a single shared list is unfair. Keyed by (organism_key, drug); falls back to RULES[drug].
#   - Enterococcus ampicillin resistance = PBP5 mutations, NOT beta-lactamases (blaTEM/SHV absent in
#     Enterococcus), so the generic beta-lactamase baseline scores ~chance; PBP5 is the fair baseline.
#   - S. pneumoniae penicillin = PBP1a/2b/2x mosaics only (drop the staph blaZ/mecA noise).
ORG_RULES = {
    ("efaecium", "ampicillin"): ["pbp5"],
    ("spneumoniae", "penicillin"): ["pbp1a", "pbp2b", "pbp2x"],
}


def rules_for(organism: str, drug: str) -> list[str]:
    """The known-gene rule tokens for this organism+drug (organism override, else the shared default)."""
    return ORG_RULES.get((organism, drug)) or RULES.get(drug, [])


ORGANISMS = {
    "ecoli": {"taxon": 562, "amrfinder": "Escherichia", "ids": "data/raw/ecoli_ids.txt",
              "len": (4.0, 6.0),
              "drugs": {"ciprofloxacin": "ciprofloxacin", "gentamicin": "gentamicin",
                        "trimethoprim_sulfamethoxazole": "trimethoprim%2Fsulfamethoxazole",
                        "ampicillin": "ampicillin", "ceftriaxone": "ceftriaxone"}},
    "saureus": {"taxon": 1280, "amrfinder": "Staphylococcus_aureus", "ids": "data/raw/saureus_ids.txt",
                "len": (2.5, 3.3),  # S. aureus ~2.8 Mbp
                "drugs": {"oxacillin": "oxacillin", "cefoxitin": "cefoxitin",
                          "ciprofloxacin": "ciprofloxacin", "erythromycin": "erythromycin",
                          "clindamycin": "clindamycin"}},
    "abaumannii": {"taxon": 470, "amrfinder": "Acinetobacter_baumannii",
                   "ids": "data/raw/abaumannii_ids.txt",
                   "len": (3.4, 4.6),  # A. baumannii ~3.9 Mbp
                   "drugs": {"meropenem": "meropenem", "imipenem": "imipenem",
                             "ciprofloxacin": "ciprofloxacin", "gentamicin": "gentamicin",
                             "amikacin": "amikacin"}},
    # --- Round-4 expansion candidates (scaffolded; acquire+train per-organism when BV-BRC is up) ---
    "paeruginosa": {"taxon": 287, "amrfinder": "Pseudomonas_aeruginosa",
                    "ids": "data/raw/paeruginosa_ids.txt",
                    "len": (5.5, 7.5),  # P. aeruginosa ~6.3 Mbp
                    "drugs": {"meropenem": "meropenem", "ceftazidime": "ceftazidime",
                              "ciprofloxacin": "ciprofloxacin", "tobramycin": "tobramycin"}},
    "senterica": {"taxon": 28901, "amrfinder": "Salmonella",
                  "ids": "data/raw/senterica_ids.txt",
                  "len": (4.4, 5.4),  # S. enterica ~4.8 Mbp
                  "drugs": {"ampicillin": "ampicillin", "ceftriaxone": "ceftriaxone",
                            "ciprofloxacin": "ciprofloxacin", "chloramphenicol": "chloramphenicol",
                            "trimethoprim_sulfamethoxazole": "trimethoprim%2Fsulfamethoxazole"}},
    "efaecium": {"taxon": 1352, "amrfinder": "Enterococcus_faecium",
                 "ids": "data/raw/efaecium_ids.txt",
                 "len": (2.4, 3.4),  # E. faecium ~2.9 Mbp
                 "drugs": {"vancomycin": "vancomycin", "ampicillin": "ampicillin",
                           "tetracycline": "tetracycline"}},   # cipro/erythromycin dropped (too few)
    "ecloacae": {"taxon": 550, "amrfinder": "Enterobacter_cloacae",
                 "ids": "data/raw/ecloacae_ids.txt",
                 "len": (4.6, 6.0),  # E. cloacae ~5.2 Mbp
                 "drugs": {"meropenem": "meropenem", "ceftazidime": "ceftazidime",
                           "ciprofloxacin": "ciprofloxacin", "gentamicin": "gentamicin",
                           "trimethoprim_sulfamethoxazole": "trimethoprim%2Fsulfamethoxazole"}},
    "spneumoniae": {"taxon": 1313, "amrfinder": "Streptococcus_pneumoniae",
                    "ids": "data/raw/spneumoniae_ids.txt",
                    "len": (1.9, 2.4),  # S. pneumoniae ~2.1 Mbp
                    "drugs": {"penicillin": "penicillin", "erythromycin": "erythromycin",
                              "tetracycline": "tetracycline",
                              "trimethoprim_sulfamethoxazole": "trimethoprim%2Fsulfamethoxazole"}},  # FQ dropped
    "cjejuni": {"taxon": 197, "amrfinder": "Campylobacter",
                "ids": "data/raw/cjejuni_ids.txt",
                "len": (1.5, 2.0),  # C. jejuni ~1.7 Mbp
                "drugs": {"ciprofloxacin": "ciprofloxacin", "tetracycline": "tetracycline",
                          "erythromycin": "erythromycin", "gentamicin": "gentamicin"}},
}


def annotate_one(gid: str, organism: str) -> str:
    tsv = AMR_CACHE / f"{gid}.tsv"
    if tsv.exists() and tsv.stat().st_size > 0:
        return gid
    fna = AMR_CACHE / f"{gid}.fna"
    fna.write_bytes(gzip.decompress((GENOMES / f"{gid}.fna.gz").read_bytes()))
    r = subprocess.run([AMR_BIN, "-n", str(fna), "--organism", organism, "--threads", "1",
                        "-o", str(tsv)], capture_output=True, text=True)
    fna.unlink(missing_ok=True)
    if r.returncode != 0:
        tsv.unlink(missing_ok=True)
    return gid


def mlst_one(gid: str) -> tuple[str, str]:
    out = MLST_CACHE / f"{gid}.tsv"
    if out.exists() and out.stat().st_size > 0:
        parts = out.read_text().split("\t")
        return gid, parts[2] if len(parts) > 2 else "-"
    fna = MLST_CACHE / f"{gid}.fna"
    fna.write_bytes(gzip.decompress((GENOMES / f"{gid}.fna.gz").read_bytes()))
    r = subprocess.run([MLST_BIN, str(fna)], capture_output=True, text=True, env=MLST_ENV)
    fna.unlink(missing_ok=True)
    if r.returncode != 0 or "\t" not in r.stdout:
        return gid, "-"
    out.write_text(r.stdout)
    return gid, r.stdout.split("\t")[2]


def parse_dets(tsv: Path) -> set[str]:
    return {(r.get("Element symbol") or "").strip()
            for r in csv.DictReader(tsv.open(), delimiter="\t")
            if (r.get("Type") or "").strip().upper() == "AMR" and (r.get("Element symbol") or "").strip()}


def fetch_labels(taxon: int, drug_enc: str) -> dict[str, str]:
    q = (f"and(eq(taxon_id,{taxon}),{LAB},eq(antibiotic,%22{drug_enc}%22))"
         f"&select(genome_id,resistant_phenotype)&limit(25000)")
    r = requests.get(f"{AMR_API}?{q}", headers={"Accept": "application/json"}, timeout=180)
    byg: dict[str, set] = defaultdict(set)
    for x in r.json():
        p = x.get("resistant_phenotype")
        if p in ("Resistant", "Susceptible"):
            byg[x["genome_id"]].add(p)
    return {g: next(iter(v)) for g, v in byg.items() if len(v) == 1}


def num(x, k, d):
    v = x.get(k)
    return float(v) if v not in (None, "") else d


def run(tag: str, cfg: dict, models_only: bool, workers: int = 8) -> int:
    lo, hi = cfg["len"]
    ids = [g for g in (REPO / cfg["ids"]).read_text().split()
           if (GENOMES / f"{g}.fna.gz").exists()]
    print(f"[{tag}] {len(ids)} genomes on disk")

    # 1. QC (MIMAG + organism length window)
    qc = {r["genome_id"]: r for r in fetch_qc(ids)}
    master = [g for g in ids if g in qc and num(qc[g], "checkm_completeness", 0) >= 90
              and num(qc[g], "checkm_contamination", 99) <= 5
              and num(qc[g], "contigs", 1e9) <= 500
              and lo <= num(qc[g], "genome_length", 0) / 1e6 <= hi]
    print(f"[{tag}] QC passed: {len(master)}")

    if not models_only:
        AMR_CACHE.mkdir(parents=True, exist_ok=True); MLST_CACHE.mkdir(parents=True, exist_ok=True)
        todo = [g for g in master if not (AMR_CACHE / f"{g}.tsv").exists()]
        print(f"[{tag}] annotating {len(todo)} (AMRFinderPlus, --organism {cfg['amrfinder']})...")
        with ThreadPoolExecutor(max_workers=workers) as ex:
            for i, _ in enumerate(as_completed([ex.submit(annotate_one, g, cfg["amrfinder"]) for g in todo]), 1):
                if i % 100 == 0 or i == len(todo):
                    print(f"    annotate {i}/{len(todo)}")
        print(f"[{tag}] MLST typing...")
        lineage = {}
        with ThreadPoolExecutor(max_workers=max(2, workers // 2)) as ex:
            for gid, st in [f.result() for f in as_completed([ex.submit(mlst_one, g) for g in master])]:
                lineage[gid] = f"ST{st}" if st not in ("-", "", None) else f"novel_{gid}"
        with (REPO / f"data/processed/{tag}_lineages.csv").open("w", newline="") as f:
            w = csv.writer(f); w.writerow(["genome_id", "lineage"])
            w.writerows([[g, lineage[g]] for g in master])
    lineage = {r["genome_id"]: r["lineage"]
               for r in csv.DictReader((REPO / f"data/processed/{tag}_lineages.csv").open())}

    # feature matrix + labels
    per_genome = {g: parse_dets(AMR_CACHE / f"{g}.tsv") for g in master if (AMR_CACHE / f"{g}.tsv").exists()}
    cols = sorted(set().union(*per_genome.values()))
    feats = pd.DataFrame({g: [1 if c in per_genome[g] else 0 for c in cols] for g in per_genome},
                         index=cols).T
    feats.to_csv(REPO / f"data/processed/{tag}_features.csv", index_label="genome_id")
    labs = {col: fetch_labels(cfg["taxon"], enc) for col, enc in cfg["drugs"].items()}

    # per-drug models
    results = {}
    # persist the lineage split (train/test STs) so the zero-leakage test can assert it per organism
    split_rows = []

    print(f"\n[{tag}] per-drug models (VME<=3% operating point, unseen-lineage test)")
    print(f"{'drug':30s}{'ROC':>7}{'PR':>7}{'VME':>7}{'ME':>7}{'rulesROC':>9}{'Brier':>8}")
    for col in cfg["drugs"]:
        g = [x for x in feats.index if labs[col].get(x) in ("Resistant", "Susceptible")]
        X = feats.loc[g].values.astype(int)
        y = np.array([1 if labs[col][x] == "Resistant" else 0 for x in g])
        groups = np.array([lineage[x] for x in g])
        if len(set(y)) < 2 or min((y == 1).sum(), (y == 0).sum()) < 20:
            print(f"{col:30s}  (too few labels: R={int((y==1).sum())} S={int((y==0).sum())})")
            continue
        sgkf = StratifiedGroupKFold(n_splits=4, shuffle=True, random_state=42)
        tr, te = next(sgkf.split(X, y, groups))
        if not split_rows:   # one representative lineage split, for the zero-leakage test
            split_rows = [(g[i], groups[i], "train") for i in tr] + \
                         [(g[i], groups[i], "test") for i in te]
        idx = [i for i, c in enumerate(cols) if any(t in c for t in rules_for(tag, col))]
        rule_s = (X[te][:, idx].sum(1) > 0).astype(float) if idx else np.zeros(len(te))
        lr = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=42).fit(X[tr], y[tr])
        oof = cross_val_predict(LogisticRegression(max_iter=2000, class_weight="balanced", random_state=42),
                                X[tr], y[tr], groups=groups[tr],
                                cv=StratifiedGroupKFold(4, shuffle=True, random_state=42),
                                method="predict_proba")[:, 1]
        thr = pick_threshold(y[tr], oof, 0.03)
        lr_s = lr.predict_proba(X[te])[:, 1]
        rep = full_report(y[te].tolist(), (lr_s >= thr).astype(int).tolist(), lr_s.tolist())
        rep_rule = full_report(y[te].tolist(), rule_s.astype(int).tolist(), rule_s.tolist())
        # calibrated XGBoost + Brier (non-negotiable #7: always calibrate + report calibration)
        spw = (y[tr] == 0).sum() / max((y[tr] == 1).sum(), 1)
        xgb = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.1, subsample=0.9,
                            colsample_bytree=0.8, eval_metric="logloss", scale_pos_weight=spw,
                            random_state=42, n_jobs=4)
        xoof = cross_val_predict(xgb, X[tr], y[tr], groups=groups[tr],
                                 cv=StratifiedGroupKFold(4, shuffle=True, random_state=42),
                                 method="predict_proba")[:, 1]
        iso = IsotonicRegression(out_of_bounds="clip").fit(xoof, y[tr])
        xthr = pick_threshold(y[tr], iso.transform(xoof), 0.03)
        xgb.fit(X[tr], y[tr])
        xs = iso.transform(xgb.predict_proba(X[te])[:, 1])
        rep_cal = full_report(y[te].tolist(), (xs >= xthr).astype(int).tolist(), xs.tolist())
        rep_cal["brier"] = float(brier_score_loss(y[te], xs))
        results[col] = {"logreg": rep, "xgboost_calibrated": rep_cal, "rules": rep_rule,
                        "n_test": len(te), "test_R": int(y[te].sum()),
                        "test_lineages": len(set(groups[te]))}
        print(f"{col:30s}{rep['roc_auc']:7.3f}{rep['pr_auc']:7.3f}"
              f"{rep['very_major_error_rate']:7.3f}{rep['major_error_rate']:7.3f}"
              f"{rep_rule['roc_auc']:9.3f}{rep_cal['brier']:8.3f}")

    import json
    (REPO / f"results/metrics/{tag}_metrics.json").write_text(json.dumps(results, indent=2))
    if split_rows:
        with (REPO / f"data/processed/{tag}_split.csv").open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["genome_id", "lineage", "fold"])
            w.writerows(split_rows)
    print(f"\nwrote data/processed/{tag}_features.csv, results/metrics/{tag}_metrics.json, "
          f"data/processed/{tag}_split.csv")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True, choices=list(ORGANISMS))
    ap.add_argument("--models-only", action="store_true")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    return run(args.organism, ORGANISMS[args.organism], args.models_only, args.workers)


if __name__ == "__main__":
    sys.exit(main())
