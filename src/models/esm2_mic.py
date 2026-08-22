"""Foundation-model experiment #2 — ESM-2 protein language model for MIC prediction.

Motivation. Determinant presence/absence discards **allele-level** variation: a genome is scored the
same whether it carries blaKPC-2 or a weak KPC point variant, or a porin with vs without a
resistance mutation. That variation is exactly where the hard MIC drugs (meropenem, cefoxitin) may
have residual signal — presence/absence saturates, but the *sequence* of the resistance protein still
differs. ESM-2 (a protein language model) embeds each resistance-gene protein into a vector that
encodes its sequence, so two alleles of the same gene get different embeddings.

Pipeline (each stage caches to disk; re-runs are incremental):
  1. extract  — per genome, pull each AMRFinderPlus AMR hit's nucleotide region from the FASTA using
                (contig, start, stop, strand), translate to protein. Cache genome -> [(sym, seq)].
  2. embed    — dedupe unique protein sequences (alleles repeat across genomes), embed each ONCE with
                ESM-2 150M (mean-pool over residues), cache sha1(seq) -> 640-d vector.
  3. features — per genome, build a PLM vector. Two aggregations: mean over ALL its AMR proteins, and
                a per-drug-class mean (carbapenem proteins for meropenem, cephalosporin for cefoxitin).
  4. compare  — MIC regression (log2 mg/L) under the phylogeny-aware split, three feature sets:
                determinants only (baseline) / PLM only / determinants + PLM. Metric: Essential
                Agreement (±1 doubling dilution, the CLSI/FDA criterion). Does ESM-2 raise EA?

    python -m src.models.esm2_mic --stage extract      # no torch needed
    python -m src.models.esm2_mic --stage embed        # needs ESM-2 (MPS)
    python -m src.models.esm2_mic --stage compare
    python -m src.models.esm2_mic --stage all
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import pickle
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np
import pandas as pd

REPEATS = 6   # 6 seeds x 5 folds = 30 paired folds

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
GENOMES = REPO / "data/raw/genomes"
AMRDIR = REPO / "data/interim/amrfinder"
CACHE = REPO / "data/interim/esm"
OUT = REPO / "results/reports/summary_21_esm2_mic.md"

# Organism registry — the pipeline is organism-agnostic; only these paths + taxon + drug list change.
# The embedding cache is keyed by sequence hash + model, so it is SHARED across organisms.
ORGANISMS = {
    "kpneu": {"name": "Klebsiella pneumoniae", "taxon": 573, "proteins": "proteins.pkl",
              "features": "data/processed/thin_slice_cipro_features.csv",
              "lineages": "data/processed/thin_slice_cipro_lineages.csv",
              "drugs": ["meropenem", "cefoxitin", "gentamicin", "ciprofloxacin"]},
    "abaum": {"name": "Acinetobacter baumannii", "taxon": 470, "proteins": "proteins_abaum.pkl",
              "features": "data/processed/abaumannii_features.csv",
              "lineages": "data/processed/abaumannii_lineages.csv",
              "drugs": ["meropenem", "imipenem"]},
}
# back-compat module defaults (K. pneumoniae) for helper scripts that import these names
FEATURES = REPO / ORGANISMS["kpneu"]["features"]
LINEAGES = REPO / ORGANISMS["kpneu"]["lineages"]
PROT_PKL = CACHE / ORGANISMS["kpneu"]["proteins"]

# ESM-2 model registry: key -> (fair-esm loader name, final layer index, embed dim)
MODELS = {"150M": ("esm2_t30_150M_UR50D", 30, 640),
          "650M": ("esm2_t33_650M_UR50D", 33, 1280)}


def emb_path(model_key: str) -> Path:
    return CACHE / f"embeddings_{model_key}.pkl"   # sha1(seq) -> np.float32[dim]

# drug -> AMRFinderPlus Subclass tokens whose proteins are mechanistically relevant to that drug
DRUG_CLASS = {"meropenem": {"CARBAPENEM"}, "imipenem": {"CARBAPENEM"},
              "cefoxitin": {"CEPHALOSPORIN", "CARBAPENEM", "BETA-LACTAM"},
              "gentamicin": {"AMINOGLYCOSIDE"}, "ciprofloxacin": {"QUINOLONE"}}


def fetch_mic_taxon(drug: str, taxon: int) -> dict[str, float]:
    """genome_id -> median log2(MIC mg/L) for any organism (taxon-parametrised)."""
    import re
    import requests
    from collections import defaultdict
    lab = 'eq(evidence,%22Laboratory%20Method%22)'
    q = (f"and(eq(taxon_id,{taxon}),{lab},eq(antibiotic,%22{drug}%22))"
         f"&select(genome_id,measurement,measurement_unit)&limit(25000)")
    r = None
    for _ in range(4):
        r = requests.get(f"https://www.bv-brc.org/api/genome_amr/?{q}",
                         headers={"Accept": "application/json"}, timeout=180)
        if r.status_code == 200 and isinstance(r.json(), list):
            break
    byg: dict[str, list] = defaultdict(list)
    for x in r.json():
        m, unit = x.get("measurement"), x.get("measurement_unit")
        if not m or unit != "mg/L":
            continue
        try:
            v = float(re.sub(r"[<>=]", "", str(m)).strip())
        except ValueError:
            continue
        if v > 0:
            byg[x["genome_id"]].append(np.log2(v))
    return {g: float(np.median(vs)) for g, vs in byg.items() if vs}

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_CODON = {  # standard bacterial table 11 (ATG/GTG/TTG start handled as M only if first; keep simple)
}


def _translate(nt: str) -> str:
    from Bio.Seq import Seq
    try:
        return str(Seq(nt).translate(table=11, to_stop=False)).rstrip("*")
    except Exception:  # noqa: BLE001
        return ""


def _read_fasta(gz: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur = None
    with gzip.open(gz, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur is not None:
                seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def stage_extract(genome_ids: list[str], prot_pkl: Path = PROT_PKL) -> dict:
    CACHE.mkdir(parents=True, exist_ok=True)
    out: dict[str, list] = {}
    if prot_pkl.exists():
        out = pickle.loads(prot_pkl.read_bytes())
    todo = [g for g in genome_ids if g not in out]
    print(f"extract: {len(out)} cached, {len(todo)} to do")
    for i, g in enumerate(todo):
        fa, tsv = GENOMES / f"{g}.fna.gz", AMRDIR / f"{g}.tsv"
        if not (fa.exists() and tsv.exists()):
            out[g] = []
            continue
        contigs = _read_fasta(fa)
        prots = []
        for ln in tsv.read_text().splitlines()[1:]:
            c = ln.split("\t")
            if len(c) < 12 or c[8] != "AMR":
                continue
            cid, start, stop, strand, sym, subclass = c[1], c[2], c[3], c[4], c[5], c[11]
            seq = contigs.get(cid)
            if not seq:
                continue
            try:
                nt = seq[int(start) - 1:int(stop)]
            except ValueError:
                continue
            if strand == "-":
                nt = nt.translate(_COMP)[::-1]
            aa = _translate(nt)
            if 20 <= len(aa) <= 1500:
                prots.append((sym, subclass, aa))
        out[g] = prots
        if (i + 1) % 250 == 0:
            print(f"  {i + 1}/{len(todo)} genomes...", flush=True)
            prot_pkl.write_bytes(pickle.dumps(out))
    prot_pkl.write_bytes(pickle.dumps(out))
    npro = sum(len(v) for v in out.values())
    uniq = len({s for v in out.values() for _, _, s in v})
    print(f"extract done: {len(out)} genomes, {npro} proteins, {uniq} unique sequences")
    return out


def stage_embed(proteins: dict, model_key: str = "150M") -> dict:
    import torch
    import esm
    path = emb_path(model_key)
    emb: dict[str, np.ndarray] = pickle.loads(path.read_bytes()) if path.exists() else {}
    uniq = {}
    for v in proteins.values():
        for _, _, s in v:
            h = hashlib.sha1(s.encode()).hexdigest()
            if h not in emb:
                uniq[h] = s
    seqs = list(uniq.items())
    name, L, _ = MODELS[model_key]
    print(f"embed [{model_key} = {name}]: {len(emb)} cached, {len(seqs)} unique sequences to embed")
    if not seqs:
        return emb
    dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    model, alphabet = getattr(esm.pretrained, name)()
    bc = alphabet.get_batch_converter()
    model = model.eval().to(dev)
    valid = set("ACDEFGHIKLMNPQRSTVWY")
    def clean(s):  # internal stops / ambiguous codons -> X (unknown), which ESM-2 handles
        return "".join(c if c in valid else "X" for c in s)
    seqs.sort(key=lambda kv: len(kv[1]))          # length-sort so batches pad efficiently
    B = 16 if model_key == "150M" else 6          # 650M is heavier — smaller batch
    for i in range(0, len(seqs), B):
        chunk = seqs[i:i + B]
        batch = [(h, clean(s[:1022])) for h, s in chunk]  # ESM-2 context cap + sanitised
        _, _, toks = bc(batch)
        toks = toks.to(dev)
        with torch.no_grad():
            rep = model(toks, repr_layers=[L])["representations"][L]
        for j, (h, s) in enumerate(chunk):
            n = min(len(s), 1022)
            emb[h] = rep[j, 1:1 + n].mean(0).float().cpu().numpy().astype(np.float32)
        if (i // B) % 25 == 0:
            print(f"  {i + len(chunk)}/{len(seqs)} embedded...", flush=True)
            path.write_bytes(pickle.dumps(emb))
    path.write_bytes(pickle.dumps(emb))
    print(f"embed done: {len(emb)} vectors cached")
    return emb


def genome_plm_vectors(proteins: dict, emb: dict, drug: str) -> tuple[dict, dict]:
    """Return (mean-over-all-AMR-proteins vec, mean-over-drug-class-proteins vec) per genome."""
    classes = DRUG_CLASS.get(drug, set())
    dim = len(next(iter(emb.values())))
    allv, drugv = {}, {}
    for g, plist in proteins.items():
        vecs = [emb[hashlib.sha1(s.encode()).hexdigest()] for _, _, s in plist
                if hashlib.sha1(s.encode()).hexdigest() in emb]
        if vecs:
            allv[g] = np.mean(vecs, axis=0)
        dv = [emb[hashlib.sha1(s.encode()).hexdigest()] for _, sub, s in plist
              if any(t in sub for t in classes) and hashlib.sha1(s.encode()).hexdigest() in emb]
        drugv[g] = np.mean(dv, axis=0) if dv else np.zeros(dim, np.float32)
    return allv, drugv


def stage_compare(proteins: dict, emb: dict, model_key: str = "150M", org_key: str = "kpneu") -> dict:
    import yaml
    from scipy.stats import pearsonr
    from sklearn.model_selection import StratifiedGroupKFold
    from xgboost import XGBRegressor

    org = ORGANISMS[org_key]
    taxon, drugs = org["taxon"], org["drugs"]
    cfg = yaml.safe_load((REPO / "config/config.yaml").read_text())
    seed = cfg["project"]["seed"]
    feats = pd.read_csv(REPO / org["features"], dtype={"genome_id": str}).set_index("genome_id")
    lin = pd.read_csv(REPO / org["lineages"], dtype={"genome_id": str}
                      ).set_index("genome_id")["lineage"].to_dict()

    def ea_cv(X, y, groups, repeats=REPEATS):
        """Repeated lineage-grouped CV. Returns (mean EA, std EA, mean r, per-fold EA array).

        The per-fold array is retained so that two feature sets evaluated on IDENTICAL folds can be
        compared with a PAIRED test. Previously only the mean/std were kept, which made the reported
        Wilcoxon p-values impossible to reproduce (audit, decision #60).
        """
        eas, rs = [], []
        yb = (y >= np.median(y)).astype(int)
        for rep in range(repeats):
            skf = StratifiedGroupKFold(5, shuffle=True, random_state=seed + rep)
            for tr, te in skf.split(X, yb, groups):
                m = XGBRegressor(n_estimators=400, max_depth=4, learning_rate=0.08, subsample=0.9,
                                 colsample_bytree=0.8, random_state=seed, n_jobs=4).fit(X[tr], y[tr])
                p = m.predict(X[te])
                eas.append(float(np.mean(np.abs(p - y[te]) <= 1.0)))
                if len(set(y[te])) > 1:
                    rs.append(pearsonr(p, y[te])[0])
        return (float(np.mean(eas)), float(np.std(eas)),
                float(np.mean(rs)) if rs else float("nan"), np.array(eas))

    results = {}
    print(f"\n[{org['name']}]  {'drug':13s}{'n':>6}  {'determ':>16}{'PLM':>16}{'determ+PLM':>16}")
    for drug in drugs:
        mic = fetch_mic_taxon(drug, taxon)
        allv, drugv = genome_plm_vectors(proteins, emb, drug)
        g = [x for x in feats.index if x in mic and x in allv]
        if len(g) < 150:
            print(f"{drug:13s} too few ({len(g)})")
            continue
        Xd = feats.loc[g].values.astype(np.float32)
        Xp = np.vstack([np.concatenate([allv[x], drugv[x]]) for x in g]).astype(np.float32)
        y = np.array([mic[x] for x in g], dtype=np.float32)
        groups = np.array([lin.get(x, x) for x in g])
        base = ea_cv(Xd, y, groups)
        plm = ea_cv(Xp, y, groups)
        both = ea_cv(np.hstack([Xd, Xp]), y, groups)
        # paired Wilcoxon signed-rank on identical folds — the claim the report makes must be
        # computed here, or it cannot be made at all.
        try:
            from scipy.stats import wilcoxon
            wstat, pval = wilcoxon(both[3], base[3])
            pval = float(pval)
        except Exception:
            pval = float("nan")
        results[drug] = {"n": len(g), "n_folds": int(len(base[3])),
                         "determ_ea": base[0], "determ_ea_std": base[1],
                         "plm_ea": plm[0], "plm_ea_std": plm[1],
                         "both_ea": both[0], "both_ea_std": both[1],
                         "delta_both_minus_determ": both[0] - base[0],
                         "wilcoxon_p_both_vs_determ": pval,
                         "fold_ea_determ": base[3].tolist(),
                         "fold_ea_both": both[3].tolist()}
        print(f"{drug:13s}{len(g):>6}  {base[0]:>7.1%}±{base[1]:<6.1%}{plm[0]:>7.1%}±{plm[1]:<6.1%}"
              f"{both[0]:>7.1%}±{both[1]:<6.1%}  Δ={both[0]-base[0]:+.1%} p={pval:.4f}")

    name, _, dim = MODELS[model_key]
    md = [f"# Summary #21 — ESM-2 {model_key} for MIC Prediction · {org['name']}\n",
          f"**Date:** 2026-07-10 · ESM-2 {model_key} (`{name}`, {dim}-d) on MPS · *{org['name']}* "
          "(taxon {}) · 5-fold lineage-grouped CV · metric = **Essential Agreement** (±1 doubling "
          "dilution, CLSI/FDA).\n".format(taxon),
          "Each resistance-gene protein is extracted from the genome (AMRFinderPlus coordinates → "
          "translate) and embedded by ESM-2, so **different alleles of the same gene get different "
          "vectors** — the allele-level signal that presence/absence throws away. Per genome we mean-"
          "pool the embeddings (all AMR proteins + the drug-class subset). MIC is regressed from three "
          "feature sets on identical folds.\n",
          "| Drug | n (MIC) | Determinants (baseline) | ESM-2 PLM only | Determinants + PLM | Δ EA |",
          "|---|---|---|---|---|---|"]
    for d, r in results.items():
        md.append(f"| {d} | {r['n']} | {r['determ_ea']:.1%} | {r['plm_ea']:.1%} | "
                  f"**{r['both_ea']:.1%}** | {r['delta_both_minus_determ']:+.1%} |")
    md.append("\n**Reading.** Δ EA = (determinants + PLM) − determinants. A positive Δ on the hard "
              "drugs (meropenem, cefoxitin) means ESM-2's allele-resolution adds MIC signal beyond "
              "presence/absence. If Δ ≈ 0, the curated determinants already carry the phenotype and "
              "allele sequence adds little — an honest, useful finding either way. PLM-only shows how "
              "much the embeddings carry alone (no gene-identity features).")
    if org_key == "kpneu":
        tag = model_key if model_key != "150M" else ""
        out = OUT if not tag else OUT.with_name(f"summary_21_esm2_mic_{tag}.md")
        mtag = model_key
    else:
        out = OUT.with_name(f"summary_21_esm2_mic_{org_key}_{model_key}.md")
        mtag = f"{org_key}_{model_key}"
    out.write_text("\n".join(md) + "\n")
    (REPO / f"results/metrics/esm2_mic_{mtag}_metrics.json").write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out}")
    return results


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["extract", "embed", "compare", "all"], default="all")
    ap.add_argument("--model", choices=list(MODELS), default="150M")
    ap.add_argument("--organism", choices=list(ORGANISMS), default="kpneu")
    args = ap.parse_args()
    org = ORGANISMS[args.organism]
    prot_pkl = CACHE / org["proteins"]
    feats = pd.read_csv(REPO / org["features"], dtype={"genome_id": str})
    genome_ids = feats["genome_id"].tolist()

    proteins = pickle.loads(prot_pkl.read_bytes()) if prot_pkl.exists() else {}
    ep = emb_path(args.model)
    emb = pickle.loads(ep.read_bytes()) if ep.exists() else {}
    if args.stage in ("extract", "all"):
        proteins = stage_extract(genome_ids, prot_pkl)
    if args.stage in ("embed", "all"):
        emb = stage_embed(proteins, args.model)
    if args.stage in ("compare", "all"):
        stage_compare(proteins, emb, args.model, args.organism)
    return 0


if __name__ == "__main__":
    sys.exit(main())
