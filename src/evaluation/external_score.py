"""External validation, Step 2b — score the FROZEN models on the independent cohort.

Executes the plan pre-registered in `results/reports/summary_32_external_validation_scoping.md`
before any genome was downloaded. Nothing here is fitted:

  * models, calibrators and VME<=3% thresholds are loaded from `results/models/kpneu/*.joblib`
    exactly as shipped — no retraining, no re-tuning, no threshold re-selection;
  * features are rebuilt from our own AMRFinderPlus run and aligned to each model's stored
    `feature_cols`, so determinants absent from the external genomes are honest zeros;
  * all five drugs are reported whatever the numbers;
  * the known-gene **rules baseline** is scored on the same isolates, so "what ML adds" is
    measured externally too.

Two cohorts are reported side by side, because the cohort was de-duplicated against training by
*accession*, which does not remove near-clonal siblings:

  ALL      — every external isolate with a usable label
  ST-NOVEL — only isolates whose MLST sequence type is absent from the training set

The gap between them is the honest headline. Reporting only ALL would reproduce, on external data,
the population-structure leakage this project exists to avoid.

Usage:
    python -m src.evaluation.external_score
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

import joblib
import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

COHORT = ROOT / "results" / "metrics" / "external_cohort_kpneu.csv"
AMR_OUT = ROOT / "data" / "interim" / "external_validation" / "amrfinder"
MLST_OUT = ROOT / "data" / "interim" / "external_validation" / "mlst"
MODELS = ROOT / "results" / "models" / "kpneu"
TRAIN_LINEAGES = ROOT / "data" / "processed" / "thin_slice_cipro_lineages.csv"
OUT_JSON = ROOT / "results" / "metrics" / "external_validation.json"

PANEL_CSV = {  # cohort-CSV column -> model filename stem
    "meropenem": "meropenem",
    "gentamicin": "gentamicin",
    "ciprofloxacin": "ciprofloxacin",
    "trimethoprim-sulfamethoxazole": "trimethoprim_sulfamethoxazole",
    "cefoxitin": "cefoxitin",
}
SEED = 42


def determinants(acc: str) -> set[str] | None:
    p = AMR_OUT / f"{acc}.tsv"
    if not p.exists():
        return None
    dets = set()
    for r in csv.DictReader(p.read_text().splitlines(), delimiter="\t"):
        if (r.get("Type") or "").strip().upper() == "AMR":
            sym = (r.get("Element symbol") or "").strip()
            if sym:
                dets.add(sym)
    return dets


def sequence_type(acc: str) -> str | None:
    p = MLST_OUT / f"{acc}.tsv"
    if not p.exists() or not p.read_text().strip():
        return None
    parts = p.read_text().splitlines()[0].split("\t")
    st = parts[2].strip() if len(parts) > 2 else ""
    return st if st and st != "-" else None


def training_sts() -> set[str]:
    """Sequence types present in training.

    Reads the raw `ST` column, which is either a bare ST number or "-" for untypeable genomes
    (those get a per-genome `novel_<id>` lineage and cannot match an external ST). Using this
    column avoids parsing the `lineage` string, where "ST258" and "novel_573.12908" coexist.

    Being *conservative* matters here: an ST wrongly omitted from this set would let a training
    lineage leak into the ST_NOVEL cohort and inflate the headline number.
    """
    if not TRAIN_LINEAGES.exists():
        raise FileNotFoundError(
            f"{TRAIN_LINEAGES} missing — cannot establish which lineages are novel, and the "
            "ST_NOVEL cohort would be meaningless. Refusing to proceed."
        )
    sts = set()
    for r in csv.DictReader(open(TRAIN_LINEAGES)):
        st = (r.get("ST") or "").strip()
        if st and st != "-":
            sts.add(st)
    return sts


def rules_score(dets: set[str], rules: list[str]) -> int:
    """Transparent baseline: resistant iff any known determinant for this drug is present."""
    return int(any(any(d == g or d.startswith(g) for d in dets) for g in rules))


def clinical(y: np.ndarray, pred: np.ndarray) -> dict:
    """VME = resistant called susceptible; ME = susceptible called resistant."""
    R, S = y == 1, y == 0
    vme = float(((pred == 0) & R).sum() / max(R.sum(), 1))
    me = float(((pred == 1) & S).sum() / max(S.sum(), 1))
    ca = float((pred == y).mean())
    return {"vme": vme, "me": me, "ca": ca,
            "vme_count": int(((pred == 0) & R).sum()), "me_count": int(((pred == 1) & S).sum())}


def boot_ci(y: np.ndarray, s: np.ndarray, groups: np.ndarray, fn, n: int = 2000) -> list[float]:
    """Lineage-clustered bootstrap: resample sequence types, not isolates."""
    rng = np.random.default_rng(SEED)
    uniq = np.unique(groups)
    vals = []
    for _ in range(n):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([np.flatnonzero(groups == g) for g in pick])
        yy = y[idx]
        if len(np.unique(yy)) < 2:
            continue
        try:
            vals.append(fn(yy, s[idx]))
        except Exception:  # noqa: BLE001
            continue
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))] if vals else [float("nan")] * 2


def evaluate(drug_csv: str, stem: str, rows: list[dict], label: str) -> dict | None:
    bundle = joblib.load(MODELS / f"{stem}.joblib")
    cols = bundle["feature_cols"]
    idx = {c: i for i, c in enumerate(cols)}

    X, y, grp, dets_list = [], [], [], []
    for r in rows:
        call = (r["calls"].get(drug_csv) or "").strip().upper()
        if call not in ("R", "S"):
            continue
        v = np.zeros(len(cols), dtype=np.float32)
        for d in r["dets"]:
            j = idx.get(d)
            if j is not None:
                v[j] = 1.0
        X.append(v)
        y.append(1 if call == "R" else 0)
        grp.append(r["st"] or f"_solo_{r['acc']}")
        dets_list.append(r["dets"])
    if len(y) < 20 or len(set(y)) < 2:
        return None

    X = np.vstack(X)
    y = np.asarray(y)
    grp = np.asarray(grp)

    raw = bundle["model"].predict_proba(X)[:, 1]
    cal = bundle["calibrator"].transform(raw)
    pred = (cal >= bundle["threshold"]).astype(int)
    rules_pred = np.array([rules_score(d, bundle["rules"]) for d in dets_list])

    out = {
        "cohort": label, "n": int(len(y)), "n_R": int(y.sum()), "n_lineages": int(len(set(grp))),
        "roc_auc": float(roc_auc_score(y, raw)),
        "roc_auc_ci": boot_ci(y, raw, grp, roc_auc_score),
        "pr_auc": float(average_precision_score(y, raw)),
        "threshold": float(bundle["threshold"]),
        "brier": float(np.mean((cal - y) ** 2)),
        **clinical(y, pred),
        "rules_roc_auc": float(roc_auc_score(y, rules_pred)) if len(set(rules_pred)) > 1 else 0.5,
        "rules": clinical(y, rules_pred),
    }
    out["vme_ci"] = boot_ci(y, cal, grp,
                            lambda yy, ss: clinical(yy, (ss >= bundle["threshold"]).astype(int))["vme"])
    return out


def main() -> int:
    if not COHORT.exists():
        print(f"missing {COHORT} — run `make extscope` first", file=sys.stderr)
        return 1

    train_st = training_sts()
    rows, missing = [], 0
    for r in csv.DictReader(open(COHORT)):
        acc = (r.get("asm_acc") or "").strip()
        dets = determinants(acc)
        if dets is None:
            missing += 1
            continue
        rows.append({
            "acc": acc, "dets": dets, "st": sequence_type(acc),
            "calls": {d: r.get(d, "") for d in PANEL_CSV},
        })
    print(f"external isolates annotated: {len(rows):,}  (missing annotation: {missing:,})")
    print(f"training sequence types: {len(train_st):,}")

    typed = [r for r in rows if r["st"]]
    novel = [r for r in typed if r["st"] not in train_st]
    st_counts = Counter(r["st"] for r in typed)
    print(f"  typed by mlst          : {len(typed):,}")
    print(f"  ST absent from training: {len(novel):,}  "
          f"({100*len(novel)/max(len(typed),1):.1f}% of typed)")
    print(f"  most common STs        : {st_counts.most_common(5)}")

    results = {"meta": {
        "n_isolates": len(rows), "n_typed": len(typed), "n_st_novel": len(novel),
        "n_training_sts": len(train_st), "missing_annotation": missing,
        "top_sts": st_counts.most_common(10),
    }, "drugs": {}}

    for drug_csv, stem in PANEL_CSV.items():
        results["drugs"][stem] = {}
        for label, subset in (("ALL", rows), ("ST_NOVEL", novel)):
            res = evaluate(drug_csv, stem, subset, label)
            if res:
                results["drugs"][stem][label] = res

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    json.dump(results, open(OUT_JSON, "w"), indent=1)

    for label in ("ALL", "ST_NOVEL"):
        print(f"\n=== {label} ===")
        print(f"{'drug':32s} {'n':>5} {'R':>5} {'ROC':>6} {'rules':>6} {'VME':>7} {'ME':>7}")
        for stem in PANEL_CSV.values():
            d = results["drugs"][stem].get(label)
            if not d:
                print(f"{stem:32s}  (insufficient data)")
                continue
            print(f"{stem:32s} {d['n']:5,} {d['n_R']:5,} {d['roc_auc']:6.3f} "
                  f"{d['rules_roc_auc']:6.3f} {d['vme']*100:6.1f}% {d['me']*100:6.1f}%")
    print(f"\nwrote {OUT_JSON.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
