"""Phase 11 — the demo: genome in -> per-drug resistant/susceptible + calibrated confidence +
the determinants behind each call ("why this strain resists drug X").

Runs AMRFinderPlus on the input genome (or reuses a cached annotation for a known genome_id),
builds the determinant feature vector, and for each drug reports the calibrated probability, the
R/S call at the clinical VME<=target operating threshold, and the SHAP-ranked determinants that
drove the call.

Usage:
    python -m src.app.predict --genome path/to/genome.fna[.gz]
    python -m src.app.predict --genome-id 573.12771     # reuse cached annotation (instant)
"""
from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
from pathlib import Path

import joblib
import numpy as np
import shap

REPO = Path(__file__).resolve().parents[2]
ENV_BIN = "/opt/homebrew/anaconda3/envs/amr-resistance-predictor/bin"
AMR = f"{ENV_BIN}/amrfinder"
MODELS = REPO / "results/models"
CACHE = REPO / "data/interim/amrfinder"
GENOMES = REPO / "data/raw/genomes"
DRUG_ORDER = ["meropenem", "gentamicin", "ciprofloxacin", "trimethoprim_sulfamethoxazole", "cefoxitin"]


def determinants_from_tsv(tsv_text: str) -> set[str]:
    dets = set()
    for r in csv.DictReader(tsv_text.splitlines(), delimiter="\t"):
        if (r.get("Type") or "").strip().upper() == "AMR":
            sym = (r.get("Element symbol") or "").strip()
            if sym:
                dets.add(sym)
    return dets


def annotate(genome_path: Path) -> set[str]:
    """Run AMRFinderPlus on a genome FASTA (decompressing .gz if needed)."""
    tmp = None
    fna = genome_path
    if genome_path.suffix == ".gz":
        tmp = CACHE / f"_demo_{genome_path.stem}"
        tmp.parent.mkdir(parents=True, exist_ok=True)
        tmp.write_bytes(gzip.decompress(genome_path.read_bytes()))
        fna = tmp
    print(f"annotating {genome_path.name} with AMRFinderPlus (~1-3 min)...", flush=True)
    r = subprocess.run([AMR, "-n", str(fna), "--organism", "Klebsiella_pneumoniae",
                        "--threads", "4"], capture_output=True, text=True)
    if tmp:
        tmp.unlink(missing_ok=True)
    if r.returncode != 0:
        sys.exit(f"AMRFinderPlus failed: {r.stderr[-300:]}")
    return determinants_from_tsv(r.stdout)


def load_determinants(genome_id: str | None, genome: str | None) -> tuple[str, set[str]]:
    if genome_id:
        cached = CACHE / f"{genome_id}.tsv"
        if cached.exists():
            return genome_id, determinants_from_tsv(cached.read_text())
        gz = GENOMES / f"{genome_id}.fna.gz"
        if gz.exists():
            return genome_id, annotate(gz)
        sys.exit(f"genome_id {genome_id} not found in cache or data/raw/genomes.")
    p = Path(genome)
    if not p.exists():
        sys.exit(f"genome file not found: {p}")
    return p.name, annotate(p)


def explain(bundle, x_row, present_dets, top_n=4):
    """SHAP-ranked determinants present in THIS genome that drove the call."""
    model, cols = bundle["model"], bundle["feature_cols"]
    sv = shap.TreeExplainer(model).shap_values(x_row.reshape(1, -1))[0]
    present_idx = [i for i, c in enumerate(cols) if c in present_dets]
    ranked = sorted(present_idx, key=lambda i: -abs(sv[i]))
    out = []
    for i in ranked[:top_n]:
        if abs(sv[i]) < 1e-6:
            continue
        out.append((cols[i], "↑resistant" if sv[i] > 0 else "↓susceptible", float(sv[i])))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Predict per-drug resistance from a K. pneumoniae genome.")
    ap.add_argument("--genome", help="path to genome FASTA (.fna or .fna.gz)")
    ap.add_argument("--genome-id", help="BV-BRC genome_id already downloaded/annotated")
    args = ap.parse_args()
    if not (args.genome or args.genome_id):
        sys.exit("provide --genome or --genome-id")

    bundles = {d: joblib.load(MODELS / f"{d}.joblib") for d in DRUG_ORDER
               if (MODELS / f"{d}.joblib").exists()}
    if not bundles:
        sys.exit("no models found — run: python -m src.models.save_models")

    name, dets = load_determinants(args.genome_id, args.genome)
    cols = next(iter(bundles.values()))["feature_cols"]
    x = np.array([1 if c in dets else 0 for c in cols], dtype=int)
    known = [d for d in dets if d in set(cols)]

    print("\n" + "=" * 68)
    print(f"  RESISTANCE PREDICTION — {name}  (Klebsiella pneumoniae)")
    print(f"  {len(dets)} resistance determinants detected ({len(known)} in model vocabulary)")
    print("=" * 68)
    for drug in DRUG_ORDER:
        if drug not in bundles:
            continue
        b = bundles[drug]
        raw = b["model"].predict_proba(x.reshape(1, -1))[0, 1]
        prob = float(b["calibrator"].transform([raw])[0])
        call = "RESISTANT" if prob >= b["threshold"] else "susceptible"
        # certainty from how decisive the calibrated probability is (not from the threshold)
        flag = "" if (prob >= 0.8 or prob <= 0.2) else "  [uncertain]"
        drivers = explain(b, x, dets)
        print(f"\n  {drug.replace('_', '/'):30s} {call:11s}  P(resistant)={prob:.2f}{flag}")
        if drivers:
            print("     drivers: " + ", ".join(f"{d} {arrow}" for d, arrow, _ in drivers))
        else:
            print("     drivers: none of this genome's determinants are model-relevant for this drug")
    print("\n" + "-" * 68)
    print("  Operating point: threshold tuned to cap very-major error at the clinical target.")
    print("  Research/decision-support only — not a substitute for phenotypic testing.")
    print("=" * 68 + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
