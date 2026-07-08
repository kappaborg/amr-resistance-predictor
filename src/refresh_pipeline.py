"""One-command refresh — rebuild the ENTIRE pipeline on the full genome set (existing + Phase-2b
top-up), so nothing is left at the 1,472-genome stage. Run this once the top-up download finishes.

Steps (each resumable / idempotent):
  1. QC every genome on disk -> master modeling list (data/processed/panel_genomes.csv)
  2. AMRFinderPlus annotate the master list           (only new genomes actually run)
  3. MLST -> lineages for the master list             (only new genomes actually run)
  4. rebuild the determinant feature matrix
  5. rebuild multi-drug panel labels
  6. re-run per-drug models, SHAP, and saved demo models

Downstream files keep their existing paths (…_features.csv, …_lineages.csv, panel_labels.csv), so
train_panel / shap_analysis / save_models automatically pick up the enlarged data.

Usage:
    python -m src.refresh_pipeline --config config/config.yaml [--workers 8]
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from src.data.build_panel_labels import DRUGS, fetch as fetch_labels  # noqa: E402
from src.data.qc_metadata import fetch_qc  # noqa: E402
from src.features.build import annotate_all, build_matrix  # noqa: E402
from src.split.make_split import assign_lineages  # noqa: E402

GENOMES = REPO / "data/raw/genomes"
MASTER = REPO / "data/processed/panel_genomes.csv"
PANEL = REPO / "data/processed/panel_labels.csv"
LINEAGES = REPO / "data/processed/thin_slice_cipro_lineages.csv"


def num(x, k, d):
    v = x.get(k)
    return float(v) if v not in (None, "") else d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    cfg = yaml.safe_load((REPO / args.config).read_text())
    q = cfg["qc"]
    lo, hi = q["genome_length_mbp"]

    # 1. QC every genome on disk -> master modeling list
    ids_disk = sorted(p.name.replace(".fna.gz", "") for p in GENOMES.glob("*.fna.gz"))
    print(f"[1/6] QC on {len(ids_disk)} genomes on disk...")
    qc = {r["genome_id"]: r for r in fetch_qc(ids_disk)}

    def passes(g):
        x = qc.get(g)
        return bool(x) and (num(x, "checkm_completeness", 0) >= q["min_completeness"]
                            and num(x, "checkm_contamination", 99) <= q["max_contamination"]
                            and num(x, "contigs", 1e9) <= q["max_contigs"]
                            and lo <= num(x, "genome_length", 0) / 1e6 <= hi)

    master = [g for g in ids_disk if passes(g)]
    MASTER.parent.mkdir(parents=True, exist_ok=True)
    with MASTER.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["genome_id"])
        w.writerows([[g] for g in master])
    print(f"      master modeling set: {len(master)} genomes (dropped {len(ids_disk)-len(master)} on QC)")

    # 2. annotate (resumable)
    print(f"[2/6] AMRFinderPlus annotation ({args.workers} workers)...")
    annotate_all(master, workers=args.workers, threads=1)

    # 3. MLST lineages (resumable). assign_lineages writes LINEAGES for exactly these ids.
    print("[3/6] MLST lineage typing...")
    assign_lineages(master, workers=max(2, args.workers // 2))

    # 4. feature matrix
    print("[4/6] building feature matrix...")
    build_matrix(master)

    # 5. panel labels for the master set
    print("[5/6] fetching multi-drug labels...")
    per_drug = {col: fetch_labels(enc) for col, enc in DRUGS.items()}
    with PANEL.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["genome_id"] + list(DRUGS))
        for g in master:
            w.writerow([g] + [per_drug[col].get(g, "") for col in DRUGS])
    for col in DRUGS:
        R = sum(1 for g in master if per_drug[col].get(g) == "Resistant")
        S = sum(1 for g in master if per_drug[col].get(g) == "Susceptible")
        print(f"      {col:32s} R={R:5d} S={S:5d}")

    # 6. models + interpretation + saved demo models
    print("[6/6] re-running models, SHAP, saved models...")
    py = sys.executable
    for mod in ("src.models.train_panel", "src.interpret.shap_analysis", "src.models.save_models"):
        print(f"      -> {mod}")
        r = subprocess.run([py, "-m", mod, "--config", args.config], cwd=str(REPO))
        if r.returncode != 0:
            print(f"      ! {mod} failed (rc={r.returncode})")
            return 1

    print("\n✅ refresh complete on the full genome set. Update results/reports/REPORT.md numbers "
          "from results/reports/summary_05_week2_panel.md + summary_06_shap.md.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
