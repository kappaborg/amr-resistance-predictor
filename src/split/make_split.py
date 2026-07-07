"""Phase 5 — lineage assignment (MLST) + phylogeny-aware train/test split.

The project's cardinal rule: NEVER split randomly. Bacteria are clonal, so a random split leaks
lineage identity and inflates every metric. We assign each genome a sequence type (ST) with `mlst`,
then hold out whole STs for the test set so the model must generalize to unseen lineages.
`tests/test_split.py` asserts no ST appears in both train and test.

Novel/untypeable genomes (mlst reports ST '-') each get a unique synthetic lineage id so they are
never pooled into one giant pseudo-lineage (which would leak).

Usage:
    python -m src.split.make_split --config config/config.yaml           # type + split
    python -m src.split.make_split --config config/config.yaml --workers 8
"""
from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import yaml
from sklearn.model_selection import GroupShuffleSplit

REPO = Path(__file__).resolve().parents[2]
MLST = "/opt/homebrew/anaconda3/envs/amr-resistance-predictor/bin/mlst"
GENOMES = REPO / "data/raw/genomes"
CACHE = REPO / "data/interim/mlst"
LABELS = REPO / "data/processed/thin_slice_cipro_labels_qc.csv"
LINEAGES = REPO / "data/processed/thin_slice_cipro_lineages.csv"
SPLIT = REPO / "data/processed/thin_slice_cipro_split.csv"


def genome_ids() -> list[str]:
    return [r["genome_id"] for r in csv.DictReader(LABELS.open())]


def type_one(gid: str) -> tuple[str, str]:
    """Run mlst on one genome; return (gid, ST). Cached. ST '-' means untypeable."""
    out = CACHE / f"{gid}.tsv"
    if out.exists() and out.stat().st_size > 0:
        st = out.read_text().split("\t")[2] if "\t" in out.read_text() else "-"
        return gid, st
    fna = CACHE / f"{gid}.fna"
    fna.write_bytes(gzip.decompress((GENOMES / f"{gid}.fna.gz").read_bytes()))
    r = subprocess.run([MLST, str(fna)], capture_output=True, text=True)
    fna.unlink(missing_ok=True)
    if r.returncode != 0:
        return gid, "-"
    out.write_text(r.stdout)
    parts = r.stdout.split("\t")
    return gid, parts[2] if len(parts) > 2 else "-"


def assign_lineages(ids: list[str], workers: int) -> dict[str, str]:
    CACHE.mkdir(parents=True, exist_ok=True)
    st: dict[str, str] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(type_one, g): g for g in ids}
        for i, fut in enumerate(as_completed(futs), 1):
            gid, s = fut.result()
            st[gid] = s
            if i % 100 == 0 or i == len(ids):
                print(f"  mlst [{i}/{len(ids)}]")
    # untypeable -> unique synthetic lineage so they don't pool into one leaky group
    lineage: dict[str, str] = {}
    for g, s in st.items():
        lineage[g] = f"ST{s}" if s not in ("-", "", None) else f"novel_{g}"
    with LINEAGES.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["genome_id", "ST", "lineage"])
        for g in ids:
            w.writerow([g, st.get(g, "-"), lineage[g]])
    return lineage


def make_split(ids: list[str], lineage: dict[str, str], labels: dict[str, str],
               test_fraction: float, seed: int) -> None:
    groups = [lineage[g] for g in ids]
    y = [labels[g] for g in ids]
    gss = GroupShuffleSplit(n_splits=1, test_size=test_fraction, random_state=seed)
    train_idx, test_idx = next(gss.split(ids, y, groups))
    fold = {ids[i]: "train" for i in train_idx}
    fold.update({ids[i]: "test" for i in test_idx})

    train_lin = {lineage[ids[i]] for i in train_idx}
    test_lin = {lineage[ids[i]] for i in test_idx}
    overlap = train_lin & test_lin
    with SPLIT.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["genome_id", "lineage", "label", "fold"])
        for g in ids:
            w.writerow([g, lineage[g], labels[g], fold[g]])

    from collections import Counter
    tr = Counter(labels[ids[i]] for i in train_idx)
    te = Counter(labels[ids[i]] for i in test_idx)
    print(f"\nLINEAGES: {len(set(groups))} distinct ({sum(1 for v in Counter(groups).values() if v==1)} singletons)")
    print(f"TRAIN: {len(train_idx)} genomes / {len(train_lin)} lineages  R/S={tr.get('Resistant')}/{tr.get('Susceptible')}")
    print(f"TEST : {len(test_idx)} genomes / {len(test_lin)} lineages  R/S={te.get('Resistant')}/{te.get('Susceptible')}")
    print(f"LINEAGE OVERLAP (must be 0): {len(overlap)}")
    if overlap:
        print("  !! LEAKAGE:", sorted(overlap)[:10])
        sys.exit(1)
    print("PASS: no lineage in both train and test.")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    cfg = yaml.safe_load((REPO / args.config).read_text())
    ids = genome_ids()
    labels = {r["genome_id"]: r["ciprofloxacin"] for r in csv.DictReader(LABELS.open())}
    print(f"typing {len(ids)} genomes with mlst...")
    lineage = assign_lineages(ids, args.workers)
    make_split(ids, lineage, labels,
               cfg["split"]["test_fraction"], cfg["project"]["seed"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
