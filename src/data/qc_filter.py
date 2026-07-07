"""Phase 3 — apply the QC filter + finalize the modeling label set.

Keeps genomes that are (a) successfully downloaded AND (b) pass all config QC thresholds
(literature-aligned: MIMAG completeness/contamination + K. pneumoniae length range). Writes the
clean per-genome R/S label table that Phase 4/5 consume. Nothing is deleted from disk.

Usage:
    python -m src.data.qc_filter --config config/config.yaml --tag thin_slice_cipro --drug ciprofloxacin
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import yaml


def num(x, k, default):
    v = x.get(k)
    return float(v) if v not in (None, "") else default


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--tag", default="thin_slice_cipro")
    ap.add_argument("--drug", default="ciprofloxacin")
    args = ap.parse_args()
    repo = Path(__file__).resolve().parents[2]
    cfg = yaml.safe_load((repo / args.config).read_text())
    q = cfg["qc"]
    comp_t, cont_t, contig_t = q["min_completeness"], q["max_contamination"], q["max_contigs"]
    len_lo, len_hi = q["genome_length_mbp"]

    qc = {r["genome_id"]: r for r in csv.DictReader((repo / f"data/interim/{args.tag}_qc.csv").open())}
    labels = {r["genome_id"]: r[args.drug]
              for r in csv.DictReader((repo / f"data/raw/{args.tag}_labels.csv").open())
              if r[args.drug] in ("Resistant", "Susceptible")}
    downloaded = {p.stem.replace(".fna", "") for p in (repo / "data/raw/genomes").glob("*.fna.gz")}

    def passes(g):
        x = qc.get(g, {})
        return (num(x, "checkm_completeness", 0) >= comp_t
                and num(x, "checkm_contamination", 99) <= cont_t
                and num(x, "contigs", 1e9) <= contig_t
                and len_lo <= num(x, "genome_length", 0) / 1e6 <= len_hi)

    kept, drop_qc, drop_missing = [], [], []
    for g, lab in labels.items():
        if g not in downloaded:
            drop_missing.append(g)
        elif passes(g):
            kept.append((g, lab))
        else:
            drop_qc.append(g)

    out = repo / f"data/processed/{args.tag}_labels_qc.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["genome_id", args.drug])
        w.writerows(sorted(kept))

    bal = Counter(l for _, l in kept)
    print(f"labeled genomes (R/S)     : {len(labels)}")
    print(f"dropped: not downloaded   : {len(drop_missing)} {drop_missing}")
    print(f"dropped: failed QC        : {len(drop_qc)}")
    print(f"KEPT for modeling         : {len(kept)}")
    print(f"  class balance           : {dict(bal)} "
          f"(minority frac {min(bal.values())/sum(bal.values()):.2f})")
    print(f"QC thresholds applied     : completeness>={comp_t}%, contamination<={cont_t}%, "
          f"contigs<={contig_t}, length {len_lo}-{len_hi} Mbp")
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
