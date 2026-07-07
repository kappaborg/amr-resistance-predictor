"""Live progress for the long-running Phase 4/5 jobs. Run anytime:

    python -m src.progress            # one snapshot
    python -m src.progress --watch    # refresh every 15s until both done
"""
from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
LABELS = REPO / "data/processed/thin_slice_cipro_labels_qc.csv"
AMR = REPO / "data/interim/amrfinder"
MLST = REPO / "data/interim/mlst"


def total() -> int:
    return sum(1 for _ in csv.DictReader(LABELS.open()))


def bar(done: int, tot: int, width: int = 40) -> str:
    frac = done / tot if tot else 0
    filled = int(frac * width)
    return f"[{'█' * filled}{'░' * (width - filled)}] {done}/{tot} ({100*frac:.1f}%)"


def snapshot(tot: int) -> tuple[int, int]:
    a = len(list(AMR.glob("*.tsv"))) if AMR.exists() else 0
    m = len(list(MLST.glob("*.tsv"))) if MLST.exists() else 0
    return a, m


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--watch", action="store_true")
    args = ap.parse_args()
    tot = total()
    prev_a = prev_m = None
    t0 = time.time()
    while True:
        a, m = snapshot(tot)
        # rate/ETA using progress since first sample this session
        line_a = bar(a, tot)
        line_m = bar(m, tot)
        print("\033[2J\033[H" if args.watch else "", end="")
        print(f"Thin-slice progress ({tot} genomes)\n")
        print(f"  AMRFinderPlus (features) {line_a}")
        print(f"  MLST (lineages)          {line_m}")
        if args.watch and prev_a is not None:
            da = a - prev_a
            if da > 0:
                eta = (tot - a) * (15 / da)
                print(f"\n  features rate ~{da/15*60:.0f}/min, ETA ~{eta/3600:.1f}h")
        if a >= tot and m >= tot:
            print("\n  ✅ both complete.")
            return 0
        if not args.watch:
            return 0
        prev_a, prev_m = a, m
        time.sleep(15)


if __name__ == "__main__":
    raise SystemExit(main())
