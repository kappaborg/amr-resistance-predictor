"""Phase 3 prep — fetch BV-BRC CheckM QC metadata for the acquired genomes and analyze
threshold sensitivity, so QC cutoffs are validated against the real data distribution rather
than guessed. Writes data/interim/<tag>_qc.csv.

Usage:
    python -m src.data.qc_metadata --tag thin_slice_cipro
"""
from __future__ import annotations

import argparse
import csv
import statistics as st
import sys
from collections import Counter
from pathlib import Path

import requests

GEN = "https://www.bv-brc.org/api/genome/"
FIELDS = ["genome_id", "genome_status", "genome_quality", "genome_length", "contigs",
          "contig_n50", "checkm_completeness", "checkm_contamination",
          "coarse_consistency", "fine_consistency", "gc_content", "cds"]


def fetch_qc(ids: list[str]) -> list[dict]:
    out: list[dict] = []
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        q = "in(genome_id,(" + ",".join(chunk) + f"))&select({','.join(FIELDS)})&limit(200)"
        r = requests.get(f"{GEN}?{q}", headers={"Accept": "application/json"}, timeout=120)
        r.raise_for_status()
        out.extend(r.json())
    return out


def pct(vals, p):
    vals = sorted(vals)
    return vals[min(len(vals) - 1, int(p / 100 * len(vals)))]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="thin_slice_cipro")
    args = ap.parse_args()
    repo = Path(__file__).resolve().parents[2]

    ids = [r["genome_id"] for r in csv.DictReader((repo / f"data/raw/{args.tag}_labels.csv").open())]
    recs = fetch_qc(ids)
    outp = repo / f"data/interim/{args.tag}_qc.csv"
    outp.parent.mkdir(parents=True, exist_ok=True)
    with outp.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for x in recs:
            w.writerow({k: x.get(k, "") for k in FIELDS})
    print(f"fetched QC for {len(recs)}/{len(ids)} genomes -> {outp}\n")

    def col(name):
        return [float(x[name]) for x in recs if x.get(name) not in (None, "")]

    comp, cont = col("checkm_completeness"), col("checkm_contamination")
    contigs, n50, length = col("contigs"), col("contig_n50"), col("genome_length")
    fine = col("fine_consistency")
    print("QUALITY LABEL (BV-BRC):", dict(Counter(x.get("genome_quality") for x in recs)))
    print("STATUS:", dict(Counter(x.get("genome_status") for x in recs)))
    print(f"\ncompleteness  : min={min(comp):.1f} p5={pct(comp,5):.1f} p50={pct(comp,50):.1f} max={max(comp):.1f}")
    print(f"contamination : max={max(cont):.1f} p95={pct(cont,95):.1f} p50={pct(cont,50):.1f} min={min(cont):.1f}")
    print(f"contigs       : min={min(contigs):.0f} p50={pct(contigs,50):.0f} p95={pct(contigs,95):.0f} max={max(contigs):.0f}")
    print(f"contig_n50    : min={min(n50):.0f} p5={pct(n50,5):.0f} p50={pct(n50,50):.0f}")
    print(f"genome_length : min={min(length)/1e6:.2f} p50={pct(length,50)/1e6:.2f} max={max(length)/1e6:.2f} Mbp")
    print(f"fine_consist. : min={min(fine):.1f} p5={pct(fine,5):.1f} p50={pct(fine,50):.1f}")

    def num(x, k, default):
        """Safe float: treats missing/'' as default but preserves a real 0 (avoids `0 or d` trap)."""
        v = x.get(k)
        return float(v) if v not in (None, "") else default

    print("\n== threshold sensitivity (genomes DROPPED / kept of {}) ==".format(len(recs)))
    for comp_t, cont_t, contig_t, len_lo, len_hi in [
        (90, 5, 500, 4.5, 7.0), (95, 5, 400, 4.8, 6.5), (98, 3, 300, 5.0, 6.2),
    ]:
        kept = [x for x in recs
                if num(x, "checkm_completeness", 0) >= comp_t
                and num(x, "checkm_contamination", 99) <= cont_t
                and num(x, "contigs", 1e9) <= contig_t
                and len_lo <= num(x, "genome_length", 0) / 1e6 <= len_hi]
        print(f"  comp>={comp_t} cont<={cont_t} contigs<={contig_t} len[{len_lo}-{len_hi}]Mbp "
              f"-> kept {len(kept)}, dropped {len(recs)-len(kept)}")

    # also report BV-BRC's own 'Good' label + species length sanity as a simple, calibrated baseline
    good_len = [x for x in recs if x.get("genome_quality") == "Good"
                and 4.8 <= num(x, "genome_length", 0) / 1e6 <= 6.5]
    print(f"  BV-BRC 'Good' AND len[4.8-6.5]Mbp -> kept {len(good_len)}, dropped {len(recs)-len(good_len)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
