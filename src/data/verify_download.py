"""Phase 2 — integrity verification for a downloaded genome set.

Checks: file count, SHA256 re-match, gzip decompression, FASTA validity, per-genome assembly size
(guards against the API-truncation bug), label completeness, and class balance. Exits non-zero on
any failure so it can gate the pipeline.

Usage:
    python -m src.data.verify_download --tag thin_slice_cipro --label-col ciprofloxacin
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import statistics as st
import sys
from collections import Counter
from pathlib import Path

MIN_ASSEMBLY_MB = 4.0   # K. pneumoniae ~5.3 Mbp; anything well below signals truncation/QC issue


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="thin_slice_cipro")
    ap.add_argument("--label-col", default="ciprofloxacin")
    args = ap.parse_args()
    repo = Path(__file__).resolve().parents[2]
    gdir = repo / "data/raw/genomes"

    files = sorted(gdir.glob("*.fna.gz"))
    ck = {r["genome_id"]: r for r in csv.DictReader((repo / f"data/raw/{args.tag}_checksums.csv").open())}
    lab = {r["genome_id"]: r[args.label_col]
           for r in csv.DictReader((repo / f"data/raw/{args.tag}_labels.csv").open())}

    fails: list[str] = []
    bad_sha = bad_gz = bad_fa = small = 0
    sizes_mb: list[float] = []
    contigs: list[int] = []
    disk_ids = set()

    for f in files:
        gid = f.stem.replace(".fna", "")
        disk_ids.add(gid)
        data = f.read_bytes()
        if gid in ck and hashlib.sha256(data).hexdigest() != ck[gid]["sha256"]:
            bad_sha += 1
        try:
            raw = gzip.decompress(data)
        except Exception:
            bad_gz += 1
            continue
        if not raw.startswith(b">"):
            bad_fa += 1
        mb = len(raw) / 1e6
        sizes_mb.append(mb)
        contigs.append(raw.count(b">"))
        if mb < MIN_ASSEMBLY_MB:
            small += 1

    print(f"file count            : {len(files)}")
    print(f"checksum rows         : {len(ck)}")
    print(f"label rows            : {len(lab)}")
    print(f"SHA256 mismatches     : {bad_sha}")
    print(f"gzip failures         : {bad_gz}")
    print(f"non-FASTA files       : {bad_fa}")
    if sizes_mb:
        print(f"assembly size (Mbp)   : mean={st.mean(sizes_mb):.2f} min={min(sizes_mb):.2f} "
              f"max={max(sizes_mb):.2f}")
        print(f"contigs/genome        : mean={st.mean(contigs):.0f} max={max(contigs)}")
    print(f"undersized (<{MIN_ASSEMBLY_MB} Mbp): {small}")

    unlabeled = disk_ids - set(lab)
    no_checksum = disk_ids - set(ck)
    missing = set(lab) - disk_ids
    dl_balance = Counter(lab[g] for g in disk_ids if g in lab)
    print(f"on disk, unlabeled    : {len(unlabeled)}")
    print(f"on disk, no checksum  : {len(no_checksum)}")
    print(f"labeled, not on disk  : {len(missing)} {sorted(missing) if missing else ''}")
    print(f"downloaded balance    : {dict(dl_balance)}")

    if bad_sha: fails.append(f"{bad_sha} checksum mismatches")
    if bad_gz: fails.append(f"{bad_gz} gzip failures")
    if bad_fa: fails.append(f"{bad_fa} non-FASTA files")
    if small: fails.append(f"{small} undersized assemblies (truncation?)")
    if unlabeled: fails.append(f"{len(unlabeled)} unlabeled genomes on disk")
    if no_checksum: fails.append(f"{len(no_checksum)} genomes without checksum")

    if fails:
        print("\nFAIL: " + "; ".join(fails))
        return 1
    print("\nPASS: all integrity checks clean.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
