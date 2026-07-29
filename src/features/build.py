"""Phase 4 — feature extraction: AMRFinderPlus -> determinant feature matrix.

Runs AMRFinderPlus on each QC-passed genome (with --organism Klebsiella_pneumoniae so gyrA/parC
and other point mutations are called), caches the per-genome TSV (resumable), then assembles a
binary genome x determinant matrix. Determinants = AMRFinderPlus "Element symbol" for Type==AMR
rows (acquired genes + point mutations); VIRULENCE/STRESS elements are excluded.

Usage:
    python -m src.features.build --workers 8 --threads 1            # full run (resumable)
    python -m src.features.build --workers 8 --threads 1 --limit 8  # timing test on 8 genomes
    python -m src.features.build --matrix-only                      # rebuild matrix from cached TSVs
"""
from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
AMR = "/opt/homebrew/anaconda3/envs/amr-resistance-predictor/bin/amrfinder"
ORGANISM = "Klebsiella_pneumoniae"
GENOMES = REPO / "data/raw/genomes"
CACHE = REPO / "data/interim/amrfinder"
LABELS = REPO / "data/processed/thin_slice_cipro_labels_qc.csv"
MATRIX = REPO / "data/processed/thin_slice_cipro_features.csv"


def genome_ids() -> list[str]:
    return [r["genome_id"] for r in csv.DictReader(LABELS.open())]


def annotate_one(gid: str, threads: int) -> tuple[str, int, float]:
    """Annotate a single genome; cache TSV. Returns (gid, returncode, seconds). Resumable."""
    tsv = CACHE / f"{gid}.tsv"
    if tsv.exists() and tsv.stat().st_size > 0:
        return gid, 0, 0.0  # already done
    gz = GENOMES / f"{gid}.fna.gz"
    fna = CACHE / f"{gid}.fna"
    fna.write_bytes(gzip.decompress(gz.read_bytes()))
    t = time.time()
    r = subprocess.run([AMR, "-n", str(fna), "--organism", ORGANISM, "--threads", str(threads),
                        "-o", str(tsv)], capture_output=True, text=True)
    fna.unlink(missing_ok=True)
    if r.returncode != 0:
        tsv.unlink(missing_ok=True)  # don't cache a failed run
        return gid, r.returncode, time.time() - t
    return gid, 0, time.time() - t


def annotate_all(ids: list[str], workers: int, threads: int) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    todo = [g for g in ids if not ((CACHE / f"{g}.tsv").exists() and (CACHE / f"{g}.tsv").stat().st_size > 0)]
    print(f"annotating {len(todo)}/{len(ids)} genomes ({len(ids)-len(todo)} cached) "
          f"@ {workers} workers x {threads} threads")
    done = fail = 0
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(annotate_one, g, threads): g for g in todo}
        for fut in as_completed(futs):
            gid, rc, el = fut.result()
            if rc == 0:
                done += 1
            else:
                fail += 1
                print(f"  ! {gid} failed rc={rc}")
            if (done + fail) % 25 == 0 or (done + fail) == len(todo):
                elapsed = time.time() - t0
                rate = elapsed / max(done + fail, 1)
                remain = rate * (len(todo) - done - fail)
                print(f"  [{done+fail}/{len(todo)}] ok={done} fail={fail} "
                      f"elapsed={elapsed/60:.1f}m rate={rate:.0f}s/genome ETA={remain/60:.0f}m")


def parse_determinants(tsv: Path) -> set[str]:
    """Element symbols for Type==AMR rows (acquired genes + point mutations)."""
    dets: set[str] = set()
    rows = list(csv.DictReader(tsv.open(), delimiter="\t"))
    for r in rows:
        if (r.get("Type") or "").strip().upper() == "AMR":
            sym = (r.get("Element symbol") or "").strip()
            if sym:
                dets.add(sym)
    return dets


def encode_matrix(per_genome: dict[str, set[str]],
                  ids: list[str]) -> tuple[list[str], list[tuple[str, list[int]]]]:
    """Pure, testable core: build the binary genome x determinant matrix.

    Columns = sorted union of all determinants seen. Returns (cols, [(genome_id, row_of_0/1), ...])
    for every id present in `per_genome`, in `ids` order. No file I/O.
    """
    all_dets: set[str] = set()
    for d in per_genome.values():
        all_dets |= d
    cols = sorted(all_dets)
    rows = [(g, [1 if c in per_genome[g] else 0 for c in cols]) for g in ids if g in per_genome]
    return cols, rows


def build_matrix(ids: list[str]) -> None:
    per_genome: dict[str, set[str]] = {}
    missing = 0
    for g in ids:
        tsv = CACHE / f"{g}.tsv"
        if not tsv.exists():
            missing += 1
            continue
        per_genome[g] = parse_determinants(tsv)
    cols, rows = encode_matrix(per_genome, ids)
    with MATRIX.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["genome_id"] + cols)
        for g, row in rows:
            w.writerow([g] + row)
    print(f"\nmatrix: {len(rows)} genomes x {len(cols)} determinants "
          f"({missing} genomes missing TSV) -> {MATRIX}")
    # quick view of the most common determinants
    freq = {c: sum(1 for d in per_genome.values() if c in d) for c in cols}
    top = sorted(freq.items(), key=lambda kv: -kv[1])[:12]
    print("most common determinants:")
    for c, n in top:
        print(f"  {c:24s} {n}/{len(per_genome)} ({100*n/len(per_genome):.0f}%)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0, help="annotate only first N (timing test)")
    ap.add_argument("--matrix-only", action="store_true", help="skip annotation, rebuild matrix")
    args = ap.parse_args()

    ids = genome_ids()
    if args.limit:
        ids = ids[:args.limit]
    if not args.matrix_only:
        annotate_all(ids, args.workers, args.threads)
    build_matrix(ids)
    return 0


if __name__ == "__main__":
    sys.exit(main())
