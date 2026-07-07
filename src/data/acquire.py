"""Phase 2 — BV-BRC data acquisition (phenotypes + assembled genomes).

Thin-slice mode (Week-1): balanced ciprofloxacin subsample of K. pneumoniae, download genomes,
label table, checksums, manifest. Full mode (Week-2, after the gate) pulls the whole panel.

Verified API contract (Phase 1/2 probes):
  - phenotypes: genome_amr, lab only (evidence="Laboratory Method").
  - genome FASTA: genome_sequence API, http_accept=application/dna+fasta, one file per genome_id.

Resumable (skips genomes already on disk), rate-limited, SHA256 per file. Never fabricates: failed
downloads are logged and skipped, subsample size is recorded — no silent capping.

Usage:
    python -m src.data.acquire --config config/config.yaml --thin-slice
    python -m src.data.acquire --config config/config.yaml --thin-slice --per-class 750
    python -m src.data.acquire --config config/config.yaml --full        # after Week-1 gate
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import random
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import requests
import yaml

AMR = "https://www.bv-brc.org/api/genome_amr/"
SEQ = "https://www.bv-brc.org/api/genome_sequence/"
LAB = 'eq(evidence,%22Laboratory%20Method%22)'
PAUSE_S = 0.15


def fetch_phenotypes(taxon_id: int, antibiotic: str) -> dict[str, str]:
    """Return {genome_id: 'Resistant'|'Susceptible'} keeping only genomes with a single consistent label."""
    ab = antibiotic.replace(" ", "%20").replace("/", "%2F")
    q = (f"and(eq(taxon_id,{taxon_id}),{LAB},eq(antibiotic,%22{ab}%22))"
         f"&select(genome_id,resistant_phenotype)&limit(25000)")
    r = requests.get(f"{AMR}?{q}", headers={"Accept": "application/json"}, timeout=180)
    r.raise_for_status()
    byg: dict[str, set] = defaultdict(set)
    for x in r.json():
        p = x.get("resistant_phenotype")
        if p in ("Resistant", "Susceptible"):
            byg[x["genome_id"]].add(p)
    return {g: next(iter(v)) for g, v in byg.items() if len(v) == 1}


def balanced_sample(labels: dict[str, str], per_class: int, seed: int) -> dict[str, str]:
    rng = random.Random(seed)
    out: dict[str, str] = {}
    for cls in ("Resistant", "Susceptible"):
        ids = sorted(g for g, l in labels.items() if l == cls)
        rng.shuffle(ids)
        take = ids[:per_class]
        if len(take) < per_class:
            print(f"  ! only {len(take)} {cls} available (< {per_class} requested) — taking all")
        out.update({g: cls for g in take})
    return out


def download_genome(genome_id: str, dest: Path) -> tuple[bool, int, str]:
    """Download one genome's contigs FASTA, gzip to dest. Returns (ok, bytes_gzipped, sha256)."""
    if dest.exists():  # resumable: skip already-downloaded
        data = dest.read_bytes()
        return True, len(data), hashlib.sha256(data).hexdigest()
    # limit(25000) is REQUIRED: the API defaults to 25 rows (contigs), silently truncating
    # multi-contig WGS assemblies to ~2 Mbp. Without it, genomes are incomplete (verified 2026-07-07).
    url = f"{SEQ}?eq(genome_id,{genome_id})&limit(25000)&http_accept=application/dna+fasta"
    # retry transient network failures (timeouts/resets) so one bad genome can't abort the run.
    content = b""
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=120)
            if r.status_code == 200 and r.content.startswith(b">"):
                content = r.content
                break
        except requests.RequestException:
            pass
        time.sleep(2 * (attempt + 1))  # backoff: 2s, 4s, 6s
    if not content:
        return False, 0, ""
    # sanity guard: a real K. pneumoniae assembly is ~5 Mbp; reject implausibly small pulls
    # (catches API truncation before it reaches the feature matrix).
    if len(content) < 2_500_000:
        print(f"  ! {genome_id}: only {len(content)/1e6:.2f} MB FASTA — likely truncated, skipping")
        return False, 0, ""
    gz = gzip.compress(content)
    dest.write_bytes(gz)
    time.sleep(PAUSE_S)
    return True, len(gz), hashlib.sha256(gz).hexdigest()


def run(cfg: dict, thin_slice: bool, per_class: int) -> int:
    repo = Path(__file__).resolve().parents[2]
    taxon = cfg["organism"]["taxon_id"]
    seed = cfg["project"]["seed"]
    max_gb = cfg["data"]["max_download_gb"]
    drugs = [cfg["thin_slice_drug"]] if thin_slice else cfg["drugs"]
    tag = "thin_slice_cipro" if thin_slice else "full_panel"

    gdir = repo / "data/raw/genomes"
    gdir.mkdir(parents=True, exist_ok=True)

    # 1. phenotypes + selection
    selected: dict[str, dict[str, str]] = {}  # genome_id -> {drug: label}
    for drug in drugs:
        labels = fetch_phenotypes(taxon, drug)
        print(f"{drug}: {len(labels)} consistent-label genomes {dict(Counter(labels.values()))}")
        chosen = balanced_sample(labels, per_class, seed) if thin_slice else labels
        for g, l in chosen.items():
            selected.setdefault(g, {})[drug] = l
    genome_ids = sorted(selected)
    print(f"\nselected {len(genome_ids)} distinct genomes for download ({tag})")

    # 2. write label table
    lab_path = repo / f"data/raw/{tag}_labels.csv"
    with lab_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["genome_id"] + drugs)
        for g in genome_ids:
            w.writerow([g] + [selected[g].get(d, "") for d in drugs])

    # 3. download genomes (resumable, checksummed, size-guarded)
    total = 0
    ok = fail = 0
    ck_path = repo / f"data/raw/{tag}_checksums.csv"
    with ck_path.open("w", newline="") as ckf:
        ckw = csv.writer(ckf)
        ckw.writerow(["genome_id", "bytes_gzipped", "sha256"])
        for i, g in enumerate(genome_ids, 1):
            success, nbytes, sha = download_genome(g, gdir / f"{g}.fna.gz")
            if success:
                ok += 1
                total += nbytes
                ckw.writerow([g, nbytes, sha])
            else:
                fail += 1
                print(f"  ! download failed: {g}")
            if i % 100 == 0 or i == len(genome_ids):
                print(f"  [{i}/{len(genome_ids)}] ok={ok} fail={fail} total={total/1e9:.2f} GB")
            if total / 1e9 > max_gb:
                print(f"\n!! hit max_download_gb={max_gb} GB at {i} genomes — stopping. "
                      f"Raise the ceiling in config to continue.")
                break

    # 4. manifest row
    print(f"\nDONE: {ok} genomes ({total/1e9:.2f} GB gzipped), {fail} failed.")
    print(f"labels -> {lab_path}\nchecksums -> {ck_path}\ngenomes -> {gdir}/")
    print("\n>>> Add to data/manifest.md: BV-BRC K. pneumoniae "
          f"{tag}, {ok} genomes, {total/1e9:.2f} GB, seed={seed}, date 2026-07-07")
    return 0 if fail == 0 else 1


def main() -> int:
    ap = argparse.ArgumentParser(description="BV-BRC acquisition (phenotypes + genomes).")
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--thin-slice", action="store_true", help="Week-1 balanced cipro subsample")
    ap.add_argument("--full", action="store_true", help="full 5-drug panel (after gate)")
    ap.add_argument("--per-class", type=int, default=750, help="genomes per class for thin slice")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[2]
    cfg_path = Path(args.config) if args.config.startswith("/") else repo / args.config
    cfg = yaml.safe_load(cfg_path.read_text())

    if not (args.thin_slice or args.full):
        print("Specify --thin-slice or --full.")
        return 2
    return run(cfg, thin_slice=args.thin_slice, per_class=args.per_class)


if __name__ == "__main__":
    sys.exit(main())
