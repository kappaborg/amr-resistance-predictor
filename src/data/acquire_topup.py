"""Phase 2b — targeted top-up: pull additional RESISTANT genomes for the data-limited drugs
(meropenem, gentamicin, cefoxitin) plus cefoxitin-susceptible (its set is smallest). Stabilises the
VME operating points without pulling the full 12 GB panel.

Downloads only genomes we don't already have. Reuses the fixed downloader (limit(25000) + retry +
size guard). Writes data/raw/topup_ids.txt for downstream steps.

Usage:
    python -m src.data.acquire_topup
"""
from __future__ import annotations

import csv
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

from src.data.acquire import download_genome

REPO = Path(__file__).resolve().parents[2]
AMR = "https://www.bv-brc.org/api/genome_amr/"
LAB = 'eq(evidence,%22Laboratory%20Method%22)'
GENOMES = REPO / "data/raw/genomes"
# meropenem/gentamicin: add Resistant only (they are S-heavy). cefoxitin: add both (smallest set).
NEED = {"meropenem": ["Resistant"], "gentamicin": ["Resistant"],
        "cefoxitin": ["Resistant", "Susceptible"]}


def consistent(drug: str) -> dict[str, str]:
    q = (f"and(eq(taxon_id,573),{LAB},eq(antibiotic,%22{drug}%22))"
         f"&select(genome_id,resistant_phenotype)&limit(25000)")
    r = requests.get(f"{AMR}?{q}", headers={"Accept": "application/json"}, timeout=180)
    r.raise_for_status()
    byg: dict[str, set] = defaultdict(set)
    for x in r.json():
        p = x.get("resistant_phenotype")
        if p in ("Resistant", "Susceptible"):
            byg[x["genome_id"]].add(p)
    return {g: next(iter(v)) for g, v in byg.items() if len(v) == 1}


def main() -> int:
    have = {r["genome_id"] for r in
            csv.DictReader((REPO / "data/processed/thin_slice_cipro_labels_qc.csv").open())}
    candidates: set[str] = set()
    for drug, classes in NEED.items():
        lab = consistent(drug)
        for cls in classes:
            new = {g for g, v in lab.items() if v == cls and g not in have}
            candidates |= new
            print(f"{drug} {cls}: +{len(new)} new")
    candidates -= {p.stem.replace(".fna", "") for p in GENOMES.glob("*.fna.gz")}  # skip already downloaded
    ids = sorted(candidates)
    print(f"\ntop-up candidates to download: {len(ids)}")
    (REPO / "data/raw/topup_ids.txt").write_text("\n".join(ids) + "\n")

    GENOMES.mkdir(parents=True, exist_ok=True)
    # BV-BRC is slow under load (~11s/genome); downloads are I/O-bound so parallelise with threads.
    ck = REPO / "data/raw/topup_checksums.csv"
    total = ok = fail = done = 0
    with ck.open("w", newline="") as f, ThreadPoolExecutor(max_workers=6) as ex:
        w = csv.writer(f)
        w.writerow(["genome_id", "bytes_gzipped", "sha256"])
        futs = {ex.submit(download_genome, g, GENOMES / f"{g}.fna.gz"): g for g in ids}
        for fut in as_completed(futs):
            success, nbytes, sha = fut.result()
            done += 1
            if success:
                ok += 1
                total += nbytes
                w.writerow([futs[fut], nbytes, sha])
            else:
                fail += 1
            if done % 100 == 0 or done == len(ids):
                print(f"  [{done}/{len(ids)}] ok={ok} fail={fail} {total/1e9:.2f} GB", flush=True)
    print(f"\nDONE: {ok} downloaded ({total/1e9:.2f} GB), {fail} failed. ids -> data/raw/topup_ids.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
