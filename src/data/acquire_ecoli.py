"""Second-organism generalization — acquire a balanced E. coli genome set.

Panel: ciprofloxacin, gentamicin, trimethoprim/sulfamethoxazole (shared with K. pneumoniae for a
direct cross-organism comparison) + ampicillin, ceftriaxone (E. coli-relevant beta-lactams).

Selects a balanced ciprofloxacin thin slice + resistant top-up for the R-poor drugs (gentamicin,
ceftriaxone), downloads to the shared data/raw/genomes dir (562.* ids don't collide with 573.*),
and writes data/raw/ecoli_ids.txt. Reuses the fixed, retry-safe downloader.

Usage:
    python -m src.data.acquire_ecoli
"""
from __future__ import annotations

import csv
import random
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

from src.data.acquire import download_genome

REPO = Path(__file__).resolve().parents[2]
AMR = "https://www.bv-brc.org/api/genome_amr/"
LAB = 'eq(evidence,%22Laboratory%20Method%22)'
GENOMES = REPO / "data/raw/genomes"
SEED = 42
PANEL = {"ciprofloxacin": "ciprofloxacin", "gentamicin": "gentamicin",
         "trimethoprim_sulfamethoxazole": "trimethoprim%2Fsulfamethoxazole",
         "ampicillin": "ampicillin", "ceftriaxone": "ceftriaxone"}


def consistent(drug_enc: str) -> dict[str, str]:
    q = (f"and(eq(taxon_id,562),{LAB},eq(antibiotic,%22{drug_enc}%22))"
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
    rng = random.Random(SEED)
    labs = {col: consistent(enc) for col, enc in PANEL.items()}
    selected: set[str] = set()

    # 1. balanced ciprofloxacin thin slice (750 R / 750 S)
    cip = labs["ciprofloxacin"]
    for cls in ("Resistant", "Susceptible"):
        ids = sorted(g for g, v in cip.items() if v == cls)
        rng.shuffle(ids)
        selected |= set(ids[:750])
    # 2. resistant top-up for the R-poor drugs; both classes for the small ceftriaxone set
    for col, r_cap, s_cap in [("gentamicin", 600, 0), ("ceftriaxone", 500, 500),
                              ("ampicillin", 400, 0)]:
        for cls, cap in (("Resistant", r_cap), ("Susceptible", s_cap)):
            if not cap:
                continue
            ids = sorted(g for g, v in labs[col].items() if v == cls and g not in selected)
            rng.shuffle(ids)
            selected |= set(ids[:cap])

    ids = sorted(selected)
    (REPO / "data/raw/ecoli_ids.txt").write_text("\n".join(ids) + "\n")
    print(f"selected {len(ids)} E. coli genomes to download")
    for col in PANEL:
        R = sum(1 for g in ids if labs[col].get(g) == "Resistant")
        S = sum(1 for g in ids if labs[col].get(g) == "Susceptible")
        print(f"  {col:32s} R={R:5d} S={S:5d}")

    GENOMES.mkdir(parents=True, exist_ok=True)
    todo = [g for g in ids if not (GENOMES / f"{g}.fna.gz").exists()]
    print(f"\ndownloading {len(todo)} (parallel)...")
    ck = REPO / "data/raw/ecoli_checksums.csv"
    ok = fail = done = 0
    total = 0
    with ck.open("w", newline="") as f, ThreadPoolExecutor(max_workers=6) as ex:
        w = csv.writer(f)
        w.writerow(["genome_id", "bytes_gzipped", "sha256"])
        futs = {ex.submit(download_genome, g, GENOMES / f"{g}.fna.gz"): g for g in todo}
        for fut in as_completed(futs):
            success, nbytes, sha = fut.result()
            done += 1
            if success:
                ok += 1
                total += nbytes
                w.writerow([futs[fut], nbytes, sha])
            else:
                fail += 1
            if done % 100 == 0 or done == len(todo):
                print(f"  [{done}/{len(todo)}] ok={ok} fail={fail} {total/1e9:.2f} GB", flush=True)
    print(f"\nDONE: {ok} downloaded ({total/1e9:.2f} GB), {fail} failed. ids -> data/raw/ecoli_ids.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
