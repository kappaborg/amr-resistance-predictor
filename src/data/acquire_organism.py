"""General organism acquisition — survey drug balance, select a balanced genome set, download.

Reuses the fixed retry-safe downloader. Config-driven per organism (S. aureus, A. baumannii).
Selects a balanced primary-drug thin slice + resistant top-up for the R-poor drugs, so every panel
drug has enough R/S. Writes data/raw/<tag>_ids.txt.

    python -m src.data.acquire_organism --organism saureus
    python -m src.data.acquire_organism --organism abaumannii
"""
from __future__ import annotations

import argparse
import csv
import random
import sys
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

# organism -> (taxon, primary drug, panel drugs).  All drug names url-safe (no '/').
ORG_ACQ = {
    "saureus": {"taxon": 1280, "primary": "ciprofloxacin",
                "panel": ["oxacillin", "cefoxitin", "ciprofloxacin", "erythromycin", "clindamycin"]},
    "abaumannii": {"taxon": 470, "primary": "ciprofloxacin",
                   "panel": ["meropenem", "imipenem", "ciprofloxacin", "gentamicin", "amikacin"]},
}


def consistent(taxon: int, drug: str) -> dict[str, str]:
    q = (f"and(eq(taxon_id,{taxon}),{LAB},eq(antibiotic,%22{drug}%22))"
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True, choices=list(ORG_ACQ))
    ap.add_argument("--per-class", type=int, default=750)
    args = ap.parse_args()
    cfg = ORG_ACQ[args.organism]
    tag = args.organism
    rng = random.Random(SEED)

    labs = {d: consistent(cfg["taxon"], d) for d in cfg["panel"]}
    print(f"[{tag}] per-drug consistent-label genomes:")
    for d in cfg["panel"]:
        c = {v: sum(1 for x in labs[d].values() if x == v) for v in ("Resistant", "Susceptible")}
        print(f"    {d:16s} R={c['Resistant']:5d} S={c['Susceptible']:5d}")

    selected: set[str] = set()
    # balanced primary thin slice
    prim = labs[cfg["primary"]]
    for cls in ("Resistant", "Susceptible"):
        ids = sorted(g for g, v in prim.items() if v == cls)
        rng.shuffle(ids)
        selected |= set(ids[:args.per_class])
    # resistant top-up for the other drugs (+ some S for small ones)
    for d in cfg["panel"]:
        if d == cfg["primary"]:
            continue
        for cls, cap in (("Resistant", 500), ("Susceptible", 300)):
            ids = sorted(g for g, v in labs[d].items() if v == cls and g not in selected)
            rng.shuffle(ids)
            selected |= set(ids[:cap])

    ids = sorted(selected)
    (REPO / f"data/raw/{tag}_ids.txt").write_text("\n".join(ids) + "\n")
    print(f"\n[{tag}] selected {len(ids)} genomes")
    for d in cfg["panel"]:
        R = sum(1 for g in ids if labs[d].get(g) == "Resistant")
        S = sum(1 for g in ids if labs[d].get(g) == "Susceptible")
        print(f"    {d:16s} R={R:5d} S={S:5d}")

    GENOMES.mkdir(parents=True, exist_ok=True)
    todo = [g for g in ids if not (GENOMES / f"{g}.fna.gz").exists()]
    print(f"\ndownloading {len(todo)} (parallel)...")
    ck = REPO / f"data/raw/{tag}_checksums.csv"
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
    print(f"\nDONE: {ok} downloaded ({total/1e9:.2f} GB), {fail} failed. ids -> data/raw/{tag}_ids.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
