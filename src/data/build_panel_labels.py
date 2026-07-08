"""Phase 2b/7 — build the multi-drug label table for the already-annotated genome set.

Reuses the 1472 genomes we already downloaded + annotated (no new compute). Fetches lab R/S for
all 5 panel drugs, keeps consistent labels, writes one row per genome with a column per drug.

Usage:
    python -m src.data.build_panel_labels
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import requests

REPO = Path(__file__).resolve().parents[2]
AMR = "https://www.bv-brc.org/api/genome_amr/"
LAB = 'eq(evidence,%22Laboratory%20Method%22)'
DRUGS = {  # column name -> BV-BRC antibiotic (url-encoded)
    "meropenem": "meropenem",
    "gentamicin": "gentamicin",
    "ciprofloxacin": "ciprofloxacin",
    "trimethoprim_sulfamethoxazole": "trimethoprim%2Fsulfamethoxazole",
    "cefoxitin": "cefoxitin",
}


def fetch(drug_enc: str) -> dict[str, str]:
    q = (f"and(eq(taxon_id,573),{LAB},eq(antibiotic,%22{drug_enc}%22))"
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
    genomes = [r["genome_id"] for r in
               csv.DictReader((REPO / "data/processed/thin_slice_cipro_labels_qc.csv").open())]
    gset = set(genomes)
    per_drug = {col: fetch(enc) for col, enc in DRUGS.items()}
    out = REPO / "data/processed/panel_labels.csv"
    with out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["genome_id"] + list(DRUGS))
        for g in genomes:
            w.writerow([g] + [per_drug[col].get(g, "") for col in DRUGS])
    # coverage report
    print(f"panel labels for {len(gset)} annotated genomes -> {out}\n")
    for col in DRUGS:
        labeled = {g: per_drug[col][g] for g in genomes if g in per_drug[col]}
        R = sum(1 for v in labeled.values() if v == "Resistant")
        S = sum(1 for v in labeled.values() if v == "Susceptible")
        print(f"  {col:32s} labeled={R+S:4d}  R={R:4d} S={S:4d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
