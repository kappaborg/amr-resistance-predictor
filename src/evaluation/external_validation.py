"""External validation, Step 1 — scope an independent K. pneumoniae cohort from NCBI Pathogen Detection.

All training data came from BV-BRC. BV-BRC ingests AMR phenotypes *from* NCBI BioSample/Antibiogram
records, so NCBI Pathogen Detection sits **upstream** of our corpus: a naive external validation
against it would score the model on isolates it was trained on.

This module therefore does three things and nothing else:
  1. parse the NCBI PD Klebsiella AMR metadata table (AST phenotypes keyed by accession),
  2. fetch BioSample/assembly accessions for our training genomes from BV-BRC,
  3. **anti-join** to keep only isolates the model has never seen, and report per-drug R/S counts.

It downloads metadata only (~179 MB) and no genome assemblies. It trains nothing and scores nothing:
producing external metrics is Step 2, against the pre-registered plan in
`results/reports/summary_32_external_validation_scoping.md`.

Usage:
    python -m src.evaluation.external_validation            # scope (re-uses cached download)
    python -m src.evaluation.external_validation --refresh  # re-download the NCBI table
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "external_validation"
OUT_CSV = ROOT / "results" / "metrics" / "external_cohort_kpneu.csv"

# NCBI Pathogen Detection release used for the scoping run. Pinned deliberately: a moving "latest"
# would make the cohort irreproducible.
PD_RELEASE = "PDG000000012.2502"
PD_URL = (
    f"https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/{PD_RELEASE}"
    f"/AMR/{PD_RELEASE}.amr.metadata.tsv"
)
PD_TSV = RAW / "kleb_amr_metadata.tsv"
ACC_JSON = RAW / "ours_accessions.json"

KPNEU_TAXID = "573"
PANEL = [
    "meropenem",
    "gentamicin",
    "ciprofloxacin",
    "trimethoprim-sulfamethoxazole",
    "cefoxitin",
]

BVBRC_API = "https://www.bv-brc.org/api/genome/"
BVBRC_HDRS = {
    "Content-Type": "application/rqlquery+x-www-form-urlencoded",
    "Accept": "application/json",
    "User-Agent": "reading-resistance-external-validation/1.0 (research)",
}


def download_pd_table(refresh: bool = False) -> Path:
    """Fetch the NCBI PD Klebsiella AMR metadata table (~179 MB), unless already cached."""
    RAW.mkdir(parents=True, exist_ok=True)
    if PD_TSV.exists() and not refresh:
        print(f"[cache] {PD_TSV.name} ({PD_TSV.stat().st_size:,} bytes)")
        return PD_TSV
    print(f"[get ] {PD_URL}")
    urllib.request.urlretrieve(PD_URL, PD_TSV)
    print(f"[ok  ] {PD_TSV.stat().st_size:,} bytes")
    return PD_TSV


def parse_pd_ast(path: Path) -> list[tuple]:
    """Return (biosample, assembly, collection_date, {drug: SIR}) for K. pneumoniae with panel AST."""
    csv.field_size_limit(10**9)
    recs = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("species_taxid") != KPNEU_TAXID:
                continue
            ast = (row.get("AST_phenotypes") or "").strip()
            if not ast or ast.upper() in {"NULL", "NA"}:
                continue
            calls = {}
            for tok in ast.split(","):
                if "=" not in tok:
                    continue
                drug, val = tok.rsplit("=", 1)
                drug = drug.strip().lower()
                if drug in PANEL:
                    calls[drug] = val.strip().upper()
            if calls:
                recs.append(
                    (row.get("biosample_acc"), row.get("asm_acc"), row.get("collection_date"), calls)
                )
    return recs


def fetch_our_accessions(genome_ids: list[str], batch: int = 200) -> dict[str, tuple]:
    """Map our BV-BRC genome_ids -> (biosample_accession, assembly_accession). Metadata only."""
    out: dict[str, tuple] = {}
    for i in range(0, len(genome_ids), batch):
        chunk = genome_ids[i : i + batch]
        query = (
            "in(genome_id,({}))"
            "&select(genome_id,biosample_accession,assembly_accession)&limit(1000)"
        ).format(",".join(chunk))
        req = urllib.request.Request(
            BVBRC_API, data=query.encode(), headers=BVBRC_HDRS, method="POST"
        )
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req, timeout=120) as resp:
                    data = json.load(resp)
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(2 * (attempt + 1))
        for rec in data:
            out[rec["genome_id"]] = (
                rec.get("biosample_accession"),
                rec.get("assembly_accession"),
            )
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refresh", action="store_true", help="re-download the NCBI PD table")
    args = ap.parse_args()

    download_pd_table(refresh=args.refresh)
    recs = parse_pd_ast(PD_TSV)
    print(f"\nNCBI PD {PD_RELEASE}: K. pneumoniae isolates with >=1 panel drug: {len(recs):,}")

    panel_csv = ROOT / "data" / "processed" / "panel_genomes.csv"
    if not panel_csv.exists():
        print(f"missing {panel_csv} — run the K. pneumoniae pipeline first", file=sys.stderr)
        return 1
    ids = [r["genome_id"] for r in csv.DictReader(open(panel_csv))]

    if ACC_JSON.exists() and not args.refresh:
        ours = json.load(open(ACC_JSON))
        print(f"[cache] accessions for {len(ours):,} training genomes")
    else:
        print(f"[get ] BV-BRC accessions for {len(ids):,} training genomes")
        ours = fetch_our_accessions(ids)
        json.dump(ours, open(ACC_JSON, "w"))

    our_bs = {v[0] for v in ours.values() if v[0]}
    our_asm = {v[1].split(".")[0] for v in ours.values() if v[1]}
    # genomes we cannot anti-join at all -> bounds the undetected-overlap risk
    unmatchable = sum(1 for v in ours.values() if not v[0] and not v[1])

    fresh = [
        r
        for r in recs
        if not ((r[0] and r[0] in our_bs) or (r[1] and r[1].split(".")[0] in our_asm))
    ]
    overlap = len(recs) - len(fresh)
    downloadable = [r for r in fresh if r[1] and r[1] not in {"NULL", ""}]

    print(f"  overlapping our training set (excluded): {overlap:,} ({100*overlap/len(recs):.1f}%)")
    print(f"  genuinely unseen                       : {len(fresh):,}")
    print(f"  ...with a downloadable assembly        : {len(downloadable):,}")

    print(f"\n{'drug':34s} {'R':>5} {'S':>5} {'R+S':>6} {'bal':>5}")
    for drug in PANEL:
        c = collections.Counter(r[3].get(drug) for r in downloadable if drug in r[3])
        r_, s_ = c["R"], c["S"]
        rs = r_ + s_
        print(f"{drug:34s} {r_:5,} {s_:5,} {rs:6,} {(min(r_,s_)/rs if rs else 0):5.2f}")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["biosample_acc", "asm_acc", "collection_date", *PANEL])
        for bs, asm, cdate, calls in downloadable:
            w.writerow([bs, asm, cdate, *[calls.get(d, "") for d in PANEL]])
    print(f"\nwrote {OUT_CSV.relative_to(ROOT)} ({len(downloadable):,} rows)")

    bound = 100 * unmatchable / len(downloadable) if downloadable else 0
    print(
        f"\nresidual risk: {unmatchable:,} training genomes ({100*unmatchable/len(ours):.1f}%) carry "
        f"no accession and could not be anti-joined\n"
        f"  -> undetected overlap bounded above by {bound:.1f}% of the cohort (worst case)"
    )
    print(
        f"\nStep 2 estimate: ~{len(downloadable)*5/1024:.1f} GB of assemblies, "
        f"~{len(downloadable)*47/3600:.1f} h of AMRFinderPlus. Not started."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
