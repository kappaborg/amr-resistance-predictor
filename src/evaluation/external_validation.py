"""External validation, Step 1 — scope an independent K. pneumoniae cohort.

All training data came from BV-BRC, which ingests AMR phenotypes *from* NCBI BioSample/Antibiogram
records. NCBI Pathogen Detection therefore sits **upstream** of our corpus: a naive external
validation against it would score the model on isolates it was trained on.

Two sources are combined, each filtered before use:
  * **EBI AMR Portal (CABBAGE)** — carries a `database` provenance column, so every row PATRIC/BV-BRC
    touched is dropped *by construction*. This is the stronger independence guarantee.
  * **NCBI Pathogen Detection** — no provenance column, so independence rests entirely on the
    accession anti-join; measured overlap with our training set is 16.1%.

Both are then anti-joined against the BioSample/assembly accessions of our training genomes, and
merged (isolate/drug pairs where the two sources disagree are dropped, not arbitrated).

Downloads metadata only (~196 MB) and no genome assemblies. Trains nothing, scores nothing:
producing external metrics is Step 2, against the pre-registered plan in
`results/reports/summary_32_external_validation_scoping.md`.

CAVEAT carried into Step 2: this anti-join is **identity-level**. It does not remove near-clonal
siblings of training isolates, so Step 2 must run `mlst` and report performance both on all external
isolates and restricted to sequence types absent from training.

Usage:
    python -m src.evaluation.external_validation            # scope (re-uses cached downloads)
    python -m src.evaluation.external_validation --refresh  # re-download both tables
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

# EMBL-EBI AMR Portal (CABBAGE). Preferred source: it carries a `database` provenance column, so
# rows PATRIC/BV-BRC touched can be excluded *by construction* rather than by accession luck.
# Dickens et al., bioRxiv 2025, doi:10.1101/2025.11.12.688105. CC BY 4.0, no login.
PORTAL_RELEASE = "2026-07"
PORTAL_URL = (
    f"https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/{PORTAL_RELEASE}/phenotype.csv.gz"
)
PORTAL_GZ = RAW / "phenotype.csv.gz"

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


def download_portal(refresh: bool = False) -> Path:
    """Fetch the EBI AMR Portal phenotype table (~17 MB), unless already cached."""
    RAW.mkdir(parents=True, exist_ok=True)
    if PORTAL_GZ.exists() and not refresh:
        print(f"[cache] {PORTAL_GZ.name} ({PORTAL_GZ.stat().st_size:,} bytes)")
        return PORTAL_GZ
    print(f"[get ] {PORTAL_URL}")
    urllib.request.urlretrieve(PORTAL_URL, PORTAL_GZ)
    print(f"[ok  ] {PORTAL_GZ.stat().st_size:,} bytes")
    return PORTAL_GZ


def parse_portal(path: Path) -> dict[str, dict]:
    """K. pneumoniae isolates with panel AST, EXCLUDING every row of PATRIC provenance.

    Returns {biosample: {"asm": str, "calls": {drug: 'R'|'S'}}}. Rows where the same isolate/drug
    carries contradictory calls are dropped rather than arbitrated.
    """
    import gzip

    out: dict[str, dict] = {}
    with gzip.open(path, "rt", newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("organism") != "Klebsiella pneumoniae":
                continue
            if "PATRIC" in (row.get("database") or "").upper():
                continue  # provenance filter — the whole point of using this source
            drug = (row.get("antibiotic_name") or "").strip().lower()
            if drug not in PANEL:
                continue
            val = {"resistant": "R", "susceptible": "S"}.get(
                (row.get("resistance_phenotype") or "").strip().lower()
            )
            if not val:
                continue
            bs = row.get("BioSample_ID")
            if not bs:
                continue
            rec = out.setdefault(bs, {"asm": row.get("assembly_ID") or "", "calls": {}})
            prev = rec["calls"].get(drug)
            rec["calls"][drug] = "CONFLICT" if (prev and prev != val) else val
    for rec in out.values():
        rec["calls"] = {d: v for d, v in rec["calls"].items() if v in ("R", "S")}
    return {b: r for b, r in out.items() if r["calls"]}


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

    def unseen(bs: str | None, asm: str | None) -> bool:
        return not ((bs and bs in our_bs) or (asm and asm.split(".")[0] in our_asm))

    fresh = [r for r in recs if unseen(r[0], r[1])]
    overlap = len(recs) - len(fresh)
    ncbi_ok = {r[0]: r for r in fresh if r[1] and r[1] not in {"NULL", ""}}

    print(f"  overlapping our training set (excluded): {overlap:,} ({100*overlap/len(recs):.1f}%)")
    print(f"  genuinely unseen                       : {len(fresh):,}")
    print(f"  ...with a downloadable assembly        : {len(ncbi_ok):,}")

    # --- Source 1: EBI AMR Portal (provenance-filtered) ---
    download_portal(refresh=args.refresh)
    portal_all = parse_portal(PORTAL_GZ)
    portal_ok = {
        b: r for b, r in portal_all.items() if unseen(b, r["asm"]) and r["asm"]
    }
    print(f"\nAMR Portal {PORTAL_RELEASE}: PATRIC-filtered + anti-joined, with assembly: {len(portal_ok):,}")

    # --- merge; drop any isolate/drug where the two sources disagree ---
    merged: dict[str, dict] = {}
    conflicts = 0
    for bs in set(portal_ok) | set(ncbi_ok):
        asm = portal_ok.get(bs, {}).get("asm") or (ncbi_ok.get(bs) or (None, "", "", {}))[1]
        cdate = (ncbi_ok.get(bs) or (None, None, "", {}))[2] or ""
        calls = {}
        for drug in PANEL:
            pv = portal_ok.get(bs, {}).get("calls", {}).get(drug)
            nv = (ncbi_ok.get(bs) or (None, None, None, {}))[3].get(drug)
            pv = pv if pv in ("R", "S") else None
            nv = nv if nv in ("R", "S") else None
            if pv and nv and pv != nv:
                conflicts += 1
                continue
            if pv or nv:
                calls[drug] = pv or nv
        if calls:
            merged[bs] = {"asm": asm, "date": cdate, "calls": calls}

    inter = set(portal_ok) & set(ncbi_ok)
    print(f"  intersection with NCBI: {len(inter):,}  | union with usable labels: {len(merged):,}")
    print(f"  cross-source label conflicts dropped: {conflicts}")

    print(f"\n{'drug':34s} {'R':>5} {'S':>5} {'R+S':>6} {'bal':>5}")
    for drug in PANEL:
        c = collections.Counter(m["calls"].get(drug) for m in merged.values() if drug in m["calls"])
        r_, s_ = c["R"], c["S"]
        rs = r_ + s_
        print(f"{drug:34s} {r_:5,} {s_:5,} {rs:6,} {(min(r_,s_)/rs if rs else 0):5.2f}")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["biosample_acc", "asm_acc", "collection_date", "sources", *PANEL])
        for bs, m in sorted(merged.items()):
            src = "+".join(
                s for s, present in (("portal", bs in portal_ok), ("ncbi", bs in ncbi_ok)) if present
            )
            w.writerow([bs, m["asm"], m["date"], src, *[m["calls"].get(d, "") for d in PANEL]])
    print(f"\nwrote {OUT_CSV.relative_to(ROOT)} ({len(merged):,} rows)")
    downloadable = merged

    bound = 100 * unmatchable / len(downloadable) if downloadable else 0
    print(
        f"\nresidual risk: {unmatchable:,} training genomes ({100*unmatchable/len(ours):.1f}%) carry "
        f"no accession and could not be anti-joined\n"
        f"  -> undetected overlap bounded above by {bound:.1f}% of the cohort (worst case)"
    )
    print(
        f"\nStep 2 estimate: ~{len(downloadable)*1.61/1024:.1f} GB of assemblies "
        f"(1.61 MB/assembly, measured), ~{len(downloadable)*47/3600:.1f} h AMRFinderPlus "
        f"+ ~{len(downloadable)*2/3600:.1f} h mlst. Not started.\n"
        "NOTE: Step 2 MUST run mlst and report performance BOTH on all external isolates AND\n"
        "restricted to sequence types absent from training — accession de-duplication does not\n"
        "remove near-clonal siblings (non-negotiable #1 at the external boundary)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
