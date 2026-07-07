"""Phase 1 — BV-BRC organism/drug survey (COUNTS ONLY, no genome downloads).

Facets the `genome_amr` table per candidate organism to list antibiotics, then counts
Resistant / Susceptible / Intermediate genomes per drug via the Content-Range header.
Applies the shortlist criteria from docs/phase1_query_plan.md and writes:
  - results/metrics/phase1_counts.csv   (machine-readable)
  - results/reports/phase1_shortlist.md (the table for microbiology review)

This is metadata only: responses are integer counts, not sequences. Nothing here downloads
a genome. Organism/drug choice remains a GATE requiring microbiology sign-off.

Usage:
    python -m src.data.survey --config config/config.yaml
    python -m src.data.survey --config config/config.yaml --dry-run   # print queries, no network
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import requests
import yaml

API = "https://www.bv-brc.org/api/genome_amr/"

# Candidate organisms (taxon_id verified on first live call; see docs/phase1_query_plan.md).
CANDIDATES = [
    ("Escherichia coli", 562),
    ("Klebsiella pneumoniae", 573),
    ("Salmonella enterica", 28901),
    ("Mycobacterium tuberculosis", 1773),
    ("Staphylococcus aureus", 1280),
    ("Acinetobacter baumannii", 470),
    ("Pseudomonas aeruginosa", 287),
    ("Neisseria gonorrhoeae", 485),
    ("Streptococcus pneumoniae", 1313),
]

# Preferred taxon filter field; falls back to taxon_id if the API rejects it.
TAXON_FIELDS = ["taxon_lineage_ids", "taxon_id"]
PHENOTYPES = ["Resistant", "Susceptible", "Intermediate"]
HEADERS_JSON = {"Accept": "application/json", "Content-Type": "application/rqlquery+x-www-form-urlencoded"}
HEADERS_SOLR = {"Accept": "application/solr+json", "Content-Type": "application/rqlquery+x-www-form-urlencoded"}
PAUSE_S = 0.34  # be polite to the public API


def _post(query: str, headers: dict, dry_run: bool):
    """Send an RQL query as POST body (avoids URL-length limits). Returns response or None."""
    if dry_run:
        print(f"  [dry-run] POST {API}  body={query!r}  accept={headers['Accept']}")
        return None
    resp = requests.post(API, data=query, headers=headers, timeout=60)
    resp.raise_for_status()
    time.sleep(PAUSE_S)
    return resp


def facet_antibiotics(taxon_field: str, taxon_id: int, dry_run: bool) -> list[str]:
    """List antibiotics that have any phenotype data for this organism."""
    query = (f"eq({taxon_field},{taxon_id})&limit(1)"
             f"&facet((field,antibiotic),(mincount,1),(limit,-1))")
    resp = _post(query, HEADERS_SOLR, dry_run)
    if resp is None:
        return []
    facets = resp.json().get("facet_counts", {}).get("facet_fields", {}).get("antibiotic", [])
    # Solr returns a flat [name, count, name, count, ...] list.
    return [facets[i] for i in range(0, len(facets), 2)]


def count(taxon_field: str, taxon_id: int, antibiotic: str, phenotype: str, dry_run: bool) -> int:
    """Exact count for (organism, antibiotic, phenotype) via the Content-Range header."""
    ab = antibiotic.replace('"', '\\"')
    query = (f"and(eq({taxon_field},{taxon_id}),"
             f'eq(antibiotic,"{ab}"),eq(resistant_phenotype,{phenotype}))&limit(1)')
    resp = _post(query, HEADERS_JSON, dry_run)
    if resp is None:
        return -1
    # Content-Range: "items 0-0/1234" -> total is after the slash.
    cr = resp.headers.get("Content-Range", "")
    return int(cr.split("/")[-1]) if "/" in cr else len(resp.json())


def resolve_taxon_field(taxon_id: int, dry_run: bool) -> str:
    """Pick the first taxon filter field the API accepts (verify, don't assume)."""
    if dry_run:
        return TAXON_FIELDS[0]
    for field in TAXON_FIELDS:
        try:
            _post(f"eq({field},{taxon_id})&limit(1)", HEADERS_JSON, dry_run)
            return field
        except requests.HTTPError:
            continue
    raise RuntimeError(f"No usable taxon field among {TAXON_FIELDS} for taxon {taxon_id}")


def survey(cfg: dict, dry_run: bool) -> list[dict]:
    min_per_class = cfg.get("data", {}).get("min_genomes_per_class", 100)
    min_balance = 0.20
    min_total = 300
    rows: list[dict] = []

    for name, taxon_id in CANDIDATES:
        print(f"\n## {name} (taxon {taxon_id})")
        try:
            field = resolve_taxon_field(taxon_id, dry_run)
            drugs = facet_antibiotics(field, taxon_id, dry_run)
        except Exception as e:  # noqa: BLE001 - report and continue, never fabricate
            print(f"  ! query failed: {e}")
            continue
        if dry_run:
            drugs = ["<antibiotic>"]
        for ab in drugs:
            counts = {p: count(field, taxon_id, ab, p, dry_run) for p in PHENOTYPES}
            r, s = counts["Resistant"], counts["Susceptible"]
            total = r + s if r >= 0 and s >= 0 else -1
            balance = (min(r, s) / total) if total > 0 else 0.0
            qualifies = (min(r, s) >= min_per_class and balance >= min_balance
                         and total >= min_total)
            rows.append({
                "organism": name, "taxon_id": taxon_id, "antibiotic": ab,
                "n_resistant": r, "n_susceptible": s,
                "n_intermediate": counts["Intermediate"],
                "total_RS": total, "balance": round(balance, 3),
                "qualifies": qualifies,
            })
            if not dry_run:
                flag = "  <== qualifies" if qualifies else ""
                print(f"  {ab:24s} R={r:>5} S={s:>5} bal={balance:0.2f}{flag}")
    return rows


def write_outputs(rows: list[dict], repo: Path) -> None:
    csv_path = repo / "results/metrics/phase1_counts.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    md = ["# Phase 1 Shortlist — BV-BRC counts\n",
          "Qualifying = min(R,S) >= threshold AND balance >= 0.20 AND total >= 300. "
          "Intermediate excluded. Lineage diversity confirmed later (Phase 5).\n"]
    by_org: dict[str, list[dict]] = {}
    for row in rows:
        by_org.setdefault(row["organism"], []).append(row)
    for org, items in by_org.items():
        q = sorted([i for i in items if i["qualifies"]],
                   key=lambda i: min(i["n_resistant"], i["n_susceptible"]), reverse=True)
        md.append(f"\n## {org}  ({len(q)} qualifying drugs)\n")
        md.append("| Antibiotic | #R | #S | #Int | balance | total | qualifies |")
        md.append("|---|---|---|---|---|---|---|")
        for i in sorted(items, key=lambda x: x["qualifies"], reverse=True):
            md.append(f"| {i['antibiotic']} | {i['n_resistant']} | {i['n_susceptible']} | "
                      f"{i['n_intermediate']} | {i['balance']} | {i['total_RS']} | "
                      f"{'yes' if i['qualifies'] else 'no'} |")
    report = repo / "results/reports/phase1_shortlist.md"
    report.write_text("\n".join(md) + "\n")
    print(f"\nWrote {csv_path}\nWrote {report}")


def main() -> int:
    ap = argparse.ArgumentParser(description="BV-BRC organism/drug survey (counts only).")
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--dry-run", action="store_true", help="print planned queries, no network")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[2]
    cfg = yaml.safe_load((repo / args.config).read_text()) if not args.config.startswith("/") \
        else yaml.safe_load(Path(args.config).read_text())

    rows = survey(cfg, args.dry_run)
    if args.dry_run:
        print("\n[dry-run] no network calls made, no outputs written.")
        return 0
    if not rows:
        print("No rows collected — check API availability/field names.")
        return 1
    write_outputs(rows, repo)
    return 0


if __name__ == "__main__":
    sys.exit(main())
