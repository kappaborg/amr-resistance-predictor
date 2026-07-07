"""Phase 1 — BV-BRC organism/drug survey (COUNTS ONLY, no genome downloads).

Verified BV-BRC Data API contract (2026-07-07):
  - core `genome_amr`; filter organism by `taxon_id` (species-level; NO `taxon_lineage_ids` field).
  - `evidence == "Laboratory Method"` marks real lab phenotypes; everything else is computational.
    We count LAB ONLY — using computational predictions as labels would be circular.
  - exact counts read from the `Content-Range` response header (GET, `Accept: application/json`).
  - antibiotic list via facet: `...&facet((field,antibiotic),(mincount,1),(limit,-1))&json(nl,map)`
    with `Accept: application/solr+json` (the `json(nl,map)` is required for parseable JSON).

Two-step so thresholds can be tuned WITHOUT re-hitting the network:
  1. survey()          -> results/metrics/phase1_counts.csv   (raw R/S/I counts, queried once)
  2. build_shortlist() -> results/reports/phase1_shortlist.md  (applies thresholds to the CSV)

Metadata only: responses are integer counts, not sequences. Organism/drug choice stays a GATE
requiring microbiology sign-off.

Usage:
    python -m src.data.survey --config config/config.yaml            # query + shortlist
    python -m src.data.survey --config config/config.yaml --shortlist-only  # re-threshold, no network
    python -m src.data.survey --config config/config.yaml --dry-run
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
LAB = 'eq(evidence,%22Laboratory%20Method%22)'  # URL-encoded: evidence="Laboratory Method"
PHENOTYPES = ["Resistant", "Susceptible", "Intermediate"]
PAUSE_S = 0.30  # polite to the public API

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

# A drug needs at least this many lab rows total before we bother counting R/S per phenotype
# (a drug below this can never clear the shortlist; skipping saves hundreds of calls).
PREFILTER_MIN_LAB_ROWS = 300


def _get(query: str, accept: str, dry_run: bool):
    if dry_run:
        print(f"  [dry-run] GET {API}?{query}  accept={accept}")
        return None
    r = requests.get(f"{API}?{query}", headers={"Accept": accept}, timeout=90)
    r.raise_for_status()
    time.sleep(PAUSE_S)
    return r


def content_range_total(query: str, dry_run: bool) -> int:
    r = _get(query + "&limit(1)", "application/json", dry_run)
    if r is None:
        return -1
    cr = r.headers.get("Content-Range", "")
    return int(cr.split("/")[-1]) if "/" in cr else 0


def facet_antibiotics(taxon_id: int, dry_run: bool) -> dict[str, int]:
    """Lab-only antibiotic -> total-lab-row-count for one organism."""
    q = (f"and(eq(taxon_id,{taxon_id}),{LAB})&limit(0)"
         f"&facet((field,antibiotic),(mincount,1),(limit,-1))&json(nl,map)")
    r = _get(q, "application/solr+json", dry_run)
    if r is None:
        return {"<antibiotic>": PREFILTER_MIN_LAB_ROWS}
    return r.json()["facet_counts"]["facet_fields"]["antibiotic"]


def count_phenotype(taxon_id: int, antibiotic: str, phenotype: str, dry_run: bool) -> int:
    ab = antibiotic.replace(" ", "%20").replace("/", "%2F")
    q = (f"and(eq(taxon_id,{taxon_id}),{LAB},"
         f"eq(antibiotic,%22{ab}%22),eq(resistant_phenotype,{phenotype}))")
    return content_range_total(q, dry_run)


def survey(dry_run: bool) -> list[dict]:
    rows: list[dict] = []
    for name, taxon_id in CANDIDATES:
        print(f"\n## {name} (taxon {taxon_id})")
        try:
            drugs = facet_antibiotics(taxon_id, dry_run)
        except Exception as e:  # noqa: BLE001 - report and continue, never fabricate
            print(f"  ! facet failed: {e}")
            continue
        # only spend R/S/I calls on drugs with enough lab data to possibly qualify
        worth = {d: n for d, n in drugs.items() if dry_run or n >= PREFILTER_MIN_LAB_ROWS}
        print(f"  {len(drugs)} lab antibiotics; {len(worth)} with >= {PREFILTER_MIN_LAB_ROWS} lab rows")
        for ab in sorted(worth, key=lambda d: -worth[d]):
            c = {p: count_phenotype(taxon_id, ab, p, dry_run) for p in PHENOTYPES}
            r, s, i = c["Resistant"], c["Susceptible"], c["Intermediate"]
            total = r + s if r >= 0 and s >= 0 else -1
            balance = round(min(r, s) / total, 3) if total > 0 else 0.0
            rows.append({"organism": name, "taxon_id": taxon_id, "antibiotic": ab,
                         "n_resistant": r, "n_susceptible": s, "n_intermediate": i,
                         "total_RS": total, "balance": balance})
            if not dry_run:
                print(f"    {ab:26s} R={r:>6} S={s:>6} I={i:>5} bal={balance:0.2f}")
    return rows


def write_counts(rows: list[dict], repo: Path) -> Path:
    p = repo / "results/metrics/phase1_counts.csv"
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return p


def build_shortlist(repo: Path, min_per_class: int, min_balance: float, min_total: int) -> Path:
    csv_path = repo / "results/metrics/phase1_counts.csv"
    rows = list(csv.DictReader(csv_path.open()))
    for r in rows:
        for k in ("n_resistant", "n_susceptible", "n_intermediate", "total_RS"):
            r[k] = int(r[k])
        r["balance"] = float(r["balance"])
        r["qualifies"] = (min(r["n_resistant"], r["n_susceptible"]) >= min_per_class
                          and r["balance"] >= min_balance and r["total_RS"] >= min_total)

    md = ["# Phase 1 Shortlist — BV-BRC lab-only counts\n",
          f"**Thresholds:** min(R,S) >= {min_per_class} · balance >= {min_balance} · "
          f"total(R+S) >= {min_total}. Lab phenotypes only (evidence=Laboratory Method). "
          "Intermediate excluded. Row counts (genome may appear >1x); dedup in Phase 2/3. "
          "Lineage diversity confirmed at Phase 5.\n"]
    by_org: dict[str, list[dict]] = {}
    for r in rows:
        by_org.setdefault(r["organism"], []).append(r)
    ranking = sorted(by_org.items(),
                     key=lambda kv: sum(x["qualifies"] for x in kv[1]), reverse=True)
    md.append("## Ranking by # qualifying drugs\n")
    md.append("| Organism | qualifying drugs |\n|---|---|")
    for org, items in ranking:
        md.append(f"| {org} | {sum(x['qualifies'] for x in items)} |")
    for org, items in ranking:
        nq = sum(x["qualifies"] for x in items)
        md.append(f"\n## {org}  ({nq} qualifying)\n")
        md.append("| Antibiotic | #R | #S | #Int | balance | total | qualifies |")
        md.append("|---|---|---|---|---|---|---|")
        for i in sorted(items, key=lambda x: (x["qualifies"], min(x["n_resistant"], x["n_susceptible"])),
                        reverse=True):
            md.append(f"| {i['antibiotic']} | {i['n_resistant']} | {i['n_susceptible']} | "
                      f"{i['n_intermediate']} | {i['balance']} | {i['total_RS']} | "
                      f"{'**yes**' if i['qualifies'] else 'no'} |")
    out = repo / "results/reports/phase1_shortlist.md"
    out.write_text("\n".join(md) + "\n")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="BV-BRC organism/drug survey (counts only).")
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--dry-run", action="store_true", help="print planned queries, no network")
    ap.add_argument("--shortlist-only", action="store_true",
                    help="re-apply thresholds to existing counts CSV (no network)")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[2]
    cfg_path = Path(args.config) if args.config.startswith("/") else repo / args.config
    cfg = yaml.safe_load(cfg_path.read_text())
    data_cfg = cfg.get("data", {})
    min_per_class = data_cfg.get("min_genomes_per_class", 100)
    min_balance = data_cfg.get("min_balance", 0.20)
    min_total = data_cfg.get("min_total_rs", 300)

    if not args.shortlist_only:
        rows = survey(args.dry_run)
        if args.dry_run:
            print("\n[dry-run] no network calls made, no outputs written.")
            return 0
        if not rows:
            print("No rows collected — check API availability/field names.")
            return 1
        print(f"\nWrote {write_counts(rows, repo)}")

    out = build_shortlist(repo, min_per_class, min_balance, min_total)
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
