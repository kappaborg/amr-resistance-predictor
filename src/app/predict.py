"""Phase 11 — the demo: genome in -> per-drug resistant/susceptible + calibrated confidence +
the determinants behind each call ("why this strain resists drug X").

Runs AMRFinderPlus on the input genome (or reuses a cached annotation for a known genome_id),
builds the determinant feature vector, and for each drug reports the calibrated probability, the
R/S call at the clinical VME<=target operating threshold, and the SHAP-ranked determinants that
drove the call.

Usage:
    python -m src.app.predict --genome path/to/genome.fna[.gz]
    python -m src.app.predict --genome-id 573.12771     # reuse cached annotation (instant)
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import subprocess
import sys
from pathlib import Path

import joblib
import numpy as np
import shap

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.app.registry import ORGANISMS, ORG_ORDER, SCHEME_ORG  # noqa: E402

ENV_BIN = "/opt/homebrew/anaconda3/envs/amr-resistance-predictor/bin"
AMR = f"{ENV_BIN}/amrfinder"
MLST = f"{ENV_BIN}/mlst"
MLST_ENV = {**os.environ, "PATH": f"{ENV_BIN}:{os.environ.get('PATH', '')}"}  # mlst needs blastn on PATH
MODELS = REPO / "results/models"
CACHE = REPO / "data/interim/amrfinder"
MLST_CACHE = REPO / "data/interim/mlst"
GENOMES = REPO / "data/raw/genomes"


def determinants_from_tsv(tsv_text: str) -> set[str]:
    dets = set()
    for r in csv.DictReader(tsv_text.splitlines(), delimiter="\t"):
        if (r.get("Type") or "").strip().upper() == "AMR":
            sym = (r.get("Element symbol") or "").strip()
            if sym:
                dets.add(sym)
    return dets


def _read_genome_bytes(genome_path: Path) -> bytes:
    """Return decompressed FASTA bytes. Detects gzip by MAGIC BYTES (0x1f8b), not extension — so a
    mislabelled .gz (actually plain) or a plain-named file that is really gzipped both work."""
    raw = genome_path.read_bytes()
    if raw[:2] == b"\x1f\x8b":                       # gzip magic
        try:
            raw = gzip.decompress(raw)
        except OSError:
            sys.exit("The uploaded .gz file is not a valid gzip archive. Re-compress it with gzip, "
                     "or upload the plain FASTA (.fna/.fasta/.fa).")
    return raw


def _looks_like_protein(data: bytes) -> bool:
    """Distinguish a protein FASTA (.faa) from a nucleotide one by alphabet: DNA is almost entirely
    A/C/G/T/U/N, whereas amino-acid sequences are dominated by the other ~14 letters."""
    seq = b"".join(l for l in data.splitlines() if not l.startswith(b">"))[:4000].upper()
    seq = bytes(c for c in seq if 65 <= c <= 90)     # letters only
    if not seq:
        return False
    nuc = sum(1 for c in seq if c in b"ACGTUN")
    return nuc / len(seq) < 0.85                       # <85% nucleotide letters => protein


def annotate(genome_path: Path, amrfinder_org: str) -> set[str]:
    """Run AMRFinderPlus on a genome FASTA. Nucleotide (.fna/.fasta/.fa/.gz) -> `-n`; protein
    (.faa, e.g. predicted proteins) -> `-p`. Molecule type is detected by content, not extension."""
    data = _read_genome_bytes(genome_path)
    if not data.lstrip()[:1] == b">":                # basic FASTA sanity check
        sys.exit("That doesn't look like a FASTA (no '>' header found). Please upload an assembled "
                 "nucleotide genome (.fna/.fasta/.fa/.gz) or a protein FASTA (.faa).")
    is_protein = _looks_like_protein(data) or genome_path.suffix.lower() == ".faa"
    flag, mol = ("-p", "protein") if is_protein else ("-n", "nucleotide")
    tmp = CACHE / f"_demo_{genome_path.stem}.fa"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    tmp.write_bytes(data)
    print(f"annotating {genome_path.name} ({mol}) with AMRFinderPlus --organism {amrfinder_org} "
          "(~1-3 min)...", flush=True)
    try:
        r = subprocess.run([AMR, flag, str(tmp), "--organism", amrfinder_org, "--threads", "4"],
                           capture_output=True, text=True)
    finally:
        tmp.unlink(missing_ok=True)
    if r.returncode != 0:
        sys.exit(f"AMRFinderPlus failed: {r.stderr[-300:]}")
    return determinants_from_tsv(r.stdout)


# ---- Species verification via MLST (a scheme is species-specific) --------------------------------
def _parse_scheme(mlst_line: str) -> str | None:
    parts = mlst_line.split("\t")
    s = parts[1].strip() if len(parts) > 1 else ""
    return s if s and s != "-" else None


def genome_scheme(genome_id: str | None, genome: str | None) -> str | None:
    """Detected MLST scheme for the input (or None). Cached genomes reuse data/interim/mlst/<id>.tsv;
    uploads run mlst live on the nucleotide FASTA. Protein (.faa) uploads return None (unverifiable)."""
    if genome_id:
        tsv = MLST_CACHE / f"{genome_id}.tsv"
        if tsv.exists() and tsv.read_text().strip():
            return _parse_scheme(tsv.read_text().splitlines()[0])
        gz = GENOMES / f"{genome_id}.fna.gz"
        return _run_mlst(gz) if gz.exists() else None
    p = Path(genome)
    data = _read_genome_bytes(p)
    if _looks_like_protein(data) or p.suffix.lower() == ".faa":
        return None                                        # mlst needs nucleotide input
    tmp = CACHE / f"_mlst_{p.stem}.fna"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    tmp.write_bytes(data)
    try:
        return _run_mlst(tmp)
    finally:
        tmp.unlink(missing_ok=True)


def _run_mlst(fna: Path) -> str | None:
    src = fna
    if fna.suffix == ".gz":                                # decompress cached .fna.gz to a temp
        src = CACHE / f"_mlst_{fna.stem}"
        src.write_bytes(gzip.decompress(fna.read_bytes()))
    try:
        r = subprocess.run([MLST, str(src)], capture_output=True, text=True, env=MLST_ENV, timeout=600)
    except Exception:  # noqa: BLE001
        return None
    finally:
        if src is not fna:
            src.unlink(missing_ok=True)
    if r.returncode != 0 or not r.stdout.strip():
        return None
    return _parse_scheme(r.stdout.splitlines()[0])


def _marker_hits(dets: set[str]) -> set[str]:
    """Organism keys whose intrinsic signature genes are present (family prefix match). Works on both
    nucleotide and PROTEIN input, so it catches wrong-species uploads MLST can't (e.g. protein .faa)."""
    hits = set()
    for k, o in ORGANISMS.items():
        for m in o.get("markers", []):
            if any(d == m or d.startswith(m) for d in dets):
                hits.add(k)
                break
    return hits


def species_verdict(org_key: str, scheme: str | None, dets: set[str] | None = None) -> tuple[str, str]:
    """('ok'|'mismatch'|'unverified', message). Two independent signals, conservatively combined:
      - MLST scheme (nucleotide only) — species-specific typing.
      - Intrinsic-marker genes (works on protein too) — e.g. fosA/oqxA => K. pneumoniae, blaADC =>
        A. baumannii. Catches a genome carrying ANOTHER panel organism's signature.
    A mismatch from EITHER signal withholds; either positive confirms; otherwise unverified."""
    sel = ORGANISMS[org_key]["display"]
    mapped = SCHEME_ORG.get(scheme) if scheme else None
    hits = _marker_hits(dets or set())
    others = hits - {org_key}

    # --- mismatch (strong): either MLST or intrinsic markers point to a different panel organism ---
    if mapped and mapped != org_key:
        return "mismatch", (f"MLST types this genome to the {ORGANISMS[mapped]['display']} scheme "
                            f"('{scheme}'), not {sel}")
    if others and org_key not in hits:
        names = ", ".join(ORGANISMS[o]["display"] for o in sorted(others))
        return "mismatch", f"carries the signature intrinsic genes of {names}, not {sel}"

    # --- confirmed: MLST scheme matches, or the selected organism's own markers are present ---
    if mapped == org_key:
        return "ok", f"MLST scheme '{scheme}' is consistent with {sel}"
    if org_key in hits:
        return "ok", f"carries {sel} intrinsic marker gene(s)"

    return "unverified", (f"species not confirmed (MLST '{scheme or '—'}' inconclusive, no signature "
                          f"markers) — predictions assume {sel}")


def load_determinants(genome_id: str | None, genome: str | None,
                      amrfinder_org: str) -> tuple[str, set[str]]:
    if genome_id:
        cached = CACHE / f"{genome_id}.tsv"
        if cached.exists():
            return genome_id, determinants_from_tsv(cached.read_text())
        gz = GENOMES / f"{genome_id}.fna.gz"
        if gz.exists():
            return genome_id, annotate(gz, amrfinder_org)
        sys.exit(f"genome_id {genome_id} not found in cache or data/raw/genomes.")
    p = Path(genome)
    if not p.exists():
        sys.exit(f"genome file not found: {p}")
    return p.name, annotate(p, amrfinder_org)


def explain(bundle, x_row, present_dets, top_n=4):
    """SHAP-ranked determinants present in THIS genome that drove the call."""
    model, cols = bundle["model"], bundle["feature_cols"]
    sv = shap.TreeExplainer(model).shap_values(x_row.reshape(1, -1))[0]
    present_idx = [i for i, c in enumerate(cols) if c in present_dets]
    ranked = sorted(present_idx, key=lambda i: -abs(sv[i]))
    out = []
    for i in ranked[:top_n]:
        if abs(sv[i]) < 1e-6:
            continue
        out.append((cols[i], "↑resistant" if sv[i] > 0 else "↓susceptible", float(sv[i])))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Predict per-drug resistance from a bacterial genome.")
    ap.add_argument("--organism", choices=ORG_ORDER, default="kpneu",
                    help="which organism's models to use (default kpneu)")
    ap.add_argument("--genome", help="path to genome FASTA (.fna or .fna.gz)")
    ap.add_argument("--genome-id", help="BV-BRC genome_id already downloaded/annotated")
    args = ap.parse_args()
    if not (args.genome or args.genome_id):
        sys.exit("provide --genome or --genome-id")

    org = ORGANISMS[args.organism]
    drug_order = list(org["drugs"])
    mdir = MODELS / args.organism
    bundles = {d: joblib.load(mdir / f"{d}.joblib") for d in drug_order
               if (mdir / f"{d}.joblib").exists()}
    if not bundles:
        sys.exit(f"no models for {args.organism} — run: python -m src.models.save_models "
                 f"--organism {args.organism}")

    name, dets = load_determinants(args.genome_id, args.genome, org["amrfinder"])
    cols = next(iter(bundles.values()))["feature_cols"]
    x = np.array([1 if c in dets else 0 for c in cols], dtype=int)
    known = [d for d in dets if d in set(cols)]

    # Species verification: MLST scheme (nucleotide) + intrinsic markers (works on protein too)
    status, msg = species_verdict(args.organism, genome_scheme(args.genome_id, args.genome), dets)

    print("\n" + "=" * 68)
    print(f"  RESISTANCE PREDICTION — {name}  ({org['display']})")
    print(f"  {len(dets)} resistance determinants detected ({len(known)} in model vocabulary)")
    print(f"  species check: {status} — {msg}")
    print("=" * 68)
    if status == "mismatch":
        print(f"\n  ⚠  Wrong species — {msg}.")
        print(f"     The {org['display']} models cannot make a valid prediction for a different "
              "species. Predictions withheld.")
        print(f"     Select the matching organism in --organism, or upload a {org['display']} genome.")
        print("=" * 68 + "\n")
        return 0
    if not known and org.get("markers"):
        # Organisms with near-universal intrinsic determinants (fosA/oqxA, blaPDC, blaADC, ...): zero
        # determinants means the input is almost certainly a different species / not an assembled genome.
        print(f"\n  ⚠  No known resistance determinants found — this does not look like a "
              f"{org['display']} genome.")
        print("     The models only see catalogued resistance genes and assume the input IS the "
              "selected\n     organism, so no meaningful prediction can be made. Predictions withheld.")
        print("     Upload an assembled genome of the selected species (or pick the right --organism).")
        print("=" * 68 + "\n")
        return 0
    if not known:
        # Markerless organism (E. coli, Salmonella): 0 determinants is a genuinely possible
        # pan-susceptible isolate, so predict rather than reject. (Low P(resistant) => likely S;
        # the VME<=3% threshold may still label low-prevalence drugs RESISTANT — read the probability.)
        print(f"\n  No known resistance determinants detected — likely pan-susceptible "
              f"{org['display']} (read P(resistant); the VME-safety threshold over-calls low-R drugs).")
    for drug in drug_order:
        if drug not in bundles:
            continue
        b = bundles[drug]
        raw = b["model"].predict_proba(x.reshape(1, -1))[0, 1]
        prob = float(b["calibrator"].transform([raw])[0])
        call = "RESISTANT" if prob >= b["threshold"] else "susceptible"
        # certainty from how decisive the calibrated probability is (not from the threshold)
        flag = "" if (prob >= 0.8 or prob <= 0.2) else "  [uncertain]"
        drivers = explain(b, x, dets)
        print(f"\n  {drug.replace('_', '/'):30s} {call:11s}  P(resistant)={prob:.2f}{flag}")
        if drivers:
            print("     drivers: " + ", ".join(f"{d} {arrow}" for d, arrow, _ in drivers))
        else:
            print("     drivers: none of this genome's determinants are model-relevant for this drug")
    print("\n" + "-" * 68)
    print("  Operating point: threshold tuned to cap very-major error at the clinical target.")
    print("  Research/decision-support only — not a substitute for phenotypic testing.")
    print("=" * 68 + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
