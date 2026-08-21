"""External validation, Step 2a — fetch and annotate the independent cohort.

Downloads the assemblies listed in `results/metrics/external_cohort_kpneu.csv`, then runs
**AMRFinderPlus 4.2.7 / DB 2026-05-15.1** and **mlst** on each — the same tools and versions used
for training, deliberately re-run rather than reusing the sources' precomputed calls (NCBI ships
refgene DB 2026-01-21.1 and the EBI Portal discloses no DB version; a feature-extraction mismatch
would depress external performance for reasons unrelated to generalization).

Every stage is **resumable**: an accession whose output already exists is skipped, so the run can be
interrupted and restarted freely.

`mlst` is not optional here. The cohort was de-duplicated against training by *accession*, which does
not remove near-clonal siblings; Step 2b reports performance both on all external isolates and
restricted to sequence types absent from training (see `summary_32`).

Usage:
    python -m src.evaluation.external_fetch_annotate                 # full run, resumable
    python -m src.evaluation.external_fetch_annotate --limit 5       # smoke test
    python -m src.evaluation.external_fetch_annotate --stage download
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import shutil
import subprocess
import sys
import time
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
COHORT = ROOT / "results" / "metrics" / "external_cohort_kpneu.csv"
GENOMES = ROOT / "data" / "raw" / "external_validation" / "genomes"
AMR_OUT = ROOT / "data" / "interim" / "external_validation" / "amrfinder"
MLST_OUT = ROOT / "data" / "interim" / "external_validation" / "mlst"

AMRFINDER_ORGANISM = "Klebsiella_pneumoniae"
NCBI_DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{acc}/download"


def _env_bin() -> str:
    cand = Path(sys.executable).parent
    if (cand / "amrfinder").exists() or (cand / "mlst").exists():
        return str(cand)
    for tool in ("amrfinder", "mlst"):
        found = shutil.which(tool)
        if found:
            return str(Path(found).parent)
    return "/opt/homebrew/anaconda3/envs/amr-resistance-predictor/bin"


ENV_BIN = _env_bin()
AMR = f"{ENV_BIN}/amrfinder"
MLST = f"{ENV_BIN}/mlst"
TOOL_ENV = {**os.environ, "PATH": f"{ENV_BIN}:{os.environ.get('PATH', '')}"}


def load_cohort(limit: int | None = None) -> list[str]:
    accs = []
    for row in csv.DictReader(open(COHORT)):
        acc = (row.get("asm_acc") or "").strip()
        if acc and acc not in {"NULL", ""}:
            accs.append(acc)
    accs = sorted(set(accs))
    return accs[:limit] if limit else accs


# ---------------------------------------------------------------- download
def fetch_one(acc: str, retries: int = 4) -> tuple[str, str]:
    """Download one assembly to <acc>.fna.gz. Returns (acc, 'ok'|'skip'|'fail: ...')."""
    dest = GENOMES / f"{acc}.fna.gz"
    if dest.exists() and dest.stat().st_size > 10_000:
        return acc, "skip"
    url = NCBI_DATASETS.format(acc=acc) + "?include_annotation_type=GENOME_FASTA"
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "reading-resistance/1.0"})
            with urllib.request.urlopen(req, timeout=300) as resp:
                blob = resp.read()
            with zipfile.ZipFile(io.BytesIO(blob)) as z:
                names = [n for n in z.namelist() if n.endswith((".fna", ".fasta"))]
                if not names:
                    return acc, "fail: no FASTA in archive"
                data = z.read(names[0])
            if not data.lstrip()[:1] == b">":
                return acc, "fail: not FASTA"
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(gzip.compress(data))
            return acc, "ok"
        except Exception as exc:  # noqa: BLE001
            if attempt == retries - 1:
                return acc, f"fail: {type(exc).__name__}"
            time.sleep(2 * (attempt + 1))
    return acc, "fail: exhausted"


def stage_download(accs: list[str], workers: int = 8) -> None:
    GENOMES.mkdir(parents=True, exist_ok=True)
    todo = [a for a in accs if not (GENOMES / f"{a}.fna.gz").exists()]
    print(f"[download] {len(todo):,} to fetch ({len(accs)-len(todo):,} already present)")
    if not todo:
        return
    done = fails = 0
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for fut in as_completed({ex.submit(fetch_one, a): a for a in todo}):
            acc, status = fut.result()
            done += 1
            if status.startswith("fail"):
                fails += 1
                print(f"  ! {acc}: {status}", flush=True)
            if done % 100 == 0 or done == len(todo):
                rate = done / max(time.time() - t0, 1)
                print(f"  {done:,}/{len(todo):,}  ({rate:.1f}/s, {fails} failed)", flush=True)
    print(f"[download] complete: {len(todo)-fails:,} ok, {fails} failed")


# ---------------------------------------------------------------- annotate
def _decompress_to(acc: str, tmpdir: Path) -> Path:
    src = GENOMES / f"{acc}.fna.gz"
    out = tmpdir / f"{acc}.fna"
    out.write_bytes(gzip.decompress(src.read_bytes()))
    return out


def amrfinder_one(acc: str, threads: int = 3) -> tuple[str, str]:
    dest = AMR_OUT / f"{acc}.tsv"
    if dest.exists():
        return acc, "skip"
    if not (GENOMES / f"{acc}.fna.gz").exists():
        return acc, "fail: no genome"
    tmpdir = AMR_OUT / "_tmp"
    tmpdir.mkdir(parents=True, exist_ok=True)
    fna = _decompress_to(acc, tmpdir)
    try:
        r = subprocess.run(
            [AMR, "-n", str(fna), "--organism", AMRFINDER_ORGANISM, "--threads", str(threads)],
            capture_output=True, text=True, env=TOOL_ENV, timeout=1800,
        )
        if r.returncode != 0:
            return acc, f"fail: {r.stderr.strip()[-160:]}"
        if not r.stdout.startswith("Protein id"):        # guard against a truncated/empty run
            return acc, "fail: unexpected AMRFinderPlus output"
        tmp_out = dest.with_suffix(".tsv.part")           # atomic: write then rename, so killing
        tmp_out.write_text(r.stdout)                      # the job never leaves a partial .tsv
        tmp_out.replace(dest)
        return acc, "ok"
    except Exception as exc:  # noqa: BLE001
        return acc, f"fail: {type(exc).__name__}"
    finally:
        fna.unlink(missing_ok=True)


def mlst_one(acc: str) -> tuple[str, str]:
    dest = MLST_OUT / f"{acc}.tsv"
    if dest.exists():
        return acc, "skip"
    if not (GENOMES / f"{acc}.fna.gz").exists():
        return acc, "fail: no genome"
    tmpdir = MLST_OUT / "_tmp"
    tmpdir.mkdir(parents=True, exist_ok=True)
    fna = _decompress_to(acc, tmpdir)
    try:
        r = subprocess.run([MLST, str(fna)], capture_output=True, text=True,
                           env=TOOL_ENV, timeout=900)
        if r.returncode != 0 or not r.stdout.strip():
            return acc, "fail: mlst empty"
        dest.write_text(r.stdout)
        return acc, "ok"
    except Exception as exc:  # noqa: BLE001
        return acc, f"fail: {type(exc).__name__}"
    finally:
        fna.unlink(missing_ok=True)


def _run_stage(name: str, fn, accs: list[str], outdir: Path, workers: int) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    todo = [a for a in accs if not (outdir / f"{a}.tsv").exists()]
    print(f"[{name}] {len(todo):,} to process ({len(accs)-len(todo):,} cached)")
    if not todo:
        return
    done = fails = 0
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for fut in as_completed({ex.submit(fn, a): a for a in todo}):
            acc, status = fut.result()
            done += 1
            if status.startswith("fail"):
                fails += 1
                print(f"  ! {acc}: {status}", flush=True)
            if done % 25 == 0 or done == len(todo):
                el = time.time() - t0
                eta = (len(todo) - done) * el / max(done, 1) / 3600
                print(f"  {done:,}/{len(todo):,}  ({el/60:.0f} min elapsed, "
                      f"~{eta:.1f} h left, {fails} failed)", flush=True)
    shutil.rmtree(outdir / "_tmp", ignore_errors=True)
    print(f"[{name}] complete: {len(todo)-fails:,} ok, {fails} failed")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=None, help="process only the first N (smoke test)")
    ap.add_argument("--stage", choices=["download", "amrfinder", "mlst", "all"], default="all")
    ap.add_argument("--amr-workers", type=int, default=3,
                    help="concurrent AMRFinderPlus processes (each uses --threads)")
    ap.add_argument("--amr-threads", type=int, default=3)
    args = ap.parse_args()

    if not COHORT.exists():
        print(f"missing {COHORT} — run `make extscope` first", file=sys.stderr)
        return 1
    accs = load_cohort(args.limit)
    print(f"cohort: {len(accs):,} assemblies\n")

    if args.stage in ("download", "all"):
        stage_download(accs)
    if args.stage in ("amrfinder", "all"):
        _run_stage("amrfinder", lambda a: amrfinder_one(a, args.amr_threads),
                   accs, AMR_OUT, args.amr_workers)
    if args.stage in ("mlst", "all"):
        _run_stage("mlst", mlst_one, accs, MLST_OUT, 4)

    have = {
        "genomes": sum(1 for a in accs if (GENOMES / f"{a}.fna.gz").exists()),
        "amrfinder": sum(1 for a in accs if (AMR_OUT / f"{a}.tsv").exists()),
        "mlst": sum(1 for a in accs if (MLST_OUT / f"{a}.tsv").exists()),
    }
    print(f"\nready: {have} of {len(accs):,}")
    print("next: python -m src.evaluation.external_score")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
