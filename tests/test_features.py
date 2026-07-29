"""Feature-builder tests (proposal §4): AMRFinderPlus TSV -> binary determinant matrix.

Covers the two pure functions in src/features/build.py:
  - parse_determinants: keep Type==AMR Element symbols, drop VIRULENCE/STRESS and blanks.
  - encode_matrix: correct 0/1 encoding, sorted+consistent columns, genome order preserved.
Plus an integrity check on the real matrix if it has been built.
"""
import csv
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from src.features.build import encode_matrix, parse_determinants  # noqa: E402

_HEADER = ["Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol",
           "Element name", "Scope", "Type", "Subtype", "Class", "Subclass"]


def _write_tsv(path: Path, rows: list[dict]) -> Path:
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(_HEADER)
        for r in rows:
            w.writerow([r.get(h, "") for h in _HEADER])
    return path


def test_parse_determinants_keeps_only_amr(tmp_path):
    tsv = _write_tsv(tmp_path / "g.tsv", [
        {"Element symbol": "blaKPC-2", "Type": "AMR"},
        {"Element symbol": "gyrA_S83I", "Type": "AMR"},
        {"Element symbol": "ybtS", "Type": "VIRULENCE"},   # excluded
        {"Element symbol": "arsB", "Type": "STRESS"},       # excluded
    ])
    assert parse_determinants(tsv) == {"blaKPC-2", "gyrA_S83I"}


def test_parse_determinants_ignores_blank_and_case(tmp_path):
    tsv = _write_tsv(tmp_path / "g.tsv", [
        {"Element symbol": "sul1", "Type": "amr"},    # lower-case Type still counts
        {"Element symbol": "", "Type": "AMR"},          # blank symbol dropped
    ])
    assert parse_determinants(tsv) == {"sul1"}


def test_parse_determinants_empty(tmp_path):
    assert parse_determinants(_write_tsv(tmp_path / "g.tsv", [])) == set()


def test_encode_matrix_binary_and_sorted():
    per_genome = {"g1": {"blaKPC-2", "sul1"}, "g2": {"sul1"}, "g3": set()}
    cols, rows = encode_matrix(per_genome, ["g1", "g2", "g3"])
    assert cols == ["blaKPC-2", "sul1"]                 # sorted union
    d = dict(rows)
    assert d["g1"] == [1, 1]
    assert d["g2"] == [0, 1]                            # absence -> 0
    assert d["g3"] == [0, 0]
    assert all(v in (0, 1) for _, r in rows for v in r)  # strictly binary


def test_encode_matrix_preserves_id_order_and_skips_missing():
    per_genome = {"g2": {"a"}, "g1": {"a"}}
    cols, rows = encode_matrix(per_genome, ["g1", "g2", "g_absent"])
    assert [g for g, _ in rows] == ["g1", "g2"]         # ids order kept, absent skipped
    assert cols == ["a"]


def test_real_matrix_is_binary_with_unique_ids():
    """If the real K. pneumoniae matrix exists, every cell is 0/1 and genome_ids are unique."""
    m = REPO / "data/processed/thin_slice_cipro_features.csv"
    if not m.exists():
        pytest.skip("feature matrix not built yet")
    rows = list(csv.DictReader(m.open()))
    ids = [r["genome_id"] for r in rows]
    assert len(ids) == len(set(ids)), "duplicate genome_id in feature matrix"
    cols = [c for c in rows[0] if c != "genome_id"]
    for r in rows[:200]:                                # sample for speed
        assert all(r[c] in ("0", "1") for c in cols), "non-binary value in feature matrix"
