"""Data-join tests (proposal §4): features ⋈ labels ⋈ lineages by genome_id.

Every model is trained on genomes that have (a) a determinant feature row, (b) a consistent R/S
label, and (c) an MLST lineage. A silent misalignment there — a duplicated id, a labelled genome
with no lineage — would corrupt the phylogeny-aware split. These tests lock the join contract with a
synthetic fixture (always runs) and check the real processed files when present.
"""
import csv
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


def inner_join_ids(features, labels, lineages):
    """Genomes usable for training = present in ALL three sources. Pure set intersection."""
    return set(features) & set(labels) & set(lineages)


def test_inner_join_keeps_only_common_ids():
    feats = {"g1", "g2", "g3", "g4"}
    labels = {"g1", "g2", "g3"}          # g4 unlabelled
    lineages = {"g1", "g2", "g4"}        # g3 untyped
    assert inner_join_ids(feats, labels, lineages) == {"g1", "g2"}


def test_inner_join_empty_when_no_overlap():
    assert inner_join_ids({"a"}, {"b"}, {"c"}) == set()


def _ids(path, col="genome_id"):
    return [r[col] for r in csv.DictReader(path.open())]


def test_real_join_integrity():
    """On the real K. pneumoniae files: unique ids, and every labelled+featured genome is typed."""
    feats = REPO / "data/processed/thin_slice_cipro_features.csv"
    lin = REPO / "data/processed/thin_slice_cipro_lineages.csv"
    panel = REPO / "data/processed/panel_labels.csv"
    if not (feats.exists() and lin.exists() and panel.exists()):
        pytest.skip("processed K. pneumoniae files not built yet")

    fids, lids, pids = _ids(feats), _ids(lin), _ids(panel)
    # no duplicate keys in any source (a dup would double-weight a genome / break the join)
    for name, ids in [("features", fids), ("lineages", lids), ("panel", pids)]:
        assert len(ids) == len(set(ids)), f"duplicate genome_id in {name}"

    fset, lset = set(fids), set(lids)
    # every genome with features + a lab R/S label in the panel must also carry a lineage,
    # otherwise it would be dropped from (or silently unbalance) the phylogeny-aware split.
    rows = list(csv.DictReader(panel.open()))
    drugs = [c for c in rows[0] if c != "genome_id"]
    labelled_featured = {r["genome_id"] for r in rows
                         if r["genome_id"] in fset
                         and any(r.get(d) in ("Resistant", "Susceptible") for d in drugs)}
    assert labelled_featured, "no labelled+featured genomes found"
    untyped = labelled_featured - lset
    assert not untyped, f"{len(untyped)} labelled+featured genomes have no lineage: {sorted(untyped)[:5]}"


def test_real_lineage_file_has_expected_columns():
    lin = REPO / "data/processed/thin_slice_cipro_lineages.csv"
    if not lin.exists():
        pytest.skip("lineage file not built yet")
    header = next(csv.reader(lin.open()))
    assert "genome_id" in header and "lineage" in header
