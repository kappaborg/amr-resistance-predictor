"""The most important test in the project: prove NO lineage appears in both train and test.

A random train/test split leaks lineage identity (bacteria are clonal) and inflates every metric.
This test is the guardrail for the phylogeny-aware split (Phase 5). It runs on a synthetic fixture
now so the invariant is locked in before real data exists; the real split plugs into the same
`assert_no_lineage_overlap` contract.
"""
import pytest


def assert_no_lineage_overlap(train_lineages, test_lineages):
    """Core invariant: the set of lineages in train and test must be disjoint."""
    overlap = set(train_lineages) & set(test_lineages)
    assert not overlap, f"Lineage leakage: {sorted(overlap)} appear in both train and test"


def test_disjoint_lineages_pass():
    assert_no_lineage_overlap(["ST11", "ST15", "ST258"], ["ST101", "ST307"])


def test_shared_lineage_fails():
    with pytest.raises(AssertionError, match="Lineage leakage"):
        assert_no_lineage_overlap(["ST11", "ST15"], ["ST15", "ST307"])


# All organisms' phylogeny-aware splits (K. pneumoniae uses the thin-slice tag; the rest are written
# by src.organism_pipeline as <tag>_split.csv). The zero-leakage invariant must hold for EVERY one.
_SPLIT_TAGS = ["thin_slice_cipro", "ecoli", "saureus", "abaumannii", "paeruginosa", "senterica",
               "efaecium", "spneumoniae"]


@pytest.mark.parametrize("tag", _SPLIT_TAGS)
def test_real_split_has_no_leakage(tag):
    """For each organism's phylogeny-aware split, assert zero shared lineage between train and test."""
    import csv
    from pathlib import Path
    split = Path(__file__).resolve().parents[1] / f"data/processed/{tag}_split.csv"
    if not split.exists():
        pytest.skip(f"{tag} split not built yet (run src.organism_pipeline --organism {tag})")
    rows = list(csv.DictReader(split.open()))
    train = [r["lineage"] for r in rows if r["fold"] == "train"]
    test = [r["lineage"] for r in rows if r["fold"] == "test"]
    assert train and test, f"{tag} split has an empty fold"
    assert_no_lineage_overlap(train, test)
