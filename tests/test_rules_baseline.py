"""Guard the honest baseline: every shipped rules list must actually discriminate.

The project's central claim is "what the model adds *over a known-gene lookup*". That claim is only
meaningful if the lookup is a fair opponent. A rule list that can never fire, or that fires on every
genome, is not a weak baseline — it is a **broken measurement**, and it silently flatters the model.

This has now bitten the project twice:

  * **P. aeruginosa (decision #42):** the provisional ceftazidime/tobramycin rules included
    `blaPDC` and `aph(3')-IIb`, which are *intrinsic* to P. aeruginosa and present in ~100% of
    genomes. Firing everywhere, they could not separate R from S -> a fake ROC of 0.500.
  * **K. pneumoniae cefoxitin (decision #57):** the shared `RULES["cefoxitin"]` entry is the
    *S. aureus* `mecA`/`mecC` screen (cefoxitin is the standard MRSA surrogate test). In Klebsiella
    those genes occur in **0 of 3,850** genomes, so the baseline never fired at all and scored a
    degenerate 0.500 by predicting a constant. That artifact reached the report, the poster, both
    slide decks and both walkthroughs before it was caught by hand.

Neither failure was detectable from the metrics alone: a broken baseline and a genuinely useless
baseline both produce ROC ~0.5. The only way to tell them apart is to ask whether the gene list
*occurs* in the organism at all. That is what this test does, for every shipped model.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

# A baseline firing outside this band cannot discriminate, whatever its ROC says.
MIN_FIRE = 0.01   # below this it is effectively never present
MAX_FIRE = 0.99   # above this it is effectively always present (intrinsic gene)

# Documented biological exceptions — a genuinely absent mechanism, not a wiring mistake.
# Each entry must cite the decision that established it, so exceptions cannot be added silently.
KNOWN_ABSENT: dict[tuple[str, str], str] = {
    ("spneumoniae", "trimethoprim_sulfamethoxazole"):
        "decision #48 — S. pneumoniae TMP-SMX resistance is driven by folA/folP point mutations, "
        "which AMRFinderPlus does not catalogue; the acquired sul/dfr genes genuinely do not occur "
        "in pneumococcus. Reported as a scope limit, not as an ML-adds headline.",
}


def _shipped_models() -> list[tuple[str, str]]:
    models = REPO / "results" / "models"
    if not models.exists():
        return []
    return sorted(
        (mf.parent.name, mf.stem)
        for mf in models.rglob("*.joblib")
    )


_MODELS = _shipped_models()


def fire_rate(feature_columns, rules: list[str], matrix) -> float:
    """Fraction of genomes carrying at least one determinant matching the rule list.

    Matching is prefix-based, exactly as the rules baseline itself scores (`blaKPC` must match
    `blaKPC-2`), so this measures what the baseline actually sees — not what we hope it sees.
    """
    matched = [c for c in feature_columns if any(c == g or c.startswith(g) for g in rules)]
    if not matched:
        return 0.0
    return float((matrix[matched].max(axis=1) > 0).mean())


@pytest.mark.skipif(not _MODELS, reason="no trained models on disk (run `make models`)")
@pytest.mark.parametrize("organism,drug", _MODELS, ids=[f"{o}/{d}" for o, d in _MODELS])
def test_rules_baseline_can_discriminate(organism: str, drug: str):
    """Every shipped rules baseline must occur in 1–99% of its own organism's genomes."""
    joblib = pytest.importorskip("joblib")
    pd = pytest.importorskip("pandas")
    from src.app.registry import ORGANISMS  # noqa: PLC0415

    meta = ORGANISMS.get(organism)
    assert meta is not None, f"{organism} has shipped models but no registry entry"

    features = REPO / meta["features"]
    if not features.exists():
        pytest.skip(f"{organism} feature matrix not built ({features.name})")

    bundle = joblib.load(REPO / "results" / "models" / organism / f"{drug}.joblib")
    rules = bundle["rules"]
    assert rules, f"{organism}/{drug} ships an EMPTY rule list — the baseline cannot be scored"

    matrix = pd.read_csv(features, index_col=0)
    rate = fire_rate(matrix.columns, rules, matrix)

    excuse = KNOWN_ABSENT.get((organism, drug))
    if excuse:
        assert rate < MIN_FIRE, (
            f"{organism}/{drug} is listed in KNOWN_ABSENT but its baseline now fires in "
            f"{rate:.1%} of genomes. The documented exception no longer applies — re-check it "
            f"and remove the entry.\n  exception on record: {excuse}"
        )
        return

    assert rate >= MIN_FIRE, (
        f"BROKEN BASELINE — {organism}/{drug} rule list {rules} occurs in {rate:.1%} of "
        f"{organism} genomes, so it can never predict 'resistant'.\n"
        f"Its ROC will be a degenerate ~0.500 that makes the model look better than it is.\n"
        f"Most likely another organism's rule list is leaking through `rules_for()` — check "
        f"ORG_RULES in src/organism_pipeline.py (see decision #57, K. pneumoniae cefoxitin).\n"
        f"If the mechanism is genuinely absent in this organism, add it to KNOWN_ABSENT with the "
        f"decision that establishes it."
    )
    assert rate <= MAX_FIRE, (
        f"BROKEN BASELINE — {organism}/{drug} rule list {rules} occurs in {rate:.1%} of "
        f"{organism} genomes, i.e. essentially all of them, so it cannot separate R from S.\n"
        f"This is the intrinsic-gene artifact from decision #42 (blaPDC / aph(3')-IIb in "
        f"P. aeruginosa). Remove the intrinsic gene from the rule list."
    )


def test_known_absent_entries_are_documented():
    """Every exception must cite the decision that justifies it — no silent escape hatches."""
    for key, reason in KNOWN_ABSENT.items():
        assert "decision #" in reason, f"{key} exception cites no decision-log entry"


def test_fire_rate_detects_both_failure_modes():
    """The detector itself: a never-present and an always-present list must both be caught."""
    pd = pytest.importorskip("pandas")
    matrix = pd.DataFrame({"blaKPC-2": [1, 0, 1, 0], "intrinsic": [1, 1, 1, 1]})
    assert fire_rate(matrix.columns, ["blaKPC"], matrix) == 0.5      # a useful baseline
    assert fire_rate(matrix.columns, ["mecA"], matrix) == 0.0        # never fires  -> caught
    assert fire_rate(matrix.columns, ["intrinsic"], matrix) == 1.0   # always fires -> caught
