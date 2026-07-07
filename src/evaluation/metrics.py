"""Clinical + discrimination metrics for AMR prediction.

Very-major error (VME) = a truly Resistant strain called Susceptible — the dangerous mistake.
Major error (ME)      = a truly Susceptible strain called Resistant.
These lead every results table (per the proposal), alongside ROC-AUC and PR-AUC.

Convention: positive class = Resistant.
"""
from __future__ import annotations

from sklearn.metrics import average_precision_score, roc_auc_score


def clinical_rates(y_true: list[int], y_pred: list[int]) -> dict:
    """y_* are 1=Resistant, 0=Susceptible. VME/ME are rates over the relevant true class."""
    n_R = sum(1 for t in y_true if t == 1)
    n_S = sum(1 for t in y_true if t == 0)
    vme = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 0)  # R called S
    me = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 1)   # S called R
    return {
        "very_major_error_rate": vme / n_R if n_R else float("nan"),  # missed resistance
        "major_error_rate": me / n_S if n_S else float("nan"),
        "vme_count": vme, "me_count": me, "n_resistant": n_R, "n_susceptible": n_S,
    }


def discrimination(y_true: list[int], y_score: list[float]) -> dict:
    return {
        "roc_auc": roc_auc_score(y_true, y_score),
        "pr_auc": average_precision_score(y_true, y_score),  # positive = Resistant
    }


def full_report(y_true, y_pred, y_score) -> dict:
    out = {}
    out.update(discrimination(y_true, y_score))
    out.update(clinical_rates(y_true, y_pred))
    return out
