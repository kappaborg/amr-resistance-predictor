"""Phase 10 bonus — cross-resistance / linked-determinant discovery.

Flags determinants the model leans on (high mean |SHAP|, presence pushing toward Resistant) that are
NOT the drug's own causal mechanism. Each is classified:
  - own mechanism      : a known determinant for THIS drug (expected; sanity check, excluded below).
  - co-selection (X)   : a known determinant for ANOTHER panel drug — co-carried in MDR lineages.
  - cross-resistance   : a known AMR gene for a drug/class OUTSIDE the panel (e.g. tet, cat, mph, arr)
                         that statistically predicts this drug's resistance → a linked-resistance
                         marker for microbiology review (⚕).

HONEST LIMITATION: features come from AMRFinderPlus, a *curated catalog of KNOWN determinants*, so no
feature here can be a genuinely novel gene — this surfaces linked/cross-resistance structure, not new
mechanisms. True novel-gene discovery would need pan-genome / k-mer features (explicitly out of scope
for the sprint). Candidates are LEADS, not conclusions.

Usage:
    python -m src.interpret.novel_determinants --config config/config.yaml
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import shap
import yaml
from sklearn.model_selection import StratifiedGroupKFold
from xgboost import XGBClassifier

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.interpret.shap_analysis import KNOWN  # noqa: E402
from src.models.train_panel import RULES  # noqa: E402

FEATURES = REPO / "data/processed/thin_slice_cipro_features.csv"
LINEAGES = REPO / "data/processed/thin_slice_cipro_lineages.csv"
PANEL = REPO / "data/processed/panel_labels.csv"
OUT = REPO / "results/reports/summary_11_novel_determinants.md"

# Union of every determinant token known to the panel (rules + curated mechanisms), per drug.
KNOWN_BY_DRUG = {d: set(RULES.get(d, [])) | set(KNOWN.get(d, {})) for d in RULES}


def classify(feature: str, drug: str) -> str:
    if any(tok in feature for tok in KNOWN_BY_DRUG[drug]):
        return "own mechanism"
    others = [d for d in RULES if d != drug and any(tok in feature for tok in KNOWN_BY_DRUG[d])]
    if others:
        return f"co-selection ({others[0].replace('_', '/')})"
    return "cross-resistance"  # known AMR gene for a non-panel drug/class


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--top", type=int, default=25, help="features per drug to inspect")
    args = ap.parse_args()
    cfg = yaml.safe_load((REPO / args.config).read_text())
    seed, test_frac = cfg["project"]["seed"], cfg["split"]["test_fraction"]

    feats = pd.read_csv(FEATURES, dtype={"genome_id": str}).set_index("genome_id")
    lineage = {r["genome_id"]: r["lineage"] for r in csv.DictReader(LINEAGES.open())}
    labels = pd.read_csv(PANEL, dtype=str)
    cols = list(feats.columns)

    md = ["# Summary #11 — Cross-Resistance / Linked-Determinant Discovery\n",
          "High-|SHAP| determinants (presence → Resistant) that are NOT the drug's own causal "
          "mechanism, per drug. **cross-resistance** = a known AMR gene for a drug/class *outside* "
          "the panel that predicts this drug's resistance → a linked-resistance lead for microbiology "
          "review (⚕).\n\n"
          "**Honest limitation:** features are AMRFinderPlus *catalogued* determinants, so nothing "
          "here is a novel gene — this reveals linked/cross-resistance structure. Genuine novel-gene "
          "discovery needs pan-genome/k-mer features (out of scope). Leads, not conclusions.\n"]
    all_novel = {}
    for drug in RULES:
        sub = labels[labels[drug].isin(["Resistant", "Susceptible"])]
        sub = sub[sub["genome_id"].isin(feats.index)]
        g = sub["genome_id"].values
        X = feats.loc[g].values.astype(int)
        y = (sub[drug].values == "Resistant").astype(int)
        groups = np.array([lineage[x] for x in g])
        tr, te = next(StratifiedGroupKFold(max(2, round(1 / test_frac)), shuffle=True,
                                           random_state=seed).split(X, y, groups))
        spw = (y[tr] == 0).sum() / max((y[tr] == 1).sum(), 1)
        model = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.1, subsample=0.9,
                              colsample_bytree=0.8, eval_metric="logloss", scale_pos_weight=spw,
                              random_state=seed, n_jobs=4).fit(X[tr], y[tr])
        sv = shap.TreeExplainer(model).shap_values(X[te])
        mean_abs = np.abs(sv).mean(0)
        order = np.argsort(mean_abs)[::-1]

        rows, novel = [], []
        for i in order[:args.top]:
            present = X[te][:, i] == 1
            direction = sv[present, i].mean() if present.any() else 0.0
            if direction <= 0:  # only presence-pushes-Resistant candidates
                continue
            cls = classify(cols[i], drug)
            if cls == "own mechanism":
                continue
            rows.append((cols[i], float(mean_abs[i]), cls))
            if cls == "cross-resistance":
                novel.append(cols[i])
        all_novel[drug] = novel
        md.append(f"\n## {drug}  ({len(novel)} cross-resistance marker(s))\n")
        md.append("| Determinant | mean\\|SHAP\\| | classification |")
        md.append("|---|---|---|")
        for feat, ms, cls in rows[:12]:
            mark = f"**{cls}**" if cls == "cross-resistance" else cls
            md.append(f"| {feat} | {ms:.3f} | {mark} |")
        print(f"{drug:32s} cross-resistance markers: {', '.join(novel) if novel else '(none)'}")

    OUT.write_text("\n".join(md) + "\n")
    (REPO / "results/metrics/novel_determinants.json").write_text(json.dumps(all_novel, indent=2))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
