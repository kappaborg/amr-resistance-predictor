"""Meropenem MIC: 30-paired-fold comparison of determinants vs +ESM-2 150M vs +ESM-2 650M.

Same paired CV protocol used to validate the 150M gain (6 repeats x 5 lineage-grouped folds,
identical folds across feature sets), so the two model sizes are compared on exactly the same splits.
"""
from __future__ import annotations
import os, sys, pickle, hashlib
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import numpy as np, pandas as pd
from pathlib import Path
from sklearn.model_selection import StratifiedGroupKFold
from xgboost import XGBRegressor
from scipy import stats

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from src.models.esm2_mic import genome_plm_vectors, FEATURES, LINEAGES, emb_path
from src.models.mic_regression import fetch_mic

proteins = pickle.loads((REPO / "data/interim/esm/proteins.pkl").read_bytes())
feats = pd.read_csv(FEATURES, dtype={"genome_id": str}).set_index("genome_id")
lin = pd.read_csv(LINEAGES, dtype={"genome_id": str}).set_index("genome_id")["lineage"].to_dict()
mic = fetch_mic("meropenem")

embs = {k: pickle.loads(emb_path(k).read_bytes()) for k in ("150M", "650M")}
# common genome set: has MIC + has all-AMR embedding under BOTH models
plm = {}
for k in ("150M", "650M"):
    allv, drugv = genome_plm_vectors(proteins, embs[k], "meropenem")
    plm[k] = (allv, drugv)
g = [x for x in feats.index if x in mic and x in plm["150M"][0] and x in plm["650M"][0]]
Xd = feats.loc[g].values.astype(np.float32)
y = np.array([mic[x] for x in g], np.float32)
gr = np.array([lin.get(x, x) for x in g])
yb = (y >= np.median(y)).astype(int)
def plmX(k):
    allv, drugv = plm[k]
    return np.vstack([np.concatenate([allv[x], drugv[x]]) for x in g]).astype(np.float32)
feat_sets = {"determ": Xd,
             "determ+150M": np.hstack([Xd, plmX("150M")]),
             "determ+650M": np.hstack([Xd, plmX("650M")])}

def ea(X, tr, te):
    m = XGBRegressor(n_estimators=400, max_depth=4, learning_rate=0.08, subsample=0.9,
                     colsample_bytree=0.8, random_state=0, n_jobs=4).fit(X[tr], y[tr])
    return np.mean(np.abs(m.predict(X[te]) - y[te]) <= 1.0)

scores = {k: [] for k in feat_sets}
for seed in range(6):
    for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(Xd, yb, gr):
        for k, X in feat_sets.items():
            scores[k].append(ea(X, tr, te))
for k in scores:
    scores[k] = np.array(scores[k])

print(f"meropenem  n={len(g)}  |  30 paired folds (6 repeats x 5), identical folds across feature sets\n")
print(f"{'feature set':16s}{'EA mean':>10}{'EA std':>9}{'vs determ Δ':>13}{'Wilcoxon p':>13}")
base = scores["determ"]
for k in feat_sets:
    d = scores[k] - base
    p = stats.wilcoxon(d).pvalue if k != "determ" and np.any(d != 0) else float("nan")
    print(f"{k:16s}{scores[k].mean():>9.1%}{scores[k].std():>9.1%}{d.mean():>+12.1%}{p:>13.4f}")
# head-to-head 650M vs 150M
d = scores["determ+650M"] - scores["determ+150M"]
p = stats.wilcoxon(d).pvalue if np.any(d != 0) else float("nan")
print(f"\n650M vs 150M (paired): mean Δ {d.mean():+.1%} | 650M wins {np.mean(d>0):.0%} of folds | Wilcoxon p={p:.4f}")
# since our test result/scores for 650M is not improving our 150M result no need to use 650M