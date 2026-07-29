"""Foundation-model experiment #1 — isolate-similarity Graph Neural Network (AMR-GNN style).

Motivation. The published AMR-GNN framework (Nat. Commun. 2026) does NOT build a gene graph; it
treats **each isolate (genome) as a node**, connects isolates by genomic similarity, and does node
classification. The graph's job is to model **population / clonal structure directly** — the exact
confounder this project is built around (lineage-held-out splits). Its documented gains were largest
on the *hardest* drugs (cefepime +28.8%, aztreonam +18.9% AUROC). That maps onto our hard cases
(meropenem, cefoxitin), so it is worth an honest test.

The honest test. For each drug, under the SAME phylogeny-aware (lineage-grouped) CV split, we compare
three models on identical folds:
  1. XGBoost   — gradient-boosted trees on the determinant matrix (our inductive baseline).
  2. MLP       — a 2-layer neural net on the SAME features, NO graph (controls for "is it just a NN?").
  3. GCN       — a 2-layer graph conv net: features + isolate kNN-similarity graph (transductive).
If GCN > MLP, the graph structure adds value *beyond* the features. If GCN ~= XGBoost, trees already
capture it. Either way it is a real, reportable result — not a claim we assume.

Honesty notes carried into the report:
  * The GCN is **transductive** (it sees every node's FEATURES and the graph at train time, but only
    TRAIN-lineage LABELS enter the loss). XGBoost/MLP are inductive. We state this explicitly.
  * Lineage grouping means every member of a lineage shares a fold, so there is no trivial
    same-lineage label leak through the graph.

    python -m src.models.gnn_resistance --config config/config.yaml
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")  # benign OpenMP dup on this env

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import yaml
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from xgboost import XGBClassifier

REPO = Path(__file__).resolve().parents[2]
FEATURES = REPO / "data/processed/thin_slice_cipro_features.csv"
LINEAGES = REPO / "data/processed/thin_slice_cipro_lineages.csv"
PANEL = REPO / "data/processed/panel_labels.csv"
OUT = REPO / "results/reports/summary_20_gnn.md"
DRUGS = ["meropenem", "gentamicin", "ciprofloxacin", "trimethoprim_sulfamethoxazole", "cefoxitin"]

DEV = torch.device("mps" if torch.backends.mps.is_available() else "cpu")


def knn_adjacency(X: np.ndarray, k: int = 10) -> torch.Tensor:
    """Symmetric-normalized adjacency of a kNN graph over isolates (cosine similarity of features).

    Returns dense  Â = D̃^{-1/2}(A+I)D̃^{-1/2}  as a float tensor on DEV.
    """
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-9)
    sim = Xn @ Xn.T                                   # cosine similarity, [n, n]
    np.fill_diagonal(sim, -1.0)                        # exclude self for neighbour selection
    nn_idx = np.argpartition(-sim, k, axis=1)[:, :k]   # top-k neighbours per node
    n = X.shape[0]
    A = np.zeros((n, n), dtype=np.float32)
    rows = np.repeat(np.arange(n), k)
    A[rows, nn_idx.ravel()] = 1.0
    A = np.maximum(A, A.T)                              # symmetrise
    A[np.diag_indices(n)] = 1.0                         # add self-loops (A + I)
    deg = A.sum(1)
    dinv = 1.0 / np.sqrt(deg)
    Ahat = (A * dinv[:, None]) * dinv[None, :]         # D^-1/2 (A+I) D^-1/2
    return torch.tensor(Ahat, dtype=torch.float32, device=DEV)


class GCN(nn.Module):
    def __init__(self, d_in: int, d_hid: int = 128, p: float = 0.4):
        super().__init__()
        self.w0 = nn.Linear(d_in, d_hid)
        self.w1 = nn.Linear(d_hid, 2)
        self.p = p

    def forward(self, ahat: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        h = F.relu(self.w0(ahat @ x))
        h = F.dropout(h, self.p, self.training)
        return self.w1(ahat @ h)


class MLP(nn.Module):
    def __init__(self, d_in: int, d_hid: int = 128, p: float = 0.4):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(d_in, d_hid), nn.ReLU(), nn.Dropout(p),
                                 nn.Linear(d_hid, 2))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


def train_torch(model, forward_fn, x, y, tr, te, epochs=300, lr=5e-3, wd=5e-4):
    """Generic trainer: loss only on `tr` mask; returns P(resistant) on `te`. Class-weighted."""
    model = model.to(DEV)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    ytr = y[tr]
    w = torch.tensor([1.0 / max(1, (ytr == 0).sum()), 1.0 / max(1, (ytr == 1).sum())],
                     dtype=torch.float32, device=DEV)
    w = w / w.sum() * 2
    yt = torch.tensor(y, dtype=torch.long, device=DEV)
    tr_t = torch.tensor(tr, dtype=torch.long, device=DEV)
    for _ in range(epochs):
        model.train()
        opt.zero_grad()
        logits = forward_fn(model, x)
        loss = F.cross_entropy(logits[tr_t], yt[tr_t], weight=w)
        loss.backward()
        opt.step()
    model.eval()
    with torch.no_grad():
        prob = F.softmax(forward_fn(model, x), dim=1)[:, 1].cpu().numpy()
    return prob[te]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--k", type=int, default=10, help="kNN graph degree")
    args = ap.parse_args()
    cfg = yaml.safe_load((REPO / args.config).read_text())
    seed = cfg["project"]["seed"]
    torch.manual_seed(seed)
    np.random.seed(seed)

    feats = pd.read_csv(FEATURES, dtype={"genome_id": str}).set_index("genome_id")
    lin = pd.read_csv(LINEAGES, dtype={"genome_id": str}).set_index("genome_id")["lineage"].to_dict()
    panel = pd.read_csv(PANEL, dtype={"genome_id": str}).set_index("genome_id")

    print(f"device: {DEV} | features: {feats.shape} | kNN k={args.k} | {args.folds}-fold lineage CV\n")
    print(f"{'drug':32s}{'n':>6}{'%R':>6}{'XGB':>8}{'MLP':>8}{'GCN':>8}{'GCN-MLP':>9}")

    results = {}
    for drug in DRUGS:
        lab = panel[drug]
        g = [x for x in feats.index if lab.get(x) in ("Resistant", "Susceptible")]
        if len(g) < 200:
            print(f"{drug:32s} too few labels ({len(g)})")
            continue
        Xnp = feats.loc[g].values.astype(np.float32)
        y = np.array([1 if lab[x] == "Resistant" else 0 for x in g])
        groups = np.array([lin.get(x, x) for x in g])
        pctR = y.mean()

        # one adjacency per drug's labelled subgraph; feature tensor shared across folds
        ahat = knn_adjacency(Xnp, k=args.k)
        Xt = torch.tensor(Xnp, device=DEV)

        skf = StratifiedGroupKFold(n_splits=args.folds, shuffle=True, random_state=seed)
        auc = {"xgb": [], "mlp": [], "gcn": []}
        for tr, te in skf.split(Xnp, y, groups):
            if len(set(y[te])) < 2 or len(set(y[tr])) < 2:
                continue
            xgb = XGBClassifier(n_estimators=400, max_depth=4, learning_rate=0.08, subsample=0.9,
                                colsample_bytree=0.8, eval_metric="logloss", n_jobs=4,
                                random_state=seed).fit(Xnp[tr], y[tr])
            auc["xgb"].append(roc_auc_score(y[te], xgb.predict_proba(Xnp[te])[:, 1]))
            auc["mlp"].append(roc_auc_score(y[te], train_torch(
                MLP(Xnp.shape[1]), lambda m, x: m(x), Xt, y, tr, te)))
            auc["gcn"].append(roc_auc_score(y[te], train_torch(
                GCN(Xnp.shape[1]), lambda m, x: m(ahat, x), Xt, y, tr, te)))

        m = {k: float(np.mean(v)) for k, v in auc.items()}
        s = {k: float(np.std(v)) for k, v in auc.items()}
        results[drug] = {"n": len(g), "pct_resistant": float(pctR),
                         "xgb_auc": m["xgb"], "mlp_auc": m["mlp"], "gcn_auc": m["gcn"],
                         "xgb_std": s["xgb"], "mlp_std": s["mlp"], "gcn_std": s["gcn"],
                         "gcn_minus_mlp": m["gcn"] - m["mlp"], "gcn_minus_xgb": m["gcn"] - m["xgb"]}
        print(f"{drug:32s}{len(g):>6}{pctR:>6.0%}{m['xgb']:>8.3f}{m['mlp']:>8.3f}"
              f"{m['gcn']:>8.3f}{m['gcn']-m['mlp']:>+9.3f}")

    # ---- report ----
    md = ["# Summary #20 — Isolate-Similarity GNN (AMR-GNN-style)\n",
          f"**Date:** 2026-07-09 · device `{DEV}` · {args.folds}-fold lineage-grouped CV · "
          f"kNN graph k={args.k} on cosine feature similarity.\n",
          "Each **isolate is a node**; edges join genomically similar isolates. A 2-layer GCN does "
          "node classification (R/S), modelling population structure directly — the AMR-GNN idea "
          "(Nat. Commun. 2026), whose largest gains were on the hardest drugs. Compared, on identical "
          "phylogeny-aware folds, against XGBoost (our baseline) and an MLP on the same features with "
          "**no graph** (the control).\n",
          "| Drug | n | %R | XGBoost | MLP (no graph) | **GCN (+graph)** | GCN − MLP |",
          "|---|---|---|---|---|---|---|"]
    for d, r in results.items():
        md.append(f"| {d} | {r['n']} | {r['pct_resistant']:.0%} | {r['xgb_auc']:.3f} | "
                  f"{r['mlp_auc']:.3f} | **{r['gcn_auc']:.3f}** | {r['gcn_minus_mlp']:+.3f} |")
    md += ["\n**Reading (mean unseen-lineage ROC-AUC).** *GCN − MLP* isolates the value of the graph: "
           "both use the identical determinant features and folds, so a positive gap means the "
           "isolate-similarity structure adds signal the features alone don't carry. *GCN vs XGBoost* "
           "asks whether the graph model beats our strong tree baseline.\n",
           "**Honesty.** The GCN is **transductive** — it sees all nodes' *features* and the graph "
           "during training, but only training-lineage *labels* enter the loss; XGBoost and the MLP "
           "are inductive. Lineage grouping puts every member of a lineage in one fold, so no trivial "
           "same-lineage label leaks across the graph. This is the same evaluation contract as the "
           "rest of the project (unseen lineages), applied fairly to all three models."]
    OUT.write_text("\n".join(md) + "\n")
    (REPO / "results/metrics/gnn_metrics.json").write_text(json.dumps(results, indent=2))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
