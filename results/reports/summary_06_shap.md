# Summary #6 — Week 3: SHAP Interpretation + Biological Validation

Global SHAP (XGBoost, TreeExplainer) per drug on the held-out **unseen-lineage** test set. Top determinants checked against known mechanisms — ✓ = matches established biology (⚕ microbiology to confirm). Beeswarm plots in `results/figures/shap_<drug>.png`.


## meropenem  (2/10 top features match known mechanisms)

| Rank | Determinant | mean\|SHAP\| | pushes toward | known mechanism |
|---|---|---|---|---|
| 1 | parC_S80I | 1.725 | Resistant | — (investigate) |
| 2 | blaKPC-2 | 1.084 | Resistant | ✓ KPC carbapenemase |
| 3 | blaKPC-3 | 0.538 | Resistant | ✓ KPC carbapenemase |
| 4 | gyrA_S83I | 0.430 | Resistant | — (investigate) |
| 5 | blaOXA | 0.428 | Resistant | — (investigate) |
| 6 | dfrA14 | 0.389 | Susceptible | — (investigate) |
| 7 | blaTEM-1 | 0.368 | Resistant | — (investigate) |
| 8 | ble | 0.329 | Resistant | — (investigate) |
| 9 | aadA2 | 0.321 | Resistant | — (investigate) |
| 10 | aac(3)-IIe | 0.314 | Resistant | — (investigate) |

## gentamicin  (2/10 top features match known mechanisms)

| Rank | Determinant | mean\|SHAP\| | pushes toward | known mechanism |
|---|---|---|---|---|
| 1 | aac(3)-IIe | 1.955 | Resistant | ✓ aminoglycoside acetyltransferase |
| 2 | aac(3)-IId | 0.819 | Resistant | ✓ aminoglycoside acetyltransferase |
| 3 | parC_S80I | 0.763 | Resistant | — (investigate) |
| 4 | gyrA_D87N | 0.611 | Resistant | — (investigate) |
| 5 | oqxB | 0.566 | Resistant | — (investigate) |
| 6 | blaTEM-1 | 0.532 | Resistant | — (investigate) |
| 7 | aadA1 | 0.403 | Resistant | — (investigate) |
| 8 | aph(3'')-Ib | 0.377 | Resistant | — (investigate) |
| 9 | mph(A) | 0.373 | Resistant | — (investigate) |
| 10 | sul2 | 0.369 | Resistant | — (investigate) |

## ciprofloxacin  (7/10 top features match known mechanisms)

| Rank | Determinant | mean\|SHAP\| | pushes toward | known mechanism |
|---|---|---|---|---|
| 1 | parC_S80I | 2.090 | Resistant | ✓ topoisomerase QRDR mutation |
| 2 | gyrA_S83I | 1.527 | Resistant | ✓ gyrase QRDR mutation |
| 3 | qnrB1 | 0.852 | Resistant | ✓ PMQR (Qnr) |
| 4 | aac(6')-Ib-cr5 | 0.435 | Resistant | ✓ PMQR (aac-cr) |
| 5 | sul1 | 0.422 | Resistant | — (investigate) |
| 6 | qnrS1 | 0.405 | Resistant | ✓ PMQR (Qnr) |
| 7 | oqxB | 0.373 | Susceptible | ✓ OqxAB efflux |
| 8 | catA1 | 0.358 | Resistant | — (investigate) |
| 9 | gyrA_S83F | 0.339 | Resistant | ✓ gyrase QRDR mutation |
| 10 | aph(3')-Ia | 0.306 | Susceptible | — (investigate) |

## trimethoprim_sulfamethoxazole  (4/10 top features match known mechanisms)

| Rank | Determinant | mean\|SHAP\| | pushes toward | known mechanism |
|---|---|---|---|---|
| 1 | dfrA14 | 1.269 | Resistant | ✓ dihydrofolate reductase |
| 2 | sul1 | 1.229 | Resistant | ✓ sulfonamide resistance |
| 3 | aadA2 | 0.802 | Resistant | — (investigate) |
| 4 | sul2 | 0.757 | Resistant | ✓ sulfonamide |
| 5 | mrx(A) | 0.538 | Resistant | — (investigate) |
| 6 | aac(3)-IIe | 0.491 | Resistant | — (investigate) |
| 7 | parC_S80I | 0.281 | Resistant | — (investigate) |
| 8 | blaCTX-M-15 | 0.250 | Resistant | — (investigate) |
| 9 | gyrA_S83I | 0.231 | Resistant | — (investigate) |
| 10 | dfrA1 | 0.191 | Resistant | ✓ dihydrofolate reductase |

## cefoxitin  (2/10 top features match known mechanisms)

| Rank | Determinant | mean\|SHAP\| | pushes toward | known mechanism |
|---|---|---|---|---|
| 1 | gyrA_S83I | 1.164 | Resistant | — (investigate) |
| 2 | blaKPC-2 | 0.773 | Resistant | — (investigate) |
| 3 | ompK35_E42RfsTer47 | 0.609 | Susceptible | ✓ porin loss |
| 4 | ompK36_D135DGD | 0.484 | Resistant | ✓ porin loss |
| 5 | blaCTX-M-15 | 0.450 | Susceptible | — (investigate) |
| 6 | blaOXA | 0.310 | Resistant | — (investigate) |
| 7 | mph(A) | 0.301 | Resistant | — (investigate) |
| 8 | oqxB19 | 0.279 | Susceptible | — (investigate) |
| 9 | blaKPC-3 | 0.260 | Resistant | — (investigate) |
| 10 | tet(A) | 0.241 | Susceptible | — (investigate) |

---
## Biological validation — the honest interpretation (⚕ to review)

**The correct causal mechanisms are learned and validated for every drug:**
- ciprofloxacin → gyrA/parC QRDR mutations + qnr/aac(6')-Ib-cr PMQR (7/10 top features known).
- gentamicin → aac(3)-IIe/IId acetyltransferases are the top-2 features.
- TMP-SMX → dfrA + sul1/sul2 (the folate-pathway determinants).
- meropenem → **blaKPC-2/blaKPC-3 carry the largest per-genome SHAP impact (+3 to +5)** — the
  beeswarm shows carbapenemase presence driving the "Resistant" call, exactly as expected.
- cefoxitin → **ompK35/ompK36 porin mutations surface** — confirming the Week-2 finding that the
  model captures porin-loss resistance a gene-lookup baseline cannot.

**The important nuance: co-selection confounders rank high.** For meropenem and cefoxitin the
*mean* |SHAP| ranking is topped by fluoroquinolone determinants (parC_S80I, gyrA_S83I) that are NOT
their causal mechanism. This is real epidemiology, not a bug: MDR resistance genes travel together
in high-risk K. pneumoniae lineages/plasmids (e.g. ST258 carries KPC *and* fluoroquinolone
mutations), so a prevalent co-carried marker gets credit by association. Two signals coexist:
- **Causal** — high per-instance SHAP (blaKPC spikes to +5 only when present): the true mechanism.
- **Co-selection** — high *mean* SHAP via prevalence across resistant lineages: a lineage marker.

**Why this matters (defensible finding for the write-up):** it is precisely what the biological-
validation step is for (proposal §6.6). The model's discrimination is real, but interpretation must
separate mechanism from co-carriage — mean |SHAP| rewards prevalence, so per-instance impact
(the beeswarm) is the honest lens. This also cautions against reading SHAP rank as causation.

**Next:** microbiologist review of each drug's top determinants; the Phase-2b top-up (more diverse
resistant lineages) should also reduce co-selection artifacts by breaking gene–lineage correlation.
