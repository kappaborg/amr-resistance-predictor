# Summary #6 — Week 3: SHAP Interpretation + Biological Validation

Global SHAP (XGBoost, TreeExplainer) per drug on the held-out **unseen-lineage** test set. Top determinants checked against known mechanisms — ✓ = matches established biology (⚕ microbiology to confirm). Beeswarm plots in `results/figures/shap_<drug>.png`.


## meropenem  (4/10 top features match known mechanisms)

| Rank | Determinant | mean\|SHAP\| | pushes toward | known mechanism |
|---|---|---|---|---|
| 1 | blaKPC-2 | 1.642 | Resistant | ✓ KPC carbapenemase |
| 2 | blaKPC-3 | 0.840 | Resistant | ✓ KPC carbapenemase |
| 3 | parC_S80I | 0.721 | Resistant | — (investigate) |
| 4 | ble | 0.547 | Resistant | — (investigate) |
| 5 | blaOXA-48 | 0.373 | Resistant | ✓ OXA-48 carbapenemase |
| 6 | ompK36_D135DGD | 0.372 | Resistant | ✓ porin loss |
| 7 | oqxB | 0.256 | Susceptible | — (investigate) |
| 8 | aac(6')-Ib | 0.226 | Resistant | — (investigate) |
| 9 | blaOXA | 0.206 | Resistant | — (investigate) |
| 10 | blaTEM-1 | 0.205 | Resistant | — (investigate) |

## gentamicin  (4/10 top features match known mechanisms)

| Rank | Determinant | mean\|SHAP\| | pushes toward | known mechanism |
|---|---|---|---|---|
| 1 | aac(3)-IIe | 1.868 | Resistant | ✓ aminoglycoside acetyltransferase |
| 2 | aac(3)-IId | 1.085 | Resistant | ✓ aminoglycoside acetyltransferase |
| 3 | parC_S80I | 0.753 | Resistant | — (investigate) |
| 4 | aadA1 | 0.273 | Resistant | — (investigate) |
| 5 | blaOXA | 0.234 | Susceptible | — (investigate) |
| 6 | ant(2'')-Ia | 0.231 | Resistant | ✓ nucleotidyltransferase |
| 7 | blaSHV-1 | 0.228 | Susceptible | — (investigate) |
| 8 | floR | 0.201 | Resistant | — (investigate) |
| 9 | sul2 | 0.198 | Resistant | — (investigate) |
| 10 | aac(3)-IVa | 0.175 | Susceptible | ✓ aminoglycoside acetyltransferase |

## ciprofloxacin  (6/10 top features match known mechanisms)

| Rank | Determinant | mean\|SHAP\| | pushes toward | known mechanism |
|---|---|---|---|---|
| 1 | parC_S80I | 2.643 | Resistant | ✓ topoisomerase QRDR mutation |
| 2 | gyrA_S83I | 1.358 | Resistant | ✓ gyrase QRDR mutation |
| 3 | qnrB1 | 0.598 | Resistant | ✓ PMQR (Qnr) |
| 4 | aac(6')-Ib-cr5 | 0.543 | Resistant | ✓ PMQR (aac-cr) |
| 5 | sul1 | 0.492 | Resistant | — (investigate) |
| 6 | sul2 | 0.381 | Susceptible | — (investigate) |
| 7 | gyrA_S83F | 0.330 | Resistant | ✓ gyrase QRDR mutation |
| 8 | dfrA14 | 0.283 | Resistant | — (investigate) |
| 9 | qnrS1 | 0.283 | Resistant | ✓ PMQR (Qnr) |
| 10 | aadA2 | 0.278 | Resistant | — (investigate) |

## trimethoprim_sulfamethoxazole  (4/10 top features match known mechanisms)

| Rank | Determinant | mean\|SHAP\| | pushes toward | known mechanism |
|---|---|---|---|---|
| 1 | sul1 | 1.334 | Resistant | ✓ sulfonamide resistance |
| 2 | dfrA14 | 1.185 | Resistant | ✓ dihydrofolate reductase |
| 3 | dfrA12 | 1.129 | Resistant | ✓ dihydrofolate reductase |
| 4 | sul2 | 0.866 | Resistant | ✓ sulfonamide |
| 5 | mrx(A) | 0.439 | Resistant | — (investigate) |
| 6 | aph(6)-Id | 0.423 | Resistant | — (investigate) |
| 7 | parC_S80I | 0.392 | Resistant | — (investigate) |
| 8 | aac(3)-IIe | 0.369 | Resistant | — (investigate) |
| 9 | blaTEM-1 | 0.349 | Susceptible | — (investigate) |
| 10 | catA1 | 0.334 | Resistant | — (investigate) |

## cefoxitin  (1/10 top features match known mechanisms)

| Rank | Determinant | mean\|SHAP\| | pushes toward | known mechanism |
|---|---|---|---|---|
| 1 | blaKPC-2 | 1.205 | Resistant | — (investigate) |
| 2 | ble | 0.596 | Resistant | — (investigate) |
| 3 | ompK36_D135DGD | 0.452 | Resistant | ✓ porin loss |
| 4 | parC_S80I | 0.408 | Resistant | — (investigate) |
| 5 | blaKPC-3 | 0.254 | Resistant | — (investigate) |
| 6 | ramR_A19V | 0.245 | Susceptible | — (investigate) |
| 7 | gyrA_S83F | 0.234 | Resistant | — (investigate) |
| 8 | blaCTX-M-15 | 0.214 | Susceptible | — (investigate) |
| 9 | mph(A) | 0.202 | Resistant | — (investigate) |
| 10 | dfrA12 | 0.170 | Susceptible | — (investigate) |
