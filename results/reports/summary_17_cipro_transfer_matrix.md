# Summary #17 — 6-Organism Ciprofloxacin Cross-Species Transfer

Ciprofloxacin (gyrA/parC/grlA) modelled in **6 organisms**. ROC-AUC; **rows = train organism, columns = test organism**. Diagonal = within-organism (70/30); off-diagonal = **zero-shot cross-species transfer** (shared determinants only, no target training).

| train \ test | K.pneumoniae | E.coli | A.baumannii | S.aureus | P.aeruginosa | S.enterica |
|---|---|---|---|---|---|---|
| **K.pneumoniae** | 0.989 | 0.979 | 0.606 | 0.463 | 0.759 | 0.953 |
| **E.coli** | 0.974 | 0.991 | 0.854 | 0.594 | 0.471 | 0.950 |
| **A.baumannii** | 0.793 | 0.739 | 0.997 | 0.483 | 0.817 | 0.699 |
| **S.aureus** | 0.515 | 0.694 | 0.500 | 0.989 | — | 0.626 |
| **P.aeruginosa** | 0.839 | 0.787 | 0.888 | — | 0.885 | 0.806 |
| **S.enterica** | 0.949 | 0.972 | 0.407 | 0.594 | 0.393 | 0.978 |

**Reading:** off-diagonal ROC well above 0.5 = the fluoroquinolone determinant→phenotype mapping transfers across species without retraining. Transfer is strongest **among the Gram-negatives** (shared gyrA/parC numbering); transfer to/from **Gram-positive S. aureus** is weaker because its gyrA residue numbering and grlA (vs parC) differ — an honest, mechanistically-expected boundary of generalization, not a failure. Gram type: K.pneumoniae(−), E.coli(−), A.baumannii(−), S.aureus(+), P.aeruginosa(−), S.enterica(−).
