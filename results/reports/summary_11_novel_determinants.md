# Summary #11 — Cross-Resistance / Linked-Determinant Discovery

High-|SHAP| determinants (presence → Resistant) that are NOT the drug's own causal mechanism, per drug. **cross-resistance** = a known AMR gene for a drug/class *outside* the panel that predicts this drug's resistance → a linked-resistance lead for microbiology review (⚕).

**Honest limitation:** features are AMRFinderPlus *catalogued* determinants, so nothing here is a novel gene — this reveals linked/cross-resistance structure. Genuine novel-gene discovery needs pan-genome/k-mer features (out of scope). Leads, not conclusions.


## meropenem  (7 cross-resistance marker(s))

| Determinant | mean\|SHAP\| | classification |
|---|---|---|
| parC_S80I | 0.721 | co-selection (ciprofloxacin) |
| ble | 0.547 | **cross-resistance** |
| aac(6')-Ib | 0.226 | **cross-resistance** |
| blaOXA | 0.206 | **cross-resistance** |
| blaTEM-1 | 0.205 | **cross-resistance** |
| blaSHV-12 | 0.190 | **cross-resistance** |
| sul1 | 0.186 | co-selection (trimethoprim/sulfamethoxazole) |
| aadA1 | 0.127 | **cross-resistance** |
| gyrA_S83F | 0.121 | co-selection (ciprofloxacin) |
| qnrB1 | 0.118 | co-selection (ciprofloxacin) |
| aph(3'')-Ib | 0.076 | **cross-resistance** |
| sul2 | 0.073 | co-selection (trimethoprim/sulfamethoxazole) |

## gentamicin  (8 cross-resistance marker(s))

| Determinant | mean\|SHAP\| | classification |
|---|---|---|
| parC_S80I | 0.753 | co-selection (ciprofloxacin) |
| aadA1 | 0.273 | **cross-resistance** |
| floR | 0.201 | **cross-resistance** |
| sul2 | 0.198 | co-selection (trimethoprim/sulfamethoxazole) |
| qnrB1 | 0.156 | co-selection (ciprofloxacin) |
| blaKPC-2 | 0.152 | co-selection (meropenem) |
| qnrS1 | 0.149 | co-selection (ciprofloxacin) |
| mph(A) | 0.140 | **cross-resistance** |
| dfrA14 | 0.122 | co-selection (trimethoprim/sulfamethoxazole) |
| mrx(A) | 0.117 | **cross-resistance** |
| msr(E) | 0.116 | **cross-resistance** |
| gyrA_D87N | 0.111 | co-selection (ciprofloxacin) |

## ciprofloxacin  (6 cross-resistance marker(s))

| Determinant | mean\|SHAP\| | classification |
|---|---|---|
| sul1 | 0.492 | co-selection (trimethoprim/sulfamethoxazole) |
| dfrA14 | 0.283 | co-selection (trimethoprim/sulfamethoxazole) |
| aadA2 | 0.278 | **cross-resistance** |
| tet(A) | 0.212 | **cross-resistance** |
| catA1 | 0.169 | **cross-resistance** |
| blaTEM-1 | 0.159 | **cross-resistance** |
| blaOXA | 0.156 | **cross-resistance** |
| catB3 | 0.131 | **cross-resistance** |

## trimethoprim_sulfamethoxazole  (7 cross-resistance marker(s))

| Determinant | mean\|SHAP\| | classification |
|---|---|---|
| mrx(A) | 0.439 | **cross-resistance** |
| aph(6)-Id | 0.423 | **cross-resistance** |
| parC_S80I | 0.392 | co-selection (ciprofloxacin) |
| aac(3)-IIe | 0.369 | co-selection (gentamicin) |
| catA1 | 0.334 | **cross-resistance** |
| oqxB19 | 0.271 | co-selection (ciprofloxacin) |
| qnrB1 | 0.227 | co-selection (ciprofloxacin) |
| blaKPC-2 | 0.129 | co-selection (meropenem) |
| blaSHV-11 | 0.128 | **cross-resistance** |
| blaSHV-12 | 0.128 | **cross-resistance** |
| aadA1 | 0.115 | **cross-resistance** |
| ompK35_E42RfsTer47 | 0.112 | co-selection (meropenem) |

## cefoxitin  (8 cross-resistance marker(s))

| Determinant | mean\|SHAP\| | classification |
|---|---|---|
| blaKPC-2 | 1.205 | co-selection (meropenem) |
| ble | 0.596 | **cross-resistance** |
| parC_S80I | 0.408 | co-selection (ciprofloxacin) |
| blaKPC-3 | 0.254 | co-selection (meropenem) |
| gyrA_S83F | 0.234 | co-selection (ciprofloxacin) |
| mph(A) | 0.202 | **cross-resistance** |
| qnrB1 | 0.161 | co-selection (ciprofloxacin) |
| blaNDM-1 | 0.153 | co-selection (meropenem) |
| blaOXA | 0.148 | **cross-resistance** |
| blaOXA-48 | 0.139 | co-selection (meropenem) |
| sul2 | 0.133 | co-selection (trimethoprim/sulfamethoxazole) |
| gyrA_S83I | 0.132 | co-selection (ciprofloxacin) |
