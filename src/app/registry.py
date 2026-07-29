"""Shared deployment registry — the eight organisms the demo serves.

One place that maps an organism key to everything the demo/model-saver needs: display name, NCBI
taxon, the AMRFinderPlus `--organism` string (note E. coli is "Escherichia", NOT "Escherichia_coli"),
the processed feature/lineage files, the label source (K. pneumoniae has a saved panel; the others are
fetched live), and the per-drug BV-BRC antibiotic encodings. `save_models`, `predict`, and the
Streamlit app all import this so they never drift out of sync.
"""
from __future__ import annotations

import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def load_dotenv(path: Path = REPO / ".env") -> None:
    """Minimal zero-dependency .env loader. Reads KEY=value lines from the project-root .env and
    populates os.environ WITHOUT overriding variables already set in the real shell (real export
    wins, then .env, matching standard precedence). Supports `export KEY=...`, comments, and quotes.
    Called once at import so every entry point (predict / report / streamlit / save_models) sees it.
    """
    if not path.exists():
        return
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[len("export "):]
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        k, v = k.strip(), v.strip().strip('"').strip("'")
        if k:
            os.environ.setdefault(k, v)   # real shell export takes precedence over .env


load_dotenv()

ORGANISMS: dict[str, dict] = {
    "kpneu": {
        "display": "Klebsiella pneumoniae", "taxon": 573, "amrfinder": "Klebsiella_pneumoniae",
        "features": "data/processed/thin_slice_cipro_features.csv",
        "lineages": "data/processed/thin_slice_cipro_lineages.csv",
        "labels_file": "data/processed/panel_labels.csv",       # saved panel (columns = drug names)
        "schemes": ["klebsiella"],                              # MLST scheme(s) that identify this species
        "markers": ["oqxA", "oqxB"],                           # intrinsic (~91%); fosA dropped — also ~99% in P. aeruginosa
        "drugs": {"meropenem": "meropenem", "gentamicin": "gentamicin",
                  "ciprofloxacin": "ciprofloxacin",
                  "trimethoprim_sulfamethoxazole": "trimethoprim%2Fsulfamethoxazole",
                  "cefoxitin": "cefoxitin"},
    },
    "ecoli": {
        "display": "Escherichia coli", "taxon": 562, "amrfinder": "Escherichia",
        "features": "data/processed/ecoli_features.csv",
        "lineages": "data/processed/ecoli_lineages.csv",
        "labels_file": None,                                     # fetched live per drug
        "schemes": ["ecoli", "ecoli_achtman_4", "escherichia"],
        "markers": ["uhpT_E350Q", "gyrA_S83L", "ptsI_V25I"],   # E. coli-specific (weak positives; 0% in others)
        "drugs": {"ciprofloxacin": "ciprofloxacin", "gentamicin": "gentamicin",
                  "trimethoprim_sulfamethoxazole": "trimethoprim%2Fsulfamethoxazole",
                  "ampicillin": "ampicillin", "ceftriaxone": "ceftriaxone"},
    },
    "abaumannii": {
        "display": "Acinetobacter baumannii", "taxon": 470, "amrfinder": "Acinetobacter_baumannii",
        "features": "data/processed/abaumannii_features.csv",
        "lineages": "data/processed/abaumannii_lineages.csv",
        "labels_file": None,
        "schemes": ["abaumannii", "abaumannii_2"],
        "markers": ["blaADC"],                                 # Acinetobacter-derived cephalosporinase (~99%, species-specific)
        "drugs": {"meropenem": "meropenem", "imipenem": "imipenem",
                  "ciprofloxacin": "ciprofloxacin", "gentamicin": "gentamicin",
                  "amikacin": "amikacin"},
    },
    "saureus": {
        "display": "Staphylococcus aureus", "taxon": 1280, "amrfinder": "Staphylococcus_aureus",
        "features": "data/processed/saureus_features.csv",
        "lineages": "data/processed/saureus_lineages.csv",
        "labels_file": None,
        "schemes": ["saureus"],
        "markers": ["blaZ", "blaI", "blaR1"],                  # staph beta-lactamase operon (~90%, 0% in Gram-negatives)
        "drugs": {"oxacillin": "oxacillin", "cefoxitin": "cefoxitin",
                  "ciprofloxacin": "ciprofloxacin", "erythromycin": "erythromycin",
                  "clindamycin": "clindamycin"},
    },
    "paeruginosa": {
        "display": "Pseudomonas aeruginosa", "taxon": 287, "amrfinder": "Pseudomonas_aeruginosa",
        "features": "data/processed/paeruginosa_features.csv",
        "lineages": "data/processed/paeruginosa_lineages.csv",
        "labels_file": None,
        "schemes": ["paeruginosa"],
        "markers": ["aph(3')-IIb", "catB7", "nalC_G71E"],      # intrinsic P. aeruginosa (~95-99%, 0% elsewhere)
        "drugs": {"meropenem": "meropenem", "ceftazidime": "ceftazidime",
                  "ciprofloxacin": "ciprofloxacin", "tobramycin": "tobramycin"},
    },
    "senterica": {
        "display": "Salmonella enterica", "taxon": 28901, "amrfinder": "Salmonella",
        "features": "data/processed/senterica_features.csv",
        "lineages": "data/processed/senterica_lineages.csv",
        "labels_file": None,
        "schemes": ["salmonella"],
        "markers": [],   # no near-universal Salmonella-specific determinant; species relies on MLST
        "drugs": {"ampicillin": "ampicillin", "ceftriaxone": "ceftriaxone",
                  "ciprofloxacin": "ciprofloxacin", "chloramphenicol": "chloramphenicol",
                  "trimethoprim_sulfamethoxazole": "trimethoprim%2Fsulfamethoxazole"},
    },
    "efaecium": {
        "display": "Enterococcus faecium", "taxon": 1352, "amrfinder": "Enterococcus_faecium",
        "features": "data/processed/efaecium_features.csv",
        "lineages": "data/processed/efaecium_lineages.csv",
        "labels_file": None,
        "schemes": ["efaecium"],
        "markers": ["msr(C)", "pbp5_N496K"],   # intrinsic E. faecium (~98/91%, 0% in other panel orgs)
        "drugs": {"vancomycin": "vancomycin", "ampicillin": "ampicillin",
                  "tetracycline": "tetracycline"},
    },
    "spneumoniae": {
        "display": "Streptococcus pneumoniae", "taxon": 1313, "amrfinder": "Streptococcus_pneumoniae",
        "features": "data/processed/spneumoniae_features.csv",
        "lineages": "data/processed/spneumoniae_lineages.csv",
        "labels_file": None,
        "schemes": ["spneumoniae"],
        "markers": [],   # no near-universal S. pneumoniae-specific determinant; species relies on MLST
        "drugs": {"penicillin": "penicillin", "erythromycin": "erythromycin",
                  "tetracycline": "tetracycline",
                  "trimethoprim_sulfamethoxazole": "trimethoprim%2Fsulfamethoxazole"},
    },
}

# Display order for the UI
ORG_ORDER = ["kpneu", "ecoli", "abaumannii", "saureus", "paeruginosa", "senterica", "efaecium",
             "spneumoniae"]

# Reverse map: MLST scheme -> organism key (for species verification of uploads)
SCHEME_ORG = {s: k for k, o in ORGANISMS.items() for s in o["schemes"]}
