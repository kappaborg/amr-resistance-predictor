# Presentation package

Everything needed to deliver a detailed talk on **Reading Resistance**.

| File | What it is | How to use |
|---|---|---|
| **`slides.html`** | Self-contained visual slide deck (51 slides, embeds the real result figures). | Double-click to open in any browser. **← →** navigate · **F** fullscreen · **S** toggle speaker notes · **P** print → save as PDF. |
| **`SPEAKER_GUIDE.md`** | The full narrative — every slide's *"what's on screen"* + *"what to say"*, all the numbers, the code behind each step, positives/negatives, and a **Q&A defense** section. | Read once end-to-end to own the story, then present from `slides.html` with this open beside you. |

## Notes
- Open `slides.html` from **inside this folder** (or from the repo) so its `../figures/*.png`
  references resolve. If a figure looks missing, regenerate figures with `make figures` from the repo root.
- The deck and guide use **only verified numbers** from `results/` and the code, as of the handoff.
- ⚕ in the guide marks slides that are the **microbiologist's** home turf (organism choice, label QC,
  the VME/ME clinical framing, and the biological validation of SHAP determinants).
- Deeper source material: `results/reports/REPORT.md` (full write-up), `docs/decisions.md` (53-entry
  lab notebook), `data/manifest.md` (provenance), `results/poster/poster.html` (A0 poster).

## To export a PDF handout
Open `slides.html`, press **P** (or Ctrl/Cmd-P). In the print dialog choose "Save as PDF" and,
for a slide-per-page layout, set margins to *None* and enable background graphics. Print mode
automatically stacks all slides and shows the speaker notes under each.
