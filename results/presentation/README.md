# Presentation package

Everything needed to deliver a detailed talk on **Reading Resistance**.

| File | What it is | How to use |
|---|---|---|
| **`slides.html`** | English visual slide deck (51 slides, embeds the real result figures). | Double-click to open in any browser. **← →** navigate · **F** fullscreen · **S** toggle speaker notes · **P** print → save as PDF. |
| **`slides_zh.html`** | 中文幻灯片（简体，51 页，结构与英文版一致）。 | 浏览器打开。**← →** 翻页 · **F** 全屏 · **S** 演讲者备注 · **P** 打印/导出 PDF。 |
| **`slides.pdf`** | English deck as a ready PDF (16:9, one slide per page, no speaker notes). | For printing / emailing / projecting without a browser. |
| **`slides_zh.pdf`** | 中文幻灯片 PDF（16:9，每页一张，不含备注）。 | 便于打印 / 发送 / 投影。 |
| **`SPEAKER_GUIDE.md`** | The full narrative — every slide's *"what's on screen"* + *"what to say"*, all the numbers, the code behind each step, positives/negatives, and a **Q&A defense** section. | Read once end-to-end to own the story, then present from the deck with this open beside you. (English; the Chinese decks carry per-slide speaker notes inline — press **S**.) |

> **Translation quality.** The Chinese deck was checked by two independent review passes — a
> terminology/number-fidelity audit and a back-translation/fluency audit — and their corrections
> applied (e.g. isotonic → 保序, causal genes → 作为耐药根源的). Gene names (`gyrA`, `mecA`,
> `blaKPC`…) and standard acronyms (ROC-AUC, VME, MIC, SHAP) are intentionally kept in the original.

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
