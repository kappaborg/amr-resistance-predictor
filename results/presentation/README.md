# Presentation package

Everything needed to deliver a detailed talk on **Reading Resistance**.

| File | What it is | How to use |
|---|---|---|
| **`slides.html`** | English visual slide deck (51 slides, embeds the real result figures). | Double-click to open in any browser. **← →** navigate · **F** fullscreen · **S** toggle speaker notes · **P** print → save as PDF. |
| **`slides_zh.html`** | 中文幻灯片（简体，51 页，结构与英文版一致）。 | 浏览器打开。**← →** 翻页 · **F** 全屏 · **S** 演讲者备注 · **P** 打印/导出 PDF。 |
| **`slides.pdf`** | English deck as a ready PDF (16:9, one slide per page, no speaker notes). | For printing / emailing / projecting without a browser. |
| **`slides_zh.pdf`** | 中文幻灯片 PDF（16:9，每页一张，不含备注）。 | 便于打印 / 发送 / 投影。 |
| **`how_it_works.html`** | **Plain-language walkthrough for a non-technical audience** — the ten project steps in order, with five hand-drawn diagrams, analogies, a two-minute summary and a glossary. Assumes no coding or ML knowledge. | Open in any browser. Written to be read aloud or handed to a colleague who needs to understand *what we did* without the maths. |
| **`how_it_works_zh.html`** | **中文通俗版讲解** — 与英文版结构完全一致（十个步骤、五张图解、两分钟速讲、术语表），面向没有编程与机器学习基础的读者。 | 浏览器打开。可直接读给同事听，或发给需要了解「我们做了什么」但不需要数学细节的同事。 |
| **`SPEAKER_GUIDE.md`** | The full narrative — every slide's *"what's on screen"* + *"what to say"*, all the numbers, the code behind each step, positives/negatives, and a **Q&A defense** section. | Read once end-to-end to own the story, then present from the deck with this open beside you. (English; the Chinese decks carry per-slide speaker notes inline — press **S**.) |

> **Translation quality.** The Chinese deck was checked by two independent review passes — a
> terminology/number-fidelity audit and a back-translation/fluency audit — and their corrections
> applied (e.g. isotonic → 保序, causal genes → 作为耐药根源的). Gene names (`gyrA`, `mecA`,
> `blaKPC`…) and standard acronyms (ROC-AUC, VME, MIC, SHAP) are intentionally kept in the original.

## Corrections applied (2026-08-22, second pass — IMPORTANT)
A documentation-integrity audit found a **real bug in the shipped cefoxitin baseline**, now fixed.
The *K. pneumoniae* gene-lookup baseline had been using `mecA`/`mecC` — the ***S. aureus*** methicillin
genes, because cefoxitin is the standard MRSA surrogate test. Those genes occur in **0 of 3,850**
*Klebsiella* genomes, so the baseline never fired and scored a degenerate **0.500 by predicting a
constant**. Corrected to plasmid AmpC (`blaCMY/DHA/ACT/FOX/MOX`), the honest baseline is **0.518**,
missing **96.3%** of resistant strains (1,469 of 1,525).
**The headline claim is unchanged and now better supported** — the gene lookup genuinely fails on
cefoxitin, rather than failing because it was given the wrong organism's genes. Every deck, PDF, the
poster and both walkthroughs now quote 0.518 / 96.3%. **Use this version.**

## Corrections applied (2026-08-22, first pass)
A numbers audit against `results/metrics/*.json` found and fixed two stale figures that had been
carried over from an earlier, smaller run. **If you downloaded an earlier copy of this pack, use this one.**
- The cefoxitin rules baseline misses **96.6%** of resistant strains (368 of 381), not "~92%". The
  corrected number makes the headline *stronger*, not weaker.
- The VME claim is now stated precisely: every drug's VME **point estimate** is ≤3%, but the
  lineage-clustered 95% CIs are wide (meropenem 1.3% [0.3–3.5%] crosses the bar at its upper bound),
  so the deck now says the target is met **in expectation, not guaranteed**. Expect a judge to probe
  this — the honest answer is now on the slide.

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
