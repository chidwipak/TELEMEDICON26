# Presentation (TELEMEDICON26 · Team Genesis)

Six-slide 16:9 pitch deck for **CliffForge**, written in LaTeX Beamer. All diagrams and charts are native TikZ vector graphics
(no raster screenshots), so they stay crisp at any projection size.

| File | Purpose |
|---|---|
| `main.tex` | deck source (6 slides: cover, problem & validation, solution & innovation, methodology, prototype/proof of concept, impact & next steps) |
| `references.bib` | verified bibliography (entries copied from the project's research bibliography; only cited works are kept) |
| `TELEMEDICON26_CliffForge_Pitch.pdf` | compiled deck (committed so it can be viewed without LaTeX) |
| `build.sh` | reproducible build, works from any directory |

## Build

Requires `pdflatex`, `latexmk`, `biber` and the `beamer`, `biblatex`, `pgf/tikz` packages (TeX Live, or TinyTeX with
`tlmgr install beamer biblatex biber`).

```bash
bash presentation/build.sh        # -> presentation/TELEMEDICON26_CliffForge_Pitch.pdf (+ build/ with intermediates, git-ignored)
```

The script finds TinyTeX's `biber` automatically, runs `latexmk` (three passes incl. biber) and prints the page count (6).

## Where the numbers come from

Every figure on slides 2 and 5 is taken from `results/summary.json` and `results/metrics/*.json` (the frozen reference run):

* within-series Spearman 0.192 (ChemProp) vs 0.581 (relative model), 12 series, paired Wilcoxon p = 0.005 - `results/metrics/activity_cliff_metrics.json`
* 46 cliff pairs (27 train-train, 14 test-train, 5 strict held-out), 44 close held-out pairs - `results/metrics/activity_cliff_detection_metrics.json`
* 2 of 5 vs 0 of 5 cliffs flagged, average precision 0.73 vs 0.35 - same file, setting A
* 324 -> 231 -> 34 -> 34 funnel - `results/metrics/funnel.json`; 82 of 316 recovered pairs - `results/metrics/generation_recall.json`
* FEP+ 0.386 -> 0.644 is labelled as **prior internship research**, not a new result.

`python scripts/validate_submission.py` checks that `main.tex` and the README contain exactly these values and that the PDF has six pages.

## Notes

* The slides present the cliff-detection result as a *preliminary* signal on five held-out positives; the lead-optimization
  setting (14 cliffs) is mentioned as mixed. See the README and `results/metrics/` for details.
* A PowerPoint version is not provided: converting Beamer/TikZ to PPTX loses layout fidelity and editability, so the verified
  PDF and editable LaTeX source are the deliverables.
