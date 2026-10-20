# Presentation (TELEMEDICON26 · Team Genesis)

Eleven-slide, 16:9 **CliffForge** pitch deck in LaTeX Beamer. It explains drug-discovery terms for a CS/AI/ML audience, separates the proposed protein-aware research pipeline from the implemented 2D prototype, and compares CliffForge with ChemProp on two Mpro datasets. Diagrams and result displays are editable TikZ vector graphics in the existing navy-and-teal theme.

| File | Purpose |
|---|---|
| main.tex | Deck source, slides 1–11 |
| references.bib | Bibliography cited by the deck |
| TELEMEDICON26_CliffForge_Pitch.pdf | Compiled presentation for viewing and projection |
| main.pdf | Identical compiled copy beside the LaTeX source |
| speaker_notes_10min.md | Ten-minute, slide-by-slide speaking script |
| build.sh | Reproducible PDF build |

## Slide narrative

1. Cover: relative-affinity analogue prioritization and activity-cliff analysis.
2. Short vocabulary: ligand, pocket, affinity, relative affinity, congeneric series, activity cliff.
3. Linked pocket-interaction diagram: a small edit can change contacts and measured activity.
4. Absolute-activity baselines: one-molecule prediction followed by subtraction for a pair.
5. Research focus: rank analogues within a series and pay attention to activity cliffs.
6. **Proposed** protein-aware generation, pose/scoring, shared pair encoder, and held-out evaluation.
7. **Implemented** RDKit matched-pair generation and 2D relative-activity model.
8. The COVID Moonshot and Dablander Mpro datasets, including their different cliff rules.
9. One actual prototype run from seed through filtering and model-ranked shortlist.
10. Retrospective ChemProp versus CliffForge results for cliff ranking, local ranking, and lead selection.
11. Impact, users, sustainability, proposed business model, and six-stage roadmap.

## Build

Requires pdflatex, latexmk, biber, and the beamer, biblatex, and pgf/tikz packages (TeX Live or TinyTeX).

    bash presentation/build.sh

The output is presentation/TELEMEDICON26_CliffForge_Pitch.pdf, with an identical presentation/main.pdf; build intermediates are written to a git-ignored directory. The deck should have **11 pages**.

## Evidence used on the slides

**COVID Moonshot noncovalent Mpro:** results/moonshot_mpro_tanimoto/ contains the fixed reference run. The standardized set has 1,102 molecules (961 train, 141 test). A cliff is a pair with Morgan-fingerprint Tanimoto similarity at least 0.80 and measured |ΔpIC50| ≥ 1. The full set has 46 such pairs; 44 close test–test pairs contain five strict held-out cliffs. The deck reports cliff average precision 0.731 (CliffForge) versus 0.349 (ChemProp), two versus zero of the five true cliffs flagged at the fixed threshold, and mean within-series Spearman 0.581 versus 0.192 across 12 held-out series. Source files: metrics/activity_cliff_detection_metrics.json and metrics/activity_cliff_metrics.json.

The slide 9 demo is a **model output**, not an assay result: 324 attempted edits, 231 unique non-parent candidates, and 34 passing and scored candidates (metrics/funnel.json, metrics/shortlist.md). The retrospective generation check recovered 82 of 316 reachable test–test single-edit pairs. Direction was correct for 53/82 with CliffForge and 47/82 with ChemProp; the paired difference was not statistically significant. The recovered set is selected by training-observed transformations (metrics/generation_recall.json).

**Dablander cleaned Mpro matched pairs:** results/dablander_mmp_5fold/ contains 1,924 molecules and 12,594 supplied matched molecular pairs, with five molecule-level folds and three seeds per neural model. A binary cliff has measured |ΔpIC50| ≥ 2; a noncliff has |Δ| ≤ 1. Among 2,132 binary test–test pairs, 106 are cliffs; 343 intermediate pairs are excluded only from binary detection. Cliff average precision is 0.391 versus 0.137; at the fixed predicted-difference threshold, 20/106 versus 0/106 true cliffs are flagged (CliffForge also has 16 false positives). Mean local-reference Spearman is 0.509 versus 0.343. Selecting the top 10% of proposals gives mean measured gains of 1.116 versus 1.002 pIC50 units. These numbers and the planned six-comparison adjusted inference are in performance_summary.json and protocol.json.

The adjusted Dablander comparison supports CliffForge’s local ranking advantage over ChemProp; its cliff average-precision and top-10% gain differences versus ChemProp remain inconclusive. Moonshot’s five strict held-out cliffs are preliminary evidence only. The two datasets have different cliff definitions and evaluation protocols, so their rates and counts are not pooled.

## Interpretation and limits

- The full protein-aware pipeline on slide 6 is the **research plan** from CliffForge_Teammate_Brief.pdf, pp. 2–3. The prototype uses training-observed RDKit substituent swaps and a 2D Siamese D-MPNN. It does not run DiffSBDD, docking, protein-pocket encoding, or stereochemistry-aware evaluation.
- The measured label in both datasets is **pIC50 assay potency**, not a direct binding-affinity measurement. Docking differences in the proposed pipeline would be noisy proxy labels. Generated analogues are candidate hypotheses, not confirmed activity cliffs.
- The prototype’s ranking subtracts three-model ensemble standard deviation from the predicted relative activity. This is a practical uncertainty signal, not a calibrated probability.
- Slide 11’s enterprise licensing or funded co-development route is a **proposal**. [CHAI Discovery’s commercial access page](https://www.chaidiscovery.com/product) and [its news page listing the Pfizer license](https://www.chaidiscovery.com/news) provide the business-model analogy. No pricing, revenue share, customer demand, or partnership is claimed for CliffForge.

The PDF and editable LaTeX source are the deliverables. A PowerPoint conversion is not included because it would reduce fidelity and editability of the TikZ diagrams.
