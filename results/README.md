# Results

| Location | Experiment |
|---|---|
| `summary.json`, `metrics/`, `figures/` | Original COVID Moonshot Mpro reference run: 1,102 molecules, one 961/141 train/test split. `summary.json` is at the root of `results/`; the figures and detailed metrics are in their respective subfolders. |
| `dablander_mmp_5fold/` | Separate five-fold evaluation on Dablander's cleaned Mpro molecule and supplied MMP tables. Its summaries, fold outputs, predictions, logs, and checkpoints are kept together. |

## Original Moonshot Mpro reference results

The root-level `summary.json`, `metrics/`, and `figures/` are the **frozen reference results** of the original
CliffForge proof-of-concept run. They were produced by the code in `src/cliffforge/`; re-running the scripts writes
new files to `outputs/` (git-ignored). `scripts/validate_submission.py` checks selected reference claims in the
project README and slide deck against the frozen files.

| File | Content |
|---|---|
| `summary.json` | machine-readable summary: dataset, models, split, seeds, thresholds, metrics, counts, warnings |
| `metrics/activity_cliff_detection_metrics.json` | direct cliff-detection evaluation (settings A/B, 3 models, uncertainty, decision rule) |
| `metrics/activity_cliff_metrics.json` | descriptive cliff counts, pair-level direction accuracy, within-series Spearman (ChemProp vs relative) |
| `metrics/generation_recall.json` | retrospective generate-and-score check on held-out measured analogue pairs |
| `metrics/validation_thresholds.json` | F1-optimal cliff thresholds chosen on out-of-fold train-partition validation pairs |
| `metrics/funnel.json` and `metrics/{ranked_candidates.csv,generated_candidates.csv,shortlist.md,seed_selection_table.csv,retrospective_heldout_check.csv}` | end-to-end demo outputs for the reference seed (324 edits -> 231 unique -> 34 pass filters -> 34 scored) |
| `figures/*.png` | slide/README figures (workflow, funnel, within-series ranking, cliff counts, confusion matrices, precision-recall) |

Per-pair prediction tables (which contain many Moonshot structures with measured activities) are intentionally **not**
committed; regenerate them with `scripts/run_evaluation.sh`.

All values are retrospective. Predicted quantities (model-predicted delta-pIC50, scores) are never measured affinities, and
no generated molecule has been synthesised, assayed or docked.

## Dablander MMP five-fold evaluation

`dablander_mmp_5fold/` is the completed Colab run on the **1,924 source molecules** and **12,594 supplied MMP pairs**.
It contains five molecule-level outer folds. Each `fold_*/` directory holds the split IDs, evaluation pairs, predictions,
training logs, and three CliffForge plus three ChemProp 1.6.1 checkpoints. For test–train pairs, the known outer-training
reference activity is used for ChemProp and 5-NN. The notebook keeps Dablander's supplied MMP rows without the older
Moonshot similarity or pair filters.

Start with `dablander_mmp_5fold/RESULTS_README.md` for the headline tables and interpretation, then
`dablander_mmp_5fold/performance_summary.json` for exact values and planned confidence intervals. The
`dablander_mmp_5fold/protocol.json` file records the dataset, fold, model, and metric settings; the
`CliffForge_Dablander_5Fold_Colab.ipynb` notebook contains the evaluation code.
`dablander_mmp_5fold/bootstrap_paired_differences.csv` preserves the draws used for those intervals and for the
exploratory unadjusted 90% cliff-detection interval shown on slide 6. The run's `sha256_manifest.csv` records the
content hashes of the original Colab export; paths in that manifest are relative to `dablander_mmp_5fold/`.

Raw source tables, molecule activities, pair-level predictions, logs, and weights remain local and git-ignored.
The earlier 1,102-molecule proof-of-concept reference results remain in `summary.json`, `metrics/`, and `figures/`.
