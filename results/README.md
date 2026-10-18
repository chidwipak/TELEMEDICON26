# Reference results

These files are the **frozen reference results** of the CliffForge proof-of-concept run (Mpro, split 961/141). They are
produced by the code in `src/cliffforge/`; re-running the scripts writes new results to `outputs/` (git-ignored) and never
overwrites this folder. `scripts/validate_submission.py` checks that the README and slide deck agree with these files.

| File | Content |
|---|---|
| `summary.json` | machine-readable summary: dataset, models, split, seeds, thresholds, metrics, counts, warnings |
| `metrics/activity_cliff_detection_metrics.json` | direct cliff-detection evaluation (settings A/B, 3 models, uncertainty, decision rule) |
| `metrics/activity_cliff_metrics.json` | descriptive cliff counts, pair-level direction accuracy, within-series Spearman (ChemProp vs relative) |
| `metrics/generation_recall.json` | retrospective generate-and-score check on held-out measured analogue pairs |
| `metrics/validation_thresholds.json` | F1-optimal cliff thresholds chosen on out-of-fold train-partition validation pairs |
| `metrics/funnel.json`, `ranked_candidates.csv`, `generated_candidates.csv`, `shortlist.md`, `seed_selection_table.csv`, `retrospective_heldout_check.csv` | end-to-end demo outputs for the reference seed (324 edits -> 231 unique -> 34 pass filters -> 34 scored) |
| `figures/*.png` | slide/README figures (workflow, funnel, within-series ranking, cliff counts, confusion matrices, precision-recall) |

Per-pair prediction tables (which contain many Moonshot structures with measured activities) are intentionally **not**
committed; regenerate them with `scripts/run_evaluation.sh`.

All values are retrospective. Predicted quantities (model-predicted delta-pIC50, scores) are never measured affinities, and
no generated molecule has been synthesised, assayed or docked.
