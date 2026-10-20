# CliffForge Dablander five-fold results

This bundle contains the exact source CSVs and hashes, molecule folds, outer train/test IDs and ChemProp internal split CSVs,
training logs, all CliffForge and ChemProp checkpoints, per-molecule and per-pair predictions,
paired bootstrap draws, and performance_summary.json.

## Primary outcomes

| Metric | CliffForge | ChemProp | 5-NN |
|---|---:|---:|---:|
| strict_AP | 0.3905 | 0.1369 | 0.0923 |
| within_series_rho | 0.5092 | 0.3429 | 0.2866 |
| lead_top10_gain | 1.1165 | 1.0021 | 0.8238 |

## Overall relative affinity (held-out MMPs)

Numerical errors use signed ΔpIC50. Test–train pairs use the measured outer-training reference for ChemProp and 5-NN.

| Setting | Metric | CliffForge | ChemProp | 5-NN |
|---|---|---:|---:|---:|
| all_heldout_pairs | mae_delta_pIC50 | 0.3442 | 0.4159 | 0.4680 |
| all_heldout_pairs | rmse_delta_pIC50 | 0.5275 | 0.5592 | 0.6423 |
| all_heldout_pairs | pearson_r | 0.7533 | 0.7158 | 0.6261 |
| all_heldout_pairs | spearman_rho | 0.7211 | 0.6637 | 0.6005 |
| all_heldout_pairs | direction_accuracy_nonzero | 0.7984 | 0.7521 | 0.7441 |
| test_test | mae_delta_pIC50 | 0.4587 | 0.5216 | 0.5693 |
| test_test | rmse_delta_pIC50 | 0.6578 | 0.7261 | 0.7949 |
| test_test | pearson_r | 0.5753 | 0.4022 | 0.2673 |
| test_test | spearman_rho | 0.5354 | 0.3613 | 0.2979 |
| test_test | direction_accuracy_nonzero | 0.7236 | 0.6398 | 0.6265 |
| test_train | mae_delta_pIC50 | 0.3302 | 0.4030 | 0.4556 |
| test_train | rmse_delta_pIC50 | 0.5093 | 0.5353 | 0.6210 |
| test_train | pearson_r | 0.7724 | 0.7448 | 0.6581 |
| test_train | spearman_rho | 0.7413 | 0.6922 | 0.6288 |
| test_train | direction_accuracy_nonzero | 0.8076 | 0.7659 | 0.7586 |
| within_series_reference_pairs | mae_delta_pIC50 | 0.3105 | 0.3950 | 0.4496 |
| within_series_reference_pairs | rmse_delta_pIC50 | 0.4703 | 0.5230 | 0.6137 |
| within_series_reference_pairs | pearson_r | 0.8073 | 0.7542 | 0.6645 |
| within_series_reference_pairs | spearman_rho | 0.7743 | 0.6999 | 0.6336 |
| within_series_reference_pairs | direction_accuracy_nonzero | 0.8241 | 0.7712 | 0.7626 |

## Lead optimization

The top 10% of scored candidate/reference proposals are selected separately in each fold.

| Metric | CliffForge | ChemProp | 5-NN |
|---|---:|---:|---:|
| n_selected_proposals | 2026 | 2026 | 2026 |
| fraction_actual_gain_gt_0 | 0.9284 | 0.8159 | 0.6900 |
| fraction_actual_gain_gte_0p5 | 0.8391 | 0.7029 | 0.6209 |

## Cliff detection

Binary metrics exclude halfAC pairs; thresholded metrics use |predicted ΔpIC50| > 1.5.

| Setting | Metric | CliffForge | ChemProp | 5-NN |
|---|---|---:|---:|---:|
| test_test | average_precision | 0.3905 | 0.1369 | 0.0923 |
| test_test | roc_auc | 0.7999 | 0.7341 | 0.6138 |
| test_test | accuracy_at_1p5 | 0.9522 | 0.9503 | 0.9479 |
| test_test | mcc_at_1p5 | 0.3050 | 0.0000 | -0.0111 |
| test_test | f1_at_1p5 | 0.2817 | 0.0000 | 0.0000 |
| test_test | precision_at_1p5 | 0.5556 | 0.0000 | 0.0000 |
| test_test | recall_at_1p5 | 0.1887 | 0.0000 | 0.0000 |
| test_test | tp | 20 | 0 | 0 |
| test_test | fp | 16 | 0 | 5 |
| test_test | tn | 2010 | 2026 | 2021 |
| test_test | fn | 86 | 106 | 106 |
| test_train | average_precision | 0.6577 | 0.5950 | 0.4011 |
| test_train | roc_auc | 0.8705 | 0.9226 | 0.8234 |
| test_train | accuracy_at_1p5 | 0.9709 | 0.9676 | 0.9570 |
| test_train | mcc_at_1p5 | 0.6208 | 0.5671 | 0.4368 |
| test_train | f1_at_1p5 | 0.6041 | 0.5463 | 0.4451 |
| test_train | precision_at_1p5 | 0.8616 | 0.8248 | 0.5792 |
| test_train | recall_at_1p5 | 0.4651 | 0.4084 | 0.3614 |
| test_train | tp | 386 | 339 | 300 |
| test_train | fp | 62 | 72 | 218 |
| test_train | tn | 16508 | 16498 | 16352 |
| test_train | fn | 444 | 491 | 530 |

## Planned superiority contrasts

Each interval is a two-sided Bonferroni-adjusted paired cluster-bootstrap interval
across the six comparisons. Superiority is supported only if the lower bound is positive.

| Contrast | Difference | Adjusted interval | Supported |
|---|---:|---:|---|
| strict_AP__vs__chemprop | 0.25362349659292555 | [-0.0984, 0.5638] | False |
| strict_AP__vs__knn_5nn | 0.2981667477987395 | [-0.0310, 0.6017] | False |
| within_series_rho__vs__chemprop | 0.16632613405979807 | [0.1342, 0.1975] | True |
| within_series_rho__vs__knn_5nn | 0.2225714549134103 | [0.1898, 0.2539] | True |
| lead_top10_gain__vs__chemprop | 0.1143793379536373 | [-0.0587, 0.2655] | False |
| lead_top10_gain__vs__knn_5nn | 0.29271411288258276 | [0.0742, 0.5084] | True |

**All six supported:** False

Strict test-test and known-reference test-train MMPs are separate tasks. The source paper used
three repeats of two-fold CV; this notebook uses one prespecified five-fold molecule CV.
Pairs can be correlated and a source MMP can appear twice as a known-reference scenario
in opposite query/reference orientations. The archive preserves pair, molecule and fold IDs.
