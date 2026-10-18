# CliffForge PoC shortlist — seed `moonshot_idx237`

**PROTOTYPE OUTPUT. Candidates are model-prioritized, not validated. Scores are 2D-model predicted ΔpIC50 vs the seed (not measured potency, not docking, not binding free energy).**

- Seed SMILES (2D): `CNC(=O)CN1CC(C(=O)Nc2cncc3ccc(NS(C)(=O)=O)cc23)c2cc(Cl)ccc2C1=O`
- Seed partition: `test`; measured pIC50: 7.18; PostEra ID(s): MAT-POS-2e8b2191-11
- Target: SARS-CoV-2 Mpro (COVID Moonshot noncovalent set). Pocket-aware docking: **not run**.
- Funnel: 324 attempted edits → 324 valid → 231 unique non-parent → 34 pass all filters → 34 model-scored

## Top 10 by lower-confidence-bound (predicted Δ − ensemble SD)

| # | candidate SMILES | edit | Tan. | pred ΔpIC50 vs seed (±SD) | cliff alert | caveat |
|--:|---|---|--:|--:|---|---|
| 1 | `CNC(=O)CN1CC(C(=O)Nc2cncc3ccc(NS(C)(=O)=O)cc23)c2cc(NS(C)(=O)=O)ccc2C1=O` | `Cl[*:1]>>CS(=O)(=O)N[*:1]` | 0.90 | +0.18 ± 0.01 | - | ensemble SD 0.01; 2 train precedent pair(s); 2D model: no pocket / stereo / docking |
| 2 | `CNC(=O)C1(N2CC(C(=O)Nc3cncc4ccc(NS(C)(=O)=O)cc34)c3cc(Cl)ccc3C2=O)CC1` | `CNC(=O)C[*:1]>>CNC(=O)C1([*:1])CC1` | 0.72 | +0.15 ± 0.02 | - | ensemble SD 0.02; 2 train precedent pair(s); 2D model: no pocket / stereo / docking |
| 3 | `CNC(=O)CN1CC(C(=O)Nc2cncc3ccc(OC)cc23)c2cc(Cl)ccc2C1=O` | `CS(=O)(=O)N[*:1]>>CO[*:1]` | 0.75 | +0.08 ± 0.05 | - | ensemble SD 0.05; 2 train precedent pair(s); 2D model: no pocket / stereo / docking |
| 4 | `CNC(=O)CN1CC(C(=O)Nc2cncc3ccc(NS(=O)(=O)N4CC(C#N)C4)cc23)c2cc(Cl)ccc2C1=O` | `CS(=O)(=O)[*:1]>>N#CC1CN(S(=O)(=O)[*:1])C1` | 0.79 | +0.06 ± 0.04 | - | ensemble SD 0.04; 2 train precedent pair(s); 2D model: no pocket / stereo / docking |
| 5 | `CNC(=O)CN1CC(C(=O)Nc2cncc3ccc(NS(=O)(=O)CC4(C#N)CC4)cc23)c2cc(Cl)ccc2C1=O` | `CS(=O)(=O)[*:1]>>N#CC1(CS(=O)(=O)[*:1])CC1` | 0.77 | +0.02 ± 0.03 | - | ensemble SD 0.03; 2 train precedent pair(s); 2D model: no pocket / stereo / docking |
| 6 | `CS(=O)(=O)Nc1ccc2cncc(NC(=O)C3CN(CC(N)=O)C(=O)c4ccc(Cl)cc43)c2c1` | `CN[*:1]>>N[*:1]` | 0.85 | -0.02 ± 0.01 | - | ensemble SD 0.01; 2 train precedent pair(s); 2D model: no pocket / stereo / docking |
| 7 | `CCNC(=O)CN1CC(C(=O)Nc2cncc3ccc(NS(C)(=O)=O)cc23)c2cc(Cl)ccc2C1=O` | `C[*:1]>>CC[*:1]` | 0.86 | -0.03 ± 0.01 | - | ensemble SD 0.01; 2 train precedent pair(s); 2D model: no pocket / stereo / docking |
| 8 | `CNC(=O)CN1CC(C(=O)Nc2cncc3ccc(NS(C)(=O)=O)cc23)c2cc(F)ccc2C1=O` | `Cl[*:1]>>F[*:1]` | 0.84 | -0.07 ± 0.01 | - | ensemble SD 0.01; 10 train precedent pair(s); 2D model: no pocket / stereo / docking |
| 9 | `CS(=O)(=O)Nc1ccc2cncc(NC(=O)C3CN(CC(=O)NCC4CC4)C(=O)c4ccc(Cl)cc43)c2c1` | `C[*:1]>>C1CC1C[*:1]` | 0.83 | +0.00 ± 0.11 | - | ensemble SD 0.11; 4 train precedent pair(s); 2D model: no pocket / stereo / docking |
| 10 | `CC(=O)NC(=O)CN1CC(C(=O)Nc2cncc3ccc(NS(C)(=O)=O)cc23)c2cc(Cl)ccc2C1=O` | `C[*:1]>>CC(=O)[*:1]` | 0.86 | -0.08 ± 0.04 | - | ensemble SD 0.04; 2 train precedent pair(s); 2D model: no pocket / stereo / docking |

## Retrospective check on measured held-out analogues

1 generated candidates (0 of them pass the filters) are real Moonshot compounds from the held-out test partition (never used to mine transformations or to train). Direction accuracy of predicted vs measured Δ(candidate − seed): 1.0. See `retrospective_heldout_check.csv` (tiny n; illustrative only).
