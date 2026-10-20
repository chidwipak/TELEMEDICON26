# Model weights

**Model weights are not committed to Git**. The original Moonshot commands expect checkpoints under `models/`;
the completed Dablander run's checkpoints are present locally under `results/dablander_mmp_5fold/` and are
git-ignored. The two runs use different checkpoint directories.

## 1. Relative-affinity model: `models/rel_siamese_s{0,1,2}.pt`

Three seeds of a whole-molecule Siamese D-MPNN trained on pIC50 *differences* between same-series **training** compounds.

| Setting | Value |
|---|---|
| Architecture | shared D-MPNN encoder (hidden 128, depth 4, dropout 0.1), head on `[g_a, g_b, abs(g_a-g_b)]` (`cliffforge/networks.py`) |
| Node / edge features | 48 / 8 (whole-molecule 2D graph; `cliffforge/features.py`) |
| Training data | `data/moonshot/all.csv`, `split == train` only (4,578 same-series pairs, up to 600 per series) |
| Optimisation | Adam, lr 5e-4, weight decay 1e-5, batch 32, SmoothL1 loss, 70 epochs, gradient clip 5.0, seeds 0/1/2 |
| Checkpoint format | `torch.save({"state_dict": ..., "seed": ..., "epochs": ..., "n_train_pairs": ..., ...})` |

Create them (about 10 minutes per seed on CPU):

```bash
bash scripts/prepare_data.sh
bash scripts/train_relative_models.sh        # writes models/rel_siamese_s0.pt, _s1.pt, _s2.pt
```

Reference checkpoints from the research run (SHA-256), to check whether you hold the *same* weights:

```
a479ad899f1632b7fd1a1836ee9be6994939972ef0af422afbbd003d262d7799  rel_siamese_s0.pt
cc75f85d156e743580cd5058a068ca457d81e672500e1239596076a316c4ec24  rel_siamese_s1.pt
4438b08d5129e9d5c85a97cac49bdd261a14e75954303167f6c40d6977dd8a90  rel_siamese_s2.pt
```

Retraining with the same seeds on a different CPU/BLAS or library build can change low-order digits; expect the same
qualitative results but not necessarily identical numbers. The reference numbers in `results/` were produced with the
checkpoints above.

## 2. ChemProp absolute-affinity baseline: `models/chemprop_moonshot/`

ChemProp 1.6.1 regression model (D-MPNN, depth 3, hidden 300, 60 epochs, defaults otherwise, ChemProp's internal random
0.8/0.1/0.1 split) trained on the Moonshot **train** pIC50 table. The directory must contain
`fold_0/model_0/model.pt`; ChemProp also writes `args.json` at the directory root.

```bash
bash scripts/train_chemprop_baseline.sh       # writes models/chemprop_moonshot/
```

Reference checkpoint SHA-256 (`fold_0/model_0/model.pt`):

```
b0efb884d863ca05d0ad449658c5556081078744beb0692e6cefc9ecaab5b18d
```

## 3. Dablander MMP five-fold evaluation models

The completed run is kept together under `results/dablander_mmp_5fold/`. For each outer fold `0` through `4`, its three
CliffForge checkpoints are `fold_<fold>/cliffforge/seed_<seed>/final.pt`; its three ChemProp 1.6.1 checkpoints are
`fold_<fold>/chemprop/seed_<seed>/fold_0/model_0/model.pt`, with seeds `0`, `1`, and `2`. This is **15 checkpoints per
neural model**, paired with the run's split IDs, training logs, predictions, and checksum manifest. The checkpoints
remain local and git-ignored. Use the checkpoint from the matching outer fold when reproducing a held-out prediction;
the five-fold aggregate in `results/dablander_mmp_5fold/performance_summary.json` is an evaluation result, not a
single deployable checkpoint.

The Dablander notebook trains CliffForge on supplied MMP rows whose two molecules are in the outer training pool.
It saves CliffForge's final checkpoint after 70 epochs (batch size 32). ChemProp trains on the same outer-training
molecules, makes its own 80/10/10 split within that pool, and saves the best internal-validation checkpoint from
up to 60 epochs. The 5-NN baseline has no trained checkpoint. The Moonshot training scripts above do not recreate
these Dablander fold checkpoints; use `CliffForge_Dablander_5Fold_Colab.ipynb` for that run.

## What needs which weights

| Command | Needs |
|---|---|
| `scripts/run_demo.sh` | `rel_siamese_s*.pt` (+ `data/moonshot/all.csv`); `--no-score` runs generation only |
| `scripts/run_evaluation.sh` | `rel_siamese_s*.pt`, `chemprop_moonshot/`, `data/moonshot/all.csv` |
| `scripts/run_validation_cv.sh` (optional, ~20 min) | `data/moonshot/all.csv`, ChemProp CLI |
| `PYTHONPATH=src python -m unittest discover -s tests -p test_generation.py` | none |

The full `unittest discover -s tests` command currently fails because legacy `tests/test_mmp_pairs.py` imports a
removed MMP-pair helper. The generation tests above pass without data or weights; data/weight-dependent tests are
skipped when those files are absent.
