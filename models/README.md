# Model weights

**No model weights are bundled in this repository** (they are git-ignored). Place or generate them locally as below.

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
`fold_0/model_0/model.pt` (and `args.json`).

```bash
bash scripts/train_chemprop_baseline.sh       # writes models/chemprop_moonshot/
```

Reference checkpoint SHA-256 (`fold_0/model_0/model.pt`):

```
b0efb884d863ca05d0ad449658c5556081078744beb0692e6cefc9ecaab5b18d
```

## What needs which weights

| Command | Needs |
|---|---|
| `scripts/run_demo.sh` | `rel_siamese_s*.pt` (+ `data/moonshot/all.csv`); `--no-score` runs generation only |
| `scripts/run_evaluation.sh` | `rel_siamese_s*.pt`, `chemprop_moonshot/`, `data/moonshot/all.csv` |
| `scripts/run_validation_cv.sh` (optional, ~20 min) | `data/moonshot/all.csv`, ChemProp CLI |
| `python -m unittest discover -s tests` | none (weight/data tests are skipped when absent) |
