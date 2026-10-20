# Local evaluator demo

This page loads the **three fold-3 CliffForge checkpoints** and **three fold-3 ChemProp 1.6.1 checkpoints** from the Dablander MMP run. It performs live CPU inference for a manually entered candidate SMILES, reference SMILES, and measured reference IC₅₀ or pIC₅₀. The five processing stages show the structures, actual atom/bond graph, model routes, each seed's output, and the final comparison. It does not train models.

## Run

From the repository root on the machine with the model environment installed:

```bash
./demo/run_local.sh
```

Open [http://127.0.0.1:8765](http://127.0.0.1:8765). Stop the server with Ctrl-C. The service binds only to `127.0.0.1`. Set `CLIFFFORGE_DEMO_PORT` if 8765 is occupied. Set `CLIFFFORGE_RESULTS_DIR` to another complete Dablander run folder if needed.
Use `CLIFFFORGE_PYTHON=/path/to/venv/bin/python` to select a Python environment. The repository root `requirements.txt` documents an older training environment; use the verified versions below for this demo.

The launcher checks ChemProp 1.6.1, NumPy 1.x, hashes of the model source files recorded in `protocol.json`, and hashes of all six weights against `sha256_manifest.csv`. In a Codex worktree, it can copy the six ignored weights from the local primary checkout at `/Volumes/Work/TELEMEDICON26/results/dablander_mmp_5fold`; it does not overwrite existing weights. The model environment used for verification was Python 3.11.0, Torch 2.10.0, PyG 2.6.1, RDKit 2024.03.5, NumPy 1.26.4, ChemProp 1.6.1, FastAPI 0.116.1, and Uvicorn 0.35.0. The existing environment also needs the dependencies of the repository's `src/cliffforge` package. Run `./demo/run_local.sh --verify-only` to check files and dependencies without opening the server.

## Built-in example

Click **Load detected cliff example**, then **Predict with both models**. This is source MMP 10765 from the fold-3 `test_train` evaluation: candidate molecule 1655 was held out and reference molecule 1614 was in training. The reference's measured IC₅₀ is 85.8740204655224 µM (pIC₅₀ 4.0661382). The candidate's measured value is hidden until **Reveal measured outcome**. Its measured candidate-minus-reference ΔpIC₅₀ is +2.6073021, a Dablander activity cliff. CliffForge's archived three-seed mean is +2.0093143 and crosses the run's prediction cutoff of `|ΔpIC₅₀| > 1.5`.

For ChemProp, the page predicts **only the candidate's** pIC₅₀, then subtracts the **measured** reference pIC₅₀. CliffForge predicts the difference directly from both molecular graphs. A manually entered pair is labelled as having unverified MMP status. The selected example illustrates model behavior; it is not evidence of overall or statistically significant superiority.

## Regression check

```bash
python3 demo/test_demo.py
```

This compares live per-seed predictions with the archived fold-3 CSV, checks measured-reference arithmetic and graph drawing counts, and verifies invalid SMILES are rejected. All numbers are CPU predictions from loaded checkpoints; differences below `1e-4` pIC₅₀ from the Colab archive are accepted for runtime precision.
