# Data setup

This project has two separate SARS-CoV-2 Mpro data workflows. The original Moonshot reference run uses
`scripts/prepare_data.sh` and stores its local inputs and derived tables under `data/`. The newer Dablander
five-fold run uses `CliffForge_Dablander_5Fold_Colab.ipynb`, which downloads its own source tables in Colab.
The raw Moonshot table and the Dablander source tables are not committed to this repository. The local
`data/raw/` and `data/moonshot/` tables, and the large files in the Dablander result bundle, are git-ignored.

## Original Moonshot reference run

The following input is required for the scripts in `scripts/` and `src/cliffforge/data_prep.py`.

| Item | Value |
|---|---|
| Expected path | `data/raw/cdd_noncovalent_dates_2023_10_18_filt.csv` |
| Source | `COVID_moonshot_data_clean`: <https://github.com/asapdiscovery/COVID_moonshot_data_clean> (cleaned export of <https://github.com/asapdiscovery/COVID_moonshot_submissions>, downloaded from the COVID Moonshot CDD Vault) |
| Upstream license | MIT (repository), open-science data of the COVID Moonshot consortium |
| Assay / citation | Boby et al., *Science* 382, eabo7201 (2023), doi:10.1126/science.abo7201 |
| Preparation columns | `suspected_SMILES`, `IC50_(µM)`; `Canonical_PostEra_ID` is used only by the optional raw-record lookup in the demo |
| Raw size | 2,146 data rows in the referenced upstream CSV |

```bash
mkdir -p data/raw
curl -L -o data/raw/cdd_noncovalent_dates_2023_10_18_filt.csv \
  https://raw.githubusercontent.com/asapdiscovery/COVID_moonshot_data_clean/main/cdd_noncovalent_dates_2023_10_18_filt.csv
```

(If the upstream file layout changes, download the repository and copy the file named above into `data/raw/`.)

## Derived tables (generated, not downloaded)

```bash
bash scripts/prepare_data.sh          # == python -m cliffforge.data_prep
```

writes `data/moonshot/all.csv`, `train.csv`, `test.csv`, and `moonshot_prep.json`:

* IC50 (µM) -> pIC50 = 6 - log10(IC50); non-numeric or non-positive IC50 values and invalid SMILES are dropped. The preparation code does **not** inspect `Semiquant` or confidence-interval columns, or explicitly reject non-finite numeric values, so a numeric semi-quantitative measurement can remain;
* structures are canonicalised **without stereochemistry** and replicate measurements are averaged
  -> **1,102 unique 2D molecules**;
* Bemis-Murcko scaffolds define series; the 14 series with >= 8 members are split within-series
  (35 % test, `SEED = 17`), everything else stays in train -> **961 train / 141 test**.

`all.csv` columns: `smiles, pIC50, series, series_size, split, n_rep, rep_sd`.
`train.csv` and `test.csv` contain only `smiles,pIC50`. The split uses a fixed random seed; exact reproduction
also depends on using the same source file and software versions.

## Dablander five-fold run

The [five-fold Colab notebook](../CliffForge_Dablander_5Fold_Colab.ipynb) downloads
`molecule_data_clean.csv` and `MMP_data_clean.csv` from the Mpro directory of the
[Dablander source repository](https://github.com/MarkusFerdinandDablander/QSAR-activity-cliff-experiments),
pinned to commit `3081a1a40e3f5105a60b6d8b41f1ce69f8a0f453`. These files are separate from the Moonshot
CSV above: the notebook does **not** use `scripts/prepare_data.sh`, the 35% series split, or the Moonshot
pair and cliff filters.

The notebook expects **1,924 molecules** and **12,594 supplied MMP rows**. It reads `SMILES` and
`f_avg_IC50 [uM]` from the molecule table and `smiles_1`, `smiles_2` from the MMP table. It requires
positive numeric IC50 values and exact matches between MMP endpoints and molecule SMILES, then computes
pIC50 as above. It preserves source molecule and MMP identities, including stereochemistry; it does not
canonicalize, deduplicate, or generate replacement pairs. Measured `|ΔpIC50| ≥ 2` is a cliff,
`|ΔpIC50| ≤ 1` is a non-cliff, and intermediate pairs are excluded only from binary cliff detection.
Five shuffled **molecule-level** outer folds (seed 17017) supply the evaluation splits.

The completed export is under [`results/dablander_mmp_5fold/`](../results/dablander_mmp_5fold/).
Its `raw/` directory contains the downloaded source files locally; `molecules.csv`,
`source_mmp_pairs.csv`, and `molecule_folds.csv` record the derived data and split. These large files
are git-ignored. The committed `raw_sha256.json` records source-file hashes, and `protocol.json`
records the run settings. The `results/README.md` explains the bundle's remaining files.

## What is *not* needed

The FEP+ public benchmark, TNKS2 and VEGFR2 data from the earlier internship research are **not** used by this
repository and are not included. The FEP+ numbers quoted in the README are previously reported results of that
research (see the README, "References and prior work").
