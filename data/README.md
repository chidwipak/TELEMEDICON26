# Data setup

No datasets are redistributed in this repository. Everything under `data/` except this file and `.gitkeep` placeholders
is git-ignored. CliffForge needs one public dataset.

## Required input: COVID Moonshot non-covalent Mpro IC50 table

| Item | Value |
|---|---|
| Expected path | `data/raw/cdd_noncovalent_dates_2023_10_18_filt.csv` |
| Source | `COVID_moonshot_data_clean`: <https://github.com/asapdiscovery/COVID_moonshot_data_clean> (cleaned export of <https://github.com/asapdiscovery/COVID_moonshot_submissions>, downloaded from the COVID Moonshot CDD Vault) |
| Upstream license | MIT (repository), open-science data of the COVID Moonshot consortium |
| Assay / citation | Boby et al., *Science* 382, eabo7201 (2023), doi:10.1126/science.abo7201 |
| Used columns | `suspected_SMILES`, `Canonical_PostEra_ID`, `IC50_(µM)` |
| Raw size | 2,146 rows (1,889 PostEra IDs); 1,450 rows have a numeric IC50 > 0 |

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

writes `data/moonshot/all.csv`, `train.csv`, `test.csv`, `moonshot_prep.json`:

* IC50 (µM) -> pIC50 = 6 - log10(IC50); censored / non-numeric / non-positive IC50 rows are dropped;
* structures are canonicalised **without stereochemistry** and replicate measurements are averaged
  -> **1,102 unique 2D molecules**;
* Bemis-Murcko scaffolds define series; the 14 series with >= 8 members are split within-series
  (35 % test, `SEED = 17`), everything else stays in train -> **961 train / 141 test**.

`all.csv` columns: `smiles, pIC50, series, series_size, split, n_rep, rep_sd`.
The split is deterministic: regenerating from the raw file reproduces the reference tables byte-for-byte.

## What is *not* needed

The FEP+ public benchmark, TNKS2 and VEGFR2 data from the earlier internship research are **not** used by this
repository and are not included. The FEP+ numbers quoted in the README are previously reported results of that
research (see the README, "References and prior work").
