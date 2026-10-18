# Examples

`seed_ligand.smi` - the 2D seed used for the reference demo: a held-out COVID Moonshot Mpro compound
(`all.csv` row 237, PostEra ID `MAT-POS-2e8b2191-11`, measured pIC50 7.18, test partition). It is selected by an
outcome-blind rule (held-out potent compounds with the most filter-passing candidates; see
`results/metrics/seed_selection_table.csv`).

Run the generation + scoring pipeline on it (needs `data/moonshot/all.csv` and the relative-model weights):

```bash
python -m cliffforge.demo --seed-smiles "$(cut -d' ' -f1 examples/seed_ligand.smi)" --seed-id moonshot_idx237
# or the rule-based auto-seed used for the reference results:
bash scripts/run_demo.sh
```

Any other seed SMILES can be supplied the same way. The generator only proposes edits seen in the training partition of
the Moonshot data, so seeds from Mpro-like chemotypes yield the most candidates.

Reference outputs for this seed are in `results/metrics/` (`ranked_candidates.csv`, `generated_candidates.csv`,
`shortlist.md`, `funnel.json`).
