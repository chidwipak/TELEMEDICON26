# CliffForge: 10-minute speaker notes

**Run of show:** Eleven slides, approximately ten minutes. The wording is for a computer science and AI/ML audience; the diagrams on slides 3, 6, and 9 supply the physical and model context. There is no separate simulation cue in this version.

## Slide 1 — CliffForge (0:00–0:35)

“CliffForge asks a local drug-discovery question: given a promising molecule, which small edit should a chemist test next? Our research vision combines constrained analogue generation, relative-affinity scoring, and uncertainty-aware ranking. Today we have a working computational prototype and retrospective evaluations. The protein-aware stages and laboratory validation are the next research steps.”

## Slide 2 — Drug discovery vocabulary (0:35–1:20)

“A ligand is a small molecule that may interact with a target protein. The pocket is the three-dimensional region it occupies. Binding affinity describes how tightly it binds; a lower dissociation constant, or K-d, means tighter binding. We compare relative affinity when we ask whether ligand B binds better than ligand A to the same target. A congeneric series shares a chemical core and varies by small edits. An activity cliff is a pair of close analogues with a large measured activity difference. Our data report pIC50, an assay potency measure, so a pIC50 difference is not itself a direct binding-affinity measurement.”

## Slide 3 — A beneficial activity cliff (1:20–2:05)

“The two pocket sketches show the same protein pocket and the same ligand core. Ligand B differs by one local edit; its new group could make a favorable contact that ligand A lacks. That interaction is a hypothesis, while the illustrated assay outcome is that B inhibits the target more strongly. A gain of at least one pIC50 unit is at least tenfold higher assay potency: a beneficial activity cliff between close analogues. The sketches are conceptual, not measured poses. Better assay potency could support a therapeutic improvement, but no clinical response has been measured here.”

## Slide 4 — Absolute activity and affinity models (2:05–2:50)

“These are representative architectures that output one scalar for a candidate or drug–target pair. ChemProp uses a ligand graph and directed bond message passing; in our benchmark it is the tested, ligand-only comparator predicting Mpro pIC50. DeepDTA uses separate one-dimensional CNNs for a drug SMILES string and protein sequence. GraphDTA uses a ligand graph network with a protein sequence encoder. Pafnucy uses a voxelized three-dimensional protein–ligand complex and a 3D CNN. The latter three illustrate alternative published absolute predictors; we have not benchmarked them here. To compare close analogues using any single-candidate predictor, score A and B separately and subtract. Our measured pIC50 endpoint is assay potency, not a direct binding-affinity measurement.”

## Slide 5 — Research focus (2:50–3:30)

“Our focus is predicting binding-affinity within series, especially activity cliffs. In practical terms, a chemist has a lead and several nearby analogue ideas. Synthesis and assays cost time, so the useful output is which edit is most likely to improve the lead and which edit might produce an unexpectedly large gain or loss. We test that idea using held-out measured activity; generated structures and predicted scores remain hypotheses until experiments measure them.”

## Slide 6 — Proposed protein-aware pipeline (3:30–4:45)

“This is the full proposed A-to-G research pipeline, not the current implementation. We start with protein structure, bound ligand pose, and ligand chemistry. A local fragment or R-group is masked or replaced. DiffSBDD inpainting is the planned generation route, with FEgrow or RDKit enumeration as a practical fallback. We then check chemical validity and pose plausibility and score the parent and edited ligand with SMINA, GNINA, or Vina under consistent conditions. The difference in those scores is only a noisy training proxy, not measured affinity. A Siamese network applies the same ligand-plus-pocket encoder to A and B; a difference head estimates their relative affinity or rank, with uncertainty. The decisive evaluation uses held-out experimentally measured pairs, especially activity cliffs. Generated analogues remain hypotheses until tested.”

## Slide 7 — Current prototype (4:45–5:35)

“The five boxes are the implemented workflow in order. Input is a seed SMILES plus 961 measured Mpro training compounds. Generate applies RDKit single-substituent swaps observed at least twice in training. Filter removes invalid or duplicate structures and checks similarity, edit size, PAINS, and synthetic accessibility. Score uses shared-weight Siamese directed message passing on seed and candidate 2D graphs. A three-model ensemble predicts candidate-minus-seed delta pIC50 and its standard deviation. Rank sorts by predicted change minus that spread and flags close analogues with a large predicted change. The spread is an uncertainty signal, not a calibrated probability. This version does not use protein coordinates, DiffSBDD, docking, or stereochemistry.”

## Slide 8 — Two datasets (5:35–6:20)

“Both datasets concern SARS-CoV-2 Mpro activity, but they use different pair definitions. The COVID Moonshot set has 1,102 standardized molecules, with 961 in training and 141 held out. A cliff requires fingerprint similarity of at least 0.80 and a measured pIC50 difference of at least one. There are 46 such pairs overall, but only five among the 44 close test–test pairs. The Dablander cleaned matched-pair set has 1,924 molecules and 12,594 supplied pairs. Its cliff threshold is two pIC50 units, or at least a hundredfold potency difference; noncliffs differ by at most one. Five molecule folds give 2,132 binary test–test pairs, including 106 cliffs. We do not pool the cliff rates across these protocols.”

## Slide 9 — One prototype run (6:20–7:10)

“Here is the data path for one held-out Moonshot seed. A training-observed substituent swap proposes an analogue. The generator attempted 324 edits, found 231 unique non-parent candidates, and kept and scored 34 after filtering. The model compares each candidate with the seed and returns a predicted delta and ensemble spread. The top ranked row is plus 0.18 ± 0.01 pIC50, with no cliff alert. The drawn chemical motif is schematic; the counts and shortlist values come from the actual run. None of these suggested structures was synthesized, assayed, or docked.”

## Slide 10 — Retrospective results (7:10–9:05)

“The table separates three decisions: ranking true cliffs, ordering nearby analogues, and selecting promising edits. For Moonshot, CliffForge’s cliff average precision is 0.73 versus ChemProp’s 0.35. At the fixed threshold, it flags two of five strict held-out cliffs versus none for ChemProp. Five positives make that preliminary. Across 12 held-out series, mean Spearman ranking is 0.581 versus 0.192. In a separate retrospective generation check, 82 of 316 reachable held-out pairs were recovered; CliffForge predicted the correct activity direction for 53 of 82 versus 47 of 82 for ChemProp. That difference is not statistically significant, and recovered pairs are a selected sample, so it is only a lead-optimization proxy.

“The five-fold Dablander test is larger. Cliff ranking average precision is 0.391 versus 0.137. At a predicted difference greater than 1.5, CliffForge flags 20 of 106 true cliffs, with 16 false positives; ChemProp flags none. The cliff average-precision advantage does not pass the planned familywise-adjusted comparison. For local reference neighborhoods, mean Spearman is 0.509 versus 0.343, and the adjusted comparison supports the improvement. In top-ten-percent lead selection, the mean measured gain is 1.116 versus 1.002 pIC50 units. Its difference from ChemProp remains inconclusive under adjusted inference. These experiments use distinct definitions and protocols; they should not be read as one combined success rate.”

## Slide 11 — Impact, business model, and roadmap (9:05–10:00)

“If the signal generalizes, the value is a better explained shortlist for medicinal chemists deciding what to synthesize and assay. Prospective users include academic groups, biotech teams, and pharmaceutical R-and-D; none is a confirmed partner. Possible sustainability routes include an open toolkit, research collaborations, or an optional hosted workflow. The separate business proposal is enterprise licensing of a lead-optimization workflow or funded co-development with pharmaceutical research partners. That route is inspired by the way CHAI Discovery offers commercial access to its platform and has announced a Pfizer license, but CliffForge has no tested demand, pricing, or commercial terms.

“The roadmap is: working prototype today; DiffSBDD analogue generation; structure-based scoring; stereochemistry and three-dimensional evaluation; broader cliff validation; then experimental testing to check whether computational priorities translate into useful chemistry.”

## Presenter reminders

- The research vision concerns relative binding affinity, while the reported labels are pIC50 assay potency. Do not equate pIC50, docking score, and binding free energy.
- Slide 6 is proposed. Slide 7 and the run on slide 9 are implemented. The present generator transfers training-observed substituent swaps; it is not DiffSBDD.
- The Moonshot cliff result has only five strict held-out positives. The 53/82 direction comparison is a selected retrospective proxy, not an experimental result for generated candidates.
- On Dablander, CliffForge’s local ranking improvement passes the planned adjusted comparison. Its cliff AP and lead-gain differences versus ChemProp do not. The exploratory unadjusted 90% cliff interval in the results bundle does not replace the planned adjusted analysis.
- “Uncertainty” on the prototype slides means three-model ensemble spread; no calibration claim is made.

## Source notes

- Deck and implementation: presentation/main.tex, CliffForge_Teammate_Brief.pdf (pp. 2–3), src/cliffforge/generate.py, src/cliffforge/relmodel.py, and src/cliffforge/demo.py.
- Moonshot reference run: results/moonshot_mpro_tanimoto/metrics/activity_cliff_metrics.json, activity_cliff_detection_metrics.json, generation_recall.json, funnel.json, and shortlist.md in the same directory.
- Dablander five-fold evaluation: results/dablander_mmp_5fold/performance_summary.json and protocol.json.
- Business-model analogy: [CHAI Discovery commercial access](https://www.chaidiscovery.com/product) and [CHAI Discovery news listing the Pfizer license](https://www.chaidiscovery.com/news). CHAI focuses on a different, more mature molecular-design offering; no commercial terms are inferred for CliffForge.
