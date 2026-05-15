---
name: druglike-pocket-refiner
description: Refine expert or docking-derived SMILES toward more realistic drug-like kinase inhibitor candidates by improving QED/properties/synthetic-risk while preserving pocket-filling logic; prepares candidates for smiles-to-vina-docking and downstream pose-analyzer structural triage.
---

# Druglike Pocket Refiner

Use this skill when a molecule docks well or fills the pocket well but looks chemically unrealistic, too large, too flexible, too polar/lipophilic, low-QED, or hard to synthesize.

This skill is downstream of:
- `smiles-to-vina-docking`: supplies `dock_history.csv`, poses, interaction counts, and analog generation helpers.
- `pose-analyzer`: supplies pose mode-family analysis, pocket/surface visualization, iterative grid-box tuning, and pocket-graph export for structural interpretation.

Use the user's long-running workspace as the default mental model:

```text
~/vina_task2/dock_history        persistent docking ledger and pose/log store
~/vina_task2/druglike_refinement candidate/refinement CSV outputs
~/vina_task2/vina_bin            receptor, config, and Vina executable
```

The command defaults match that layout: `--history-dir` defaults to
`~/vina_task2/dock_history`, and `--outdir` defaults to
`~/vina_task2/druglike_refinement`. Keep explicit paths in examples when clarity
matters, but do not require the user to repeat them for this stable workspace.

## Core Idea

Do not optimize docking alone. Start from one or more high-confidence anchors:
- user expert SMILES,
- best contact-efficiency molecules from `dock_history.csv`,
- optional kinase-like ChEMBL reference SMILES.

When `--expert-smiles` is supplied, refinement is expert-centered by default: the expert
molecule is the parent, and historical molecules are used mainly for traceable pose/contact
context. Add `--include-history-anchors` only when you explicitly want top historical anchors
to compete as additional parents.

Then generate cleaner analogs that:
- improve QED,
- move MW/LogP/TPSA/rotatable bonds/HBD/HBA toward oral-drug-like windows,
- lower simple synthetic-risk proxies,
- preserve parent similarity and likely kinase-like features,
- avoid blindly growing huge pocket-filling molecules.

Default refinement is trimming-oriented. It uses delete/shrink/drastic edits and does not use growth/add edits unless `--allow-add` is explicitly passed. This skill is meant to rescue overgrown pocket-filling molecules, not make them larger.

Scaffold-level edits are enabled by default:
- conservative bioisostere/linker replacements such as amide to ketone/amine/ether/thioether linkers,
- flexible side-chain compression such as propyl alcohol to methoxy/hydroxy,
- methoxy to fluoro/hydroxy and ethyl to methyl,
- phenyl C-H to compact 6-membered heteroaryl swaps such as pyridyl and pyrimidyl,
- conservative manual mono-substituted phenyl to 5-membered heteroaryl templates such as furyl, thienyl, oxazolyl, and thiazolyl,
- cautious intramolecular cyclization for short flexible paths,
- phenol-to-quinone-like exploration is labeled as an alert and penalized because quinones can create redox/reactivity liabilities.

By default, these edits must improve the local QED/window/synthetic proxy value (`qve_delta > 0`). Use `--allow-qve-loss` only when exploration is more important than drug-likeness rescue.

## Main Scripts

Generate refinement candidates from history and/or expert SMILES:

```bash
python /path/to/druglike-pocket-refiner/scripts/refine_druglike_candidates.py \
  --history-dir ~/vina_task2/dock_history \
  --expert-smiles "YOUR_EXPERT_SMILES" \
  --outdir ~/vina_task2/druglike_refinement \
  --receptor ~/vina_task2/vina_bin/target.pdbqt \
  --max-rounds 10 \
  --target-refinement-score 5.0 \
  --batch-size 300 \
  --dock-top-candidates 20
```

Useful switches:
- `--include-history-anchors`: also use top historical molecules as parents when expert SMILES are provided.
- `--max-rounds N`: maximum refinement generations; default 10. `--max-iterations` and `--rounds` are accepted as old aliases.
- `--target-refinement-score X`: stop early when any candidate reaches this 0-10 score; default 5.0.
- `--batch-size N`: maximum new candidates generated per refinement round; default 300. `--max-candidates` is accepted as an old alias.
- `--progress-interval N`: print progress every N generated candidates inside each round; default 25.
- `--dock-top-candidates N`: dock the top N refined candidates through `smiles-to-vina-docking`, then write `druglike_refinement_docked_ranked.csv`.
- `--no-scaffold-edits`: use only simple delete/shrink/drastic edits.
- `--allow-qve-loss`: keep scaffold edits even when the local drug-like proxy does not improve.
- `--allow-add`: permit growth edits; normally avoid this during rescue.

`--expert-smiles` must be a valid RDKit SMILES. Invalid expert strings fail fast instead of
silently falling back to unrelated historical anchors.

This skill is a rescue/refinement path, not a final oral-drug filter. Defaults allow large,
low-QED intermediates (`--min-qed 0.0`, `--max-mw 1200`, `--max-rot-bonds 40`) and apply only
a very light MW penalty so an expert
pocket-filling molecule can move gradually toward better QED instead of being discarded in
round one. Use stricter `--min-qed`, `--max-mw`, and `--max-rot-bonds` only for a final
polished shortlist.

Quinone/catechol/Michael-acceptor/thiol-like motifs receive `structural_alert_penalty`. They are not automatically forbidden because they may be useful for exploration, but ranked output should treat them cautiously.

Outputs:
- `anchor_candidates.csv`: best historical anchors by contact efficiency and drug-likeness.
- `druglike_refinement_candidates.csv`: all generated candidates with QED/properties/risk.
- `druglike_refinement_ranked.csv`: ranked drug-like candidates.
- `druglike_refinement_to_dock.smi`: top candidates ready for `dock_smiles.py`.
- `druglike_refinement_docked_ranked.csv`: optional post-docking re-rank when `--dock-top-candidates` is used.

When `--dock-top-candidates` is used, docked refined molecules are also written back into
`dock_history.csv`. Their history rows preserve `ancestor_smiles`, `parent_smiles`, and a
`druglike_refine_g*_...` edit label, so the docking ledger can trace each refined molecule
back to the expert anchor and refinement operation.
The same history rows also receive a refinement score snapshot: `druglike_refinement_score`,
component scores, reference-similarity scores, `qve_delta`, alert penalty, and
`refinement_generation`. If a molecule already exists in history, these fields are updated
when the incoming refinement score is equal or better, while the existing docked pose/log
record is preserved unless `--redock-refined` is used.

Optionally build a kinase-biased ChEMBL reference SMILES table:

```bash
python /path/to/druglike-pocket-refiner/scripts/collect_chembl_kinase_smiles.py \
  --out ~/vina_task2/references/chembl_kinase_smiles.csv \
  --max-targets 100 \
  --max-activities 5000
```

Then use it:

```bash
python /path/to/druglike-pocket-refiner/scripts/refine_druglike_candidates.py \
  --history-dir ~/vina_task2/dock_history \
  --expert-smiles "YOUR_EXPERT_SMILES" \
  --reference-smiles-csv ~/vina_task2/references/chembl_kinase_smiles.csv \
  --outdir ~/vina_task2/druglike_refinement
```

Mine ChEMBL for literature molecules similar to current docked/expert molecules without downloading all ChEMBL:

```bash
python /path/to/druglike-pocket-refiner/scripts/chembl_similarity_mine.py \
  --history-dir ~/vina_task2/dock_history \
  --expert-smiles "YOUR_EXPERT_SMILES" \
  --outdir ~/vina_task2/chembl_similarity_mining \
  --similarity-threshold 70 \
  --target-keyword kinase
```

Outputs:
- `dock_history_structure_clusters.csv`: fingerprint/Tanimoto clusters of current docked molecules.
- `chembl_similarity_hits.csv`: ChEMBL compounds structurally similar to cluster representatives.
- `chembl_similarity_activities.csv`: activity rows for those similar compounds, including pChEMBL and ligand-efficiency fields when ChEMBL provides them.

## Interpretation

`druglike_refinement_score` is a 0-10 score. Inspect the component columns, not only the total:
`qed_component_score`, `synthetic_component_score`, `property_component_score`,
`reference_similarity_score`, `reference_partial_similarity_score`, and
`inherited_structure_score` explain why a candidate ranked well or poorly.

`contact_efficiency` is better than raw contact count because huge molecules naturally touch more atoms. The refiner also carries a structural footprint score from historical pose context when `--receptor` is provided. Use these to identify pocket-filling anchors, but final selection should depend on actual docking and `pose-analyzer` structural triage.

Anchor selection intentionally gives QED more weight than raw contact efficiency. A low-QED pocket-filling molecule can be an expert anchor, but it should not dominate automatic parent selection just because it touches many receptor atoms.

ChEMBL similarity search is chemically analogous to BLAST: it does not align proteins; it compares molecular fingerprints, usually interpreted with Tanimoto similarity. When `--reference-smiles-csv` is supplied, the refiner rewards both whole-molecule Morgan similarity and broader PatternFingerprint partial-structure similarity. Use it to ask: "Which real, literature-backed molecules does my designed SMILES resemble?" This is often cheaper and more informative than docking a huge external database.

Ligand efficiency terms in ChEMBL are useful triage signals:
- LE: potency per heavy atom, from binding free-energy style scaling.
- LLE: lipophilic ligand efficiency, usually pChEMBL minus calculated LogP.
- BEI: binding efficiency index, potency normalized by molecular weight.
- SEI: surface efficiency index, potency normalized by polar surface area.

`synthetic_score_proxy` is not a true retrosynthesis score. Treat it as an early penalty for very large, very flexible, complex, or awkward molecules. Later versions can connect to AiZynthFinder, ASKCOS, Manifold, or commercial retrosynthesis tools.

After generating `druglike_refinement_to_dock.smi`, dock promising candidates with `smiles-to-vina-docking`, then analyze their poses with `pose-analyzer`.
