---
name: smiles-to-vina-docking
description: Dock SMILES ligands with AutoDock Vina, keep a simple docking history, generate small analogs, rank historical molecules against a chosen reference, and export Task2 result.zip files.
---

# SMILES To Vina Docking

## Purpose

Use this skill when the user wants a practical SMILES -> Vina -> history workflow.

The current design is intentionally split into small tools:

- `scripts/run_project.py`: easiest local entry point. Given one receptor, one
  Vina-style config, and one starting SMILES, it runs environment checks, docks
  the input and optional analogs, calls `pose-analyzer`, and writes a static HTML
  pose report.
- `scripts/check_local_env.py`: verify that the active conda/Python environment
  has RDKit, Meeko, a docking executable, and optional PyMOL/OpenBabel support.
- `scripts/dock_smiles.py`: daily-use command. Dock one input SMILES, optionally generate analogs, and append results to `dock_history.csv`.
- `scripts/sync_binding_scores.py`: merge a local Task2 scoring xlsx into `dock_history.csv` and analyze score-vs-docking correlations.
- `scripts/rank_history.py`: compare all historical molecules against one reference SMILES and write `ranked_history.csv`.
- `scripts/export_result_zip.py`: package one selected molecule into a Task2-style `result.zip`.
- `scripts/lead_optimizer_legacy.py`: old large script kept only as backup. Do not use it for new work unless the user explicitly asks for legacy behavior.

Avoid odd project names. User-facing files should use simple English that is easy for Chinese users to understand: `dock_history.csv`, `history_config.json`, `ranked_history.csv`, `poses/`, `logs/`, `sdf/`, and `pdbqt/`.

## Mental Model

Think of the history directory as a persistent docking ledger, not as a temporary run folder. Keep using the same `--history-dir` over time. Do not create timestamped history folders unless the user explicitly asks for that.

History writes must be live and durable. `dock_smiles.py` appends each completed
molecule to `dock_history.csv`, flushes/fsyncs it, and only then prints the
`[dock] ... affinity=...` line. If a long run is terminated, every already
printed `[dock]` row should be readable from the CSV; only the molecule currently
inside ligand preparation/docking may be absent or incomplete.

The history directory contains:

- `dock_history.csv`: one chronological row per attempted molecule.
- `history_config.json`: saved receptor/config/Vina paths after the first run, so later commands can be short.
- `poses/`: Vina output PDBQT poses.
- `logs/`: Vina logs.
- `sdf/` and `pdbqt/`: prepared ligand files.
- `ranked_history.csv`: optional ranking output created only when ranking is requested.

For the user's long-running Task2 workspace, assume this stable layout unless told otherwise:

```text
~/vina_task2/
  dock_history/
    dock_history.csv
    history_config.json
    binding_score_analysis.csv
    binding_score_correlations.csv
    binding_score_missing_smiles.csv
    sdf/
    pdbqt/
    poses/
    logs/
    prep_logs/
  druglike_refinement/
  binding_score_mode_analysis/
  vina_bin/
```

Both `dock_smiles.py` and the shared docking utilities default to
`~/vina_task2/dock_history`. Users may still pass `--history-dir` explicitly for a
different project.

There is no permanent global lead in the docking history. The ledger records ancestry directly with SMILES strings, not with a pile of separate ids.

Keep the public history table simple: write one ligand identity column named
`smiles`. Internally, helper code may still read old histories containing
`input_smiles` and `canonical_smiles`, but new `dock_history.csv` files should
not keep both columns. Use canonical SMILES for duplicate detection.

For GNINA/GLINA runs, `cnn_pose_score` records the best mode's `REMARK CNNscore`
from the output PDBQT/log. This replaces the old public `whole_rmsd` column in
new compact histories. Keep `inner_rmsd` as the local multi-pose consistency
metric; use `cnn_pose_score` as the pose plausibility signal.

## Requirements

The user must provide their own environment-specific paths.

Required tools:

- Python with RDKit.
- AutoDock Vina.
- Meeko `mk_prepare_ligand.py` or the Meeko Python module.
- OpenBabel `obabel` as a fallback ligand-preparation path.
- A receptor PDBQT.
- A Vina config file with `center_x/y/z` and `size_x/y/z`.

Recommended WSL environment check:

```bash
conda activate md
python -c "import rdkit; print('rdkit ok')"
python -m meeko.cli.mk_prepare_ligand --help
which obabel
```

Recommended local conda checks on Windows/Linux/macOS:

```bash
conda activate md
python scripts/check_local_env.py --docking-engine vina --strict

conda activate glina
python scripts/check_local_env.py --docking-engine gnina --strict
```

The `glina` environment name is user-specific. If the executable is actually
called `glina`, use `--docking-engine glina` and `run_project.py --engine glina`.
If it is called `gnina`, use `--engine gnina`.

## One-Command Local Project

For a general arbitrary protein and one starting ligand, prefer:

```bash
python /path/to/smiles-to-vina-docking/scripts/run_project.py \
  --smiles "<STARTING_SMILES>" \
  --receptor /path/to/receptor.pdbqt \
  --config /path/to/vina_config.txt \
  --analogs
```

This writes:

- `~/vina_task2/dock_history/dock_history.csv`
- `~/vina_task2/dock_history/poses/`
- `~/vina_task2/pose_analysis/`
- `~/vina_task2/pose_report/index.html`

Defaults are deliberately conservative: no analogs unless `--analogs` is passed,
one analog round when `--analogs` is used, `--batch-size 24`, `--exhaustiveness 8`,
and `--num-modes 9`.

`run_project.py` runs `check_local_env.py --strict` before docking unless
`--skip-env-check` is supplied.

Use Vina explicitly:

```bash
python /path/to/smiles-to-vina-docking/scripts/run_project.py \
  --smiles "<STARTING_SMILES>" \
  --receptor /path/to/receptor.pdbqt \
  --config /path/to/vina_config.txt \
  --engine vina \
  --analogs
```

Use a GNINA/GLINA-style executable if it accepts Vina-compatible CLI flags:

```bash
python /path/to/smiles-to-vina-docking/scripts/run_project.py \
  --smiles "<STARTING_SMILES>" \
  --receptor /path/to/receptor.pdbqt \
  --config /path/to/vina_config.txt \
  --engine gnina \
  --analogs
```

Extra GNINA/Vina-compatible options can be passed through with repeated
`--engine-arg` tokens. Use the `--engine-arg=value` form for options that begin
with `--`:

```bash
python /path/to/smiles-to-vina-docking/scripts/dock_smiles.py \
  --smiles "<SMILES>" \
  --history-dir ~/vina_task2/dock_history \
  --vina gnina \
  --engine-arg=--cnn_scoring \
  --engine-arg=refinement
```

Add `--save-pse` only when PyMOL is installed and you want automatic `.pse`
session files. The workflow always writes `.pml` view scripts when pose analysis
succeeds.

## Daily Docking
First run needs the paths:

```bash
python /path/to/smiles-to-vina-docking/scripts/dock_smiles.py \
  --smiles "<SMILES>" \
  --nickname "lead_001" \
  --history-dir ~/vina_task2/dock_history \
  --receptor /path/to/target.pdbqt \
  --config /path/to/vina_config.txt \
  --vina /path/to/vina \
  --exhaustiveness 8 \
  --num-modes 9
```

This writes `history_config.json`. Later runs can usually be short:

```bash
python /path/to/smiles-to-vina-docking/scripts/dock_smiles.py \
  --smiles "<SMILES>" \
  --history-dir ~/vina_task2/dock_history
```

Because `~/vina_task2/dock_history` is the default, this can also be:

```bash
python /path/to/smiles-to-vina-docking/scripts/dock_smiles.py \
  --smiles "<SMILES>"
```

If the molecule already exists, the tool reports the existing `seq_id` and does not redock by default. Use `--redock-existing` only when the user truly wants another docking event for the same molecule.
`--nickname` is optional and stored as the user-facing label for that row when provided.

If the existing row has no affinity or has a failure `reason`, the tool treats it as an
incomplete docking record and retries that same `seq_id` instead of silently skipping it.
This protects rows created by an earlier bad Vina path or failed ligand preparation.

If many rows suddenly show `vina_executable_not_found` or `vina_failed_or_no_affinity`
with logs like `No such file or directory: 'vina'`, pass the full Vina path once again,
for example `--vina ~/vina_task2/vina_bin/vina`. The tool also tries the common
`<history parent>/vina_bin/vina` location when an old config only contains bare `vina`.
If `--vina` is accidentally passed as the directory `~/vina_task2/vina_bin`, the code now
normalizes it to `~/vina_task2/vina_bin/vina` before saving `history_config.json`.

## Generating Analogs

`dock_smiles.py` can also generate nearby molecules before docking them:

```bash
python /path/to/smiles-to-vina-docking/scripts/dock_smiles.py \
  --smiles "<SMILES>" \
  --history-dir ~/vina_task2/dock_history \
  --max-rounds 1 \
  --batch-size 30
```

Cost controls:

- `--max-rounds 0`: dock only the input SMILES. `--batch-size` is ignored.
- `--max-rounds 1 --batch-size 30`: dock the input molecule, then up to 30 new nonduplicate analogs.
- `--max-rounds 2+`: descendants can become parents for another wave, so runtime and chemical distance increase.
- `--batch-size`: target number of new analogs per round, after skipping duplicates.
- `--cpu`: Vina worker count. By default the script uses about the physical-core count and caps at 16; on a 16-core/32-thread Ryzen 9 9950X this defaults to 16.
- `cpu = 32` in `config.txt` is also honored if present, unless overridden by `--cpu`.
- `--exhaustiveness`, `--num-modes`, `--energy-range`: Vina runtime and pose-diversity controls.

Analog edit controls:

- Default: consider all edit families.
- `--add`: only add small substituents.
- `--delete`: only remove substituents.
- `--shrink`: delete or replace bulky substituents with smaller groups.
- `--replace`: only apply conservative replacement/bioisostere edits.
- `--drastic`: include broader edits such as aromatic C-to-N swaps, replacement
  edits, and scaffold contraction-style moves.
- `--edit-mode shrink,delete,replace`: explicit comma-separated selection.

The default batch is deliberately interleaved by edit family, so a molecule with
valid replacements should not produce a long add-only run. Replacement edits
include amide/linker swaps, methoxy/hydroxy/fluoro simplifications, aromatic
C-to-N walks, and phenyl-to-five-member heteroaryl replacements when the parent
structure allows them. Some starting SMILES legitimately have no removable
substituent or shrinkable branch; in that case delete/shrink may contribute no
candidates, but replace/add/drastic should still be mixed.

Batch order is randomized by default for exploration. Use `--deterministic-batch` only for debugging.

Generated analogs pass a lightweight property gate before docking. This prevents
multi-round analog generation from spending time on molecules that are already
far outside a practical chemical space. The gate applies only to generated
analogs, not to the manually supplied input SMILES, and these gate values are not
written into `dock_history.csv`.

Default analog property gate:

- `--min-analog-qed 0.20`
- `--max-analog-mw 650`
- `--max-analog-rot-bonds 14`
- `--max-analog-tpsa 180`

For stricter analog exploration, tighten these values, for example
`--min-analog-qed 0.40 --max-analog-mw 550 --max-analog-rot-bonds 10 --max-analog-tpsa 140`.

## History Columns

`dock_history.csv` should stay readable and chronological. Prefer one table over many redundant CSVs.

Recommended columns in the current minimal ledger:

- `seq_id`: simple sequence id such as `S000001`.
- `timestamp`: when the row was created.
- `nickname`: user-facing label, preserved if the user edits it.
- `ancestor_smiles`: the original ancestor for the family.
- `parent_smiles`: the direct parent used to generate this row, empty for the first input row.
- `edit_label`: the transformation label for generated rows, or `input_smiles` for the first dock.
- `smiles`: the normalized molecule string used by the code. Old
  `input_smiles`/`canonical_smiles` columns are read for compatibility only.

Important docking and chemistry columns:

- `affinity_kcal_mol`
- `inner_rmsd`: average RMSD of retained modes against the best-affinity mode.
- `inner_cluster_fraction`: fraction of modes within `--internal-cluster-rmsd-cutoff` of the best-affinity mode.
- `cnn_pose_score`: GNINA/GLINA CNN pose score when available.
- `hbond_count`, `hydrophobic_count`, `vdw_contact_count`, `pi_contact_count`, `ch_pi_count`
- `pose_path`, `log_path`, `sdf_path`, `pdbqt_path`, `prep_log_path`
- `reason`

Interaction counts are fast geometric heuristics. They are useful for triage, but they are not a full PLIP-style interaction analysis.

## Binding Score Sync

When the user downloads the official Task2 score table as xlsx, merge it into the docking history:

```bash
python /path/to/smiles-to-vina-docking/scripts/sync_binding_scores.py \
  --xlsx ~/vina_task2/dock_history/Task2评分统计表.xlsx \
  --history-dir ~/vina_task2/dock_history \
  --receptor ~/vina_task2/vina_bin/target.pdbqt \
  --update-history
```

The xlsx must contain columns equivalent to `mol_smiles` and `binding_scores`. The script canonicalizes SMILES for matching, preserves original history rows, writes `official_binding_score` back to `dock_history.csv`, and creates:

- `binding_score_analysis.csv`: matched rows plus docking, interaction, and simple pocket-distance features.
- `binding_score_correlations.csv`: Pearson and Spearman correlations against the official score.
- `binding_score_missing_smiles.csv`: score-table SMILES not yet present in the dock history.

To automatically dock missing score-table molecules before regenerating analysis:

```bash
python /path/to/smiles-to-vina-docking/scripts/sync_binding_scores.py \
  --xlsx ~/vina_task2/dock_history/Task2评分统计表.xlsx \
  --history-dir ~/vina_task2/dock_history \
  --receptor ~/vina_task2/vina_bin/target.pdbqt \
  --config ~/vina_task2/vina_bin/idz6F_config.txt \
  --vina ~/vina_task2/vina_bin/vina \
  --dock-missing \
  --update-history
```

Use `--dock-missing-limit 5` for a small test batch.

## Ranking

Ranking is separate from docking. Use it only when the user wants to compare the whole history against a reference molecule:

```bash
python /path/to/smiles-to-vina-docking/scripts/rank_history.py \
  --smiles "<REFERENCE_SMILES>" \
  --history-dir ~/vina_task2/dock_history
```

If the reference SMILES is missing from `dock_history.csv`, `rank_history.py` docks it first by calling `dock_smiles.py`, then ranks.

Ranking output:

- `ranked_history.csv`

Ranking uses normal-distribution style scores over the current history:

- affinity: more negative is better.
- reference RMSD: lower is more similar to the selected reference.
- inner RMSD: lower means poses cluster better around the best-affinity mode.
- whole RMSD: lower means all retained modes are more internally consistent.
- interaction score: more simple contact counts is better.

Do not treat Vina affinity alone as proof. Small kcal/mol differences are noisy.

## Result.zip Export

After a molecule exists in the history, export a Task2-style package:

```bash
python /path/to/smiles-to-vina-docking/scripts/export_result_zip.py \
  --history-dir ~/vina_task2/dock_history \
  --smiles "<SELECTED_SMILES>" \
  --route "<REACTANT1.REACTANT2>>PRODUCT>" \
  --outdir ~/vina_task2/result_export
```

The script writes:

- `result.zip`
- `result.csv`
- `result.log`

`--history-dir` is preferred. `--ledger-dir` remains as a backward-compatible alias only.

If `--route` is omitted, the script may write a placeholder or simple guessed route. For official submission, the route still needs chemical review.

## Chemistry Guidance

Start with small, interpretable edits:

- Aromatic C-H substitution: F, Cl, Br, Me, OH, NH2, CN, CHO, OMe.
- Substituent deletion and shrinkage: remove or reduce bulky non-ring branches.
- Conservative heteroatom walks and replacements: aromatic C-to-N swaps,
  amide/linker bioisosteres, methoxy simplification, and phenyl-to-furyl,
  thienyl, oxazolyl, or thiazolyl replacements where chemically valid.
- Drastic edits only after small edits are not enough.

Use docking results as error exclusion and prioritization, not experimental truth. Common reasons a "scientifically reasonable" polar addition fails in Vina:

- desolvation penalty,
- wrong H-bond geometry,
- extra torsional cost,
- steric strain,
- pose flipping,
- protonation or PDBQT atom-type issues.

## Packaging This Skill

To share the skill, zip the folder:

```text
smiles-to-vina-docking/
  SKILL.md
  scripts/dock_smiles.py
  scripts/rank_history.py
  scripts/dock_utils.py
  scripts/export_result_zip.py
  references/docking_heuristics.md
```

Do not include local receptor files, ligand results, Vina binaries, or machine-specific paths unless the user intentionally wants to share those files.

Install by extracting to a Codex skills directory, for example:

```text
~/.codex/skills/smiles-to-vina-docking
```

or:

```text
C:\Users\<user>\.codex\skills\smiles-to-vina-docking
```

## WSL Path Note

On WSL, Windows `E:\...` paths become `/mnt/e/...`.

For speed, prefer keeping the active history under Linux home, for example:

```bash
~/vina_task2/dock_history
```

Use the Windows/E drive copy as backup when possible. Large histories with hundreds of small PDBQT files can be slow under `/mnt/e` because WSL file access becomes the bottleneck.
