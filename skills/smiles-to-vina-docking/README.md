# Local SMILES-to-Docking-to-Pose Report Workflow

This is a no-Docker local workflow for starting from one receptor and one ligand
SMILES, making conservative ligand analogs, docking them with Vina/GNINA/GLINA,
analyzing pose families, and writing a static HTML report.

The simplest project command is:

```bash
python scripts/run_project.py \
  --smiles "CCO" \
  --receptor /path/to/receptor.pdbqt \
  --config /path/to/vina_config.txt \
  --analogs
```

Outputs are written under `~/vina_task2` by default:

```text
~/vina_task2/
  dock_history/dock_history.csv
  dock_history/poses/
  pose_analysis/mode_family_summary.csv
  pose_analysis/visualizations/*.pml
  pose_report/index.html
```

Open `pose_report/index.html` in a browser. If PyMOL is installed, open the linked
`.pml` files from the report. Add `--save-pse` when you want the workflow to ask
PyMOL to save `.pse` session scenes too.

## Inputs

You provide:

- `--smiles`: starting ligand SMILES.
- `--receptor`: prepared receptor file, preferably PDBQT for Vina-style docking.
- `--config`: Vina-style box config containing `center_x`, `center_y`,
  `center_z`, `size_x`, `size_y`, and `size_z`.

The workflow chooses practical defaults for:

- project directory: `~/vina_task2`
- analog rounds: `0` unless `--analogs` is passed
- analog batch size: `24`
- Vina exhaustiveness: `8`
- Vina modes: `9`
- CPU: set by `dock_smiles.py`, about physical-core count and capped at `16`

`run_project.py` performs a strict environment check before docking. Use
`--skip-env-check` only when you have already verified the environment and want
to bypass that startup check.

## Conda Route A: md Environment With AutoDock Vina

Use this route when `conda activate md` has RDKit, Meeko, OpenBabel, and Vina.

```bash
conda activate md
python scripts/check_local_env.py --docking-engine vina --strict
python scripts/run_project.py \
  --smiles "YOUR_STARTING_SMILES" \
  --receptor /path/to/receptor.pdbqt \
  --config /path/to/vina_config.txt \
  --engine vina \
  --analogs
```

If Vina is not on `PATH`, pass it once:

```bash
python scripts/run_project.py \
  --smiles "YOUR_STARTING_SMILES" \
  --receptor /path/to/receptor.pdbqt \
  --config /path/to/vina_config.txt \
  --engine-exe /path/to/vina \
  --analogs
```

## Conda Route B: glina/gnina Environment

Use this route when `conda activate glina` contains GNINA/GLINA-like docking.
Some installations name the executable `gnina`; others may use `glina`.

```bash
conda activate glina
python scripts/check_local_env.py --docking-engine gnina --strict
python scripts/run_project.py \
  --smiles "YOUR_STARTING_SMILES" \
  --receptor /path/to/receptor.pdbqt \
  --config /path/to/vina_config.txt \
  --engine gnina \
  --analogs
```

If your executable is really named `glina`:

```bash
python scripts/check_local_env.py --docking-engine glina --strict
python scripts/run_project.py \
  --smiles "YOUR_STARTING_SMILES" \
  --receptor /path/to/receptor.pdbqt \
  --config /path/to/vina_config.txt \
  --engine glina \
  --analogs
```

`run_project.py` still uses Vina-style command-line arguments. If a local
GNINA/GLINA build does not accept Vina-compatible flags, use AutoDock Vina or
wrap that executable with a small compatible script.

## Minimal Package Install

No Docker image is required. A normal recent Python conda environment is enough.
The most convenient install is usually:

```bash
conda create -n docking-local -c conda-forge python=3.11 rdkit meeko openbabel numpy scikit-learn pymol-open-source vina
conda activate docking-local
python scripts/check_local_env.py --docking-engine vina --strict
```

If you prefer pip for pure-Python packages, keep RDKit/OpenBabel/Vina from
conda-forge when possible:

```bash
python -m pip install numpy scikit-learn openpyxl
```

## One-Command Examples

Dock only the starting ligand:

```bash
python scripts/run_project.py \
  --smiles "YOUR_STARTING_SMILES" \
  --receptor /path/to/receptor.pdbqt \
  --config /path/to/vina_config.txt
```

Dock the starting ligand plus one conservative analog batch:

```bash
python scripts/run_project.py \
  --smiles "YOUR_STARTING_SMILES" \
  --receptor /path/to/receptor.pdbqt \
  --config /path/to/vina_config.txt \
  --analogs
```

Use a separate output directory:

```bash
python scripts/run_project.py \
  --smiles "YOUR_STARTING_SMILES" \
  --receptor /path/to/receptor.pdbqt \
  --config /path/to/vina_config.txt \
  --project-dir ./my_target_project \
  --analogs
```

Ask PyMOL to save PSE scenes:

```bash
python scripts/run_project.py \
  --smiles "YOUR_STARTING_SMILES" \
  --receptor /path/to/receptor.pdbqt \
  --config /path/to/vina_config.txt \
  --analogs \
  --save-pse
```

## What The Report Means

`pose-analyzer` uses `official_binding_score` when that column exists. For a
general arbitrary-protein project with no external score, it automatically uses
`affinity_kcal_mol` with lower affinity treated as better, then records this as
`analysis_score_source=affinity_kcal_mol`.

That fallback is for structural triage only. It is useful for:

- choosing a reference pose for visualization,
- seeing which residues are repeatedly contacted,
- comparing pose stability and contact coverage,
- deciding which analogs deserve manual inspection.

It is not experimental validation and should not be treated as proof of binding.

## Practical Notes

- Keep receptor/config files outside the skill directory and pass their paths.
- Keep one project directory per protein/box when possible.
- First runs can be slow because ligand PDBQT preparation and docking are real
  compute jobs.
- Generated analogs are conservative by default and pass light property gates.
- PyMOL is optional for CSV/PML output. It is required only for SASA-based surface
  detection and automatic `.pse` saving; otherwise the script uses a fallback
  surface estimate and still writes `.pml`.
