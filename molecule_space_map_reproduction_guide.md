# Molecule Space Map Reproduction Guide

This guide explains how the current `Molecule Space Map` was built and how to reproduce it on a new computer.

## Goal

The map is a chemical-space coverage view, not a binding prediction model.

It answers:

- where submitted/scored molecules sit relative to the local generated registry;
- whether high, medium, and low binding molecules cluster;
- whether new submissions occupy new space or sit inside crowded basins.

## Current Output

Main HTML:

`04_results/current/phase1_molecule_system_restructure_20260513/molecule_space_umap.html`

Backup before the latest color-bin update:

`04_results/current/phase1_molecule_system_restructure_20260513/molecule_space_umap.before_binding_bins_20260515.html`

Generator script:

`src/phase1_finish_reports_from_registry.py`

Main function:
`write_submitted_highlight_map(OUT / "molecule_space_umap.html")、`
## Required Inputs
Minimum inputs:

1. Registry of generated/local molecules:

`04_results/current/phase1_molecule_system_restructure_20260513/candidate_registry.csv`

Required columns:
- `canonical_smiles`
- `candidate_id`
- optional: `binding_score`
- optional: `hypothesis_tags`
- optional: `failure_tags`

2. Shared binding sheet export:

`05_output/shared_binding_sheet_20260515_with_our_rows_local.csv`

Required columns:

- `mol_smiles`
- `binding_score`
- `成绩提交日`
- optional: `分数`
- optional: `mol_score`
- optional: `route_score`
- optional: `sa_score`
- optional: `route`

If the shared binding sheet is missing, the script falls back to local `submitted/` folders.

## Dependencies

Use the research conda Python if available:

```bash
/opt/homebrew/Caskroom/miniconda/base/envs/research/bin/python
```

Python packages:

- `rdkit`
- `numpy`
- `scikit-learn`

The map does not need Plotly or a web server. It is a standalone SVG/HTML file.
## Rebuild Command
From repo root:
```bash
/opt/homebrew/Caskroom/miniconda/base/envs/research/bin/python -c "from src.phase1_finish_reports_from_registry import OUT, write_submitted_highlight_map; write_submitted_highlight_map(OUT / 'molecule_space_umap.html')"
```

Open the resulting HTML directly in a browser:

`04_results/current/phase1_molecule_system_restructure_20260513/molecule_space_umap.html`

## Algorithm

1. Load all molecules from `candidate_registry.csv`.
2. Canonicalize SMILES with RDKit.
3. Load scored molecules from the shared binding sheet.
4. Force all scored shared-sheet molecules into the plotted set.
5. Add registry-only molecules as background points.
6. Compute Morgan fingerprints:

```python
AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
```

7. Project fingerprints to 2D with PCA:

```python
PCA(n_components=2, random_state=13)
```

8. Draw registry-only molecules as small background points.
9. Draw scored/shared-sheet molecules as large outlined markers.
10. Add tooltips using SVG `<title>`.
11. Add table of submitted/scored records below the map.

## Color Rules

Scored molecules use fixed binding-score bins:

| binding_score | color | role |
|---:|---|---|
| `<0.30` | `#f28482` coral pastel | low |
| `0.30-0.34` | `#f6bd60` apricot pastel | weak-medium |
| `0.35-0.37` | `#84a59d` sage pastel | medium |
| `0.38-0.40` | `#5dade2` blue pastel | upper-medium |
| `>=0.40` | `#7bc96f` mint green | current high |
| missing score | `#2563eb` blue | unscored |
| registry background | `#9aa0a6` gray | generated/local pool |

These colors were chosen because the previous orange/green palette was too hard to distinguish.

## Marker Rules

Background registry molecules:

- small circles;
- gray or binding-colored if a registry binding value exists;
- lower opacity.

Scored/shared-sheet molecules:

- larger filled circle;
- dark outline;
- extra faint outer ring;
- always forced into the plot even if not in the registry.

## Drag / Zoom Behavior

The HTML uses plain JavaScript, not Plotly.

Mouse wheel:

- zooms in/out around the cursor.

Drag:

- pans the map.

Reset button:

- restores original view.

Important implementation detail:

The script does not apply an SVG transform to the whole plot. Instead, it recomputes each point's `cx` and `cy` during zoom/pan while leaving `r` unchanged.

This keeps marker sizes constant during zoom, so dense overlapping regions separate without circles becoming huge.

The relevant classes/IDs are:

- `#molecule-map`
- `.mol-point`
- `#map-reset`

## Current 83-row Check

The current map was built from the 83-row shared binding sheet:

- shared rows: `83`
- valid binding rows: `83`
- `<0.30`: `41`
- `0.30-0.34`: `17`
- `0.35-0.37`: `14`
- `0.38-0.40`: `5`
- `>=0.40`: `6`

## Updating With New Scores

1. Append or refresh rows in:

`05_output/shared_binding_sheet_20260515_with_our_rows_local.csv`

2. Make sure every row has:

- `mol_smiles`
- `binding_score`
- `成绩提交日`

3. Re-run the rebuild command.

4. Verify row counts:

```bash
/opt/homebrew/Caskroom/miniconda/base/envs/research/bin/python - <<'PY'
from src.phase1_finish_reports_from_registry import load_shared_binding_sheet_rows
import math
rows = load_shared_binding_sheet_rows()
vals = []
for r in rows:
    try:
        vals.append(float(r["binding_score"]))
    except Exception:
        vals.append(math.nan)
print("rows", len(rows))
print("valid binding", sum(not math.isnan(x) for x in vals))
print("<0.30", sum(x < 0.30 for x in vals if not math.isnan(x)))
print("0.30-0.34", sum(0.30 <= x < 0.35 for x in vals if not math.isnan(x)))
print("0.35-0.37", sum(0.35 <= x < 0.38 for x in vals if not math.isnan(x)))
print("0.38-0.40", sum(0.38 <= x < 0.40 for x in vals if not math.isnan(x)))
print(">=0.40", sum(x >= 0.40 for x in vals if not math.isnan(x)))
PY
```

## Porting To A New Computer

Copy these files/directories:

- `src/phase1_finish_reports_from_registry.py`
- `04_results/current/phase1_molecule_system_restructure_20260513/candidate_registry.csv`
- `05_output/shared_binding_sheet_20260515_with_our_rows_local.csv`

Then install dependencies and run the rebuild command.

If the new machine does not have the same conda path, use any Python with RDKit and scikit-learn:

```bash
python -c "from src.phase1_finish_reports_from_registry import OUT, write_submitted_highlight_map; write_submitted_highlight_map(OUT / 'molecule_space_umap.html')"
```

## Common Failure Modes

`ModuleNotFoundError: rdkit`

Use a conda environment with RDKit installed.

`No molecules found`

Check that `candidate_registry.csv` exists and contains `canonical_smiles`.

Submitted points missing

Check that the shared sheet CSV exists and contains valid `mol_smiles`.

Wrong row count

The map uses the local CSV copy, not Google Sheets live. Refresh the local CSV before rebuilding.

## Interpretation Cautions

- PCA axes are not chemically causal features.
- Distance in the plot is approximate fingerprint distance, not binding distance.
- Similarity clusters are memory/coverage aids, not binding predictors.
- High binding color in a region means “worth inspecting,” not “neighbor must bind.”
- This map should be paired with route sanity, descriptor windows, and observed submission outcomes before selecting molecules.
