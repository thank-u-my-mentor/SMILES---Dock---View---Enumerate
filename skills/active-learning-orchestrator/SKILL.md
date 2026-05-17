---
name: active-learning-orchestrator
description: Coordinate an active-learning drug-discovery loop across smiles-to-vina-docking, druglike-pocket-refiner, and pose-analyzer. Use when the user wants an automated or semi-automated multi-iteration workflow, frontier/diversity sampling, exploration-exploitation planning, or project-level orchestration around ~/vina_task2.
---

# Active Learning Orchestrator

Use this skill when the task is no longer a single docking, refinement, or pose-analysis run, but a multi-step project loop. This is a workflow controller, not a replacement for the three domain skills.

It coordinates:

- `pose-analyzer`: train the score-space model, rebuild the dashboard, compute RF/SHAP diagnostics.
- `druglike-pocket-refiner`: use frontier seeds as parents for later candidate generation.
- `smiles-to-vina-docking`: dock selected shortlists and append to history in later stages.

The core project directory remains:

```text
~/vina_task2
```

The orchestrator script lives in the Codex skill directory:

```text
/mnt/e/Codex/skills/active-learning-orchestrator/scripts/active_learning_orchestrator.py
```

Examples assume the working environment is already active:

```bash
conda activate md
```

If `python` is not the `md` environment Python, use the full interpreter path instead:

```bash
/mnt/l/WSL/softwares/conda_envs/md/bin/python
```

## Mental Model

Think of the three chemistry/structure skills as tools and this orchestrator as the project manager.

```text
dock_history.csv + scored rows
        ↓
pose-analyzer score model + dashboard
        ↓
frontier/gap selection
        ↓
human review or druglike-pocket-refiner generation
        ↓
smiles-to-vina-docking shortlist docking
        ↓
new official scores / new history rows
        ↓
next iteration
```

The script has two modes:

- Planning mode: train/analyze, build dashboard, and write frontier seed CSVs.
- Expansion mode: additionally call `druglike-pocket-refiner` on frontier seeds to generate new SMILES; optionally dock the generated top candidates and write them back to `dock_history.csv`.

## Terminal Or IDE

Use the terminal for normal runs. The workflow calls multiple scripts and writes project outputs, so terminal logs are the clearest source of truth.

Use VS Code or PyCharm only when editing/debugging the orchestrator itself. The script is a plain Python CLI; it does not require an IDE. VS Code is usually enough. PyCharm is useful for breakpoints and variable inspection.

## Basic Commands

Show help:

```bash
python /mnt/e/Codex/skills/active-learning-orchestrator/scripts/active_learning_orchestrator.py --help
```

Run one full planning iteration: score model, dashboard, frontier seeds.

```bash
python /mnt/e/Codex/skills/active-learning-orchestrator/scripts/active_learning_orchestrator.py \
  --project ~/vina_task2 \
  --train-size 60 \
  --test-size 20 \
  --require-shap \
  --top-frontier-seeds 30 \
  --top-shortlist 100
```

`run` is the default subcommand. This also works, but the shorter form above is preferred for daily use:

```bash
python /mnt/e/Codex/skills/active-learning-orchestrator/scripts/active_learning_orchestrator.py run \
  --project ~/vina_task2 \
  --require-shap
```

Run only frontier/gap analysis from an existing feature matrix:

```bash
python /mnt/e/Codex/skills/active-learning-orchestrator/scripts/active_learning_orchestrator.py \
  frontier \
  --feature-matrix ~/vina_task2/score_space_model/binding_score_feature_matrix.csv \
  --outdir ~/vina_task2/active_learning/iteration_01 \
  --top-frontier-seeds 30 \
  --top-shortlist 100
```

Plan a 10-iteration campaign, one reviewed iteration at a time:

```bash
python /mnt/e/Codex/skills/active-learning-orchestrator/scripts/active_learning_orchestrator.py \
  --project ~/vina_task2 \
  --iterations 10 \
  --require-shap \
  --top-frontier-seeds 30 \
  --top-shortlist 100
```

Run a real expansion iteration that creates new SMILES from frontier seeds, docks the top refined candidates, and updates `dock_history.csv`:

```bash
python /mnt/e/Codex/skills/active-learning-orchestrator/scripts/active_learning_orchestrator.py \
  --project ~/vina_task2 \
  --iterations 1 \
  --require-shap \
  --generate-new-smiles \
  --dock-generated \
  --frontier-seeds-to-refine 5 \
  --refiner-rounds 2 \
  --refiner-batch-size 200 \
  --refiner-beam-size 3 \
  --dock-top-candidates 20
```

Use `--generate-new-smiles` without `--dock-generated` when you want to inspect the generated candidate CSV before spending Vina time.

Dry-run the commands before running them:

```bash
python /mnt/e/Codex/skills/active-learning-orchestrator/scripts/active_learning_orchestrator.py \
  --project ~/vina_task2 \
  --require-shap \
  --dry-run
```

## Outputs

Each planning iteration writes to:

```text
~/vina_task2/active_learning/iteration_XX/
```

Important files:

- `frontier_seeds.csv`: best parent seeds for exploration/exploitation.
- `official_score_recommendations.csv`: docked, unscored molecules recommended for the next official binding-score cycle.
- `active_learning_shortlist.csv`: wider ranked set for human inspection.
- `active_learning_shortlist.smi`: SMILES list with `seq_id`, useful for downstream tools.
- `gap_regions_pure_smiles.csv`: 2-D Pure-SMILES-space grid summary for visualization.
- `gap_regions_structural_interaction.csv`: 2-D structural-interaction-space grid summary.
- `iteration_manifest.json`: settings and human next-step note for that iteration.

`frontier_score` combines novelty, distance from scored examples, predicted/known quality, drug-likeness, ChEMBL scaffold support, and available pocket-contact signal. It is a selection heuristic, not proof that a molecule will score well.

When `~/vina_task2/high_low_classifier/high_low_classifier_predictions.csv` is
available, `frontier_score` also includes `structure_high_probability`, a
high-vs-low classifier signal from `pose-analyzer/train_high_low_classifier.py`.
This shifts generation away from superficial descriptor optimization and toward
molecules that structurally resemble the known high-score group.

`official_score_recommendations.csv` is narrower than `frontier_seeds.csv`: it prefers molecules that already have docking/space information but do not yet have `official_binding_score`. Use it as the first place to look when choosing the next molecules to score experimentally or through the official scoring source.

By default, official recommendations and refiner parent seeds require
`structure_high_probability >= 0.60` when classifier predictions exist. Adjust
with:

```bash
--min-structure-high-probability-for-official 0.60 \
--min-structure-high-probability-for-refiner 0.60
```

## Iteration Strategy

Ten iterations is a reasonable project-level campaign target, but the loop is meaningful only when new SMILES are generated and/or new docking/official-score rows are added between model updates. If no new data is added, repeated runs mostly analyze the same history again.

`--iterations 10` now runs ten planning passes. Add `--stop-after-frontier` only when you intentionally want to stop after the first frontier seed set for manual review.

A practical rhythm is:

1. Run iteration N planning.
2. Inspect `frontier_seeds.csv` and the dashboard.
3. Run expansion mode on selected/frontier seeds to generate new SMILES.
4. Dock top generated candidates so `dock_history.csv` grows.
5. Run iteration N+1.

After official scores are obtained, sync/update `binding_score_analysis.csv`; then the next iteration can learn from the new teacher signal.

Exploration/exploitation balance can be adjusted with weights:

```bash
--weight-novelty 0.35 --weight-uncertainty 0.30 --weight-quality 0.10
```

for more exploration, or:

```bash
--weight-quality 0.30 --weight-druglike 0.20 --weight-contact 0.15
```

for more conservative exploitation.

## Dashboard

After a run, start the dashboard server from the terminal:

```bash
cd ~/vina_task2/score_space_dashboard
python -m http.server 8765 --bind 127.0.0.1
```

Open:

```text
http://127.0.0.1:8765
```

Do not open `0.0.0.0` in the browser; that is only a server listening address.

When active-learning recommendations are available, the dashboard loads:

```text
~/vina_task2/score_space_dashboard/data/active_learning_recommendations.json
```

Point highlighting:

- Gold glowing outer ring: recommended for the next official binding-score cycle.
- Teal outer ring: frontier seed for exploration/refinement.
- The detail panel shows `active_learning_rank`, `frontier_score`, novelty, uncertainty, and both space coordinates.
