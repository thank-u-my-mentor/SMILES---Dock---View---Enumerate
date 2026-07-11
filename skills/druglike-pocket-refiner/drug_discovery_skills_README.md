# Vina + Pose Analyzer + Druglike Refiner + Active Learning Workflow

This README explains the local drug-design workflow built around four Codex
skills:

- `smiles-to-vina-docking`
- `druglike-pocket-refiner`
- `pose-analyzer`
- `active-learning-orchestrator`

The workflow is designed for a project directory like:

```text
~/vina_task2
├── dock_history/
│   ├── dock_history.csv
│   ├── binding_score_analysis.csv
│   ├── poses/
│   ├── pdbqt/
│   ├── sdf/
│   └── logs/
├── references/
│   └── chembl_kinase_smiles.csv
├── vina_bin/
│   ├── target.pdbqt
│   └── idz6F_config.txt
├── score_space_model/
├── high_low_classifier/
├── score_space_dashboard/
└── active_learning/
```

## Mental Model

The system has two loops.

### Data-Building Loop

```text
SMILES
  -> smiles-to-vina-docking
  -> Vina pose + affinity + descriptors
  -> dock_history.csv
  -> official binding scores added later
```

This loop builds enough history data to learn from. `dock_history.csv` is the
main memory. `binding_score_analysis.csv` or `official_binding_score` columns
provide the teacher signal.

### Learning + Generation Loop

```text
dock_history.csv + official scores
  -> pose-analyzer score/space model
  -> high-vs-low structural classifier
  -> active-learning-orchestrator
  -> druglike-pocket-refiner
  -> smiles-to-vina-docking
  -> new docked molecules added to history
```

`pose-analyzer` is no longer just visualization. It now also builds chemical
space, pocket-interaction space, high/low structure classifiers, and activity
cliff reports.

`active-learning-orchestrator` should not blindly chase novelty. It should use
`structure_high_probability` from the high/low classifier to pick refiner seeds
and official-score recommendations.

## Skill Roles

### 1. smiles-to-vina-docking

Purpose:

- Convert SMILES to 3D ligand structures.
- Run AutoDock Vina.
- Save affinity, pose files, descriptors, and contact counts into
  `dock_history.csv`.
- Append new molecules without renaming existing `Sxxxxx` IDs.

Typical command:

```bash
python /mnt/e/Codex/skills/smiles-to-vina-docking/scripts/dock_smiles.py \
  --smiles "YOUR_SMILES" \
  --nickname "YOUR_NAME" \
  --history-dir ~/vina_task2/dock_history \
  --receptor ~/vina_task2/vina_bin/target.pdbqt \
  --config ~/vina_task2/vina_bin/idz6F_config.txt \
  --vina ~/vina_task2/vina_bin \
  --exhaustiveness 8 \
  --num-modes 9
```

### 2. druglike-pocket-refiner

Purpose:

- Take expert or historical SMILES as parent molecules.
- Generate drug-like edits: trimming, shrinking, heteroaryl replacement,
  cyclization, and related transformations.
- Optionally dock top refined candidates and write them back to history.
- Use ChEMBL kinase references to avoid chemically isolated fantasy molecules.

Typical command:

```bash
python /mnt/e/Codex/skills/druglike-pocket-refiner/scripts/refine_druglike_candidates.py \
  --history-dir ~/vina_task2/dock_history \
  --expert-smiles "YOUR_PARENT_SMILES" \
  --reference-smiles-csv ~/vina_task2/references/chembl_kinase_smiles.csv \
  --outdir ~/vina_task2/druglike_refinement \
  --target-refinement-score 8.0 \
  --batch-size 300 \
  --dock-top-candidates 50 \
  --max-rounds 10
```

### 3. pose-analyzer

Purpose:

- Analyze Vina poses and ligand-receptor contacts.
- Build Pure-SMILES-space and structural-interaction-space.
- Build local 5 A pocket embeddings from receptor atoms/residues around each
  ligand pose.
- Train diagnostic Random Forest regression.
- Train high-vs-low structural classifiers.
- Generate dashboard data and activity cliff reports.

Run score-space model:

```bash
python /mnt/e/Codex/skills/pose-analyzer/scripts/score_space_model.py \
  --history-csv ~/vina_task2/dock_history/dock_history.csv \
  --analysis-csv ~/vina_task2/dock_history/binding_score_analysis.csv \
  --reference-smiles-csv ~/vina_task2/references/chembl_kinase_smiles.csv \
  --receptor ~/vina_task2/vina_bin/target.pdbqt \
  --outdir ~/vina_task2/score_space_model \
  --train-size 60 \
  --test-size 20
```

Train high/low classifier:

```bash
python /mnt/e/Codex/skills/pose-analyzer/scripts/train_high_low_classifier.py \
  --feature-matrix ~/vina_task2/score_space_model/binding_score_feature_matrix.csv \
  --outdir ~/vina_task2/high_low_classifier
```

Build dashboard:

```bash
python /mnt/e/Codex/skills/pose-analyzer/scripts/build_score_space_dashboard.py \
  --feature-matrix ~/vina_task2/score_space_model/binding_score_feature_matrix.csv \
  --outdir ~/vina_task2/score_space_dashboard
```

Serve dashboard:

```bash
cd ~/vina_task2/score_space_dashboard
python -m http.server 8765 --bind 127.0.0.1
```

Open:

```text
http://127.0.0.1:8765
```

Do not open `0.0.0.0` in the browser. It is a listening address, not the client
URL.

### 4. active-learning-orchestrator

Purpose:

- Coordinate the multi-step loop.
- Rebuild score-space model and dashboard.
- Read high/low classifier predictions.
- Select high-priority seeds and official-score recommendations.
- Optionally call refiner and dock generated candidates.

Recommended command:

```bash
python /mnt/e/Codex/skills/active-learning-orchestrator/scripts/active_learning_orchestrator.py \
  --project ~/vina_task2 \
  --iterations 10 \
  --generate-new-smiles \
  --dock-generated \
  --classifier-predictions ~/vina_task2/high_low_classifier/high_low_classifier_predictions.csv \
  --min-structure-high-probability-for-refiner 0.60 \
  --min-structure-high-probability-for-official 0.60 \
  --top-official-recommendations 3 \
  --top-frontier-seeds 30 \
  --top-shortlist 100
```

## Recommended End-To-End Run Order

1. Build or update docking history.

```bash
python /mnt/e/Codex/skills/smiles-to-vina-docking/scripts/dock_smiles.py ...
```

2. Add or update official binding scores in:

```text
~/vina_task2/dock_history/binding_score_analysis.csv
```

3. Build score-space features.

```bash
python /mnt/e/Codex/skills/pose-analyzer/scripts/score_space_model.py \
  --history-csv ~/vina_task2/dock_history/dock_history.csv \
  --analysis-csv ~/vina_task2/dock_history/binding_score_analysis.csv \
  --reference-smiles-csv ~/vina_task2/references/chembl_kinase_smiles.csv \
  --receptor ~/vina_task2/vina_bin/target.pdbqt \
  --outdir ~/vina_task2/score_space_model
```

4. Train high/low structural classifier.

```bash
python /mnt/e/Codex/skills/pose-analyzer/scripts/train_high_low_classifier.py \
  --feature-matrix ~/vina_task2/score_space_model/binding_score_feature_matrix.csv \
  --outdir ~/vina_task2/high_low_classifier
```

5. Build dashboard.

```bash
python /mnt/e/Codex/skills/pose-analyzer/scripts/build_score_space_dashboard.py \
  --feature-matrix ~/vina_task2/score_space_model/binding_score_feature_matrix.csv \
  --outdir ~/vina_task2/score_space_dashboard
```

6. Run active learning.

```bash
python /mnt/e/Codex/skills/active-learning-orchestrator/scripts/active_learning_orchestrator.py \
  --project ~/vina_task2 \
  --iterations 10 \
  --generate-new-smiles \
  --dock-generated \
  --classifier-predictions ~/vina_task2/high_low_classifier/high_low_classifier_predictions.csv
```

7. Review:

```text
~/vina_task2/active_learning/iteration_XX/frontier_seeds.csv
~/vina_task2/active_learning/iteration_XX/official_score_recommendations.csv
~/vina_task2/score_space_dashboard
```

## Environment

The current setup assumes WSL/Linux and a conda environment named `md`.

Recommended Python:

```text
Python 3.10
```

Core Python packages:

```text
rdkit
numpy
pandas
scikit-learn
scipy
torch
umap-learn
shap
chembl_webresource_client
networkx
matplotlib
seaborn
joblib
```

Docking / chemistry command-line tools:

```text
AutoDock Vina
Open Babel
Meeko
```

Optional visualization:

```text
PyMOL
```

### Conda Environment Example

```bash
conda create -n md python=3.10 -y
conda activate md
conda install -c conda-forge rdkit numpy pandas scikit-learn scipy matplotlib seaborn networkx openbabel -y
pip install torch umap-learn shap chembl_webresource_client meeko
```

AutoDock Vina can be installed via conda-forge or provided as a local binary:

```bash
conda install -c conda-forge vina -y
```

In this project, Vina/config files are usually kept under:

```text
~/vina_task2/vina_bin
```

## Docker Option

Docker is useful if collaborators need a reproducible environment. A practical
image should include:

- Python 3.10
- RDKit
- scikit-learn
- PyTorch CPU
- Open Babel
- Meeko
- AutoDock Vina
- the four skill directories mounted at `/mnt/e/Codex/skills`
- the project mounted at `/home/user/vina_task2`

Minimal Docker strategy:

```text
base image: mambaorg/micromamba or continuumio/miniconda3
install conda packages from conda-forge
pip install torch umap-learn shap chembl_webresource_client meeko
mount local data instead of baking ~/vina_task2 into the image
```

Do not bake private receptor files, official scores, or unpublished molecules
into a public image.

## Data Hygiene Rules

1. `dock_history.csv` should remain append-oriented.
2. Do not renumber existing `Sxxxxx` IDs.
3. Failed docking should be represented by blank `affinity_kcal_mol`, not by
   trusting status text alone.
4. Chemically suspicious SMILES should be reported before deletion.
5. Active-learning official recommendations should stay small: default top 3.
6. Generated/refined molecules should be deduplicated by canonical SMILES before
   docking.

Sanity report:

```bash
python /mnt/e/Codex/skills/pose-analyzer/scripts/clean_history_sanity.py \
  --history-csv ~/vina_task2/dock_history/dock_history.csv \
  --report-csv ~/vina_task2/dock_history/smiles_sanity_rejects.csv
```

Apply cleanup only after reviewing the report:

```bash
python /mnt/e/Codex/skills/pose-analyzer/scripts/clean_history_sanity.py \
  --history-csv ~/vina_task2/dock_history/dock_history.csv \
  --report-csv ~/vina_task2/dock_history/smiles_sanity_rejects.csv \
  --apply
```

The cleanup script creates a timestamped backup before rewriting history.

## Model Interpretation

Do not treat any model output as proof of activity.

- Vina affinity is a docking heuristic.
- Random Forest regression is diagnostic.
- SHAP explains the RF model, not physical truth.
- The high/low classifier is currently more useful than regression for choosing
  candidates because it asks a cleaner question: does this molecule structurally
  resemble the high-score group?
- The tiny graph-GCN is a small-data neural baseline. It is not a transformer and
  should not be trusted without more official scores.

When official scored molecules reach the hundreds, consider larger neural models
such as ChemBERTa, Graphormer, directed MPNN, or a pretrained molecular
transformer. Before that, Morgan fingerprints plus graph baselines are usually
more honest.

## Files To Share With A Collaborator

Minimum:

```text
/mnt/e/Codex/skills/smiles-to-vina-docking
/mnt/e/Codex/skills/druglike-pocket-refiner
/mnt/e/Codex/skills/pose-analyzer
/mnt/e/Codex/skills/active-learning-orchestrator
E:/Codex/drug_discovery_skills_README.md
```

Project data, if shareable:

```text
~/vina_task2/dock_history/dock_history.csv
~/vina_task2/dock_history/binding_score_analysis.csv
~/vina_task2/references/chembl_kinase_smiles.csv
~/vina_task2/vina_bin/target.pdbqt
~/vina_task2/vina_bin/idz6F_config.txt
```

Be careful with receptor files, official scores, unpublished molecules, and any
competition/private data.
