---
name: Phylogenetic-Tree
description: Build a general protein FASTA-to-homolog/SSN/phylogenetic-tree/iTOL workflow from curated seed FASTA or project CSVs, without assuming hydrolase biology. Use for HMMER/phmmer homolog retrieval, sequence QC, representative selection, tree construction, offline kingdom-level annotation from existing metadata, and post-tree manual SoluProt/NetSolP CSV import into iTOL annotation files. Defaults support /mnt/e/Tree-Metal-F and Metal-F fluorinated-ligand seed sets.
---

# Phylogenetic-Tree

General FASTA-to-tree workflow, not hydrolase-specific. Prefer this skill when the user wants phylogenetic tree, HMMER homologs, SSN, iTOL annotations, kingdom labels, or SoluProt/NetSolP import for any enzyme/protein family.

## Defaults

- Working directory: `/mnt/e/Tree-Metal-F`
- Source project for current Metal-F seeds: `/home/qin/Metal-F_project`
- Reuse mature scripts from `/mnt/e/Codex/skills/hydrolase-homolog-tree/scripts` when appropriate.
- Treat `Phylogenetic-Tree` as the new general entry point. `hydrolase-homolog-tree` remains a legacy script library until its generic scripts are migrated or wrapped here.
- Do not run per-organism taxonomy API loops by default. For seed records, infer kingdom from existing organism/source fields and simple organism keywords. For broad homolog sets, prefer UniProt/HMMER metadata already present in FASTA headers or existing SSN nodes. Only query UniProt taxonomy if the user explicitly asks.
- SoluProt/NetSolP is post-tree manual work. First build homologs, SSN representatives, and a tree; then export the 100-200 representative sequences for manual SoluProt.
- Phylum/class annotation is optional and off by default; iTOL tree visualisation only requires kingdom-level grouping unless requested otherwise.
- For homolog retrieval, prefer aligned seed FASTA or HMM profile -> EBI HMMER `hmmsearch` -> UniProt. Default strict E-value is `1e-10`; use `1e-20` for stricter runs. Cap combined FASTA near 1000 sequences before SSN/tree.

## Recommended Metal-F Setup

Prepare seed inputs:

```bash
python /mnt/e/Codex/skills/Phylogenetic-Tree/scripts/prepare_tree_metal_f_inputs.py \
  --metal-project /home/qin/Metal-F_project \
  --outdir /mnt/e/Tree-Metal-F \
  --sets CF3 C-F fluoroaryl \
  --min-length 150
```

This writes:

- `core.fasta`: deduplicated full-length seed FASTA for HMMER.
- `seed_metadata.csv`: PDB/UniProt/organism/kingdom/provenance table.
- `itol_seed_kingdom_colorstrip.txt`: seed-level kingdom iTOL annotation.
The setup script no longer treats seed sequences as the final SoluProt target. SoluProt input should be regenerated after SSN representative selection.

## Homolog/SSN/Tree

For the current Metal-F project, the fastest reproducible entry point is:

```bash
python /mnt/e/Codex/skills/Phylogenetic-Tree/scripts/run_metal_f_tree_workflow.py \
  --metal-project /home/qin/Metal-F_project \
  --outdir /mnt/e/Tree-Metal-F \
  --min-length 120 \
  --mode core-tree \
  --threads 8
```

This broadens the core beyond the strict geometry-passing set by merging local `Metal-F_protein_entities.csv`, `Metal-F_protein_entities.fasta`, `Metal-*/*_protein_entities_kept.csv`, and existing `Metal-*/*_hmmer_full_length_seed.fasta` files. It then builds a quick ClustalO/IQ-TREE tree, writes offline kingdom iTOL annotation, and prepares post-tree SoluProt input.

For a fuller homolog-expanded tree, use the same wrapper with:

```bash
python /mnt/e/Codex/skills/Phylogenetic-Tree/scripts/run_metal_f_tree_workflow.py \
  --metal-project /home/qin/Metal-F_project \
  --outdir /mnt/e/Tree-Metal-F \
  --mode hmmer-tree \
  --hmmer-mode ebi-hmmsearch \
  --hmmer-database uniprot \
  --hmmer-evalue 1e-10 \
  --max-hits 800 \
  --max-total 1000 \
  --max-rep 200
```

The HMMER step uses `fetch_hmmer_hmmsearch.py`: it aligns `tree_core.fasta` with ClustalO unless an aligned input is supplied, submits EBI HMMER `hmmsearch`, downloads UniProt FASTA hits, and writes `hmmer/homologs_plus_core.fasta`. The user-mentioned `1E10` to `1E20` cutoff should be interpreted as the strict E-value range `1e-10` to `1e-20`.

When HMMER homolog FASTA is available, or after fetching homologs, reuse the existing general scripts:

```bash
python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/construct_sequence_qc.py \
  /mnt/e/Tree-Metal-F/core.fasta \
  --sequence-type protein \
  --output /mnt/e/Tree-Metal-F/core.construct_qc.csv

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/ssn_pipeline_mmseq2.py \
  /mnt/e/Tree-Metal-F/all.fasta \
  --core-fasta /mnt/e/Tree-Metal-F/core.fasta \
  --output /mnt/e/Tree-Metal-F/ssn \
  --max-rep 2500 \
  --cdhit-identity 0.35 \
  --cdhit-coverage 0.75 \
  --max-recon 4000

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/tree_pipeline.py \
  /mnt/e/Tree-Metal-F/ssn/representatives.fasta \
  --output /mnt/e/Tree-Metal-F/tree_analysis
```

If the user wants no taxonomy API calls, avoid taxonomy enrichment flags and keep kingdom annotation from `seed_metadata.csv` or existing node fields.

The legacy SSN/tree scripts can contain Kimi/API fallback code for unknown taxonomy. For a strict offline workflow, do not set `KIMI_API_KEY` in the command environment and run the offline kingdom annotator after SSN/tree:

```bash
python /mnt/e/Codex/skills/Phylogenetic-Tree/scripts/offline_kingdom_itol.py \
  --nodes /mnt/e/Tree-Metal-F/ssn/ssn_nodes.csv \
  --fasta /mnt/e/Tree-Metal-F/ssn/representatives.fasta \
  --seed-metadata /mnt/e/Tree-Metal-F/seed_metadata.csv \
  --outdir /mnt/e/Tree-Metal-F/tree_analysis
```

This writes `nodes_with_offline_kingdom.csv` and `itol_kingdom_color_strip_offline.txt`. Unknown kingdoms are acceptable and should be left for manual review or a local taxdump/metadata table, not fixed with hidden API calls.

## Post-Tree SoluProt Input

After SSN representative selection and tree construction, generate the manual SoluProt FASTA from the representative sequences:

```bash
python /mnt/e/Codex/skills/Phylogenetic-Tree/scripts/prepare_soluprot_tree_inputs.py \
  /mnt/e/Tree-Metal-F/ssn/representatives.fasta \
  --outdir /mnt/e/Tree-Metal-F/soluprot_inputs
```

This writes the 100-200 tree representative sequences to `soluprot_inputs/representatives_for_soluprot.fasta` plus `soluprot_id_mapping.tsv`.

## SoluProt Manual Import

The user manually runs SoluProt/NetSolP on the post-tree representative FASTA and saves the result as:

```text
/mnt/e/Tree-Metal-F/soluprot_inputs/input.csv
```

Then generate iTOL annotation files:

```bash
python /mnt/e/Codex/skills/Phylogenetic-Tree/scripts/make_soluprot_annotation.py \
  --input /mnt/e/Tree-Metal-F/soluprot_inputs/input.csv \
  --outdir /mnt/e/Tree-Metal-F/tree_analysis \
  --nodes /mnt/e/Tree-Metal-F/ssn/ssn_nodes.csv \
  --id-mapping /mnt/e/Tree-Metal-F/soluprot_inputs/soluprot_id_mapping.tsv
```

This wraps the hydrolase skill's flexible `annotate_solubility_itol.py` and writes upload-ready iTOL annotation files.

## Output Discipline

Keep generated files under `/mnt/e/Tree-Metal-F`. Do not scatter intermediate FASTA, CSV, or iTOL files in the WSL home directory.
