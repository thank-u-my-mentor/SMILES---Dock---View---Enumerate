---
name: Phylogenetic-Tree
description: Build a general protein FASTA-to-homolog/SSN/phylogenetic-tree/iTOL workflow from curated seed FASTA, project CSVs, or accession ID lists, without assuming hydrolase biology. Use for UniProt/NCBI ID-to-FASTA/PDB enzyme-library dashboards, identity heatmaps/matrices, HMMER/phmmer homolog retrieval, sequence QC, representative selection, tree construction, offline kingdom-level annotation from existing metadata, and post-tree manual SoluProt/NetSolP CSV import into iTOL annotation files. Defaults support /mnt/e/Tree-Metal-F and Metal-F fluorinated-ligand seed sets.
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

## iTOL Visual Style

Kingdom annotations are a stable visual convention inherited from the Metal-F project. Treat this as default behavior for every tree unless the user explicitly requests a different palette:

- `Bacteria`: `#E8A0B0`
- `Plant`: `#A8D5BA`
- `Animal`: `#8FB8E6`
- `Fungi`: `#59A14F`
- `Archaea`: `#F4C2A1`
- `Protist`: `#F9E79F`
- `Metagenome`: `#D5D5D5`
- `Eukaryota`: `#EDC948`
- `Unknown`: `#AAAAAA`

Default kingdom strip styling is `STRIP_WIDTH 25`, `MARGIN 5`, `BORDER_WIDTH 1`, `BORDER_COLOR #000000`. The scripts expose `--kingdom-strip-width`, `--kingdom-margin`, `--kingdom-border-width`, and `--kingdom-border-color` for small visual adjustments, but the defaults should remain Metal-F-compatible.

When adding project-specific categorical strips such as EC class, enzyme subtype, or substrate class, do not reuse the kingdom visual language. Make those strips visually distinct from kingdom by default: narrower strip width, smaller margin, light or white border, and a high-contrast publication-style categorical palette. Avoid assigning adjacent warm hues to rare or conceptually similar categories, reserve neutral gray for `other`, and give compound labels such as `4.1.3.39/43` a visually distant hue from either parent category. iTOL color-strip datasets do not reliably support dashed borders, so prefer width, margin, and border contrast over unsupported dash-like options.

## Identity Heatmap HTML Style

For candidate identity matrices and enzyme-panel heatmaps, follow the HOA heatmap style established by the user's `HOA_2型醛缩酶热图.html` example unless the user explicitly asks for a different design. This is now the default visual grammar for `identity_matrix.html`, candidate comparison heatmaps, and dashboard heatmap tabs:

- Use a quiet paper-like page: `body` background `#f6f6f2`, white panels, ink `#202326`, muted text `#687076`, line `#d8ddd6`, radius `8px`, and subtle shadow `0 12px 30px rgba(33,37,41,.08)`.
- Use `Arial,"Microsoft YaHei",sans-serif`; keep letter spacing at `0`.
- Put controls in a sticky top header: search box, order selector, label selector, identity-threshold selector, pair-mode selector, and precision selector when data are available.
- Use a two-column layout on desktop: heatmap table at left and a `380px` right-side panel for `Selected Pair`, `Color Scale`, and top similar/distant pairs. Collapse to one column on narrow screens.
- Heatmap cells should be compact and stable: `20px` square cells, sticky row labels, vertical sticky column labels, diagonal cells outlined, hover outline, and hidden cell text by default unless precision display is enabled.
- Always use the full 0-100% identity color ramp, not an observed-range or 60-100% rescale:
  - `0% #C94C4C`
  - `25% #ECA76A`
  - `50% #F2E7A6`
  - `75% #9BCB9C`
  - `90% #63B6C2`
  - `100% #3F77B5`
- Show color-scale ticks at exactly `0% / 25% / 50% / 75% / 90% / 100%`.
- Compute cell colors by linear interpolation between those stops. Never map both low and mid-low identity to similar pinks; 10%, 60%, 80%, and 95% must be visually distinct.
- Include pair detail fields when metadata are present: accession pair, identity, EC/subtype/source label, organisms, PDB, DOI, and free-text description.
- Include top similar and most distant pair lists with mini bars colored by the same identity ramp.
- Keep the heatmap self-contained in a single HTML file with embedded JSON data so it can be opened directly from disk.

## Accession List Dashboard

When the user has only UniProt accessions or NCBI protein accessions and wants a small enzyme library HTML, use `scripts/build_enzyme_library_dashboard.py`. It resolves FASTA from local CSV/FASTA first, optionally fetches missing UniProt/NCBI public records, pulls direct UniProt PDB cross-references, assigns nearest in-library PDB templates by pairwise identity, and writes a dashboard with overview, FASTA, and identity heatmap tabs.

For expression-candidate dashboards modeled after the user's Metal-F, HOA, HpcH, aldolase, and VHPO workflows, the default dashboard content should include:

- Protein FASTA for every selected candidate.
- Source/original DNA FASTA when it is available from DOCX, GenBank/EMBL, UniProt cross-references, or local sidecar files.
- A clearly labeled E. coli-optimized draft CDS for every protein when the user is preparing plasmid-ordering tables. If this is generated by a simple codon heuristic, label it as heuristic and advise supplier/manual review.
- PDB IDs, DOI/PMID/literature notes, organism, kingdom/source group, construct decision, and expression caution notes when metadata exist.
- Identity context for triage: nearest core sequence, nearest-core identity, closest selected neighbor, and closest-selected identity.
- Links to the standalone `identity_matrix.html` and tree/iTOL files from the dashboard when practical.

For a fluorinase/SAM halogenase candidate set:

```bash
python /mnt/e/Codex/skills/Phylogenetic-Tree/scripts/build_enzyme_library_dashboard.py build \
  --ids "Q70GK9 W0W999 WP_014985135.1" \
  --outdir /mnt/e/Tree-Metal-F/enzyme_library \
  --title "Natural Fluorine Enzyme Library" \
  --motifs "GTTDDS APNNGLL FADAG" \
  --fetch
```

Run the fast structural check instead of browser-based QA unless layout debugging is needed:

```bash
python /mnt/e/Codex/skills/Phylogenetic-Tree/scripts/build_enzyme_library_dashboard.py check \
  --outdir /mnt/e/Tree-Metal-F/enzyme_library
```

Use `--fetch` only when network access is acceptable. Without it, the script still builds from local metadata and FASTA, marking unresolved records in the output.

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
