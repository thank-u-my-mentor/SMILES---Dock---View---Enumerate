---
name: hydrolase-homolog-tree
description: Build a homolog discovery, sequence similarity network, phylogenetic tree, seed literature provenance, expression-vector/tag review, optional solubility annotation, and iTOL workflow from a small curated protein FASTA or base.xlsx seed table, especially 5-25 literature-validated alpha/beta hydrolase or related enzyme sequences. Use when Codex needs to run or adapt HMMER/phmmer or local hmmsearch homolog retrieval, keep the homolog set bounded to hundreds or a few thousand proteins, run MMseqs2-based SSN/representative selection, extract DOI/plasmid/His-tag/E. coli BL21(DE3) evidence for seed enzymes from papers/supplements, optionally import user-provided NetSolP/SoluProt-like predictions, and generate protected phylogenetic tree visualizations plus iTOL annotation files from the bundled tree_pipeline.py.
---

# Hydrolase Homolog Tree

Local skill path on this machine:

```bash
/mnt/e/Codex/skills/hydrolase-homolog-tree
```

Windows equivalent:

```text
E:\Codex\skills\hydrolase-homolog-tree
```

Preferred runtime:

```bash
conda activate md
```

## Workflow

Use this skill for small seed FASTA sets of experimentally characterized enzymes, then expand, reduce, and visualize homologs:

1. Start from the curated seed FASTA, usually 5-25 proteins.
2. If `base.xlsx` is provided, extract a seed provenance table and enrich it with paper/vector/tag evidence.
3. If DOI is missing, query UniProt for publication and PDB candidates before using broader literature discovery.
4. Use `literature-agent` style expansion from the seed DOI/UniProt/PDB candidates to find missed enzyme papers.
   For flavin/photoenzymatic discovery projects, start from `literature-agent` output `enzyme_seed_candidates.csv` and, when present, `enzyme_seed_candidates.fasta`: manually confirm enzyme name, UniProt/PDB/sequence evidence, and DOI before adding rows to a seed FASTA or `base.xlsx`. Use `homolog_candidate_seeds.csv` as the paper-level evidence table.
5. Fetch bounded homologs with `scripts/fetch_hmmer_homologs.py`, or use the user's HMMER web output FASTA directly.
6. Run the inherited SSN/representative-selection script `scripts/ssn_pipeline_mmseq2.py`.
7. Run the inherited tree/iTOL script `scripts/tree_pipeline.py`.
8. If the curated `base.xlsx` contains manual curation flags, generate the consolidated annotation layer.
9. Optionally import the user's manually downloaded NetSolP/SoluProt CSV into that same annotation layer; do not run this by default for literature-validated seed enzymes.
10. Before cloning or ordering constructs, run sequence QC for unknown residues, ambiguous DNA bases, internal stops, and terminal secretion/anchor-like regions.
11. Return the output paths for `representatives.fasta`, `ssn_nodes.csv`, `.xgmml`, `.newick`, iTOL text files, QC CSV, and seed review dashboard.

## Daily Commands

Prepare seed provenance review from the user's current workbook:

```bash
conda activate md

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/extract_seed_metadata_from_xlsx.py \
  /mnt/e/alpha_beta_hydrolase/base.xlsx \
  --output /mnt/e/alpha_beta_hydrolase/seed_provenance_template_v2.csv

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/build_seed_review_dashboard.py \
  /mnt/e/alpha_beta_hydrolase/seed_provenance_template_v2.csv \
  --output /mnt/e/alpha_beta_hydrolase/seed_review_dashboard_v2.html
```

Fill missing DOI/PubMed/PDB candidates from UniProt before manual curation:

```bash
conda activate md

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/uniprot_seed_literature.py \
  /mnt/e/alpha_beta_hydrolase/seed_provenance_template_v2.csv \
  --outdir /mnt/e/alpha_beta_hydrolase/seed_literature_uniprot
```

Convert the enriched seed table into a literature-agent input CSV:

```bash
python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/base_to_literature_seed_csv.py \
  /mnt/e/alpha_beta_hydrolase/seed_literature_uniprot/seed_provenance_uniprot_enriched.csv \
  --output /mnt/e/alpha_beta_hydrolase/hydrolase_literature_seed.csv
```

Then use `$literature-agent` on that CSV to perform broader DOI/citation discovery before finalizing `base.xlsx`.

For flavin/photoenzymatic candidate discovery coming from `$literature-agent`, review:

```bash
/mnt/e/literature_flavin_photoenzyme/enzyme_seed_candidates.csv
/mnt/e/literature_flavin_photoenzyme/enzyme_seed_candidates.fasta
/mnt/e/literature_flavin_photoenzyme/homolog_candidate_seeds.csv
```

Use the CSV as the manual bridge table and the FASTA only as a starting point, not as an automatically trusted core set. Confirm a real enzyme accession, sequence, or PDB entry for each selected row, then add the chosen enzymes to `base.xlsx` and the curated seed FASTA before running SSN/tree expansion.

For reproducing the user's current known local result from existing HMMER output, prefer the two-step local path:

```bash
conda activate md

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/ssn_pipeline_mmseq2.py \
  /mnt/e/alpha_beta_hydrolase/all.fasta \
  --core-fasta /mnt/e/alpha_beta_hydrolase/core.fasta \
  --output /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn \
  --max-rep 2500 \
  --cdhit-identity 0.35 \
  --cdhit-coverage 0.75 \
  --max-recon 4000

cd /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/tree_pipeline.py \
  representatives.fasta \
  --output tree_analysis
```

If the user has not already produced homolog FASTA from HMMER, `scripts/run_hydrolase_pipeline.py` can run HMMER retrieval first. Use the real seed file path, usually `/mnt/e/alpha_beta_hydrolase/core.fasta`, not a placeholder `seed.fasta`.

## Inputs

- `seed.fasta`: curated literature-supported enzyme FASTA. Treat it as the core set.
- `base.xlsx`: optional seed provenance workbook. Current known columns include family, enzyme name, source organism, known promiscuity, UniProt ID, notes, and key DOI.
- Optional local protein database: use `--db-fasta proteins.fa` to run local `hmmbuild`/`hmmsearch` instead of web HMMER.
- Optional solubility CSV from NetSolP/SoluProt or a manually normalized table. Treat this as optional downstream annotation, not a required seed-quality step.
- Optional `KIMI_API_KEY`: enables Kimi taxonomy rescue in the inherited scripts. If unset, UniProt lineage still works and unresolved taxonomy remains `Unknown`.

Keep homolog retrieval bounded. Defaults are intentionally conservative: `--max-hits 2500`, `--max-total 3000`, and `--max-rep 2500`.

## Seed Provenance

Follow the `literature-agent` style: evidence table first, dashboard for review, no silent guesses.

Use `scripts/extract_seed_metadata_from_xlsx.py` to create `seed_provenance_template_v2.csv` from `base.xlsx`. It defaults to stable English headers (`family`, `enzyme_name`, `source_organism`, `known_promiscuity`, `uniprot_id`, `seed_notes`, `key_doi`) while preserving Chinese text in values. Enrich the added columns by checking the seed DOI, paper methods, supporting information, UniProt publication links, PubMed/Crossref/OpenAlex metadata, and plasmid repositories when relevant.

Prioritize these fields:

- expression host and strain, especially `E. coli BL21(DE3)`
- plasmid vector, especially pET-23 or pET-28 family
- His-tag type and N/C-terminal position
- whether tag was cleaved
- soluble expression/purification evidence
- concise evidence quote/paraphrase, source, and confidence

Use `references/seed-provenance-policy.md` before doing paper/vector extraction. It defines `high`, `medium`, and `low` confidence and the user's wet-lab ranking preference.

When key DOI is blank, do not leave the row there. Use `scripts/uniprot_seed_literature.py` to query UniProtKB by accession. Review:

- `seed_literature_uniprot/uniprot_publication_candidates.csv`
- `seed_literature_uniprot/uniprot_pdb_candidates.csv`
- `seed_literature_uniprot/seed_provenance_uniprot_enriched.csv`

If UniProt points to PDB structures but not a clear expression paper, follow the PDB entry to its primary citation and methods/supplement.

After the tree is generated, check whether non-seed tree leaves include missed reported enzymes:

```bash
conda activate md

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/triage_tree_uniprot_literature.py \
  --nodes /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/ssn_nodes.csv \
  --base /mnt/e/alpha_beta_hydrolase/seed_provenance_template_v2.csv \
  --outdir /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/tree_literature_triage \
  --workers 6 \
  --sleep 0.05
```

Review `tree_uniprot_literature_triage.csv` first. `high` means UniProt publication titles contain enzyme/activity/structure-like terms; `skip` usually means only genome/proteome-like publications were detected. This is a fast triage, not a substitute for reading methods/results.

## Consolidated Annotations

After the tree exists, use one command to build the optional iTOL annotations from `base.xlsx` and, when available, the manually downloaded SoluProt/NetSolP CSV:

```bash
conda activate md

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/build_tree_annotations.py \
  --tree-dir /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/tree_analysis \
  --nodes /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/ssn_nodes.csv \
  --base-xlsx /mnt/e/alpha_beta_hydrolase/base.xlsx \
  --soluprot-csv /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/soluprot_inputs/input.csv \
  --id-mapping /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/soluprot_inputs/soluprot_id_mapping.tsv \
  --low-threshold 0.2 \
  --tool-label SoluProt
```

This keeps user-facing iTOL upload files at the same level as `itol_core_highlight.txt`:

- `itol_extra_unvalidated_stars.txt`: yellow-star track for unvalidated manual candidates from `base.xlsx`.
- `itol_core_short_name_text.txt`: short enzyme labels from `base.xlsx`.
- `itol_soluprot_gradient_symbols.txt`: optional SoluProt red-yellow-green gradient dots.

It keeps audit/intermediate files under `tree_analysis/annotations/`:

- `annotations/base_xlsx/base_itol_annotation_summary.csv`
- `annotations/soluprot/solubility_normalized.csv`
- `annotations/soluprot/ssn_nodes_with_solubility.csv`

## SoluProt Manual Step

NetSolP and SoluProt should be treated as optional downstream annotation sources unless the user chooses one. For seed enzymes already reported as expressed/purified in literature, do not assume a solubility problem and do not run solubility prediction by default. Use this layer for representative homologs, new candidates, or when the user provides SoluProt/NetSolP output.

SoluProt is useful as a quick sequence-level screen for E. coli expression/solubility, but it is modestly predictive; treat it as ranking/annotation, not experimental truth.

Preferred pattern:

1. Prepare a clean FASTA for every representative sequence used in the tree:

```bash
conda activate md

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/prepare_soluprot_tree_inputs.py \
  /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/representatives.fasta \
  --outdir /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/soluprot_inputs
```

This writes `representatives_for_soluprot.fasta` for all tree leaves and `soluprot_id_mapping.tsv`.

2. The user manually runs the SoluProt or NetSolP website with `representatives_for_soluprot.fasta`, then downloads the result CSV as:

```bash
/mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/soluprot_inputs/input.csv
```

The CSV must contain an ID column and a score/probability column. Accepted column names include `id`, `fa_id`, `sequence_id`, `accession`, `soluprot_probability`, `probability`, `score`, `soluble`, `prediction`, or `call`. Missing rows are allowed; unmatched rows are ignored and written to an audit file when present.

3. Re-run `scripts/build_tree_annotations.py` with `--soluprot-csv .../input.csv`. This is the only non-automatic step in the standard skill flow: Codex prepares the FASTA and imports the CSV, while the user obtains the website prediction CSV.

## Construct QC

Before moving tree candidates into Benchling/plasmid design, run a lightweight FASTA QC pass:

```bash
conda activate md

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/construct_sequence_qc.py \
  /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/representatives.fasta \
  --sequence-type protein \
  --output /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/construct_sequence_qc.csv
```

For protein FASTA, the script flags:

- `X`, `B`, `Z`, `J`, `U`, `O`, or `*` residues.
- internal stop symbols.
- N-terminal signal-peptide-like hydrophobic segments.
- C-terminal hydrophobic-anchor-like segments.

For DNA exported from Benchling, run:

```bash
python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/construct_sequence_qc.py \
  benchling_orfs.fasta \
  --sequence-type dna \
  --dna-coding \
  --output benchling_orfs.construct_qc.csv
```

For DNA, the script flags invalid characters such as `X`, ambiguous IUPAC bases such as `N/R/Y`, coding lengths not divisible by 3, and internal stop codons. This check does not replace SnapGene/Benchling visual inspection or SignalP/DeepTMHMM, but it catches common design-blocking mistakes early.

## Label Simplification

Many UniProt-derived FASTA headers contain duplicated IDs such as `A0ABV6D1V5|A0ABV6D1V5_9SPHN`. Use the simplifier after tree generation to create a shorter Newick and matching iTOL files:

```bash
conda activate md

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/simplify_itol_tree_labels.py \
  /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/tree_analysis/phylogenetic_tree.newick \
  --itol \
  /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/tree_analysis/itol_kingdom_color_strip.txt \
  /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/tree_analysis/itol_phylum_label.txt \
  /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/tree_analysis/itol_core_highlight.txt
```

Upload files from `tree_analysis/simplified_labels/` together. Do not mix simplified Newick with unsimplified iTOL annotations.

## Base.xlsx iTOL Curation Annotation

When the user updates `base.xlsx` with manual candidate flags, convert those rows into separate iTOL annotations instead of changing the protected tree pipeline.

Current rules:

- `备注` containing `漏补候选`: merged directly into `itol_core_highlight.txt` and treated like the other core sequences.
- `备注` containing `额外候选`: large yellow star with a black border in `itol_extra_unvalidated_stars.txt`.
- All non-extra core/curated rows get a text label from `酶名称`.
- If `酶名称` contains half-width or full-width parentheses, use only the last parenthesized short name. For example, `Ylehd Epoxide hydrolase (YlEH)` becomes `YlEH`.

Run after SSN/tree generation:

```bash
conda activate md

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/annotate_base_xlsx_itol.py \
  /mnt/e/alpha_beta_hydrolase/base.xlsx \
  --nodes /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/ssn_nodes.csv \
  --outdir /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/tree_analysis/annotations/base_xlsx
```

This writes:

- `itol_extra_unvalidated_stars.txt`: extra unvalidated candidate yellow-star symbols.
- `itol_core_short_name_text.txt`: core enzyme short-name labels.
- `base_itol_annotation_summary.csv`: matched/unmatched audit table.

## Running Pieces Manually

Fetch homologs:

```bash
conda activate md
python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/fetch_hmmer_homologs.py seed.fasta --outdir 01_hmmer --mode web-phmmer --max-hits 2500 --max-total 3000
```

or with a local FASTA database:

```bash
conda activate md
python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/fetch_hmmer_homologs.py seed.fasta --outdir 01_hmmer --mode local-hmmsearch --db-fasta proteins.fa --max-hits 2500
```

Run SSN/representative selection:

```bash
conda activate md
python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/ssn_pipeline_mmseq2.py 01_hmmer/homologs_plus_core.fasta --core-fasta 01_hmmer/core.fasta --output 02_ssn --max-rep 2500 --cdhit-identity 0.35 --cdhit-coverage 0.75 --max-recon 3000
```

Run tree/iTOL generation:

```bash
conda activate md
cd 02_ssn
python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/tree_pipeline.py representatives.fasta --output tree_analysis
```

Important: run `tree_pipeline.py` with current working directory set to the SSN output directory, because it intentionally reads `ssn_nodes.csv` from the current directory.

## Protected Code

Treat `scripts/tree_pipeline.py` as protected user-optimized visualization code. Do not casually rewrite:

- `export_itol_files`
- Kingdom pastel color strip settings
- phylum/class label output
- core sequence binary highlight output
- the `ssn_nodes.csv` metadata handoff

Only patch these sections when the user explicitly asks for a visual/annotation change or when a minimal compatibility fix is required.

The bundled `scripts/ssn_pipeline_mmseq2.py` is also inherited project logic. Prefer parameter changes over algorithm rewrites unless the user asks to alter clustering, taxonomy rescue, or SSN edge construction.

## Outputs

Expected end-to-end output layout:

- `01_hmmer/core.fasta`
- `01_hmmer/homologs_plus_core.fasta`
- `01_hmmer/hmmer_hits.csv`
- `02_ssn/representatives.fasta`
- `02_ssn/ssn_nodes.csv`
- `02_ssn/ssn_edges.csv`
- `02_ssn/ssn_network.xgmml`
- optional `02_ssn/construct_sequence_qc.csv`
- `02_ssn/tree_analysis/phylogenetic_tree.newick`
- `02_ssn/tree_analysis/itol_kingdom_color_strip.txt`
- `02_ssn/tree_analysis/itol_phylum_label.txt`
- `02_ssn/tree_analysis/itol_core_highlight.txt`
- optional `02_ssn/tree_analysis/itol_extra_unvalidated_stars.txt`
- optional `02_ssn/tree_analysis/itol_core_short_name_text.txt`
- optional `02_ssn/tree_analysis/itol_soluprot_gradient_symbols.txt`
- optional `02_ssn/tree_analysis/annotations/`
- optional `seed_review_dashboard_v2.html`
- optional `seed_provenance_template_v2.csv`

Upload the `.newick` plus the iTOL text files to iTOL for the annotated interactive tree.

## HMMER Web Versus Terminal

The user often uses HMMER through the website. This skill supports both routes.

- If a web HMMER homolog FASTA is already available, use it directly as the SSN input and preserve the seed FASTA as `--core-fasta`.
- If running from terminal, `scripts/fetch_hmmer_homologs.py` uses the EMBL-EBI HMMER web service from Python for web-like phmmer searches, or local HMMER tools with `--db-fasta`.
- Network calls may require approval in restricted environments. Local `conda activate md` is the preferred runtime once dependencies are installed there.

## Troubleshooting

- If web HMMER fails, retry later or switch to `--mode local-hmmsearch --db-fasta proteins.fa`.
- If tree annotations are all `Unknown`, inspect `02_ssn/ssn_nodes.csv` and make sure tree generation ran from `02_ssn` as the current working directory.
- If no tree is produced, check `iqtree` or `iqtree2` availability.
- If SSN construction fails, check `mmseqs` availability.
- If MSA fails, check `clustalo`; local HMM mode can also use `mafft` for the seed alignment.
- If solubility IDs do not match the tree, inspect `unmatched_solubility_rows.csv` and normalize IDs to node IDs or UniProt IDs.

For detailed invariants and expected files, read `references/workflow-contract.md`. For literature/vector/tag evidence handling, read `references/seed-provenance-policy.md`.
