# Hydrolase Homolog Tree Workflow Contract

Use this reference when adapting or troubleshooting the hydrolase homolog-tree skill.

## Protected Script Behavior

- Treat `scripts/tree_pipeline.py` as protected visualization code. Do not rewrite `export_itol_files`, kingdom color strip generation, phylum label generation, or core highlight output unless the user explicitly asks to change the visual style.
- `tree_pipeline.py` expects `ssn_nodes.csv` in the current working directory. Run it with `cwd` set to the SSN output directory, while passing `representatives.fasta` as the input.
- Preserve the curated seed FASTA as `core.fasta` and pass it to the SSN script with `--core-fasta`; this is what drives `is_core` and the iTOL core highlight.
- Keep homolog expansion bounded. Default target is 200-3000 sequences; avoid unbounded UniProt downloads.

## Expected Outputs

HMMER stage:

- `01_hmmer/core.fasta`
- `01_hmmer/homologs.fasta`
- `01_hmmer/homologs_plus_core.fasta`
- `01_hmmer/hmmer_hits.csv`

SSN stage:

- `02_ssn/representatives.fasta`
- `02_ssn/representatives_aln.fasta`
- `02_ssn/ssn_nodes.csv`
- `02_ssn/ssn_edges.csv`
- `02_ssn/ssn_network.xgmml`
- optional `02_ssn/construct_sequence_qc.csv`

Tree/iTOL stage:

- `02_ssn/tree_analysis/alignment.fasta`
- `02_ssn/tree_analysis/phylogenetic_tree.newick`
- `02_ssn/tree_analysis/itol_kingdom_color_strip.txt`
- `02_ssn/tree_analysis/itol_phylum_label.txt`
- `02_ssn/tree_analysis/itol_core_highlight.txt`
- optional `02_ssn/tree_analysis/itol_extra_unvalidated_stars.txt`
- optional `02_ssn/tree_analysis/itol_core_short_name_text.txt`
- optional `02_ssn/tree_analysis/itol_soluprot_gradient_symbols.txt`
- optional `02_ssn/tree_analysis/annotations/base_xlsx/base_itol_annotation_summary.csv`
- optional `02_ssn/tree_analysis/annotations/soluprot/solubility_normalized.csv`
- optional `02_ssn/tree_analysis/annotations/soluprot/ssn_nodes_with_solubility.csv`
- optional `02_ssn/tree_analysis/simplified_labels/simplified_phylogenetic_tree.newick` and matching simplified iTOL files

Seed provenance stage:

- `seed_provenance_template_v2.csv` may be generated from `base.xlsx`; default headers are stable English identifiers to avoid Chinese-header encoding issues.
- `seed_literature_uniprot/uniprot_publication_candidates.csv` and `uniprot_pdb_candidates.csv` may be generated for rows with missing DOI or weak paper evidence.
- Preserve user-supplied DOI/UniProt fields.
- Add expression metadata with explicit evidence, not guesses.

Manual base.xlsx iTOL curation:

- Keep candidate flags outside the protected tree visualizer; use `scripts/build_tree_annotations.py` or `scripts/annotate_base_xlsx_itol.py`.
- `备注` containing `漏补候选` is merged directly into `itol_core_highlight.txt`.
- `备注` containing `额外候选` becomes a large yellow star with a black border in `itol_extra_unvalidated_stars.txt`.
- Core short-name text labels come from `酶名称`; if parentheses are present, use the last half-width or full-width parenthesized term only.
- Upload-ready optional iTOL files should sit directly in `tree_analysis/`; audit/intermediate files should sit under `tree_analysis/annotations/`.

## Dependency Notes

- Web HMMER mode requires network access to EMBL-EBI HMMER and UniProt REST.
- Local HMM search mode requires `clustalo` or `mafft`, plus `hmmbuild` and `hmmsearch`.
- SSN and tree stages require `mmseqs`, `clustalo`, and `iqtree`/`iqtree2` on PATH or at the paths already embedded in the scripts.
- Python packages used by the inherited scripts include `requests`, `networkx`, and optionally `tqdm`.
- NetSolP/SoluProt prediction can be imported as a user-downloaded CSV named `input.csv`. Do not hard-depend on one web service, and do not run solubility prediction by default for literature-validated seed enzymes unless the user asks or the paper text raises solubility/insolubility.
- Prefer iTOL `DATASET_SYMBOL` gradient dots for SoluProt/NetSolP scores in this workflow. The base.xlsx curation track uses `DATASET_BINARY` so it does not collide visually with solubility dots.
- BioLM SoluProt API can be run with `scripts/run_biolm_soluprot.py` when `BIOLM_API_KEY` is available; otherwise use standalone SoluProt or any CSV export with IDs and scores.
- `scripts/construct_sequence_qc.py` has no external dependency. It is a lightweight pre-cloning screen for ambiguous residues/bases, internal stops, signal-peptide-like N termini, and hydrophobic C termini. Confirm secretion/TM predictions with SignalP, DeepTMHMM, Phobius, or a similar tool before final construct design.

## Expression Provenance Fields

Use these fields when enriching seed metadata from papers:

- `expression_host`
- `expression_strain`
- `plasmid_vector`
- `promoter`
- `antibiotic`
- `tag_type`
- `tag_position`
- `tag_cleaved`
- `expression_temperature_c`
- `inducer`
- `soluble_expression_reported`
- `purification_method`
- `evidence_quote`
- `evidence_source`
- `provenance_confidence`

Record `provenance_confidence` as `high` only when the paper or supporting data explicitly states the expression vector/tag/strain. Use `medium` for strongly implied vector maps or methods text, and `low` when only inferred from repository metadata.
