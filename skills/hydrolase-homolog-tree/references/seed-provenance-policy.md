# Seed Provenance Policy

Use this reference when enriching `base.xlsx` rows or curated seed FASTA records with expression evidence.

## Dashboard-First Review

Follow the `literature-agent` pattern:

1. Start from a seed table, not from free-form notes.
2. Keep DOI, UniProt, organism, enzyme name, and user notes intact.
3. Add evidence fields in separate columns.
4. Build an HTML dashboard for human review.
5. Do not hide uncertainty. Missing plasmid/tag data is useful information.

Generate a review table from the user's workbook:

```bash
conda activate md
python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/extract_seed_metadata_from_xlsx.py \
  /mnt/e/alpha_beta_hydrolase/base.xlsx \
  --output /mnt/e/alpha_beta_hydrolase/seed_provenance_template_v2.csv
```

After enrichment, generate the dashboard:

```bash
python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/build_seed_review_dashboard.py \
  /mnt/e/alpha_beta_hydrolase/seed_provenance_template_v2.csv \
  --output /mnt/e/alpha_beta_hydrolase/seed_review_dashboard_v2.html
```

The extractor defaults to stable English column names to avoid Chinese-header encoding problems. It preserves Chinese values from the workbook. Use `--header-style original` only when the user explicitly wants original headers.

## Evidence Rules

Record expression metadata only when there is traceable evidence.

- `high`: explicit methods/supplementary statement gives vector, host/strain, and tag or purification construct.
- `medium`: methods imply the construct from a named vector system or deposited plasmid map, but one field is not explicitly stated.
- `low`: inferred from product/manual/common vector behavior, secondary summaries, or incomplete supporting data.
- blank/unknown: no reliable evidence.

Always preserve an `evidence_quote` or concise paraphrase and an `evidence_source` such as DOI, supplement URL, plasmid repository URL, or paper section.

## Fields To Extract

- `paper_title`
- `pubmed_id`
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
- `notes_for_cloning`

## Search Strategy

Use the DOI from `base.xlsx` first. Check, in order:

1. Paper methods and supporting information.
2. UniProt publication links and cross-references.
3. PubMed/Crossref/OpenAlex metadata to resolve title and PMID.
4. Plasmid repositories or supplementary plasmid maps when the paper cites a construct.
5. Full text/PDF only when accessible and needed for methods details.

When `key_doi` is blank, query UniProt first:

```bash
conda activate md
python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/uniprot_seed_literature.py \
  /mnt/e/alpha_beta_hydrolase/seed_provenance_template_v2.csv \
  --outdir /mnt/e/alpha_beta_hydrolase/seed_literature_uniprot
```

Review UniProt publication candidates and PDB cross-references before broader web search. If PDB is available, use the RCSB primary citation and structure methods as a route to expression construct evidence.

Do not assume pET-28 means N-terminal His tag unless the construct or vector map makes that clear. pET vectors often support multiple tag configurations, and the cloned insert can change tag exposure.

## Practical Ranking For Wet Lab Seeds

The user's ideal expression pattern is:

- `E. coli BL21(DE3)`
- `pET-23` or `pET-28` family
- explicit N-terminal or C-terminal His-tag information
- soluble expression and purification reported

Rank seeds higher when those fields are explicitly supported. Do not discard otherwise strong enzymes just because vector/tag is unknown; mark them for manual verification.
