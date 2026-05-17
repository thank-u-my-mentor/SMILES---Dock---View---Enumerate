---
name: literature-agent
description: Build a literature knowledge network from Zotero CSV exports, especially for photoenzymatic radical chemistry, photoredox radical reactions, broader radical chemical reactions, and optional lipase/alpha-beta-hydrolase relevance screening. Use when the user wants to import Zotero CSV files, tune LLM extraction schema/temperature/max_tokens, handle Zotero PDF attachment paths, or generate Markdown reports and network visualizations from literature metadata.
---

# Literature Agent

Use this skill when working from Zotero CSV exports to create a literature knowledge base: structured extraction, concept network, clustered paper cards, gap analysis, and plots.

Primary script:

```bash
/mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py
```

Preferred runtime:

```bash
conda activate md
```

## Default Workflow

Prefer WSL/bash commands with line continuations. Put outputs under `/mnt/e`, for example:

```bash
mkdir -p /mnt/e/literature_knowledge

python /mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py \
  --input /mnt/c/Users/Administrator/Downloads/导出的条目.csv \
  --output /mnt/e/literature_knowledge/zotero_run \
  --schema-profile photo-radical
```

This no-API mode uses local keyword fallback. It is useful for quick smoke tests, checking CSV encoding, and confirming plots/reports can be generated.

## LLM Extraction

For real knowledge-network extraction, use the LLM path:

```bash
export KIMI_API_KEY="..."

python /mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py \
  --input /mnt/c/Users/Administrator/Downloads/导出的条目.csv \
  --output /mnt/e/literature_knowledge/zotero_llm_run \
  --use-llm \
  --schema-profile photo-radical \
  --temperature 0.1 \
  --max-tokens 2500
```

Do not hard-code API keys into the script. Use `KIMI_API_KEY`, `MOONSHOT_API_KEY`, or `--api-key`.

If Kimi returns HTTP 400, first try:

```bash
python /mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py \
  --input /mnt/c/Users/Administrator/Downloads/导出的条目.csv \
  --output /mnt/e/literature_knowledge/zotero_llm_run \
  --use-llm \
  --schema-profile photo-radical \
  --model kimi-k2.6 \
  --no-json-mode
```

The script defaults to `https://api.moonshot.cn/v1`, which matches many China-region Kimi keys. Use `--base-url https://api.moonshot.ai/v1` only if the key is from that endpoint. The script prints the first part of the API error body, which is usually needed to distinguish an unsupported model, unsupported `response_format`, wrong base URL, or token-limit issue.

Recommended extraction settings:

- `temperature 0.0-0.2` for stable structured extraction.
- `max_tokens 2500` for the built-in `photo-radical` schema.
- Increase to `3000-4000` if adding PDF-derived text chunks or longer evidence fields.

## Inputs

The script expects a Zotero CSV export. It reads common Zotero columns such as:

- `Title`
- `Author`
- `Publication Year`
- `DOI`
- `Abstract Note`
- `Publication Title`
- `Url`
- `File Attachments`
- `Link Attachments`

CSV encoding is auto-detected across `utf-8-sig`, `utf-8`, `utf-16`, `utf-16-le`, and `gb18030`. If auto-detection fails, pass:

```bash
--csv-encoding utf-16
```

Chinese paths in Zotero attachment fields are acceptable in Python 3. In WSL, Windows paths may need conversion if they are later used for reading PDFs:

```text
D:\文献库\storage\xxx.pdf
/mnt/d/文献库/storage/xxx.pdf
```

Currently the script preserves attachment paths but does not yet extract full PDF text.

## Schema Profiles

Use `--schema-profile photo-radical` by default for literature discovery. It emphasizes:

- photoenzymatic radical chemistry
- photoredox radical reactions
- broader radical chemical reactions that may not yet be enzymatic
- enzyme-integration opportunity as a manual review aid
- lipase/esterase and alpha/beta-hydrolase relevance only as an optional screening field
- photoenzymatic catalysis
- photoredox/radical mechanisms
- substrate/product classes
- mechanistic evidence
- confidence score

Use `--schema-profile lipase` only when the paper set is already centered on lipase/esterase or alpha/beta-hydrolase biology. Use `--schema-profile photoenzyme` for the original narrower photoenzymatic schema.

Manual intervention is expected: after the first run, inspect `manual_review_priority`, `enzyme_integration_opportunity`, and `lipase_or_hydrolase_relevance` rather than trusting the network as a final classification.

For custom extraction fields, provide a JSON object:

```bash
python /mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py \
  --input /mnt/c/Users/Administrator/Downloads/导出的条目.csv \
  --output /mnt/e/literature_knowledge/custom_schema_run \
  --use-llm \
  --schema-file /mnt/e/literature_knowledge/my_schema.json
```

The schema file should be:

```json
{
  "field_name": "field description or allowed values",
  "confidence": "0.0-1.0"
}
```

## Outputs

The output directory contains:

- `knowledge_report.md`: main literature synthesis.
- `gap_analysis.md`: mechanism/reaction/enzyme gap notes.
- `cluster_*.md`: paper cards grouped by reaction type.
- `timeline.png`: paper timeline.
- `reaction_enzyme_matrix.png`: reaction vs enzyme heatmap.
- `mechanism_flow.png`: mechanism distribution plot.
- `knowledge_graph.html`: generated only when `pyvis` is installed.

If `pyvis` is missing, the script skips the interactive HTML graph and still creates reports/plots.

## Validation

After editing the script, run:

```bash
python -m py_compile /mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py
```

Then run a no-API smoke test:

```bash
python /mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py \
  --input /mnt/c/Users/Administrator/Downloads/导出的条目.csv \
  --output /mnt/e/literature_knowledge/smoke_test \
  --schema-profile photo-radical
```

If using `knowledge_graph.html`, install `pyvis` in the active `md` environment:

```bash
pip install pyvis
```
