---
name: literature-agent
description: Build and browse a literature knowledge network from Zotero CSV exports or any seed table that can be normalized to DOI/PMID/title/URL records, especially for photoenzymatic radical chemistry, photoredox/radical organic chemistry, enzyme/non-enzyme boundary triage, enzyme seed-table provenance, Semantic Scholar citation discovery, persistent blacklist memory, and HTML dashboards.
---

# Literature Agent

Use this skill when working from Zotero CSV exports or any seed table with DOI/PMID/title/URL columns to build a literature knowledge base, expand it through citation discovery, classify papers with LLM extraction, and browse results in an HTML dashboard.

Local skill path on this machine:

```bash
/mnt/e/Codex/skills/Literature-agent
```

Main scripts:

```bash
/mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py
/mnt/e/Codex/skills/Literature-agent/discovery.py
```

Preferred runtime:

```bash
conda activate md
```

Known local working files:

```bash
/mnt/e/zotero.csv
/mnt/d/文献库/zotero.sqlite
/mnt/e/literature_expanded_llm
/mnt/e/literature_expanded_llm/literature_blacklist_memory.json
/mnt/e/literature_expanded_llm/literature_metadata_cache.json
```

Windows equivalents:

```text
E:\Codex\skills\Literature-agent
E:\zotero.csv
D:\文献库\zotero.sqlite
E:\literature_expanded_llm
```

## Daily Command

For the user's current Zotero workflow, prefer this command:

```bash
conda activate md

python /mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py \
  --input /mnt/e/zotero.csv \
  --output /mnt/e/literature_expanded_llm \
  --schema-profile photo-radical \
  --use-llm \
  --discover \
  --expand-mode both \
  --max-depth 1 \
  --breadth-limit 30 \
  --zotero-sqlite "/mnt/d/文献库/zotero.sqlite" \
  --target-non-other 100 \
  --max-other 20
```

For complex goal-driven searches, do not rely only on the generic `photo-radical` profile. Pass the user's research objective as a task prompt so the LLM classifies papers against that goal instead of merely clustering keywords.

For the flavin/photoenzymatic discovery task, prefer:

```bash
conda activate md

python /mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py \
  --input /mnt/e/zotero.csv \
  --output /mnt/e/literature_flavin_photoenzyme \
  --schema-profile flavin-photoenzyme \
  --field "photoenzymatic field, including photoenzyme, Photoenzymology, photobiocatalysis, enzyme photocatalysis" \
  --enzyme-feature "flavin cofactor; flavin-dependent enzymes; FAD; FMN; flavoenzyme; OYE/ERED/FAP if relevant" \
  --reaction "new-to-nature reactions: rare or absent in biological metabolism but common or useful in chemical synthesis" \
  --goal "Use evidence-backed papers to find reported and unreported flavin-dependent enzymes with potential for new-to-nature photoenzymatic reactions; export candidates for later homolog/tree analysis." \
  --use-llm \
  --llm-scope criteria-ambiguous \
  --discover \
  --expand-mode both \
  --max-depth 1 \
  --breadth-limit 30 \
  --target-non-other 120 \
  --max-other 10
```

This criteria-first mode writes `criteria_spec.md`, `task_spec.md`, annotates each paper with `criteria_status` (`criteria_pass`, `criteria_borderline`, or `criteria_fail`), and exports `homolog_candidate_seeds.csv` for downstream enzyme/homolog curation. It expands obvious synonyms for the three hard anchors, including photoenzymatic/photoenzyme/Photoenzymology/photobiocatalysis and flavin/FAD/FMN/flavoenzyme/OYE/ERED/FAP. Use `--criteria-file criteria.txt` or `.json` when the reaction sentence or goal is long.

Use full `--llm-scope all` only for small curated batches. For broad Zotero/discovery runs, `--llm-scope criteria-ambiguous` is the preferred low-token setting because it sends only criterion-matching ambiguous records to the LLM.

Legacy full-prompt command:

```bash
conda activate md

python /mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py \
  --input /mnt/e/zotero.csv \
  --output /mnt/e/literature_flavin_photoenzyme \
  --schema-profile flavin-photoenzyme \
  --task-prompt "寻找做photoenzymatic领域需要用到flavin作为辅酶因子的酶：先搜罗这些flavin酶已经被应用去合成new-to-nature反应的报道，再搜索结构或功能相似、属于flavin酶且有正规酶学表征报道、但尚未在近年photoenzymatic领域被应用的候选酶文献。" \
  --use-llm \
  --discover \
  --expand-mode both \
  --max-depth 1 \
  --breadth-limit 30 \
  --target-non-other 120 \
  --max-other 10
```

This writes `task_spec.md` and automatically adds flavin/photoenzyme search queries when the task prompt mentions flavin/FAD/FMN. Use `--task-file task.txt` when the objective is long or contains quoting-sensitive characters.

If Bash reports `--use-llm: command not found` or argparse reports an empty/unrecognized argument, the line-continuation backslash was not parsed. In Bash, `\` must be the final character on the line with no trailing spaces or copied invisible characters. Use this single-line version when in doubt:

```bash
conda activate md

python /mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py --input /mnt/e/zotero.csv --output /mnt/e/literature_expanded_llm --schema-profile photo-radical --use-llm --discover --expand-mode both --max-depth 1 --breadth-limit 30 --zotero-sqlite "/mnt/d/文献库/zotero.sqlite" --target-non-other 100 --max-other 20
```

For the hydrolase/base.xlsx normalized input, use the same single-line pattern:

```bash
python /mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py --input /mnt/e/alpha_beta_hydrolase/hydrolase_literature_seed.csv --output /mnt/e/alpha_beta_hydrolase/literature_expanded_hydrolase --schema-profile lipase --use-llm --discover --expand-mode both --max-depth 1 --breadth-limit 30 --target-non-other 100 --max-other 20
```

Outputs to open first:

```text
/mnt/e/literature_expanded_llm/literature_dashboard.html
/mnt/e/literature_expanded_llm/citation_network.html
```

Do not default to reading every `cluster_*.md`; the dashboard is now the preferred browsing surface.

## Input Normalization

The pipeline script is named `Start-from-zotero-csv.py`, but the conceptual input is not limited to Zotero. Any table is acceptable once it is normalized to a Zotero-like CSV with these columns:

```text
Title, Author, Publication Year, DOI, Abstract Note, Publication Title, Url, Tags
```

The most important seed field is DOI. If DOI is present, the agent can use Crossref/OpenAlex/Semantic Scholar discovery even when other columns are messy. If DOI is absent but UniProt/PDB/PMID/title is present, first run a domain-specific normalizer or metadata lookup to recover DOI candidates, then feed the normalized CSV here.

For hydrolase seed tables, use:

```bash
conda activate md

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/extract_seed_metadata_from_xlsx.py \
  /mnt/e/alpha_beta_hydrolase/base.xlsx \
  --output /mnt/e/alpha_beta_hydrolase/seed_provenance_template_v2.csv

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/uniprot_seed_literature.py \
  /mnt/e/alpha_beta_hydrolase/seed_provenance_template_v2.csv \
  --outdir /mnt/e/alpha_beta_hydrolase/seed_literature_uniprot

python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/base_to_literature_seed_csv.py \
  /mnt/e/alpha_beta_hydrolase/seed_literature_uniprot/seed_provenance_uniprot_enriched.csv \
  --output /mnt/e/alpha_beta_hydrolase/hydrolase_literature_seed.csv
```

Then run this skill on `/mnt/e/alpha_beta_hydrolase/hydrolase_literature_seed.csv`.

## What The Pipeline Does

1. Reads Zotero-like CSV seed papers. A real Zotero export is one valid input, but normalized DOI/title tables are also valid.
2. Optionally expands by Semantic Scholar citations/references/related papers.
3. Requires DOI for discovered papers; missing-DOI discovered papers give their slot to later candidates.
4. Enriches records with Crossref/OpenAlex subjects, concepts, keywords, and abstracts when available.
5. Runs LLM extraction using `KIMI_API_KEY`, `MOONSHOT_API_KEY`, or `--api-key`.
6. Classifies both enzymatic and non-enzymatic photochemistry:
   - `photoenzymatic/enzymatic`
   - `enzymatic non-photo`
   - `organic photoredox/radical no-enzyme`
   - `organic synthesis no-enzyme`
   - `unrelated/blacklist`
   For `--schema-profile flavin-photoenzyme`, the primary classification becomes:
   - `A_used_flavin_photoenzyme_new_to_nature`
   - `B_characterized_flavin_enzyme_not_photoenzymatic`
   - `C_photoenzyme_non_flavin`
   - `D_organic_photochemistry_no_enzyme`
   - `E_background_or_unrelated`
7. Resolves `Other_enzyme`, `unclear`, and noisy `*ase` guesses with a second LLM pass.
8. For unresolved Other/unclear papers, optionally tries to fetch an introduction snippet from the DOI landing page or URL.
9. Applies a quality filter and writes permanent blacklist memory.
10. Generates compact reports plus interactive HTML dashboards.

## Dashboard-First Outputs

Main dashboard:

```text
literature_dashboard.html
```

Use it to search and filter by:

- title, DOI, abstract, metadata keywords
- paper domain
- enzyme family
- reaction type
- manual review priority

Citation network dashboard:

```text
citation_network.html
```

It shows seed vs discovered papers and citation edges. This is the preferred replacement for a pile of markdown clusters.

Still generated:

- `knowledge_report.md`
- `gap_analysis.md`
- `ambiguous_enzyme_review.md`
- `blacklisted_papers.md`
- `discovery_summary.md`
- `network_analysis.md`
- `blindspot_analysis.md`
- `discovered_papers.md`
- `citation_distribution.png`
- `knowledge_graph.html` if `pyvis` is installed

Only generate all cluster markdown files when explicitly requested:

```bash
--export-cluster-md
```

## Blacklist Memory

Permanent memory lives at:

```bash
/mnt/e/literature_expanded_llm/literature_blacklist_memory.json
```

This JSON records both kept and blacklisted papers by DOI/title hash. On later runs:

- blacklisted records are skipped before LLM extraction
- kept records provide prior `paper_domain` and `enzyme_family`
- quality-filter decisions are saved back to the same file

By default it is written inside the output directory as `literature_blacklist_memory.json`. Use a different memory file only when intentionally starting a separate literature universe:

```bash
--blacklist-memory /mnt/e/my_other_project_blacklist.json
```

## Metadata Cache

Crossref/OpenAlex metadata cache lives at:

```bash
/mnt/e/literature_expanded_llm/literature_metadata_cache.json
```

By default it is written inside the output directory as `literature_metadata_cache.json`. It stores subjects, concepts, keywords, and abstracts. These are fed into fallback and LLM prompts as context. Disable this only for debugging or offline runs:

```bash
--no-external-metadata
```

## Other Enzyme Policy

The workflow should not casually accept `Other_enzyme`.

Rules:

- If the paper has no enzyme/protein/biocatalysis signal, classify as `no enzyme`, not Other.
- Pure photoredox/radical chemistry can still be kept when it has enzyme-pocket migration value.
- Pure organic synthesis with no photoredox/radical/enzyme relevance should be blacklisted.
- False `*ase` words such as `release`, `base`, `showcase`, `disease`, and `phase` must not become enzyme families.
- Default final quota: at least 100 non-Other useful papers if available; at most 20 Other papers.

Relevant flags:

```bash
--target-non-other 100
--max-other 20
--resolve-other-with-llm
--no-resolve-other-with-llm
--fetch-introduction-for-other
--no-fetch-introduction-for-other
--disable-quality-filter
```

`--resolve-other-with-llm` and `--fetch-introduction-for-other` are enabled by default.

## Complex Task Handling

When the user gives a long, nuanced research goal, first preserve it as a task prompt instead of reducing it to one keyword query. The task prompt should define:

- the positive anchor set, such as flavin/FAD/FMN enzymes already used for photoenzymatic new-to-nature synthesis
- the candidate expansion set, such as characterized flavin enzymes with structures, kinetics, activity assays, or substrate scope but no recent photoenzymatic application
- exclusion rules, such as genome-only annotation, non-enzyme photoredox chemistry unless it suggests a transferable reaction, or papers without DOI/title evidence
- the desired review fields, such as cofactor, enzyme family, reaction type, new-to-nature status, structure/PDB evidence, and manual priority

Use `--task-prompt` for one-line objectives or `--task-file` for longer text. The prompt is injected into every LLM extraction call and is also saved as `task_spec.md` in the output directory.

For three-anchor screening, prefer explicit criteria over a single long LLM prompt:

- `--field`: domain noun or phrase. The script expands fuzzy field synonyms such as photoenzymatic, photoenzyme, Photoenzymology, photobiocatalysis, and enzyme photocatalysis.
- `--enzyme-feature`: enzyme/cofactor criterion. For flavin work it expands flavin, FAD, FMN, flavoenzyme, old yellow enzyme, OYE, ERED, and FAP.
- `--reaction`: the longer reaction sentence, such as new-to-nature reactions that are rare in biology but common/useful in chemical synthesis.
- `--goal`: final research objective. This is preserved for ranking and manual review, but it is not treated as a strict keyword that every title/abstract must contain.

With `--criteria-action llm-gate` and `--llm-scope criteria-ambiguous`, the pipeline first runs deterministic triage, then spends LLM calls only on criterion-matching ambiguous records. `--criteria-action filter` removes `criteria_fail` records from the current run and writes `criteria_filtered_papers.csv`; use it only after checking the dashboard behavior on a small run.

For the user's flavin/photoenzyme task, use `--schema-profile flavin-photoenzyme`. This profile is designed to separate:

- papers that already report flavin-dependent photoenzymatic new-to-nature reactions
- characterized flavin enzymes not yet used in photoenzymatic synthesis
- photoenzymatic papers that do not clearly involve flavin
- organic photochemistry papers without enzymes
- weak background or unrelated records

For handoff to homolog/tree work, review `homolog_candidate_seeds.csv` first. It contains DOI/title/year/journal, criteria hits, flavin cofactor evidence, enzyme family/name hints, reaction type, characterization evidence, and structure/PDB-like hints. Treat it as a curation table for choosing enzyme seeds; do not feed every row directly into homolog expansion without confirming UniProt/PDB/sequence evidence.

## Discovery Parameters

Use Semantic Scholar discovery through the main script:

```bash
--discover
--expand-mode both          # citations | references | both | related
--max-depth 1
--breadth-limit 30
--relevance-threshold 0.25
--search-query "lipase photoredox"
--ss-api-key "$SEMANTIC_SCHOLAR_API_KEY"
```

Semantic Scholar free tier is rate-limited. Increase `--breadth-limit` before increasing `--max-depth`; depth 2 can grow quickly.

## Zotero Local Linkage

The dashboard marks papers as:

- `zotero_seed`: present in the Zotero CSV input.
- `already_in_zotero`: not in the CSV seed, but found in a supplied Zotero local database.
- `new_discovery`: discovered by API expansion and not matched to Zotero.

For normal runs, the input CSV is enough to mark seeds. To compare against the full Zotero library, pass a read-only local SQLite path:

```bash
--zotero-sqlite "/mnt/d/文献库/zotero.sqlite"
```

The script only reads bibliographic fields from SQLite: DOI, title, and publication/journal title. It does not read the `storage` directory, whose folder names are Zotero attachment keys and can look like random text. Because the path contains Chinese characters, always quote it in WSL commands. If Zotero is open and the database is locked, close Zotero or copy `zotero.sqlite` to a temporary location and pass the copy.

The run also writes:

```text
zotero_overlap_report.md
```

Use this report or the dashboard's Zotero filter to distinguish CSV seeds, items already present in the full Zotero library, and genuinely new discovery results.

This local linkage is the recommended first step before building a full Zotero plugin. A future plugin could consume `literature_dashboard.html` or the JSON memory/cache outputs to tag items inside Zotero.

Standalone discovery is still available:

```bash
python /mnt/e/Codex/skills/Literature-agent/discovery.py \
  --input /mnt/e/zotero.csv \
  --output /mnt/e/literature_discovery_only \
  --mode both \
  --max-depth 1 \
  --breadth-limit 30 \
  --export-gexf \
  --export-edgelist
```

## LLM Settings

Do not hard-code API keys. Use:

```bash
export KIMI_API_KEY="..."
```

or:

```bash
export MOONSHOT_API_KEY="..."
```

Defaults:

- base URL: `https://api.moonshot.cn/v1`
- model: `kimi-k2.6`
- temperature: `1`
- max tokens: `3600`
- read timeout: `360`

If JSON mode fails:

```bash
--no-json-mode
```

If Kimi returns HTTP 200 but `message.content` is empty and only `reasoning_content` appears, the request reached the API but the extraction did not succeed. The script now treats that as a failed attempt and retries with alternate token/json-mode settings. For long reasoning models, prefer:

```bash
--max-tokens 6000
```

To reduce load on large batches, use:

```bash
--llm-scope ambiguous
```

This lets the pipeline skip obvious papers and send only ambiguous cases to the LLM.

For goal-driven criteria searches, prefer:

```bash
--llm-scope criteria-ambiguous
```

This skips `criteria_fail` records and sends only `criteria_pass`/`criteria_borderline` records whose enzyme family, flavin status, new-to-nature status, characterization evidence, or abstract is still ambiguous.

For `kimi-k2.6`, the script disables Kimi thinking by default because literature extraction needs short JSON, not long reasoning. Only enable thinking intentionally:

```bash
--enable-kimi-thinking
```

Also, `kimi-k2.6` only accepts `temperature=0.6`; the script now auto-overrides other values.

## Input CSV

The script expects a Zotero CSV. Common columns:

- `Title`
- `Author`
- `Publication Year`
- `DOI`
- `Abstract Note`
- `Publication Title`
- `Url`
- `File Attachments`
- `Link Attachments`

CSV encoding is auto-detected across UTF-8, UTF-16, and GB18030. Use this if needed:

```bash
--csv-encoding utf-16
```

## Introduction Fetching

Semantic Scholar, Crossref, and OpenAlex usually provide abstracts and metadata, not stable full introductions.

For unresolved Other/unclear papers, the main script can make a best-effort attempt to fetch a DOI landing page or URL and extract an introduction-like HTML section:

```bash
--fetch-introduction-for-other
```

This is not a full PDF parser. If introduction extraction becomes central, add a future PDF/HTML full-text parser layer.

## Validation

After editing either script:

```bash
conda activate md
python -m py_compile \
  /mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py \
  /mnt/e/Codex/skills/Literature-agent/discovery.py
```

Quick help check:

```bash
python /mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py --help
```

Small no-API smoke test:

```bash
python /mnt/e/Codex/skills/Literature-agent/Start-from-zotero-csv.py \
  --input /mnt/e/zotero.csv \
  --output /mnt/e/literature_smoke_test \
  --schema-profile photo-radical \
  --no-external-metadata \
  --no-resolve-other-with-llm
```

## Notes For Codex

- Use `/mnt/e/Codex/...`, not `/mnt/e/CodeX/...`; WSL paths can be case-sensitive.
- Prefer WSL commands with `conda activate md`.
- Keep outputs under `/mnt/e`.
- The dashboard is the primary review interface, similar in spirit to active-learning dashboards.
- Only export cluster markdown when the user explicitly asks for bulk cluster files.
