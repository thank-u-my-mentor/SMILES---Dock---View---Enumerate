---
name: metal-f-element-rescan
description: Rescan RCSB/PDB for fluorinated ligands and transition-metal components using wwPDB CCD element formulas and mmCIF coordinates. Use when working on Metal-F projects, fluorinated ligand/metal geometry, F-to-metal distance filtering, complex-subunit warnings, RCSB CCD ligand enumeration, protein/CDS FASTA sidecar export, or auditing which metal-containing CCD components are likely real cofactors versus salts or synthetic complexes.
---

# Metal-F Element Rescan

Use this skill for element-level Metal-F discovery when ligand names or comp_ids are unpredictable. Do not rely on keywords such as `CF3`, `fluoroaryl`, or `F-` alone; many true hits use arbitrary CCD IDs such as `SVF`, `B5N`, `JKB`, or `K3U`.

## Core Rules

- Treat RCSB Search API results as candidates only. Geometry must be confirmed from downloaded mmCIF coordinates.
- Define fluorinated ligand candidates from the wwPDB Chemical Component Dictionary (CCD): component formula contains element `F`.
- Define metal candidates from CCD formula containing target elements. Default target metals in the bundled script are `FE, CO, NI, MN, CR, CU`.
- Confirm hits by calculating non-solvent `HETATM` fluorine atom to non-solvent `HETATM` target-metal atom distance. Default cutoff is `3.5 Å`.
- Preserve user-curated source files. Write sidecar outputs with a new prefix instead of overwriting existing xlsx/FASTA/PSE.
- Mark `protein_entity_count > 1` as `complex_warning=YES`; keep such rows in the xlsx evidence table, but do not include them in seed FASTA by default.
- Include protein FASTA and CDS/DNA FASTA only for single-protein-entity hits unless the user explicitly asks to include complexes.

## Scripts

Run scripts with the WSL conda `md` Python when working in the user's Metal-F environment:

```bash
/mnt/l/WSL/softwares/conda_envs/md/bin/python <script> [args]
```

### Rescan Current Curated XLSX

Use this when the user has manually added PDB IDs to an existing xlsx and wants a clean sidecar table plus FASTA:

```bash
/mnt/l/WSL/softwares/conda_envs/md/bin/python \
  /mnt/e/Codex/skills/metal-f-element-rescan/scripts/metal_f_element_ligand_rescan.py \
  --pdb-ids-from-xlsx /home/qin/Metal-F_project/Metal-F_interactions/Metal-F_geometry_fast_hits.xlsx \
  --prefix /home/qin/Metal-F_project/Metal-F_interactions/Metal-F_element_ligand_rescan_from_current_xlsx
```

Expected outputs:

- `<prefix>.xlsx`
- `<prefix>_protein.fasta`
- `<prefix>_cds_dna.fasta`

The xlsx contains `geometry_hits`, `protein_fasta`, `dna_fasta`, `ccd_f_ligands`, `ccd_metal_ligands`, and `run_summary` sheets.

### Rescan Explicit PDB IDs

Use this for a focused patch set:

```bash
/mnt/l/WSL/softwares/conda_envs/md/bin/python \
  /mnt/e/Codex/skills/metal-f-element-rescan/scripts/metal_f_element_ligand_rescan.py \
  --pdb-ids "5OAE 6EI4 6JFR 6QXD" \
  --prefix /home/qin/Metal-F_project/Metal-F_interactions/Metal-F_element_ligand_rescan_patch
```

### Full CCD-Based Candidate Rescan

Full CCD enumeration can yield thousands of candidate PDB entries. Run it in batches:

```bash
/mnt/l/WSL/softwares/conda_envs/md/bin/python \
  /mnt/e/Codex/skills/metal-f-element-rescan/scripts/metal_f_element_ligand_rescan.py \
  --prefix /home/qin/Metal-F_project/Metal-F_interactions/Metal-F_element_ligand_rescan_full \
  --start 0 \
  --max-candidates 500
```

Then continue with `--start 500`, `--start 1000`, etc. The first run caches CCD and candidate lists under the project cache.

For resumable terminal runs, use the candidate cache and done ledger:

```bash
# First run: build candidate cache and process the first batch.
/mnt/l/WSL/softwares/conda_envs/md/bin/python \
  /home/qin/Metal-F_project/metal_f_element_ligand_rescan.py \
  --prefix /home/qin/Metal-F_project/Metal-F_interactions/Metal-F_element_ligand_rescan_full_batch0000 \
  --start 0 \
  --max-candidates 500 \
  --skip-dna \
  --delete-cif-nonhits \
  --quiet-candidates \
  --resume

# Later runs: reuse the candidate cache and skip PDB IDs already in the done ledger.
/mnt/l/WSL/softwares/conda_envs/md/bin/python \
  /home/qin/Metal-F_project/metal_f_element_ligand_rescan.py \
  --prefix /home/qin/Metal-F_project/Metal-F_interactions/Metal-F_element_ligand_rescan_full_batch0500 \
  --reuse-candidates \
  --start 500 \
  --max-candidates 500 \
  --skip-dna \
  --delete-cif-nonhits \
  --quiet-candidates \
  --resume

# Check how many remain in the selected candidate set.
/mnt/l/WSL/softwares/conda_envs/md/bin/python \
  /home/qin/Metal-F_project/metal_f_element_ligand_rescan.py \
  --reuse-candidates \
  --status-only
```

Caching layers:

- `components.cif.gz`: CCD cache.
- `.cache/metal_f_element_ligand_rescan/cif/*.cif`: mmCIF cache.
- `.cache/metal_f_element_ligand_rescan/Metal-F_element_ligand_candidates.csv`: candidate PDB IDs.
- `.cache/metal_f_element_ligand_rescan/Metal-F_element_ligand_done.csv`: processed PDB ledger.

`--delete-cif-nonhits` deletes CIF files for candidates that fail the distance rule only when that CIF was newly downloaded in the current run. Add `--delete-cif-nonhits-even-if-cached` only when the user explicitly wants to purge previously cached non-hit CIF files too.

Every processed candidate, including `no_geometry_hit` and `no_f_or_metal_atoms`, is recorded in both the done ledger and the output xlsx `processed_candidates` sheet. Only passing rows are written to `geometry_hits`.

## Metal Component Audit

Use the audit script when the user asks why the CCD metal component count is large, or wants examples of accidental/synthetic/non-enzymatic metal components:

```bash
/mnt/l/WSL/softwares/conda_envs/md/bin/python \
  /mnt/e/Codex/skills/metal-f-element-rescan/scripts/audit_metal_components.py \
  --out /home/qin/Metal-F_project/Metal-F_interactions/Metal-F_ccd_metal_component_audit.xlsx
```

Interpret audit classes as heuristic triage, not a chemical ontology:

- `known_biological_metal_cofactor`: examples like hemes or Fe-S clusters.
- `simple_metal_ion_or_small_salt`: free ions or very small salts.
- `metal_salt_or_counterion`: crystallization/additive/counterion-like components.
- `synthetic_coordination_complex_or_probe`: metal complexes that may not be enzymatic cofactors.
- `large_metal_organic_component`: large metal-organic ligands.
- `ambiguous_metal_component`: inspect manually.

If excluding `ambiguous_metal_component` barely changes candidate entry counts, do not keep tightening only the metal side. That usually means the large candidate count is driven by common metal ions/cofactors (`FE`, `CU`, `NI`, `MN`, `HEM`, Fe-S clusters) appearing across many protein structures. The next useful filter is the fluorinated-ligand side: remove free fluoride-only searches, inorganic fluorinated counterions, PF6/BF4-like salts, crystallization additives, and nucleic-acid complexes before spending time on full mmCIF geometry.

## Output Interpretation

Prioritize rows with:

- `status=geometry_hit`
- `nearest_F_metal_A <= 3.5`
- `complex_warning=NO`
- `fasta_seed_included=YES`
- `ecoli_expression_flag=YES` when heterologous E. coli expression is the user's priority.

Rows with `complex_warning=YES` are still useful structural evidence, but avoid using them as HMMER seed FASTA unless the user accepts multi-subunit or native-complex targets.

## Plasmid Candidate Review

Use this after the user manually curates a DOCX/list of PDB IDs for possible BL21(DE3) heterologous expression and plasmid ordering. This workflow corrects source organism versus expression host, de-duplicates by UniProt and `>90%` sequence identity, exports an interactive HTML dashboard plus CSV report, BUY protein FASTA, BUY CDS FASTA, and an all-nonduplicate PSE.

```bash
/mnt/l/WSL/softwares/conda_envs/md/bin/python \
  /mnt/e/Codex/skills/metal-f-element-rescan/scripts/build_plasmid_candidate_review.py \
  --docx /mnt/e/QZHAO-LAB/M-F酶.docx \
  --xlsx /home/qin/Metal-F_project/Metal-F_interactions/Metal-F_element_ligand_rescan_full_merged.xlsx \
  --protein-fasta /home/qin/Metal-F_project/Metal-F_interactions/Metal-F_element_ligand_rescan_full_merged_protein.fasta \
  --outdir /mnt/e/QZHAO-LAB/Metal-F_plasmid_candidate_review
```

Decision labels:

- `BUY`: bacterial source, single protein entity, reported E. coli expression, non-redundant representative.
- `REVIEW`: potentially useful but non-bacterial source, unknown expression, or other caution.
- `AVOID_DUPLICATE`: redundant UniProt or `>90%` sequence-identity cluster member.
- `AVOID_COMPLEX`: multi-protein entity/complex; avoid for simple single-gene plasmid purchase.
- `AVOID_RISK`: native/non-bacterial/unknown expression risk that is poor for first-pass plasmid ordering.

For plasmid-order triage, add literature/SI expression evidence as a curated overlay when publisher SI files cannot be fetched automatically. Keep source organism, expression host/strain, vector/helper plasmid, and decision override in separate columns. Examples: BL21(DE3) expression can upgrade non-bacterial sources to `BUY`; Rosetta-only, Fe-S helper plasmid, native-purified protein, or special vector/cofactor systems should usually become `AVOID_RISK` unless the user explicitly accepts that burden.

DNA FASTA is not automatically equivalent to an orderable construct. Treat EMBL/ENA CDS as source evidence and run translation QC against the exact PDB protein sequence. The workflow writes full BUY CDS to `*_BUY_cds_dna.fasta`, QC metrics to `*_BUY_dna_index.csv`, and only exact/PDB-extra-residue matches to `*_BUY_pdb_matched_cds_dna.fasta`. Near matches, signal-peptide/propeptide cases, tags, truncations, and mutant constructs stay in the QC table for manual construct design.

Generate the PSE for all non-redundant candidates. Include `BUY`, `REVIEW`, `AVOID_RISK`, and `AVOID_COMPLEX`; exclude only `AVOID_DUPLICATE` rows created by UniProt duplication or `>90%` pairwise identity.

```bash
cd /tmp
/mnt/l/WSL/softwares/conda_envs/md/bin/pymol -cq \
  -r /mnt/e/Codex/skills/metal-f-element-rescan/scripts/build_element_rescan_pse.py -- \
  --xlsx /home/qin/Metal-F_project/Metal-F_interactions/Metal-F_element_ligand_rescan_full_merged.xlsx \
  --output-pse /mnt/e/QZHAO-LAB/Metal-F_plasmid_candidate_review/Metal-F_plasmid_all_nonduplicate_candidates.pse \
  --summary-csv /mnt/e/QZHAO-LAB/Metal-F_plasmid_candidate_review/Metal-F_plasmid_all_nonduplicate_candidates_pse_summary.csv \
  --pdb-ids-from-csv /mnt/e/QZHAO-LAB/Metal-F_plasmid_candidate_review/Metal-F_plasmid_candidate_review.csv \
  --csv-exclude-column final_decision \
  --csv-exclude-value AVOID_DUPLICATE
```

Render every PSE scene to a per-PDB ray-traced PNG folder and write a dashboard-friendly image index:

```bash
cd /tmp
/mnt/l/WSL/softwares/conda_envs/md/bin/pymol -cq \
  -r /mnt/e/Codex/skills/metal-f-element-rescan/scripts/render_pse_scene_pngs.py -- \
  --pse /mnt/e/QZHAO-LAB/Metal-F_plasmid_candidate_review/Metal-F_plasmid_all_nonduplicate_candidates.pse \
  --summary-csv /mnt/e/QZHAO-LAB/Metal-F_plasmid_candidate_review/Metal-F_plasmid_all_nonduplicate_candidates_pse_summary.csv \
  --review-csv /mnt/e/QZHAO-LAB/Metal-F_plasmid_candidate_review/Metal-F_plasmid_candidate_review.csv \
  --outdir /mnt/e/QZHAO-LAB/Metal-F_plasmid_candidate_review/scene_png \
  --index-csv /mnt/e/QZHAO-LAB/Metal-F_plasmid_candidate_review/scene_image_index.csv \
  --width 1600 \
  --height 1200
```

After rendering, rebuild the dashboard so `scene_image_index.csv` is embedded and each candidate detail page shows a `Scene` tab.

Important: source organism determines source kingdom. Expression host only supports heterologous-expression confidence; never use E. coli expression host to relabel a human, plant, fungal, or animal source as bacterial.
