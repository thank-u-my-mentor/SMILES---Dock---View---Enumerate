# TJ Canonical Layout And Cleanup Notes

This file records the intended minimal TJ project layout for the Fe(III)-radical
MCPB -> long MD workflow. It should be treated as a practical project map, not a
chemistry paper.

## Canonical Root

After the M7 MCPB workflow is reproducible, the TJ root should be easy for a new
student to read:

```text
E:/TJ/
  1T47.pdb
  1T47_chainA_Hpp_try_keepFe_NTD_original_names.pkout.txt
  protonation_map.csv
  protonation_map_metal_contacts.csv
  mcpb_m7_400ns/
  archive_legacy_YYYYMMDD/
  TJ_WORKFLOW_README.md
```

The active workflow should live in:

```text
E:/TJ/mcpb_m7_400ns/
  01_inputs/
  02_qm_hessian/
  03_mcpb/
  04_amber/
  05_gromacs/
  06_analysis/
  logs/
  README.md
  MCPB_BUG_AVOIDANCE_ZH.md
```

## Current Top-Level Audit

Snapshot from 2026-06-18:

```text
616_ligand_inspect/             14 files,    0.03 MB  archive candidate
fe_ligand_building_blocks/      18 files,    0.02 MB  archive candidate
feiii_radical_preflight/        79 files,   66.63 MB  keep as provenance for now
hplusplus_inputs/                4 files,    0.64 MB  keep as provenance
mcpb_m7_400ns/                  49 files,    1.58 MB  active canonical workflow
md_analysis/                    12 files,   30.83 MB  old non-MCPB/harmonic analysis
md_handoff/                     50 files, 2548.61 MB  old non-MCPB/harmonic trajectory
prepared/                        2 files,    0.21 MB  archive candidate
pymol/                           3 files,    0.00 MB  archive candidate
qm_ligand_b3lyp_631gdp/          7 files,    0.34 MB  old ligand-only QM attempt
qm_ligand_b3lyp_631gdp_ACT/     15 files,    4.08 MB  old ligand-only QM attempt
```

Root files:

```text
0616_initial_guess*.pdb          keep as provenance until final MD starts
1T47.pdb                         keep raw input
1T47_*pkout.txt                   keep H++ output
616.pse / 616.png / 616_2.png     archive candidate
acetate.pdb / nbutamide.pdb       archive candidate
dock_result.pdb                   archive candidate
protonation_map*.csv              keep
```

## Keep / Archive / Delete Policy

Keep active:

```text
mcpb_m7_400ns/
1T47.pdb
1T47_chainA_Hpp_try_keepFe_NTD_original_names.pkout.txt
protonation_map.csv
protonation_map_metal_contacts.csv
```

Keep as provenance until the 400 ns MD workflow has launched successfully:

```text
feiii_radical_preflight/
hplusplus_inputs/
0616_initial_guess*.pdb
```

Archive candidates:

```text
616_ligand_inspect/
fe_ligand_building_blocks/
md_analysis/
md_handoff/
prepared/
pymol/
qm_ligand_b3lyp_631gdp/
qm_ligand_b3lyp_631gdp_ACT/
616.pse
616.png
616_2.png
acetate.pdb
nbutamide.pdb
dock_result.pdb
```

Do not delete `md_handoff/` immediately even though it is large. It contains old
trajectory evidence and should first be moved into an archive folder. Delete it
only after the MCPB-based trajectory and analysis are reproducible.

## Suggested Archive Command

Use this only after the current GAMESS/MCPB job is not using the files being
moved. Do not move `mcpb_m7_400ns/`.

```powershell
$archive = 'E:\TJ\archive_legacy_20260618'
New-Item -ItemType Directory -Force -Path $archive
Move-Item -LiteralPath `
  'E:\TJ\616_ligand_inspect',`
  'E:\TJ\fe_ligand_building_blocks',`
  'E:\TJ\md_analysis',`
  'E:\TJ\md_handoff',`
  'E:\TJ\prepared',`
  'E:\TJ\pymol',`
  'E:\TJ\qm_ligand_b3lyp_631gdp',`
  'E:\TJ\qm_ligand_b3lyp_631gdp_ACT',`
  'E:\TJ\616.pse',`
  'E:\TJ\616.png',`
  'E:\TJ\616_2.png',`
  'E:\TJ\acetate.pdb',`
  'E:\TJ\nbutamide.pdb',`
  'E:\TJ\dock_result.pdb' `
  -Destination $archive
```

If the goal is a teaching template, archive rather than delete. A clean root
helps students reproduce the workflow, while archived files still preserve the
history of failed or exploratory branches.
