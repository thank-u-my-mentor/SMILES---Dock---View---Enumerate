# Local Paths

Use these paths for this user's WSL/Windows setup unless the project says otherwise.

## Environments

AmberTools:

```bash
source /mnt/l/WSL/conda_envs/AmberTools25/amber.sh
```

GROMACS:

```bash
/home/qin/softwares/gromacs-2026.0-gpu/bin/gmx_mpi
```

GAMESS:

```bash
/home/qin/softwares/gamess/rungms
```

GAMESS scratch/restart:

```bash
/home/qin/softwares/gamess/scratch
/home/qin/softwares/gamess/restart
```

Python for MD utilities, when needed:

```bash
/mnt/l/WSL/softwares/conda_envs/md/bin/python
```

## Current TJ Project Roots

Current organized project:

```bash
/mnt/e/TJ/260706_mcpb_m2_200ns
```

Final MCPB/tleap files:

```bash
/mnt/e/TJ/260706_mcpb_m2_200ns/mcpb_final
```

Completed 200 ns MD:

```bash
/mnt/e/TJ/260706_mcpb_m2_200ns/md_200ns
```

PBC-fixed visualization files:

```bash
/mnt/e/TJ/260706_mcpb_m2_200ns/md_200ns/pbc
```

Old raw long-run project:

```bash
/mnt/e/TJ/tj2_unl_NH_C4radical_M6_20260629_v1/mcpb
```

Current handbook:

```text
E:\QZHAO-LAB\MCPB.py 金属中心流程代码本📚.md
```

Working handbook copy:

```text
E:\Codex\MCPB_manual_working.md
```

Project quick PyMOL loader:

```pymol
@E:/TJ/260706_mcpb_m2_200ns/md_200ns/pbc/pymol_view_mcpb_md.pml
```

## Monitoring Commands

Hessian bad-point monitor:

```bash
tail -300 gamess_logs/JOB.log | grep -aE "IVIB=|FINAL U-B3LYP ENERGY|SCF IS UNCONVERGED|TOO MANY ITERATIONS|EXECUTION OF GAMESS"
```

OPT monitor:

```bash
tail -300 gamess_logs/JOB.log | grep -aE "NSERCH:|MAXIMUM GRADIENT|RMS GRADIENT|FINAL U-B3LYP ENERGY|SCF IS UNCONVERGED|EXECUTION OF GAMESS"
```

Large MK monitor:

```bash
tail -300 gamess_logs/JOB.log | grep -aE "DENSITY CONVERGED|SCF IS UNCONVERGED|FINAL U-B3LYP ENERGY|S-SQUARED|ELECTROSTATIC POTENTIAL|NUMBER OF POINTS SELECTED|EXECUTION OF GAMESS"
```
