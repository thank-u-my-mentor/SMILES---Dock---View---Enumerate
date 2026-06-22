#!/usr/bin/env python3
"""Create a clean TJ M7 MCPB -> Amber/GROMACS 400 ns project directory."""

from __future__ import annotations

import argparse
import math
import shutil
from pathlib import Path


ATOMIC_NUMBERS = {
    "H": 1.0,
    "C": 6.0,
    "N": 7.0,
    "O": 8.0,
    "S": 16.0,
    "Fe": 26.0,
}


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def copy_required(src: Path, dst: Path) -> None:
    if not src.exists():
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def parse_gamess_final_coords(log: Path) -> tuple[list[str], dict[str, str]]:
    """Extract the last evaluated/converged coordinate block from a GAMESS log."""
    lines = read_text(log).splitlines()
    marker = "COORDINATES OF ALL ATOMS ARE (ANGS)"
    indices = [idx for idx, line in enumerate(lines) if marker in line]
    if not indices:
        raise RuntimeError(f"no GAMESS coordinate block in {log}")
    coords: list[str] = []
    for line in lines[indices[-1] + 1 :]:
        stripped = line.strip()
        if not stripped:
            if coords:
                break
            continue
        parts = stripped.split()
        if len(parts) == 5 and parts[0][0].isalpha():
            try:
                elem = parts[0]
                anum = float(parts[1])
                x, y, z = (float(parts[2]), float(parts[3]), float(parts[4]))
            except ValueError:
                if coords:
                    break
                continue
            coords.append(f"{elem:<2s} {anum:8.1f} {x:14.6f} {y:14.6f} {z:14.6f}")
        elif coords:
            break
    if not coords:
        raise RuntimeError(f"failed to parse final coordinates from {log}")

    meta = {
        "energy": "",
        "s_squared": "",
        "grad_max": "",
        "grad_rms": "",
    }
    for line in reversed(lines):
        if not meta["energy"] and "TOTAL ENERGY" in line and "=" in line:
            meta["energy"] = line.split("=")[-1].strip().split()[0]
        if not meta["s_squared"] and "S-SQUARED" in line and "=" in line:
            meta["s_squared"] = line.split("=")[-1].strip().split()[0]
        if "GRAD. MAX=" in line and "R.M.S.=" in line:
            parts = line.replace("=", " ").split()
            try:
                meta["grad_max"] = parts[parts.index("MAX") + 1]
                meta["grad_rms"] = parts[parts.index("R.M.S.") + 1]
            except Exception:
                pass
            break
    return coords, meta


def write_gamess_hessian_input(out: Path, coords: list[str], *, charge: int, mult: int) -> None:
    lines = [
        " $SYSTEM MEMDDI=400 MWORDS=200 $END",
        " $CONTRL",
        "  SCFTYP=UHF DFTTYP=B3LYP RUNTYP=HESSIAN",
        f"  ICHARG={charge} MULT={mult} MAXIT=200 COORD=CART UNITS=ANGS",
        " $END",
        " $SCF DIRSCF=.T. DIIS=.T. DAMP=.T. ETHRSH=2.0 MAXDII=20 $END",
        " $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=1 $END",
        " $DATA",
        "TJ FeIII radical M7 converged small model B3LYP/6-31G(d,p) Hessian",
        "C1",
        *coords,
        " $END",
        "",
    ]
    out.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def write_gamess_runner(qm_dir: Path, job_name: str, *, cores: int, gamess: Path) -> None:
    script = qm_dir / "run_hessian_gamess.sh"
    script.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
mkdir -p logs
JOB={job_name}
GAMESS={gamess}
VERSION=${{VERSION:-00}}
NCORES=${{NCORES:-{cores}}}
echo "[GAMESS] $JOB Hessian start $(date) NCORES=$NCORES" | tee logs/${{JOB}}.status
"$GAMESS" "$JOB" "$VERSION" "$NCORES" > logs/${{JOB}}.log 2>&1
rc=$?
echo "[GAMESS] $JOB exit_code=$rc end $(date)" | tee -a logs/${{JOB}}.status
if grep -q 'EXECUTION OF GAMESS TERMINATED NORMALLY' logs/${{JOB}}.log; then
  echo "[GAMESS] $JOB status=normal" | tee -a logs/${{JOB}}.status
else
  echo "[GAMESS] $JOB status=not_normal" | tee -a logs/${{JOB}}.status
fi
exit "$rc"
""",
        encoding="utf-8",
        newline="\n",
    )
    script.chmod(0o755)

    monitor = qm_dir / "monitor_hessian.sh"
    monitor.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
tail -n 80 logs/{job_name}.log 2>/dev/null || true
echo
grep -E 'RUNTYP=|HESSIAN|VIB|FREQUENC|TOTAL WALL CLOCK|TOTAL ENERGY|TERMINATED|ERROR|FAILURE|S-SQUARED' logs/{job_name}.log 2>/dev/null | tail -n 60 || true
""",
        encoding="utf-8",
        newline="\n",
    )
    monitor.chmod(0o755)


def write_mcpb_files(mcpb_dir: Path, *, group_name: str, sm_charge: int, mult: int) -> None:
    original_pdb = mcpb_dir / "mcpb_original.pdb"
    metal_serial = find_metal_serial(original_pdb, chain="A", resseq=431)
    mcpb = f"""original_pdb mcpb_original.pdb
group_name {group_name}
cut_off 2.8
ion_ids {metal_serial}
ion_mol2files FE.mol2
water_model TIP3P
force_field ff19SB
gaff 2
frcmod_files ACT.frcmod UNL.frcmod
software_version gms
large_opt 1
add_redcrd 1
scale_factor 1.0
smmodel_chg {sm_charge}
smmodel_spin {mult}
lgmodel_chg {sm_charge}
lgmodel_spin {mult}
naa_mol2files ACT.mol2 UNL.mol2 HOH.mol2
"""
    (mcpb_dir / "mcpb.in").write_text(mcpb, encoding="utf-8", newline="\n")
    (mcpb_dir / "FE.mol2").write_text(
        """@<TRIPOS>MOLECULE
FE
    1     0     1     0     0
SMALL
USER_CHARGES


@<TRIPOS>ATOM
      1 FE           0.0000     0.0000     0.0000 Fe        1 FE         3.000000
@<TRIPOS>BOND
@<TRIPOS>SUBSTRUCTURE
     1 FE          1 TEMP              0 ****  ****    0 ROOT
""",
        encoding="utf-8",
        newline="\n",
    )
    (mcpb_dir / "HOH.mol2").write_text(
        """@<TRIPOS>MOLECULE
HOH
    3     2     1     0     0
SMALL
USER_CHARGES


@<TRIPOS>ATOM
      1 O            0.0000     0.0000     0.0000 o         1 HOH       -0.834000
      2 H1           0.9572     0.0000     0.0000 h         1 HOH        0.417000
      3 H2          -0.2390     0.9266     0.0000 h         1 HOH        0.417000
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
@<TRIPOS>SUBSTRUCTURE
     1 HOH         1 TEMP              0 ****  ****    0 ROOT
""",
        encoding="utf-8",
        newline="\n",
    )
    step1 = mcpb_dir / "01_run_mcpb_step1.sh"
    step1.write_text(
        """#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
set +u
source /mnt/l/WSL/conda_envs/AmberTools25/amber.sh
set -u
MCPB.py -i mcpb.in -s 1 | tee logs/mcpb_step1.log
""",
        encoding="utf-8",
        newline="\n",
    )
    step1.chmod(0o755)

    post = mcpb_dir / "03_after_qm_run_mcpb_steps.sh"
    post.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
set +u
source /mnt/l/WSL/conda_envs/AmberTools25/amber.sh
set -u
MCPB.py -i mcpb.in -s 2 --logf ../02_qm_hessian/logs/FEIII_radical_M7_small_fc.log | tee logs/mcpb_step2.log
MCPB.py -i mcpb.in -s 3 --logf gamess_logs/{group_name}_large_mk.log | tee logs/mcpb_step3.log
MCPB.py -i mcpb.in -s 4 | tee logs/mcpb_step4.log
""",
        encoding="utf-8",
        newline="\n",
    )
    post.chmod(0o755)


def find_metal_serial(pdb: Path, *, chain: str, resseq: int) -> int:
    matches: list[int] = []
    for line in read_text(pdb).splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom = line[12:16].strip().upper()
        resn = line[17:20].strip().upper()
        atom_chain = line[21:22].strip()
        try:
            atom_resseq = int(line[22:26])
            serial = int(line[6:11])
        except ValueError:
            continue
        element = (line[76:78].strip() or atom[:2]).upper()
        if atom_chain == chain and atom_resseq == resseq and (atom == "FE" or element == "FE" or resn.startswith("FE")):
            matches.append(serial)
    if len(matches) != 1:
        raise RuntimeError(f"expected exactly one Fe at chain {chain} resseq {resseq} in {pdb}, found {matches}")
    return matches[0]


def add_m7_water_hydrogens_to_full_pdb(*, full_pdb: Path, small_model_pdb: Path, out_pdb: Path) -> None:
    """Copy M7 HOH B875 H1/H2 coordinates into the full MCPB PDB."""
    water_h_lines: list[str] = []
    for line in read_text(small_model_pdb).splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resn = line[17:20].strip().upper()
        chain = line[21:22].strip()
        resseq = line[22:26].strip()
        atom = line[12:16].strip().upper()
        if resn == "HOH" and chain == "B" and resseq == "875" and atom in {"H1", "H2"}:
            water_h_lines.append(line)
    if len(water_h_lines) != 2:
        raise RuntimeError(f"expected HOH B875 H1/H2 in {small_model_pdb}, found {len(water_h_lines)}")

    lines = read_text(full_pdb).splitlines()
    max_serial = max(int(line[6:11]) for line in lines if line.startswith(("ATOM", "HETATM")))
    output: list[str] = []
    inserted = False
    max_serial_current = max_serial
    for line in lines:
        if line.startswith("END"):
            continue
        if line.startswith(("ATOM", "HETATM")):
            chain0 = line[21:22].strip()
            resseq0 = line[22:26].strip()
            resn0 = line[17:20].strip().upper()
            if chain0 == "B" and resseq0 == "875" and resn0 == "HOH":
                line = line[:21] + "A" + line[22:]
            elif chain0 not in {"A", ""}:
                continue
            if line[21:22].strip() == "A" and line[22:26].strip() in {"187", "270"} and line[17:20].strip().upper() == "HIS":
                line = line[:17] + "HID" + line[20:]
            if line[21:22].strip() == "A" and line[22:26].strip() == "431" and line[12:16].strip().upper() == "FE":
                line = line[:17] + "FE " + line[20:]
        output.append(line)
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resn = line[17:20].strip().upper()
        chain = line[21:22].strip()
        resseq = line[22:26].strip()
        atom = line[12:16].strip().upper()
        if resn == "HOH" and chain == "A" and resseq == "875" and atom == "O":
            for h_line in water_h_lines:
                max_serial += 1
                atom_name = h_line[12:16]
                xyz = h_line[30:54]
                output.append(f"HETATM{max_serial:5d} {atom_name} HOH A 875    {xyz}  1.00  0.00           H  ")
            inserted = True
    if not inserted:
        raise RuntimeError(f"could not find remapped HOH A875 O in {full_pdb}")
    max_serial_current = max(max_serial_current, max_serial)
    output, max_serial_current = add_missing_ha_atoms(output, start_serial=max_serial_current + 1)
    output.append("END")
    out_pdb.write_text("\n".join(output) + "\n", encoding="utf-8", newline="\n")


def add_missing_ha_atoms(lines: list[str], *, start_serial: int) -> tuple[list[str], int]:
    """Add approximate HA atoms for MCPB sidechain capping when the PDB is heavy-only."""
    target_residues = {("A", "187"), ("A", "270")}
    residue_atoms: dict[tuple[str, str], dict[str, tuple[float, float, float]]] = {}
    existing: set[tuple[str, str, str]] = set()
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = line[21:22].strip()
        resseq = line[22:26].strip()
        atom = line[12:16].strip()
        if (chain, resseq) in target_residues:
            existing.add((chain, resseq, atom))
            residue_atoms.setdefault((chain, resseq), {})[atom] = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )

    additions: dict[tuple[str, str], str] = {}
    serial = start_serial
    for key in sorted(target_residues):
        chain, resseq = key
        if (chain, resseq, "HA") in existing:
            continue
        atoms = residue_atoms.get(key, {})
        if not {"CA", "N", "C", "CB"}.issubset(atoms):
            continue
        ca = atoms["CA"]
        direction = (0.0, 0.0, 0.0)
        for neighbor in ("N", "C", "CB"):
            vec = tuple(ca[i] - atoms[neighbor][i] for i in range(3))
            norm = math.sqrt(sum(v * v for v in vec)) or 1.0
            direction = tuple(direction[i] + vec[i] / norm for i in range(3))
        norm = math.sqrt(sum(v * v for v in direction)) or 1.0
        xyz = tuple(ca[i] + 1.09 * direction[i] / norm for i in range(3))
        additions[key] = (
            f"ATOM  {serial:5d}  HA  HID {chain}{int(resseq):4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00           H  "
        )
        serial += 1

    if not additions:
        return lines, serial - 1

    output: list[str] = []
    for line in lines:
        output.append(line)
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = line[21:22].strip()
        resseq = line[22:26].strip()
        atom = line[12:16].strip()
        key = (chain, resseq)
        if atom == "CA" and key in additions:
            output.append(additions[key])
    return output, serial - 1


def write_workflow_readme(out: Path, *, meta: dict[str, str], source_root: Path, mcpb_charge: int, mult: int) -> None:
    out.write_text(
        f"""# TJ M7 Fe(III)-Radical MCPB -> 400 ns MD

This is the clean handoff directory for the current TJ mechanism-oriented model.
Old exploratory folders remain in `../feiii_radical_preflight/` as provenance;
this directory keeps only inputs and scripts needed to continue.

## Chemical Model

- Fe model: Fe(III), high-spin Fe plus radical, selected multiplicity `M7`.
- Coordination sphere: His187 NE2, His270 NE2, Glu349 OE1, ACT501 O2, UNL N1, HOH875 O.
- Protein force field target: Amber `ff19SB`.
- ACT/UNL ordinary bonded/nonbonded ligand terms: `GAFF2` from AmberTools preflight.
- HOH875 is treated as the sixth coordinating ligand; its H1/H2 coordinates are copied from the M7 small model into the MCPB PDB.
- Fe-donor bonded force field: MCPB.py from the M7-converged QM model.
- QM level used so far: GAMESS B3LYP/6-31G(d,p)-style input.
- Final MCPB model charge/multiplicity: `charge={mcpb_charge}`, `MULT={mult}`.

Note: the earlier M7 spin screen used a custom 70-atom truncated model with
`charge=0`. MCPB.py builds its own 62-atom small model and 92-atom large model;
for those MCPB-generated models, `charge=+1` gives even electron counts and is
therefore the consistent setting for `MULT=7`.

## M7 Converged Geometry

- Source: `{source_root.as_posix()}`
- Energy: `{meta.get('energy', '')}` Hartree
- S-squared: `{meta.get('s_squared', '')}`
- Gradient max/RMS: `{meta.get('grad_max', '')}` / `{meta.get('grad_rms', '')}`

## Directory Map

```text
01_inputs/
  hplusplus/       H++ input/output provenance
  ligands/         ACT/UNL GAFF2 preflight mol2/frcmod
  structures/      initial, M7-patched, and visual comparison files
02_qm_hessian/     M7 small-model Hessian GAMESS input and runner
03_mcpb/           MCPB.py input, Fe mol2, ACT/UNL mol2/frcmod, step scripts
04_amber/          reserved for tleap dry/solvated Amber topology
05_gromacs/        reserved for ACPYPE/GROMACS conversion and 400 ns run
06_analysis/       reserved for PBC-fixed trajectory, RDC dashboard, PyMOL views
logs/              project-level logs
```

## Next Commands

Run the expensive M7 Hessian:

```bash
cd /mnt/e/TJ/mcpb_m7_400ns/02_qm_hessian
env NCORES=16 bash run_hessian_gamess.sh
tail -f logs/FEIII_radical_M7_small_fc.log
```

After Hessian finishes normally, generate MCPB step 1/large-model input:

```bash
cd /mnt/e/TJ/mcpb_m7_400ns/03_mcpb
bash 01_run_mcpb_step1.sh
```

Then run the MCPB large-model RESP/MK GAMESS job generated by step 1. When it
finishes, run:

```bash
cd /mnt/e/TJ/mcpb_m7_400ns/03_mcpb
bash 03_after_qm_run_mcpb_steps.sh
```

Only after MCPB step 4 succeeds should the 400 ns MD be launched.

## Important Caveat

This directory intentionally does not delete old TJ files. For a 400 ns run,
provenance matters: if the metal-site parameters are later questioned, we need
the original H++/M7/GAMESS logs.
""",
        encoding="utf-8",
        newline="\n",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare clean TJ M7 MCPB/400 ns project.")
    parser.add_argument("--root", type=Path, default=Path("/mnt/e/TJ"))
    parser.add_argument("--outdir", type=Path, default=Path("/mnt/e/TJ/mcpb_m7_400ns"))
    parser.add_argument("--cores", type=int, default=16)
    parser.add_argument("--gamess", type=Path, default=Path("/home/qin/softwares/gamess/rungms"))
    parser.add_argument("--mult", type=int, default=7)
    parser.add_argument(
        "--charge",
        type=int,
        default=1,
        help="Charge for MCPB-generated small/large models. Default +1 gives even electron counts for M7 in the MCPB model.",
    )
    args = parser.parse_args()

    root = args.root
    out = args.outdir
    src = root / "feiii_radical_preflight"
    m7_export = src / "m7_cont60_export"
    m7_log = src / "gamess_m7_cont60" / "logs" / "FEIII_radical_M7_cont60.log"

    dirs = [
        out / "01_inputs" / "hplusplus",
        out / "01_inputs" / "ligands",
        out / "01_inputs" / "structures",
        out / "02_qm_hessian" / "logs",
        out / "03_mcpb" / "logs",
        out / "04_amber",
        out / "05_gromacs",
        out / "06_analysis",
        out / "logs",
    ]
    for directory in dirs:
        directory.mkdir(parents=True, exist_ok=True)

    copies = [
        (root / "1T47_chainA_Hpp_try_keepFe_NTD_original_names.pkout.txt", out / "01_inputs" / "hplusplus" / "hplusplus.pkout.txt"),
        (root / "hplusplus_inputs" / "1T47_chainA_Hpp_try_keepFe_NTD_original_names.pdb", out / "01_inputs" / "hplusplus" / "hplusplus_input_keepFe_ligands.pdb"),
        (src / "complex_AFeIII_UNL_C4radical_dropH4A.pdb", out / "01_inputs" / "structures" / "initial_radical_complex.pdb"),
        (m7_export / "complex_AFeIII_UNL_C4radical_M7_metalcenter_patch.pdb", out / "01_inputs" / "structures" / "m7_metalcenter_patch_full_complex.pdb"),
        (m7_export / "M7_optimized_small_model_aligned_to_initial.pdb", out / "01_inputs" / "structures" / "m7_optimized_small_model.pdb"),
        (m7_export / "M7_initial_vs_optimized_compare.tsv", out / "01_inputs" / "structures" / "M7_initial_vs_optimized_compare.tsv"),
        (m7_export / "M7_export_report.md", out / "01_inputs" / "structures" / "M7_export_report.md"),
        (m7_export / "compare_initial_vs_M7.pml", out / "01_inputs" / "structures" / "compare_initial_vs_M7.pml"),
        (src / "amber_preflight" / "ACT.mol2", out / "01_inputs" / "ligands" / "ACT.mol2"),
        (src / "amber_preflight" / "ACT.frcmod", out / "01_inputs" / "ligands" / "ACT.frcmod"),
        (src / "amber_preflight" / "UNL_C4radical.mol2", out / "01_inputs" / "ligands" / "UNL.mol2"),
        (src / "amber_preflight" / "UNL_C4radical.frcmod", out / "01_inputs" / "ligands" / "UNL.frcmod"),
    ]
    for source, destination in copies:
        copy_required(source, destination)

    # MCPB working copies use the exact names referenced by mcpb.in.
    for name in ("ACT.mol2", "ACT.frcmod", "UNL.mol2", "UNL.frcmod"):
        copy_required(out / "01_inputs" / "ligands" / name, out / "03_mcpb" / name)
    add_m7_water_hydrogens_to_full_pdb(
        full_pdb=out / "01_inputs" / "structures" / "m7_metalcenter_patch_full_complex.pdb",
        small_model_pdb=out / "01_inputs" / "structures" / "m7_optimized_small_model.pdb",
        out_pdb=out / "03_mcpb" / "mcpb_original.pdb",
    )

    coords, meta = parse_gamess_final_coords(m7_log)
    job = "FEIII_radical_M7_small_fc"
    write_gamess_hessian_input(out / "02_qm_hessian" / f"{job}.inp", coords, charge=args.charge, mult=args.mult)
    write_gamess_runner(out / "02_qm_hessian", job, cores=args.cores, gamess=args.gamess)
    write_mcpb_files(out / "03_mcpb", group_name="TJ_M7_FE", sm_charge=args.charge, mult=args.mult)
    write_workflow_readme(out / "README.md", meta=meta, source_root=src, mcpb_charge=args.charge, mult=args.mult)

    print(f"project={out}")
    print(f"hessian_input={out / '02_qm_hessian' / (job + '.inp')}")
    print(f"hessian_runner={out / '02_qm_hessian' / 'run_hessian_gamess.sh'}")
    print(f"mcpb_input={out / '03_mcpb' / 'mcpb.in'}")
    print(f"readme={out / 'README.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
