#!/usr/bin/env python3
"""Prepare a compact MCPB.py small-optimization project for Fe-radical models."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
from dataclasses import dataclass
from pathlib import Path


DONOR_ELEMENTS = {"N", "O", "S"}

HOH_MOL2 = """@<TRIPOS>MOLECULE
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
"""


@dataclass(frozen=True)
class Atom:
    line: str
    serial: int
    name: str
    resn: str
    chain: str
    resseq: str
    elem: str
    x: float
    y: float
    z: float


def parse_atom(line: str) -> Atom:
    raw_elem = line[76:78].strip() if len(line) >= 78 else ""
    elem = "".join(ch for ch in raw_elem if ch.isalpha())
    if not elem:
        elem = "".join(ch for ch in line[12:16].strip() if ch.isalpha())[:1]
    elem = "Fe" if elem.upper() == "FE" else elem.capitalize()
    return Atom(
        line=(line + " " * 80)[:80],
        serial=int(line[6:11]),
        name=line[12:16].strip(),
        resn=line[17:20].strip(),
        chain=line[21:22].strip(),
        resseq=line[22:26].strip(),
        elem=elem,
        x=float(line[30:38]),
        y=float(line[38:46]),
        z=float(line[46:54]),
    )


def parse_pdb(path: Path) -> tuple[list[str], list[Atom]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    atoms = [parse_atom(line) for line in lines if line.startswith(("ATOM", "HETATM"))]
    return lines, atoms


def distance(a: Atom, b: Atom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def find_metal(atoms: list[Atom], *, chain: str, resseq: str, atom_name: str) -> Atom:
    matches = [
        atom
        for atom in atoms
        if atom.chain == chain
        and atom.resseq == str(resseq)
        and (atom.name.upper() == atom_name.upper() or atom.elem.upper() == "FE")
    ]
    if len(matches) != 1:
        labels = [f"{atom.serial}:{atom.resn}:{atom.chain}:{atom.resseq}:{atom.name}" for atom in matches]
        raise SystemExit(f"Expected exactly one metal at {chain}:{resseq}:{atom_name}, found {labels}")
    return matches[0]


def donor_rows(atoms: list[Atom], metal: Atom, cutoff: float) -> list[tuple[float, Atom]]:
    rows = [
        (distance(metal, atom), atom)
        for atom in atoms
        if atom.serial != metal.serial and atom.elem.upper() in DONOR_ELEMENTS and distance(metal, atom) <= cutoff
    ]
    rows.sort(key=lambda row: row[0])
    return rows


def copy_required(src: Path, dst: Path) -> None:
    if not src.exists():
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def write_fe_mol2(path: Path, charge: float) -> None:
    path.write_text(
        f"""@<TRIPOS>MOLECULE
FE
    1     0     1     0     0
SMALL
USER_CHARGES


@<TRIPOS>ATOM
      1 FE           0.0000     0.0000     0.0000 Fe        1 FE         {charge:.6f}
@<TRIPOS>BOND
@<TRIPOS>SUBSTRUCTURE
     1 FE          1 TEMP              0 ****  ****    0 ROOT
""",
        encoding="utf-8",
        newline="\n",
    )


def write_hoh_mol2(path: Path) -> None:
    path.write_text(HOH_MOL2, encoding="utf-8", newline="\n")


def write_mcpb_in(
    path: Path,
    *,
    group: str,
    metal_serial: int,
    cutoff: float,
    charge: int,
    mult: int,
    ligand_names: list[str],
    naa_names: list[str],
) -> None:
    frcmods = " ".join(f"{name}.frcmod" for name in ligand_names)
    mol2s = " ".join(f"{name}.mol2" for name in naa_names)
    path.write_text(
        f"""original_pdb mcpb_original.pdb
group_name {group}
cut_off {cutoff}
ion_ids {metal_serial}
ion_mol2files FE.mol2
water_model TIP3P
force_field ff19SB
gaff 2
frcmod_files {frcmods}
software_version gms
large_opt 1
add_redcrd 1
scale_factor 1.0
smmodel_chg {charge}
smmodel_spin {mult}
lgmodel_chg {charge}
lgmodel_spin {mult}
naa_mol2files {mol2s}
""",
        encoding="utf-8",
        newline="\n",
    )


def write_shell_scripts(mcpb_dir: Path, *, group: str, cores: int, session: str, gamess: str) -> None:
    (mcpb_dir / "run_01_mcpb_step1.sh").write_text(
        """#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p logs gamess_logs
set +u
source /mnt/l/WSL/conda_envs/AmberTools25/amber.sh
set -u
MCPB.py -i mcpb.in -s 1 | tee logs/mcpb_step1.log
""",
        encoding="utf-8",
        newline="\n",
    )
    (mcpb_dir / "run_01_mcpb_step1.sh").chmod(0o755)

    (mcpb_dir / "run_02_prepare_gamess_inputs.sh").write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
/mnt/l/WSL/softwares/conda_envs/md/bin/python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/write_mcpb_gamess_pipeline.py \\
  --mcpb-dir "$PWD" \\
  --group {group} \\
  --gamess {gamess} \\
  --cores {cores}
echo prepared={group}_small_opt.inp
""",
        encoding="utf-8",
        newline="\n",
    )
    (mcpb_dir / "run_02_prepare_gamess_inputs.sh").chmod(0o755)

    (mcpb_dir / "run_03_small_opt_foreground.sh").write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={group}_small_opt
NCORES=${{NCORES:-{cores}}}
GAMESS=${{GAMESS:-{gamess}}}
VERSION=${{VERSION:-00}}
test -s "${{JOB}}.inp" || {{ echo "missing=${{JOB}}.inp"; exit 2; }}
mkdir -p gamess_logs logs
echo "[GAMESS] $JOB start $(date) NCORES=$NCORES" | tee "gamess_logs/${{JOB}}.status"
set +e
"$GAMESS" "$JOB" "$VERSION" "$NCORES" > "gamess_logs/${{JOB}}.log" 2>&1
rc=$?
set -e
echo "[GAMESS] $JOB exit_code=$rc end $(date)" | tee -a "gamess_logs/${{JOB}}.status"
if grep -q 'EXECUTION OF GAMESS TERMINATED NORMALLY' "gamess_logs/${{JOB}}.log"; then
  echo "[GAMESS] $JOB status=normal" | tee -a "gamess_logs/${{JOB}}.status"
else
  echo "[GAMESS] $JOB status=not_normal" | tee -a "gamess_logs/${{JOB}}.status"
fi
exit "$rc"
""",
        encoding="utf-8",
        newline="\n",
    )
    (mcpb_dir / "run_03_small_opt_foreground.sh").chmod(0o755)

    (mcpb_dir / "run_03_small_opt_tmux.sh").write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={group}_small_opt
SESSION=${{SESSION:-{session}}}
NCORES=${{NCORES:-{cores}}}
GAMESS=${{GAMESS:-{gamess}}}
test -s "${{JOB}}.inp" || {{ echo "missing=${{JOB}}.inp"; exit 2; }}
mkdir -p gamess_logs logs
if tmux has-session -t "$SESSION" 2>/dev/null; then
  echo "tmux_session_exists=$SESSION"
  echo "monitor=cd $PWD && bash monitor_small_opt.sh"
  exit 0
fi
tmux new-session -d -s "$SESSION" "cd '$PWD'; env NCORES='$NCORES' GAMESS='$GAMESS' bash run_03_small_opt_foreground.sh; echo; echo '[tmux] finished'; exec bash"
echo "tmux_session=$SESSION"
echo "monitor=cd $PWD && bash monitor_small_opt.sh"
echo "tail=tail -f $PWD/gamess_logs/${{JOB}}.log"
""",
        encoding="utf-8",
        newline="\n",
    )
    (mcpb_dir / "run_03_small_opt_tmux.sh").chmod(0o755)

    (mcpb_dir / "monitor_small_opt.sh").write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={group}_small_opt
echo "===== status ====="
cat "gamess_logs/${{JOB}}.status" 2>/dev/null || true
echo
echo "===== recent log ====="
tail -n 80 "gamess_logs/${{JOB}}.log" 2>/dev/null || true
echo
echo "===== search/energy markers ====="
grep -E 'RUNTYP=|NSERCH:|TOTAL ENERGY|S-SQUARED|END OF GEOMETRY SEARCH|TERMINATED|ERROR|FAILURE' "gamess_logs/${{JOB}}.log" 2>/dev/null | tail -n 80 || true
""",
        encoding="utf-8",
        newline="\n",
    )
    (mcpb_dir / "monitor_small_opt.sh").chmod(0o755)


def write_readme(outdir: Path, *, pdb: Path, group: str, metal: Atom, donors: list[tuple[float, Atom]], charge: int, mult: int) -> None:
    donor_text = "\n".join(
        f"- {dist_value:.3f} A  {atom.resn} {atom.chain}{atom.resseq} {atom.name} serial={atom.serial}"
        for dist_value, atom in donors
    )
    (outdir / "README.md").write_text(
        f"""# ZFY MCPB Small-Model Optimization

This directory prepares the Hessian-before step for the hand-built 2R5V-derived
ZFY Fe-radical intermediate.

## Input

- Clean scaffold: `{pdb.as_posix()}`
- Metal: `{metal.resn} {metal.chain}{metal.resseq} {metal.name}`, serial `{metal.serial}`
- MCPB group: `{group}`
- Model assumption: Fe(III)-radical, charge `{charge}`, multiplicity `{mult}`

## Current Fe donor candidates within cutoff

{donor_text}

Inspect this donor list before launching a long GAMESS job. Coordinating HOH
is automatically added to `naa_mol2files` as TIP3P water when it appears within
the cutoff; it is not added to `frcmod_files`.

## Commands

```bash
cd {outdir.as_posix()}/mcpb
bash run_01_mcpb_step1.sh
bash run_02_prepare_gamess_inputs.sh
env NCORES=16 bash run_03_small_opt_tmux.sh
tail -f gamess_logs/{group}_small_opt.log
```

This run performs only the GAMESS geometry optimization of MCPB's small model.
Do not launch Hessian until this log ends normally and the Fe-donor distances
are inspected.
""",
        encoding="utf-8",
        newline="\n",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare MCPB small-opt project for Fe-radical intermediates.")
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--group", default="ZFY_M7_FE")
    parser.add_argument("--metal-chain", default="A")
    parser.add_argument("--metal-resseq", default="4113")
    parser.add_argument("--metal-atom", default="FE")
    parser.add_argument("--cutoff", type=float, default=2.8)
    parser.add_argument("--charge", type=int, default=1)
    parser.add_argument("--mult", type=int, default=7)
    parser.add_argument("--cores", type=int, default=16)
    parser.add_argument("--session", default="zfy_mcpb_m7_opt")
    parser.add_argument("--gamess", default="/home/qin/softwares/gamess/rungms")
    parser.add_argument("--fe-charge", type=float, default=3.0)
    parser.add_argument(
        "--ligand",
        action="append",
        nargs=3,
        metavar=("NAME", "MOL2", "FRCMOD"),
        required=True,
        help="Nonstandard residue files to copy into MCPB dir, e.g. --ligand UNK UNK.mol2 UNK.frcmod.",
    )
    args = parser.parse_args()

    _, atoms = parse_pdb(args.pdb)
    metal = find_metal(atoms, chain=args.metal_chain, resseq=args.metal_resseq, atom_name=args.metal_atom)
    donors = donor_rows(atoms, metal, args.cutoff)

    outdir = args.outdir
    mcpb_dir = outdir / "mcpb"
    input_dir = outdir / "input"
    mcpb_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(exist_ok=True)

    copy_required(args.pdb, input_dir / args.pdb.name)
    copy_required(args.pdb, mcpb_dir / "mcpb_original.pdb")

    ligand_names: list[str] = []
    for name, mol2, frcmod in args.ligand:
        ligand_names.append(name)
        copy_required(Path(mol2), input_dir / f"{name}.mol2")
        copy_required(Path(frcmod), input_dir / f"{name}.frcmod")
        copy_required(Path(mol2), mcpb_dir / f"{name}.mol2")
        copy_required(Path(frcmod), mcpb_dir / f"{name}.frcmod")

    naa_names = list(ligand_names)
    donor_resnames = {atom.resn.upper() for _, atom in donors}
    if "HOH" in donor_resnames and "HOH" not in {name.upper() for name in naa_names}:
        naa_names.append("HOH")
        write_hoh_mol2(input_dir / "HOH.mol2")
        write_hoh_mol2(mcpb_dir / "HOH.mol2")

    write_fe_mol2(mcpb_dir / "FE.mol2", charge=args.fe_charge)
    write_mcpb_in(
        mcpb_dir / "mcpb.in",
        group=args.group,
        metal_serial=metal.serial,
        cutoff=args.cutoff,
        charge=args.charge,
        mult=args.mult,
        ligand_names=ligand_names,
        naa_names=naa_names,
    )
    write_shell_scripts(mcpb_dir, group=args.group, cores=args.cores, session=args.session, gamess=args.gamess)
    write_readme(outdir, pdb=args.pdb, group=args.group, metal=metal, donors=donors, charge=args.charge, mult=args.mult)

    with (outdir / "fe_donor_candidates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["distance_A", "serial", "resname", "chain", "resseq", "atom", "element"])
        for dist_value, atom in donors:
            writer.writerow([f"{dist_value:.3f}", atom.serial, atom.resn, atom.chain, atom.resseq, atom.name, atom.elem])

    print(f"project={outdir}")
    print(f"mcpb_dir={mcpb_dir}")
    print(f"mcpb_input={mcpb_dir / 'mcpb.in'}")
    print(f"metal_serial={metal.serial}")
    print(f"donor_table={outdir / 'fe_donor_candidates.csv'}")
    print(f"step1={mcpb_dir / 'run_01_mcpb_step1.sh'}")
    print(f"small_opt_starter={mcpb_dir / 'run_03_small_opt_tmux.sh'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
