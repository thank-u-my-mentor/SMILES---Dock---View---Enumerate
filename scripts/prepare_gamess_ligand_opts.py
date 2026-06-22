from __future__ import annotations

import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


ATOMIC_NUMBERS = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Fe": 26,
}


@dataclass
class Atom:
    serial: int
    name: str
    resn: str
    chain: str
    resi: str
    elem: str
    x: float
    y: float
    z: float


def parse_pdb(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        raw_elem = line[76:80].strip() if len(line) >= 78 else ""
        elem = "".join(ch for ch in raw_elem if ch.isalpha())
        if not elem:
            elem = "".join(ch for ch in line[12:16].strip() if ch.isalpha())[:1]
        elem = elem.capitalize()
        atoms.append(
            Atom(
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                resn=line[17:20].strip(),
                chain=line[21].strip(),
                resi=line[22:26].strip(),
                elem=elem,
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
            )
        )
    return atoms


def write_gamess_opt(path: Path, title: str, atoms: list[Atom], charge: int, mult: int) -> None:
    lines = [
        " $SYSTEM MWORDS=200 MEMDDI=50 $END",
        " $CONTRL",
        f"  SCFTYP={'RHF' if mult == 1 else 'UHF'} DFTTYP=B3LYP RUNTYP=OPTIMIZE",
        f"  ICHARG={charge} MULT={mult} MAXIT=200 COORD=CART UNITS=ANGS",
        " $END",
        " $STATPT NSTEP=300 OPTTOL=0.0005 $END",
        " $SCF DIRSCF=.T. DIIS=.T. DAMP=.T. ETHRSH=2.0 MAXDII=20 $END",
        # GAMESS notation for 6-31G(d,p): N31/6 plus one d on heavy atoms and one p on hydrogens.
        " $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=1 $END",
        " $DATA",
        title,
        "C1",
    ]
    for atom in atoms:
        atomic_number = ATOMIC_NUMBERS.get(atom.elem)
        if atomic_number is None:
            raise ValueError(f"unsupported element {atom.elem!r} for atom {atom}")
        lines.append(
            f"{atom.elem:<2s} {float(atomic_number):8.1f} "
            f"{atom.x:14.6f} {atom.y:14.6f} {atom.z:14.6f}"
        )
    lines += [" $END", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: prepare_gamess_ligand_opts.py input.pdb outdir")
        return 2

    pdb = Path(sys.argv[1])
    outdir = Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)

    residues: dict[tuple[str, str, str], list[Atom]] = defaultdict(list)
    for atom in parse_pdb(pdb):
        residues[(atom.resn, atom.chain or "-", atom.resi)].append(atom)

    jobs = []
    for keys, charge, mult, title in [
        (
            [
                ("ACT", "A", "501"),
                ("ACT", "-", "501"),
                ("ACT", "-", "1"),
                ("ACT", "A", "1"),
                ("ACE", "-", "1"),
            ],
            -1,
            1,
            "ACT_current_acetate_B3LYP_631Gdp",
        ),
        ([("UNL", "A", "1")], 0, 2, "UNL_current_radical_B3LYP_631Gdp"),
    ]:
        key = next((candidate for candidate in keys if candidate in residues), keys[0])
        atoms = residues.get(key)
        if not atoms:
            print(f"missing residue candidates {keys}")
            continue
        job_name = title
        inp = outdir / f"{job_name}.inp"
        write_gamess_opt(inp, title, atoms, charge, mult)
        jobs.append(job_name)
        print(f"wrote={inp} atoms={len(atoms)} charge={charge} mult={mult}")

    run = outdir / "run_ligand_gamess_opts.sh"
    with run.open("w", encoding="utf-8") as handle:
        handle.write("#!/usr/bin/env bash\n")
        handle.write("set -euo pipefail\n")
        handle.write("cd \"$(dirname \"$0\")\"\n")
        handle.write("mkdir -p logs\n")
        handle.write("GAMESS=${GAMESS:-/home/qin/softwares/gamess/rungms}\n")
        handle.write("VERSION=${VERSION:-00}\n")
        handle.write("NCORES=${NCORES:-4}\n")
        for job in jobs:
            handle.write(f"echo '[GAMESS] {job}'\n")
            handle.write(f"\"$GAMESS\" {job} \"$VERSION\" \"$NCORES\" > logs/{job}.log 2>&1\n")
    run.chmod(0o755)
    print(f"run_script={run}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
