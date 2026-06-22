from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


DONOR_ELEMENTS = {"N", "O", "S"}


@dataclass(frozen=True)
class Atom:
    line: str
    serial: int
    name: str
    resn: str
    chain: str
    resi: str
    elem: str
    x: float
    y: float
    z: float


def parse_atom(line: str) -> Atom:
    raw_elem = line[76:80].strip() if len(line) >= 78 else ""
    elem = "".join(ch for ch in raw_elem if ch.isalpha())
    if not elem:
        elem = "".join(ch for ch in line[12:16].strip() if ch.isalpha())[:1]
    if elem.upper() == "FE":
        elem = "Fe"
    else:
        elem = elem.capitalize()
    return Atom(
        line=(line + " " * 80)[:80],
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


def distance(a: Atom, b: Atom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def atom_id(atom: Atom) -> str:
    chain = atom.chain or "-"
    return f"{atom.resn}:{chain}:{atom.resi}:{atom.name}"


def read_pdb(path: Path) -> tuple[list[str], dict[int, Atom], list[tuple[int, int]]]:
    lines = path.read_text(errors="replace").splitlines()
    atoms: dict[int, Atom] = {}
    conect: set[tuple[int, int]] = set()
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            atom = parse_atom(line)
            atoms[atom.serial] = atom
        elif line.startswith("CONECT"):
            parts = line.split()
            if len(parts) >= 3:
                try:
                    src = int(parts[1])
                except ValueError:
                    continue
                for token in parts[2:]:
                    try:
                        dst = int(token)
                    except ValueError:
                        continue
                    if src != dst:
                        conect.add(tuple(sorted((src, dst))))
    return lines, atoms, sorted(conect)


def write_filtered_pdb(
    out: Path,
    lines: list[str],
    keep_serials: set[int],
    remove_serials: set[int] | None = None,
) -> None:
    remove_serials = remove_serials or set()
    kept = set(keep_serials) - set(remove_serials)
    output: list[str] = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            serial = int(line[6:11])
            if serial in kept:
                output.append((line + " " * 80)[:80])
        elif line.startswith("CONECT"):
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                src = int(parts[1])
            except ValueError:
                continue
            partners: list[int] = []
            for token in parts[2:]:
                try:
                    dst = int(token)
                except ValueError:
                    continue
                if src in kept and dst in kept:
                    partners.append(dst)
            if partners:
                output.append("CONECT" + f"{src:5d}" + "".join(f"{dst:5d}" for dst in partners))
        elif keep_serials == set():
            output.append(line)
    if output and not output[-1].startswith("END"):
        output.append("END")
    out.write_text("\n".join(output) + "\n", encoding="utf-8")


def write_complex_without_atoms(out: Path, lines: list[str], remove_serials: set[int]) -> None:
    output: list[str] = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            serial = int(line[6:11])
            if serial not in remove_serials:
                output.append((line + " " * 80)[:80])
        elif line.startswith("CONECT"):
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                src = int(parts[1])
            except ValueError:
                continue
            if src in remove_serials:
                continue
            partners: list[int] = []
            for token in parts[2:]:
                try:
                    dst = int(token)
                except ValueError:
                    continue
                if dst not in remove_serials:
                    partners.append(dst)
            if partners:
                output.append("CONECT" + f"{src:5d}" + "".join(f"{dst:5d}" for dst in partners))
        else:
            output.append(line)
    out.write_text("\n".join(output).rstrip() + "\n", encoding="utf-8")


def find_residue(atoms: dict[int, Atom], resn: str, chain: str, resi: str) -> list[Atom]:
    return [
        atom
        for atom in atoms.values()
        if atom.resn == resn and atom.chain == chain and atom.resi == resi
    ]


def write_amber_preflight(outdir: Path) -> None:
    script = outdir / "run_amber_ligand_preflight.sh"
    script.write_text(
        """#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
set +u
source /mnt/l/WSL/conda_envs/AmberTools25/amber.sh
set -u

echo "[antechamber] ACT acetate, charge -1, singlet"
antechamber -i ../ligands/ACT_A501.pdb -fi pdb -o ACT.mol2 -fo mol2 \
  -rn ACT -c bcc -s 2 -at gaff2 -nc -1 -m 1 \
  > logs/antechamber_ACT.log 2>&1
parmchk2 -i ACT.mol2 -f mol2 -o ACT.frcmod -s gaff2 \
  > logs/parmchk2_ACT.log 2>&1

echo "[antechamber] UNL C4 radical, charge 0, doublet"
antechamber -i ../ligands/UNL_A1_C4radical.pdb -fi pdb -o UNL_C4radical.mol2 -fo mol2 \
  -rn UNL -c bcc -s 2 -at gaff2 -nc 0 -m 2 \
  > logs/antechamber_UNL_C4radical.log 2>&1
parmchk2 -i UNL_C4radical.mol2 -f mol2 -o UNL_C4radical.frcmod -s gaff2 \
  > logs/parmchk2_UNL_C4radical.log 2>&1

cat > tleap_ligand_check.in <<'LEAP'
source leaprc.gaff2
loadamberparams ACT.frcmod
loadamberparams UNL_C4radical.frcmod
ACT = loadmol2 ACT.mol2
UNL = loadmol2 UNL_C4radical.mol2
check ACT
check UNL
saveamberparm ACT ACT.prmtop ACT.inpcrd
saveamberparm UNL UNL_C4radical.prmtop UNL_C4radical.inpcrd
quit
LEAP

tleap -f tleap_ligand_check.in > logs/tleap_ligand_check.log 2>&1
echo "[ok] Amber ligand preflight completed"
""",
        encoding="utf-8",
    )
    script.chmod(0o755)


def write_readme(outdir: Path, removed_atom: Atom, metal: Atom, donor_rows: list[tuple[float, Atom]]) -> None:
    donor_text = "\n".join(
        f"- {dist_value:.3f} A  {atom_id(atom)}  element={atom.elem}" for dist_value, atom in donor_rows[:8]
    )
    (outdir / "README_next_steps.md").write_text(
        f"""# TJ Fe(III)-Radical Preflight

This directory is a cheap preflight before the expensive Fe MCPB/GAMESS workflow.

Generated changes:

- Removed `{atom_id(removed_atom)}` to make the current UNL model an odd-electron C4/benzyl-side radical candidate.
- Kept the original protein/Fe/ligand coordinates otherwise unchanged.
- Focused on `{atom_id(metal)}` as the chain-A Fe site.

Nearest donor atoms around the selected Fe:

{donor_text}

Run ligand-format preflight:

```bash
cd {outdir.as_posix()}/amber_preflight
bash run_amber_ligand_preflight.sh
```

Interpretation:

- Success means ACT/UNL names, atom graph, BCC charges, and ligand-only tleap input are at least commandable.
- It does not validate Fe coordination chemistry.
- For production metal MD, build a 6-coordinate MCPB.py model and screen Fe(III)+radical spin states, usually MULT=5 and MULT=7, before paying for the final Hessian.
""",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare a TJ Fe(III)-radical preflight directory from the hand-built active-site PDB."
    )
    parser.add_argument("--pdb", default="/mnt/e/TJ/0616_initial_guess_ACT.pdb")
    parser.add_argument("--outdir", default="/mnt/e/TJ/feiii_radical_preflight")
    parser.add_argument("--metal-chain", default="A")
    parser.add_argument("--metal-resi", default="431")
    parser.add_argument("--radical-resn", default="UNL")
    parser.add_argument("--radical-chain", default="A")
    parser.add_argument("--radical-resi", default="1")
    parser.add_argument("--drop-hydrogen", default="H4A")
    args = parser.parse_args()

    pdb = Path(args.pdb)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "ligands").mkdir(exist_ok=True)
    (outdir / "amber_preflight").mkdir(exist_ok=True)
    (outdir / "metal_site").mkdir(exist_ok=True)

    lines, atoms, _ = read_pdb(pdb)
    metal_candidates = [
        atom
        for atom in atoms.values()
        if atom.chain == args.metal_chain
        and atom.resi == args.metal_resi
        and (atom.elem == "Fe" or atom.resn.upper().startswith("FE"))
    ]
    if len(metal_candidates) != 1:
        raise SystemExit(f"Expected one selected Fe, found {len(metal_candidates)}")
    metal = metal_candidates[0]

    radical_atoms = find_residue(atoms, args.radical_resn, args.radical_chain, args.radical_resi)
    drop_atoms = [atom for atom in radical_atoms if atom.name == args.drop_hydrogen]
    if len(drop_atoms) != 1:
        raise SystemExit(f"Expected one radical hydrogen {args.drop_hydrogen}, found {len(drop_atoms)}")
    drop_atom = drop_atoms[0]

    complex_out = outdir / "complex_AFeIII_UNL_C4radical_dropH4A.pdb"
    write_complex_without_atoms(complex_out, lines, {drop_atom.serial})

    act_atoms = find_residue(atoms, "ACT", "A", "501")
    if not act_atoms:
        raise SystemExit("ACT A501 not found")
    write_filtered_pdb(outdir / "ligands" / "ACT_A501.pdb", lines, {a.serial for a in act_atoms})
    write_filtered_pdb(outdir / "ligands" / "UNL_A1_original.pdb", lines, {a.serial for a in radical_atoms})
    write_filtered_pdb(
        outdir / "ligands" / "UNL_A1_C4radical.pdb",
        lines,
        {a.serial for a in radical_atoms},
        {drop_atom.serial},
    )

    donor_rows = [
        (distance(metal, atom), atom)
        for atom in atoms.values()
        if atom.serial != metal.serial and atom.elem in DONOR_ELEMENTS and distance(metal, atom) <= 3.2
    ]
    donor_rows.sort(key=lambda row: row[0])
    with (outdir / "metal_site" / "chainA_fe_donors.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["distance_A", "serial", "resn", "chain", "resi", "atom", "element"])
        for dist_value, atom in donor_rows:
            writer.writerow([f"{dist_value:.3f}", atom.serial, atom.resn, atom.chain, atom.resi, atom.name, atom.elem])

    write_amber_preflight(outdir / "amber_preflight")
    write_readme(outdir, drop_atom, metal, donor_rows)

    print(f"outdir={outdir}")
    print(f"complex={complex_out}")
    print(f"removed={atom_id(drop_atom)} serial={drop_atom.serial}")
    print(f"metal={atom_id(metal)} serial={metal.serial}")
    print(f"donor_table={outdir / 'metal_site' / 'chainA_fe_donors.tsv'}")
    print(f"amber_preflight={outdir / 'amber_preflight' / 'run_amber_ligand_preflight.sh'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
