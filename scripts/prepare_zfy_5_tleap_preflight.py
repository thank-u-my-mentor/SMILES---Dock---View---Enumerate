from __future__ import annotations

import argparse
import shutil
from pathlib import Path


UNK_RENAMES = {
    "H01": "H02",
    "H02": "H03",
    "H04": "H04",
    "H06": "H08",
    "H08": "H09",
    "H09": "H10",
    "H10": "H11",
    "H11": "H12",
    "H12": "H13",
}


def parse_pdb_atoms(pdb: Path):
    atoms = []
    for line in pdb.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atoms.append(
            {
                "line": line,
                "record": line[:6],
                "serial": int(line[6:11]),
                "name": line[12:16].strip(),
                "resn": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resi": line[22:26].strip(),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            }
        )
    return atoms


def replace_atom_name(line: str, name: str) -> str:
    return line[:12] + f"{name:>4}" + line[16:]


def write_pdb_no_conect_with_renames(src: Path, out: Path) -> dict[tuple[str, str, str], tuple[float, float, float]]:
    coords: dict[tuple[str, str, str], tuple[float, float, float]] = {}
    lines = []
    for line in src.read_text(errors="ignore").splitlines():
        if line.startswith("CONECT"):
            continue
        if line.startswith(("ATOM  ", "HETATM")):
            name = line[12:16].strip()
            resn = line[17:20].strip()
            chain = line[21:22].strip()
            resi = line[22:26].strip()
            if resn == "UNK" and chain == "A" and resi == "501":
                name = UNK_RENAMES.get(name, name)
                line = replace_atom_name(line, name)
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords[(resn, resi, name)] = (x, y, z)
        lines.append(line)
    out.write_text("\n".join(lines) + "\n")
    return coords


def rewrite_mol2_coords(template: Path, out: Path, coords: dict[tuple[str, str, str], tuple[float, float, float]], resn: str, resi: str) -> None:
    output = []
    in_atom = False
    for line in template.read_text(errors="ignore").splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom = True
            output.append(line)
            continue
        if line.startswith("@<TRIPOS>"):
            in_atom = False
            output.append(line)
            continue
        if in_atom and line.strip():
            parts = line.split()
            if len(parts) >= 9:
                atom_name = parts[1]
                key = (resn, resi, atom_name)
                if key in coords:
                    x, y, z = coords[key]
                    line = (
                        f"{int(parts[0]):7d} {atom_name:<8s}"
                        f"{x:10.4f}{y:10.4f}{z:10.4f} "
                        f"{parts[5]:<8s}{int(parts[6]):4d} {parts[7]:<8s}"
                        f"{float(parts[8]):11.6f}"
                    )
        output.append(line)
    out.write_text("\n".join(output) + "\n")


def write_fe_mol2(out: Path, coords: dict[tuple[str, str, str], tuple[float, float, float]]) -> None:
    x, y, z = coords[("FE", "4113", "FE")]
    out.write_text(
        "@<TRIPOS>MOLECULE\n"
        "FE\n"
        " 1 0 0 0 0\n"
        "SMALL\n"
        "USER_CHARGES\n\n"
        "@<TRIPOS>ATOM\n"
        f"      1 FE        {x:10.4f}{y:10.4f}{z:10.4f} Fe3+     4113 FE         3.000000\n"
        "@<TRIPOS>BOND\n"
    )


def write_tleap(out: Path) -> None:
    out.write_text(
        "source leaprc.protein.ff19SB\n"
        "source leaprc.gaff2\n"
        "source leaprc.water.tip3p\n"
        "loadamberparams frcmod.ions234lm_126_tip3p\n"
        "loadamberparams UNK.frcmod\n"
        "loadamberparams ACT.frcmod\n"
        "UNK = loadmol2 UNK.mol2\n"
        "ACT = loadmol2 ACT.mol2\n"
        "FE = loadmol2 FE.mol2\n"
        "x = loadpdb 5_for_tleap.pdb\n"
        "check x\n"
        "charge x\n"
        "saveamberparm x preflight.prmtop preflight.inpcrd\n"
        "quit\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--unk-template", required=True)
    parser.add_argument("--act-template", required=True)
    parser.add_argument("--unk-frcmod", required=True)
    parser.add_argument("--act-frcmod", required=True)
    args = parser.parse_args()

    pdb = Path(args.pdb)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    coords = write_pdb_no_conect_with_renames(pdb, outdir / "5_for_tleap.pdb")
    rewrite_mol2_coords(Path(args.unk_template), outdir / "UNK.mol2", coords, "UNK", "501")
    rewrite_mol2_coords(Path(args.act_template), outdir / "ACT.mol2", coords, "ACT", "502")
    write_fe_mol2(outdir / "FE.mol2", coords)
    shutil.copy2(args.unk_frcmod, outdir / "UNK.frcmod")
    shutil.copy2(args.act_frcmod, outdir / "ACT.frcmod")
    write_tleap(outdir / "tleap.in")
    print(f"outdir={outdir}")
    print(f"pdb={outdir / '5_for_tleap.pdb'}")
    print(f"unk_mol2={outdir / 'UNK.mol2'}")
    print(f"act_mol2={outdir / 'ACT.mol2'}")
    print(f"tleap={outdir / 'tleap.in'}")


if __name__ == "__main__":
    main()
