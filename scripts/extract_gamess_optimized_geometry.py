from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path


BOHR_TO_ANG = 0.529177210903


@dataclass
class TemplateAtom:
    serial: int
    name: str
    resn: str
    chain: str
    resi: str
    elem: str


def parse_template_atoms(path: Path, resn: str) -> list[TemplateAtom]:
    atoms: list[TemplateAtom] = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip() != resn:
            continue
        elem = (line[76:78].strip() or re.sub(r"[^A-Za-z]", "", line[12:16].strip())[:1]).capitalize()
        atoms.append(
            TemplateAtom(
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                resn=line[17:20].strip(),
                chain=line[21].strip() or "A",
                resi=line[22:26].strip() or "1",
                elem=elem,
            )
        )
    return atoms


def extract_last_standard_orientation(log_path: Path) -> list[tuple[str, float, float, float]]:
    lines = log_path.read_text(errors="replace").splitlines()
    blocks: list[list[tuple[str, float, float, float]]] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if "COORDINATES OF ALL ATOMS ARE" not in line and "ATOM      ATOMIC                      COORDINATES" not in line:
            i += 1
            continue

        # GAMESS blocks usually have a few header lines, then rows:
        # C           6.0  x y z
        block: list[tuple[str, float, float, float]] = []
        i += 1
        while i < len(lines):
            row = lines[i].strip()
            if not row:
                if block:
                    break
                i += 1
                continue
            parts = row.split()
            if len(parts) >= 5 and re.fullmatch(r"[A-Za-z]{1,2}", parts[0]):
                try:
                    float(parts[1])
                    x, y, z = map(float, parts[2:5])
                except ValueError:
                    i += 1
                    continue
                block.append((parts[0].capitalize(), x, y, z))
            elif block:
                break
            i += 1
        if block:
            blocks.append(block)
    if not blocks:
        raise ValueError(f"no coordinate blocks found in {log_path}")

    coords = blocks[-1]
    # GAMESS may print optimized coords in Bohr in some sections. The current
    # input uses UNITS=ANGS; detect obviously-Bohr protein-scale coordinates.
    return coords


def write_xyz(path: Path, coords: list[tuple[str, float, float, float]], title: str) -> None:
    lines = [str(len(coords)), title]
    for elem, x, y, z in coords:
        lines.append(f"{elem:<2s} {x:14.6f} {y:14.6f} {z:14.6f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_pdb(path: Path, template: list[TemplateAtom], coords: list[tuple[str, float, float, float]]) -> None:
    if len(template) != len(coords):
        raise ValueError(f"template atom count {len(template)} != GAMESS atom count {len(coords)}")
    lines = []
    for i, (tmpl, (elem, x, y, z)) in enumerate(zip(template, coords), start=1):
        chain = (tmpl.chain or "A")[:1]
        try:
            resi = int(tmpl.resi)
        except ValueError:
            resi = 1
        lines.append(
            f"HETATM{i:5d} {tmpl.name:>4s} {tmpl.resn:>3s} {chain}{resi:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    if len(sys.argv) != 5:
        print("usage: extract_gamess_optimized_geometry.py log template_pdb RESN out_prefix")
        return 2

    log_path = Path(sys.argv[1])
    template_pdb = Path(sys.argv[2])
    resn = sys.argv[3]
    out_prefix = Path(sys.argv[4])
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    coords = extract_last_standard_orientation(log_path)
    template = parse_template_atoms(template_pdb, resn)
    write_xyz(out_prefix.with_suffix(".xyz"), coords, log_path.stem)
    write_pdb(out_prefix.with_suffix(".pdb"), template, coords)
    print(f"atoms={len(coords)}")
    print(f"xyz={out_prefix.with_suffix('.xyz')}")
    print(f"pdb={out_prefix.with_suffix('.pdb')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
