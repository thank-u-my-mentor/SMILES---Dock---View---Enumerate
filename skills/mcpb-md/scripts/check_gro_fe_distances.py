#!/usr/bin/env python3
from __future__ import annotations

import math
import sys
from pathlib import Path


def read_gro(path: Path):
    lines = path.read_text(errors="replace").splitlines()
    atoms = []
    for line in lines[2:-1]:
        atoms.append(
            (
                int(line[0:5]),
                line[5:10].strip(),
                line[10:15].strip(),
                int(line[15:20]),
                (float(line[20:28]), float(line[28:36]), float(line[36:44])),
            )
        )
    return atoms


def only(atoms, resname: str, atomname: str):
    matches = [a for a in atoms if a[1] == resname and a[2] == atomname]
    if len(matches) != 1:
        raise SystemExit(f"{resname} {atomname}: expected 1 atom, found {len(matches)}")
    return matches[0]


def dist(a, b) -> float:
    return math.sqrt(sum((a[4][i] - b[4][i]) ** 2 for i in range(3))) * 10.0


def main() -> int:
    gro = Path(sys.argv[1])
    atoms = read_gro(gro)
    fe = only(atoms, "FE1", "FE")
    for label, res, atom in [
        ("UL1 N1", "UL1", "N1"),
        ("HD1 NE2", "HD1", "NE2"),
        ("HD2 NE2", "HD2", "NE2"),
        ("GU1 OE1", "GU1", "OE1"),
        ("AT1 O2", "AT1", "O2"),
        ("HH1 O", "HH1", "O"),
    ]:
        print(f"Fe-{label}: {dist(fe, only(atoms, res, atom)):.3f} A")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
