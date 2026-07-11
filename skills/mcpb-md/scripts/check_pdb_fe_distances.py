#!/usr/bin/env python3
from __future__ import annotations

import math
import sys
from pathlib import Path


DONORS = [
    ("UL1", "N1", "UL1 N1"),
    ("HD1", "NE2", "HD1 NE2"),
    ("HD2", "NE2", "HD2 NE2"),
    ("GU1", "OE1", "GU1 OE1"),
    ("AT1", "O2", "AT1 O2"),
    ("HH1", "O", "HH1 O"),
]


def parse_pdb(path: Path):
    atoms = []
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atoms.append(
            {
                "name": line[12:16].strip(),
                "res": line[17:20].strip(),
                "idx": int(line[6:11]),
                "xyz": (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ),
            }
        )
    return atoms


def find_atom(atoms, res: str, name: str):
    matches = [a for a in atoms if a["res"] == res and a["name"] == name]
    if not matches:
        raise SystemExit(f"missing atom {res} {name}")
    return matches[0]


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: check_pdb_fe_distances.py structure.pdb", file=sys.stderr)
        return 2

    atoms = parse_pdb(Path(sys.argv[1]))
    fe = find_atom(atoms, "FE1", "FE")
    for res, name, label in DONORS:
        atom = find_atom(atoms, res, name)
        d = math.dist(fe["xyz"], atom["xyz"])
        print(f"Fe-{label}: {d:.3f} A")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
