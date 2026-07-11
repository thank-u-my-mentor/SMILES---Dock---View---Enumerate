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
            {
                "resid": int(line[0:5]),
                "resname": line[5:10].strip(),
                "atom": line[10:15].strip(),
                "atomid": int(line[15:20]),
                "xyz": (float(line[20:28]), float(line[28:36]), float(line[36:44])),
                "line": line,
            }
        )
    return atoms


def dist(a, b):
    return math.sqrt(sum((a["xyz"][i] - b["xyz"][i]) ** 2 for i in range(3)))


def main() -> int:
    gro = Path(sys.argv[1])
    atomid = int(sys.argv[2])
    atoms = read_gro(gro)
    target = next(a for a in atoms if a["atomid"] == atomid)
    print("target:", target["line"])
    near = []
    for atom in atoms:
        if atom is target:
            continue
        d = dist(target, atom)
        near.append((d, atom))
    near.sort(key=lambda x: x[0])
    print("nearest atoms:")
    for d, atom in near[:30]:
        print(
            f"{d*10:8.4f} A  id={atom['atomid']:5d} "
            f"{atom['resid']:5d} {atom['resname']:<5s} {atom['atom']:<5s} {atom['line']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
