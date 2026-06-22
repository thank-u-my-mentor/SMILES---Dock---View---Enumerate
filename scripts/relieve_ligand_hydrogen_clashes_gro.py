#!/usr/bin/env python3
"""Relieve severe ligand-protein hydrogen clashes in a GRO coordinate file.

This is a pre-minimization repair for handoffs where a docked heavy-atom pose is
reasonable, but added hydrogens overlap badly. It only moves hydrogens, keeping
their parent heavy-atom bond length.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Atom:
    line_index: int
    resid: int
    resname: str
    name: str
    index: int
    x: float
    y: float
    z: float

    @property
    def xyz(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)

    def set_xyz(self, xyz: tuple[float, float, float]) -> None:
        self.x, self.y, self.z = xyz


def parse_gro(path: Path) -> tuple[list[str], list[Atom]]:
    lines = path.read_text().splitlines()
    atoms: list[Atom] = []
    for i, line in enumerate(lines[2:-1], start=2):
        if len(line) < 44:
            continue
        atoms.append(
            Atom(
                line_index=i,
                resid=int(line[0:5]),
                resname=line[5:10].strip(),
                name=line[10:15].strip(),
                index=int(line[15:20]),
                x=float(line[20:28]),
                y=float(line[28:36]),
                z=float(line[36:44]),
            )
        )
    return lines, atoms


def format_gro_line(original: str, atom: Atom) -> str:
    prefix = original[:20]
    suffix = original[44:]
    return f"{prefix}{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}{suffix}"


def parse_top_bonds(path: Path) -> dict[int, set[int]]:
    bonds: dict[int, set[int]] = {}
    section = ""
    for raw in path.read_text().splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith(";") or stripped.startswith("#"):
            continue
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped.strip("[]").strip()
            continue
        if section != "bonds":
            continue
        parts = stripped.split(";", 1)[0].split()
        if len(parts) < 2:
            continue
        try:
            a = int(parts[0])
            b = int(parts[1])
        except ValueError:
            continue
        bonds.setdefault(a, set()).add(b)
        bonds.setdefault(b, set()).add(a)
    return bonds


def distance(a: Atom, b: Atom) -> float:
    return math.dist(a.xyz, b.xyz)


def is_h(atom: Atom) -> bool:
    return atom.name.upper().startswith("H")


def vec_sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vec_add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vec_scale(a: tuple[float, float, float], scale: float) -> tuple[float, float, float]:
    return (a[0] * scale, a[1] * scale, a[2] * scale)


def norm(a: tuple[float, float, float]) -> float:
    return math.sqrt(a[0] ** 2 + a[1] ** 2 + a[2] ** 2)


def unit(a: tuple[float, float, float]) -> tuple[float, float, float]:
    n = norm(a)
    if n < 1.0e-8:
        return (1.0, 0.0, 0.0)
    return (a[0] / n, a[1] / n, a[2] / n)


def choose_parent(atom: Atom, by_index: dict[int, Atom], bonds: dict[int, set[int]]) -> Atom | None:
    for bonded_index in bonds.get(atom.index, set()):
        other = by_index.get(bonded_index)
        if other and not is_h(other):
            return other
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gro", type=Path, required=True)
    parser.add_argument("--top", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--ligand-resname", default="LIG")
    parser.add_argument("--hh-threshold-nm", type=float, default=0.07)
    parser.add_argument("--h-heavy-threshold-nm", type=float, default=0.105)
    parser.add_argument("--passes", type=int, default=8)
    args = parser.parse_args()

    lines, atoms = parse_gro(args.gro)
    by_index = {atom.index: atom for atom in atoms}
    bonds = parse_top_bonds(args.top)
    ligand = [a for a in atoms if a.resname == args.ligand_resname]
    non_ligand = [a for a in atoms if a.resname != args.ligand_resname and a.resname not in {"WAT", "NA+", "CL-"}]

    moved: list[str] = []
    for _ in range(args.passes):
        clashes: list[tuple[float, Atom, Atom, float]] = []
        for a in ligand:
            for b in non_ligand:
                if not (is_h(a) or is_h(b)):
                    continue
                d = distance(a, b)
                threshold = args.hh_threshold_nm if is_h(a) and is_h(b) else args.h_heavy_threshold_nm
                if d < threshold:
                    clashes.append((d, a, b, threshold))
        if not clashes:
            break
        d, a, b, threshold = sorted(clashes, key=lambda item: item[0])[0]
        mobile = a if is_h(a) else b
        fixed = b if mobile is a else a
        parent = choose_parent(mobile, by_index, bonds)
        if parent is None:
            moved.append(f"skip atom={mobile.index} no_parent clash={d:.4f}")
            break
        bond_len = distance(mobile, parent)
        direction = unit(vec_sub(parent.xyz, fixed.xyz))
        if norm(direction) < 1.0e-8:
            direction = unit(vec_sub(mobile.xyz, parent.xyz))
        old = mobile.xyz
        mobile.set_xyz(vec_add(parent.xyz, vec_scale(direction, bond_len)))
        moved.append(
            f"moved atom={mobile.index} {mobile.resname}{mobile.resid}:{mobile.name} "
            f"parent={parent.index}:{parent.name} away_from={fixed.index}:{fixed.resname}{fixed.resid}:{fixed.name} "
            f"clash_nm={d:.4f} target_threshold_nm={threshold:.4f} "
            f"old=({old[0]:.3f},{old[1]:.3f},{old[2]:.3f}) new=({mobile.x:.3f},{mobile.y:.3f},{mobile.z:.3f})"
        )

    for atom in atoms:
        lines[atom.line_index] = format_gro_line(lines[atom.line_index], atom)
    args.out.write_text("\n".join(lines) + "\n")

    print(f"input={args.gro}")
    print(f"output={args.out}")
    for item in moved:
        print(item)
    print(f"moved_count={sum(1 for item in moved if item.startswith('moved '))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
