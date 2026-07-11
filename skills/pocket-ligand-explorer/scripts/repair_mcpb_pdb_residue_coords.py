#!/usr/bin/env python3
"""Repair selected MCPB PDB residue coordinates from a reference PDB.

Use this when a QM small-model coordinate patch distorted a full protein
residue enough to break covalent geometry, but the MCPB-derived force constants
and custom residue names should be retained.  The script keeps the target PDB
records and residue names, and only replaces coordinates for matching atom names
in user-selected residues.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_rule(text: str) -> tuple[str, int]:
    """Parse CHAIN:RESSEQ or RESSEQ."""
    parts = text.split(":")
    if len(parts) == 1:
        return "", int(parts[0])
    if len(parts) == 2:
        return parts[0], int(parts[1])
    raise argparse.ArgumentTypeError(f"Invalid residue rule: {text}")


def atom_key(line: str) -> tuple[str, int, str]:
    chain = line[21:22].strip()
    resseq = int(line[22:26])
    atom = line[12:16].strip()
    return chain, resseq, atom


def collect_reference_coords(path: Path) -> dict[tuple[str, int, str], tuple[float, float, float]]:
    coords: dict[tuple[str, int, str], tuple[float, float, float]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            key = atom_key(line)
            coords[key] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
    return coords


def with_coords(line: str, xyz: tuple[float, float, float]) -> str:
    return f"{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}"


def vec_sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vec_add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vec_scale(a: tuple[float, float, float], scale: float) -> tuple[float, float, float]:
    return (a[0] * scale, a[1] * scale, a[2] * scale)


def norm(a: tuple[float, float, float]) -> float:
    return (a[0] ** 2 + a[1] ** 2 + a[2] ** 2) ** 0.5


def residue_coords(
    path: Path,
) -> dict[tuple[str, int], dict[str, tuple[float, float, float]]]:
    residues: dict[tuple[str, int], dict[str, tuple[float, float, float]]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            chain, resseq, atom = atom_key(line)
            residues.setdefault((chain, resseq), {})[atom] = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
        except ValueError:
            continue
    return residues


def build_clamp_shifts(
    target_pdb: Path,
    reference_pdb: Path,
    replace: set[tuple[str, int]],
) -> dict[tuple[str, int], tuple[float, float, float]]:
    target = residue_coords(target_pdb)
    reference = residue_coords(reference_pdb)
    shifts: dict[tuple[str, int], tuple[float, float, float]] = {}
    for key in replace:
        chain, resseq = key
        if not chain:
            matches = [candidate for candidate in target if candidate[1] == resseq]
            if len(matches) != 1:
                raise SystemExit(f"Cannot infer chain for residue {resseq}; matches={matches}")
            key = matches[0]
        t = target.get(key, {})
        r = reference.get(key, {})
        if not {"CA", "CB"}.issubset(t) or not {"CA", "CB"}.issubset(r):
            raise SystemExit(f"Missing CA/CB in target or reference for {key}")
        target_vec = vec_sub(t["CB"], t["CA"])
        target_len = norm(target_vec)
        ref_len = norm(vec_sub(r["CB"], r["CA"]))
        if target_len <= 0.01:
            raise SystemExit(f"Bad CA-CB vector in target for {key}")
        new_cb = vec_add(t["CA"], vec_scale(target_vec, ref_len / target_len))
        shifts[key] = vec_sub(new_cb, t["CB"])
        print(
            f"clamp={key[0]}:{key[1]} target_CA_CB={target_len:.3f} "
            f"reference_CA_CB={ref_len:.3f} shift=({shifts[key][0]:.3f},"
            f"{shifts[key][1]:.3f},{shifts[key][2]:.3f})"
        )
    return shifts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-pdb", type=Path, required=True)
    parser.add_argument("--reference-pdb", type=Path, required=True)
    parser.add_argument("--out-pdb", type=Path, required=True)
    parser.add_argument(
        "--replace-residue",
        action="append",
        type=parse_rule,
        required=True,
        help="Residue to repair as CHAIN:RESSEQ, e.g. A:187. Repeatable.",
    )
    parser.add_argument(
        "--clamp-ca-cb",
        action="store_true",
        help=(
            "Instead of replacing all selected atom coordinates, preserve the "
            "target side-chain direction and translate CB/ring atoms so CA-CB "
            "matches the reference length."
        ),
    )
    args = parser.parse_args()

    reference = collect_reference_coords(args.reference_pdb)
    replace = set(args.replace_residue)
    clamp_shifts = build_clamp_shifts(args.target_pdb, args.reference_pdb, replace) if args.clamp_ca_cb else {}
    sidechain_atoms = {"CB", "CG", "CD", "CD1", "CD2", "CE", "CE1", "CE2", "ND1", "NE2", "OE1", "OE2"}
    changed = 0
    missing: list[str] = []
    out_lines: list[str] = []
    for line in args.target_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                chain, resseq, atom = atom_key(line)
            except ValueError:
                out_lines.append(line)
                continue
            if (chain, resseq) in replace or ("", resseq) in replace:
                if args.clamp_ca_cb:
                    shift = clamp_shifts.get((chain, resseq))
                    if shift and atom in sidechain_atoms:
                        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                        line = with_coords(line, vec_add(xyz, shift))
                        changed += 1
                else:
                    ref_key = (chain, resseq, atom)
                    if ref_key in reference:
                        line = with_coords(line, reference[ref_key])
                        changed += 1
                    else:
                        missing.append(f"{chain}:{resseq}:{atom}")
        out_lines.append(line)

    args.out_pdb.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    print(f"target={args.target_pdb}")
    print(f"reference={args.reference_pdb}")
    print(f"out={args.out_pdb}")
    print(f"changed_atoms={changed}")
    if missing:
        print("missing_reference_atoms=" + ",".join(missing))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
