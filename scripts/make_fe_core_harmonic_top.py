#!/usr/bin/env python
"""Create a GROMACS topology with harmonic Fe-core coordination bonds."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-top", required=True, type=Path)
    parser.add_argument("--output-top", required=True, type=Path)
    parser.add_argument(
        "--bond",
        action="append",
        required=True,
        help="Bond as ai,aj,r_nm,k_kj_mol_nm2,comment",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    bonds: list[str] = []
    for item in args.bond:
        ai, aj, r_nm, k, *comment = item.split(",", 4)
        note = comment[0] if comment else ""
        bonds.append(f"{int(ai):6d} {int(aj):6d}   1    {float(r_nm):.4e}    {float(k):.4e} ; {note}".rstrip())

    lines = args.input_top.read_text(encoding="utf-8").splitlines()
    out: list[str] = []
    inserted = False
    in_bonds = False
    for line in lines:
        out.append(line)
        if not inserted and line.strip() == "[ bonds ]":
            in_bonds = True
            continue
        if in_bonds and not inserted and line.strip().startswith(";"):
            out.append("; Fe-core harmonic coordination bonds, positive-control only; not full QM-derived MCPB")
            out.extend(bonds)
            inserted = True
            in_bonds = False
    if not inserted:
        raise SystemExit("Did not find the first [ bonds ] header/comment location")
    args.output_top.write_text("\n".join(out) + "\n", encoding="utf-8")
    print(args.output_top)


if __name__ == "__main__":
    main()
