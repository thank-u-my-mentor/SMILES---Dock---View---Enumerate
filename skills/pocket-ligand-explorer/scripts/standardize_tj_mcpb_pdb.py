#!/usr/bin/env python3
"""Normalize TJ MCPB intermediate PDB residue names for a fresh MCPB run.

This keeps coordinates unchanged, but maps prior MCPB custom names back to
human-readable input names:

HD1/HD2 -> HID
GU1 -> GLU
UL1 -> UNL
AT1 -> ACT
HH1 -> HOH
FE1 -> FE

It also fixes the element column for the Fe atom, which can be accidentally
written as F when the residue name is FE1.
"""

from __future__ import annotations

import argparse
from pathlib import Path


RESNAME_MAP = {
    "HD1": "HID",
    "HD2": "HID",
    "GU1": "GLU",
    "UL1": "UNL",
    "AT1": "ACT",
    "HH1": "HOH",
    "FE1": "FE",
}


def patch_atom_line(line: str) -> str:
    padded = (line.rstrip("\n") + " " * 80)[:80]
    resname = padded[17:20].strip()
    atom_name = padded[12:16].strip()
    new_resname = RESNAME_MAP.get(resname, resname)
    out = padded[:17] + f"{new_resname:>3s}" + padded[20:]
    if atom_name.upper() == "FE" or new_resname.upper() == "FE":
        out = out[:76] + "FE" + out[78:]
    return out.rstrip()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    lines = []
    for line in args.pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            lines.append(patch_atom_line(line))
        else:
            lines.append(line.rstrip())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8", newline="\n")
    print(f"out={args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
