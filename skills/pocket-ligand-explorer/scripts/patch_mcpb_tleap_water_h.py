#!/usr/bin/env python3
"""Patch MCPB tleap inputs when a coordinating water uses lowercase h.

MCPB.py can emit a non-standard coordinating-water residue such as HH1 where
water hydrogens are atom type ``h``.  Amber ff19SB/GAFF2 usually has no vdW entry
for that exact lowercase type, so tleap stops at saveamberparm.  This script adds
a tiny frcmod entry for that water-hydrogen type and inserts it into a tleap
input without changing coordinates, RESP charges, or MCPB Fe force constants.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


PATCH_TEXT = """Remark: TIP3P-style internal patch for MCPB coordinating-water hydrogen type h
MASS
h    1.008    0.000    MCPB coordinating-water hydrogen

BOND
Y5-h   553.0    0.9572      TIP3P O-H for MCPB coordinating water
h -h   553.0    1.5136      TIP3P H-H distance for MCPB coordinating water

ANGLE
h -Y5-h   100.0    104.52   TIP3P H-O-H for MCPB coordinating water

DIHE

IMPROPER

NONBON
h    0.0000   0.0000
"""


def patch_tleap_input(
    source: Path,
    dest: Path,
    frcmod_name: str,
    loadpdb_name: str | None,
    output_prefix: str | None,
) -> None:
    text = source.read_text()
    load_line = f"loadamberparams {frcmod_name}"

    lines = text.splitlines()
    loadpdb_at = None
    for idx, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("mol = loadpdb"):
            loadpdb_at = idx
            if loadpdb_name:
                lines[idx] = f"mol = loadpdb {loadpdb_name}"
            break
    if loadpdb_at is None:
        raise SystemExit(f"Could not find 'mol = loadpdb' in {source}")

    if load_line not in lines:
        lines.insert(loadpdb_at, load_line)

    if output_prefix:
        replaced: list[str] = []
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("savepdb mol "):
                suffix = "dry" if "_dry" in stripped else "solv"
                replaced.append(f"savepdb mol {output_prefix}_{suffix}.pdb")
            elif stripped.startswith("saveamberparm mol "):
                suffix = "dry" if "_dry" in stripped else "solv"
                replaced.append(
                    f"saveamberparm mol {output_prefix}_{suffix}.prmtop "
                    f"{output_prefix}_{suffix}.inpcrd"
                )
            else:
                replaced.append(line)
        lines = replaced

    dest.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amber-dir", required=True, help="Directory containing the tleap input.")
    parser.add_argument(
        "--tleap-in",
        default="leap_tj_mcpb_fixed.in",
        help="Existing tleap input to patch.",
    )
    parser.add_argument(
        "--tleap-out",
        default="leap_tj_mcpb_fixed_hvdw.in",
        help="Patched tleap input to write.",
    )
    parser.add_argument(
        "--frcmod-name",
        default="HH1_h_vdw.frcmod",
        help="Local frcmod patch filename to write and load.",
    )
    parser.add_argument(
        "--loadpdb",
        default=None,
        help="Optional PDB filename to use in the 'mol = loadpdb' line.",
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Optional prefix for savepdb/saveamberparm outputs.",
    )
    args = parser.parse_args()

    amber_dir = Path(args.amber_dir)
    tleap_in = amber_dir / args.tleap_in
    tleap_out = amber_dir / args.tleap_out
    frcmod = amber_dir / args.frcmod_name

    if not tleap_in.exists():
        raise SystemExit(f"Missing tleap input: {tleap_in}")

    frcmod.write_text(PATCH_TEXT)
    patch_tleap_input(tleap_in, tleap_out, args.frcmod_name, args.loadpdb, args.output_prefix)
    print(f"wrote={frcmod}")
    print(f"wrote={tleap_out}")


if __name__ == "__main__":
    main()
