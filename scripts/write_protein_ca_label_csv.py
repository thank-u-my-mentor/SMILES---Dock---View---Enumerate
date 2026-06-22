#!/usr/bin/env python3
"""Write residue labels from a reference PDB C-alpha order."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


AA3 = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "ASH",
    "CYS",
    "CYX",
    "GLN",
    "GLU",
    "GLH",
    "GLY",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "HD1",
    "HD2",
    "HD3",
    "ILE",
    "LEU",
    "LYS",
    "LYN",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    rows = []
    seen = set()
    for line in args.pdb.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[12:16].strip() != "CA":
            continue
        resname = line[17:20].strip()
        if resname not in AA3:
            continue
        chain = line[21].strip() or "A"
        resid = int(line[22:26])
        icode = line[26].strip()
        key = (chain, resid, icode, resname)
        if key in seen:
            continue
        seen.add(key)
        label = f"{resname}:{chain}:{resid}{icode}"
        rows.append(
            {
                "series": len(rows),
                "residue": label,
                "resname": resname,
                "segid": chain,
                "resid": resid,
            }
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["series", "residue", "resname", "segid", "resid"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"labels={len(rows)} out={args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
