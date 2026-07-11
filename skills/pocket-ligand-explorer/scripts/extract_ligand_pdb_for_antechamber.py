#!/usr/bin/env python3
"""Extract one ligand residue from a complex PDB with ligand-only CONECT records."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--resname", required=True)
    parser.add_argument("--chain", default=None)
    parser.add_argument("--resseq", default=None)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    keep_serials: set[int] = set()
    atom_lines: list[str] = []
    for line in args.pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip()
        chain = line[21:22].strip()
        resseq = line[22:26].strip()
        if resname != args.resname:
            continue
        if args.chain is not None and chain != args.chain:
            continue
        if args.resseq is not None and resseq != str(args.resseq):
            continue
        serial = int(line[6:11])
        keep_serials.add(serial)
        atom_lines.append("HETATM" + line[6:].rstrip())

    conect_lines: list[str] = []
    for line in args.pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("CONECT"):
            continue
        parts = line.split()
        try:
            serials = [int(part) for part in parts[1:]]
        except ValueError:
            continue
        if not serials or serials[0] not in keep_serials:
            continue
        internal = [serials[0]] + [serial for serial in serials[1:] if serial in keep_serials]
        if len(internal) >= 2:
            conect_lines.append("CONECT" + "".join(f"{serial:5d}" for serial in internal))

    if not atom_lines:
        raise SystemExit(f"No atoms found for {args.resname} {args.chain or ''} {args.resseq or ''}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(atom_lines + conect_lines + ["END"]) + "\n", encoding="utf-8", newline="\n")
    print(f"out={args.out}")
    print(f"atoms={len(atom_lines)}")
    print(f"conect={len(conect_lines)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
