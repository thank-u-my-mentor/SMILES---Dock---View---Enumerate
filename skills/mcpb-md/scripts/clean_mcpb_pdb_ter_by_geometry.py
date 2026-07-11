#!/usr/bin/env python3
from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path


PROTEIN_NAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLN", "GLU", "GLY", "HID", "HIE", "HIP", "HIS",
    "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "ASH", "GLH", "CYM", "HD1", "HD2", "GU1",
}


def parse_atom(line: str):
    return {
        "name": line[12:16].strip(),
        "resname": line[17:20].strip(),
        "chain": line[21].strip(),
        "resid": int(line[22:26]),
        "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
        "line": line,
    }


def dist(a, b) -> float:
    return math.sqrt(sum((a["xyz"][i] - b["xyz"][i]) ** 2 for i in range(3)))


def residue_key(atom):
    return atom["chain"], atom["resid"], atom["resname"]


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: clean_mcpb_pdb_ter_by_geometry.py input.pdb output.pdb")
    inp = Path(sys.argv[1])
    out = Path(sys.argv[2])
    lines = inp.read_text().splitlines()

    atoms_by_res = defaultdict(list)
    atom_records = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            atom = parse_atom(line)
            atom_records.append(atom)
            atoms_by_res[residue_key(atom)].append(atom)

    removed = 0
    kept = 0
    decisions = []
    new_lines = []
    for i, line in enumerate(lines):
        if not line.startswith("TER"):
            new_lines.append(line)
            continue

        prev_atom = None
        next_atom = None
        for j in range(i - 1, -1, -1):
            if lines[j].startswith(("ATOM", "HETATM")):
                prev_atom = parse_atom(lines[j])
                break
        for j in range(i + 1, len(lines)):
            if lines[j].startswith(("ATOM", "HETATM")):
                next_atom = parse_atom(lines[j])
                break

        remove = False
        reason = "kept"
        if prev_atom and next_atom:
            pk = residue_key(prev_atom)
            nk = residue_key(next_atom)
            prev_res = atoms_by_res[pk]
            next_res = atoms_by_res[nk]
            prev_c = [a for a in prev_res if a["name"] == "C"]
            next_n = [a for a in next_res if a["name"] == "N"]
            if (
                prev_atom["resname"] in PROTEIN_NAMES
                and next_atom["resname"] in PROTEIN_NAMES
                and prev_atom["chain"] == next_atom["chain"]
                and prev_c
                and next_n
            ):
                cn = dist(prev_c[0], next_n[0])
                if cn <= 1.75:
                    remove = True
                    reason = f"removed peptide-like C-N {cn:.3f} A"
                else:
                    reason = f"kept chain break C-N {cn:.3f} A"

        if remove:
            removed += 1
        else:
            kept += 1
            new_lines.append(line)
        if prev_atom and next_atom:
            decisions.append(
                f"{reason}: {prev_atom['resname']}{prev_atom['resid']} -> {next_atom['resname']}{next_atom['resid']}"
            )
        else:
            decisions.append(reason)

    out.write_text("\n".join(new_lines) + "\n")
    report = out.with_suffix(out.suffix + ".ter_report.txt")
    report.write_text(
        f"input={inp}\noutput={out}\nremoved={removed}\nkept={kept}\n\n"
        + "\n".join(decisions)
        + "\n"
    )
    print(f"removed TER={removed} kept TER={kept}")
    print(f"report={report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
