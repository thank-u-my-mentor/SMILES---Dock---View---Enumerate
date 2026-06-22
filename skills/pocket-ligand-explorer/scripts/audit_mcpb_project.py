#!/usr/bin/env python3
"""Audit an MCPB.py/GAMESS project for charge, donor, and SCF status."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


ATOMIC_NUMBERS = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "S": 16,
    "FE": 26,
}


def infer_element(atom_name: str, element_field: str = "") -> str:
    raw = element_field.strip()
    if raw:
        return raw.upper().replace("+", "").replace("-", "")
    atom = atom_name.strip()
    if atom.upper().startswith("FE"):
        return "FE"
    letters = "".join(ch for ch in atom if ch.isalpha())
    if not letters:
        return ""
    if len(letters) >= 2 and letters[:2].upper() == "FE":
        return "FE"
    return letters[0].upper()


def parse_pdb(path: Path) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        atom_name = line[12:16].strip()
        element = infer_element(atom_name, line[76:78] if len(line) >= 78 else "")
        atoms.append(
            {
                "record": line[:6].strip(),
                "serial": int(line[6:11]) if line[6:11].strip().isdigit() else None,
                "atom": atom_name,
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resseq": int(line[22:26]) if line[22:26].strip().lstrip("-").isdigit() else None,
                "element": element,
                "xyz": (x, y, z),
                "line": line,
            }
        )
    return atoms


def dist(a: dict[str, object], b: dict[str, object]) -> float:
    ax, ay, az = a["xyz"]  # type: ignore[misc]
    bx, by, bz = b["xyz"]  # type: ignore[misc]
    return math.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)


def mol2_charge(path: Path) -> tuple[int, float]:
    in_atom = False
    count = 0
    charge = 0.0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom = True
            continue
        if line.startswith("@<TRIPOS>"):
            in_atom = False
            continue
        if not in_atom or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 9:
            count += 1
            try:
                charge += float(parts[8])
            except ValueError:
                pass
    return count, charge


def electron_count(atoms: list[dict[str, object]], charge: int) -> int:
    total = 0
    unknown: list[str] = []
    for atom in atoms:
        elem = str(atom["element"]).upper()
        if elem not in ATOMIC_NUMBERS:
            unknown.append(f"{atom['atom']}:{elem}")
            continue
        total += ATOMIC_NUMBERS[elem]
    if unknown:
        print("unknown_elements=" + ",".join(sorted(set(unknown))))
    return total - charge


def read_mcpb_charge_mult(path: Path) -> tuple[int | None, int | None, int | None, int | None]:
    text = path.read_text(encoding="utf-8", errors="replace")
    sm = re.search(r"ion_info\s+\S+\s+(\S+)\s+(\S+)", text)
    sm_chg = re.search(r"smmodel_chg\s+(-?\d+)", text)
    sm_spin = re.search(r"smmodel_spin\s+(\d+)", text)
    lg = re.search(r"lgmodel_chg\s+(-?\d+).*?\n.*?lgmodel_spin\s+(\d+)", text, re.S)
    small_charge = int(sm.group(1)) if sm else None
    small_spin = int(sm.group(2)) if sm else None
    if small_charge is None and sm_chg:
        small_charge = int(sm_chg.group(1))
    if small_spin is None and sm_spin:
        small_spin = int(sm_spin.group(1))
    large_charge = int(lg.group(1)) if lg else None
    large_spin = int(lg.group(2)) if lg else None
    if large_charge is None:
        m = re.search(r"lgmodel_chg\s+(-?\d+)", text)
        large_charge = int(m.group(1)) if m else None
    if large_spin is None:
        m = re.search(r"lgmodel_spin\s+(\d+)", text)
        large_spin = int(m.group(1)) if m else None
    return small_charge, small_spin, large_charge, large_spin


def log_status(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")
    pieces: list[str] = []
    if "EXECUTION OF GAMESS TERMINATED NORMALLY" in text:
        pieces.append("terminated_normally")
    if "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-" in text:
        pieces.append("terminated_abnormally")
    if "SCF IS UNCONVERGED" in text or "SCF HAS NOT CONVERGED" in text:
        pieces.append("scf_unconverged")
    if "END OF GEOMETRY SEARCH" in text:
        pieces.append("geometry_search_end")
    ss = re.findall(r"S-SQUARED\s*=\s*([0-9.]+)", text)
    if ss:
        pieces.append(f"last_s2={ss[-1]}")
    return ",".join(pieces) if pieces else "unknown"


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit MCPB/GAMESS project inputs and logs.")
    parser.add_argument("--mcpb-dir", type=Path, required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--metal-resname", default="FE")
    parser.add_argument("--metal-atom", default="FE")
    parser.add_argument("--cutoff", type=float, default=3.0)
    args = parser.parse_args()

    root = args.mcpb_dir
    mcpb_in = root / "mcpb.in"
    small_pdb = root / f"{args.group}_small.pdb"
    large_pdb = root / f"{args.group}_large.pdb"
    original_pdb = root / "mcpb_original.pdb"

    if mcpb_in.exists():
        small_charge, small_spin, large_charge, large_spin = read_mcpb_charge_mult(mcpb_in)
        print(f"mcpb_small_charge={small_charge} mcpb_small_spin={small_spin}")
        print(f"mcpb_large_charge={large_charge} mcpb_large_spin={large_spin}")
    else:
        small_charge = large_charge = None

    for pdb, label, charge in [
        (small_pdb, "small", small_charge),
        (large_pdb, "large", large_charge),
    ]:
        if pdb.exists() and charge is not None:
            atoms = parse_pdb(pdb)
            ne = electron_count(atoms, charge)
            parity = "even" if ne % 2 == 0 else "odd"
            print(f"{label}_atoms={len(atoms)} {label}_electrons={ne} {label}_parity={parity}")

    for mol2 in sorted(root.glob("*.mol2")):
        count, charge = mol2_charge(mol2)
        print(f"mol2={mol2.name} atoms={count} charge_sum={charge:.6f} nearest_integer={round(charge)}")

    if original_pdb.exists():
        atoms = parse_pdb(original_pdb)
        metals = [
            a
            for a in atoms
            if str(a["resname"]).upper() == args.metal_resname.upper()
            and str(a["atom"]).upper() == args.metal_atom.upper()
        ]
        if not metals:
            metals = [a for a in atoms if str(a["element"]).upper() == "FE"]
        print(f"metal_count={len(metals)}")
        for metal in metals[:5]:
            print(
                "metal="
                f"{metal['atom']} {metal['resname']} {metal['chain']}:{metal['resseq']} "
                f"xyz={metal['xyz']}"
            )
            nearby = []
            for atom in atoms:
                if atom is metal:
                    continue
                if str(atom["element"]).upper() not in {"N", "O", "S"}:
                    continue
                d = dist(metal, atom)
                if d <= args.cutoff:
                    nearby.append((d, atom))
            for d, atom in sorted(nearby, key=lambda x: x[0]):
                print(
                    f"donor {d:5.3f} A "
                    f"{atom['atom']} {atom['resname']} {atom['chain']}:{atom['resseq']} "
                    f"element={atom['element']}"
                )

    logs = sorted((root / "gamess_logs").glob(f"{args.group}*.log"))
    for log in logs:
        print(f"log={log.name} status={log_status(log)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
