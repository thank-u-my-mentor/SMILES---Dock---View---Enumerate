#!/usr/bin/env python3
"""Build MCPB mol2 files using small-M2 RESP charges with ChgModD-like fixing."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


RESIDUE_FILE_ORDER = [
    ("1-UNL", "UL1.mol2"),
    ("187-HID", "HD1.mol2"),
    ("270-HID", "HD2.mol2"),
    ("349-GLU", "GU1.mol2"),
    ("431-FE", "FE1.mol2"),
    ("501-ACT", "AT1.mol2"),
    ("875-HOH", "HH1.mol2"),
]

FIXED_AA_ATOMS = {"N", "H", "CA", "HA", "C", "O", "CB"}
STANDARD_AA = {"HID", "HIE", "HIP", "HIS", "GLU", "GLH", "ASP", "ASH", "CYS", "CYM"}


def read_charge_map(path: Path) -> dict[str, tuple[str, float, str]]:
    rows: dict[str, tuple[str, float, str]] = {}
    for line in path.read_text().splitlines()[1:]:
        if not line.strip():
            continue
        key, source, charge, atomtype = line.split("\t")
        rows[key] = (source, float(charge), atomtype)
    return rows


def standard_keys(path: Path) -> list[str]:
    keys: list[str] = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("LINK"):
            continue
        keys.append(line.split()[0])
    return keys


def residue_id(key: str) -> str:
    resid, resname, _atom = key.split("-", 2)
    return f"{resid}-{resname}"


def atom_name(key: str) -> str:
    return key.split("-", 2)[2]


def resname(key: str) -> str:
    return key.split("-", 2)[1]


def is_fixed_chgmod3(key: str) -> bool:
    return resname(key) in STANDARD_AA and atom_name(key) in FIXED_AA_ATOMS


def read_mol2_lines(path: Path) -> tuple[list[str], list[tuple[int, str, float]]]:
    lines = path.read_text().splitlines()
    atom_rows: list[tuple[int, str, float]] = []
    in_atoms = False
    for idx, line in enumerate(lines):
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if in_atoms and line.startswith("@<TRIPOS>"):
            break
        if in_atoms and line.strip():
            parts = line.split()
            atom_rows.append((idx, parts[1], float(parts[-1])))
    return lines, atom_rows


def rewrite_mol2_charges(infile: Path, outfile: Path, key_by_atom: dict[str, str], charges: dict[str, float]) -> None:
    lines, atom_rows = read_mol2_lines(infile)
    for idx, atom, _old_charge in atom_rows:
        key = key_by_atom[atom]
        charge = charges[key]
        prefix = lines[idx][:-12]
        lines[idx] = f"{prefix}{charge:12.6f}"
    outfile.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--standard-fp", required=True)
    parser.add_argument("--small-map", required=True)
    parser.add_argument("--template-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--target-charge", type=float, default=1.0)
    args = parser.parse_args()

    standard_fp = Path(args.standard_fp).resolve()
    small_map = Path(args.small_map).resolve()
    template_dir = Path(args.template_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    cmap = read_charge_map(small_map)
    keys = standard_keys(standard_fp)

    template_charge: dict[str, float] = {}
    key_by_file_atom: dict[str, dict[str, str]] = {}
    for resid_resname, fname in RESIDUE_FILE_ORDER:
        mol2 = template_dir / fname
        _lines, atom_rows = read_mol2_lines(mol2)
        residue_keys = [k for k in keys if residue_id(k) == resid_resname]
        by_atom = {atom_name(k): k for k in residue_keys}
        key_by_file_atom[fname] = by_atom
        for _idx, atom, charge in atom_rows:
            template_charge[by_atom[atom]] = charge

    charges: dict[str, float] = {}
    sources: dict[str, str] = {}
    adjustable: list[str] = []
    for key in keys:
        source, small_charge, _atomtype = cmap[key]
        if is_fixed_chgmod3(key):
            charges[key] = template_charge[key]
            sources[key] = "fixed_ChgModD_template"
        elif source.startswith("small_M2_RESP"):
            charges[key] = small_charge
            sources[key] = "small_M2_RESP"
            adjustable.append(key)
        else:
            charges[key] = template_charge[key]
            sources[key] = "template_fallback"

    total_before = sum(charges.values())
    correction = args.target_charge - total_before
    shift = correction / len(adjustable)
    for key in adjustable:
        charges[key] += shift
        sources[key] += "_plus_free_atom_normalization"

    total_after = sum(charges.values())

    for _resid_resname, fname in RESIDUE_FILE_ORDER:
        rewrite_mol2_charges(template_dir / fname, outdir / fname, key_by_file_atom[fname], charges)

    shutil.copy2(standard_fp, outdir / standard_fp.name)
    shutil.copy2(small_map, outdir / small_map.name)

    map_path = outdir / "smallM2_chgmod3like_standard_charge_map.tsv"
    with map_path.open("w") as handle:
        handle.write("key\tsource\tcharge\n")
        for key in keys:
            handle.write(f"{key}\t{sources[key]}\t{charges[key]:.8f}\n")

    report = outdir / "smallM2_chgmod3like_report.txt"
    fixed = [k for k in keys if is_fixed_chgmod3(k)]
    report.write_text(
        "\n".join(
            [
                "small-M2 RESP + ChgModD-like charge model",
                f"standard atoms: {len(keys)}",
                f"fixed backbone/CB atoms: {len(fixed)}",
                f"adjustable small-M2 RESP atoms: {len(adjustable)}",
                f"total before normalization: {total_before:.8f}",
                f"target total charge: {args.target_charge:.8f}",
                f"normalization correction: {correction:.8f}",
                f"per adjustable atom shift: {shift:.8f}",
                f"total after normalization: {total_after:.8f}",
                "",
                "fixed atoms:",
                *fixed,
            ]
        )
        + "\n"
    )
    print(report.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
