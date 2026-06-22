#!/usr/bin/env python
"""Audit residue numbering across PDB/GRO/RDC outputs."""

from __future__ import annotations

import csv
from pathlib import Path


def pdb_ca(path: Path, start: int, end: int) -> list[tuple[int, str, str]]:
    rows = []
    if not path.exists():
        return rows
    seen = set()
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom = line[12:16].strip()
        if atom != "CA":
            continue
        chain = line[21:22].strip()
        try:
            resid = int(line[22:26])
        except ValueError:
            continue
        if start <= resid <= end and (chain, resid) not in seen:
            seen.add((chain, resid))
            rows.append((resid, line[17:20].strip(), chain))
    return rows


def gro_ca(path: Path, start: int, end: int) -> list[tuple[int, str]]:
    rows = []
    if not path.exists():
        return rows
    seen = set()
    for line in path.read_text(errors="replace").splitlines()[2:-1]:
        if len(line) < 20:
            continue
        try:
            resid = int(line[:5])
        except ValueError:
            continue
        atom = line[10:15].strip()
        if start <= resid <= end and atom == "CA" and resid not in seen:
            seen.add(resid)
            rows.append((resid, line[5:10].strip()))
    return rows


def rdc_rows(path: Path, start: int, end: int) -> list[tuple[int, str, str]]:
    if not path.exists():
        return []
    out = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                resid = int(row["resid"])
            except Exception:
                continue
            if start <= resid <= end:
                out.append((resid, row["resname"], row["residue"]))
    return sorted(out)


def main() -> None:
    start, end = 300, 310
    paths = {
        "curated_pdb": Path("/mnt/e/Pose_create/pdb2r5v_chainA_GLU108_Hpp_heavy_plus_AFe_noHHH.pdb"),
        "md_receptor": Path("/mnt/e/Pose_create/product_R_harmonic20ns/md_handoff/R_product_mode1_hppAFe_harmonic/input/receptor_md.pdb"),
        "amber_complex": Path("/mnt/e/Pose_create/product_R_harmonic20ns/md_handoff/R_product_mode1_hppAFe_harmonic/amber/complex_amber.pdb"),
        "gromacs_gro": Path("/mnt/e/Pose_create/product_R_harmonic20ns/md_handoff/R_product_mode1_hppAFe_harmonic/gromacs/product_R_harmonic_20ns.gro"),
        "rdc_all": Path("/mnt/e/Pose_create/product_R_harmonic20ns/md_analysis/product_R_harmonic_20ns_rdc/residue_ligand_rdc_all.csv"),
    }
    for name in ["curated_pdb", "md_receptor", "amber_complex"]:
        print(f"== {name} ==")
        for resid, resname, chain in pdb_ca(paths[name], start, end):
            print(f"{resid:4d} {resname:>3s} chain={chain}")
    print("== gromacs_gro ==")
    for resid, resname in gro_ca(paths["gromacs_gro"], start, end):
        print(f"{resid:4d} {resname:>3s}")
    print("== rdc_all ==")
    for resid, resname, label in rdc_rows(paths["rdc_all"], start, end):
        print(f"{resid:4d} {resname:>3s} {label}")


if __name__ == "__main__":
    main()
