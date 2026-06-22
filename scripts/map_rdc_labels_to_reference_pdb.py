#!/usr/bin/env python
"""Map MD/GROMACS residue labels in RDC outputs back to reference PDB numbering."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


EXCLUDE_RESNAMES = {"LIG", "FE", "WAT", "SOL", "NA", "CL", "Na+", "Cl-"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rdc-dir", required=True, type=Path)
    parser.add_argument("--reference-pdb", required=True, type=Path, help="PDB with the desired residue numbering")
    parser.add_argument("--gro", required=True, type=Path, help="GROMACS structure used by MDAnalysis/RDC")
    return parser.parse_args()


def reference_residues(path: Path) -> list[dict[str, object]]:
    residues: list[dict[str, object]] = []
    seen: set[tuple[str, int, str]] = set()
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM  "):
            continue
        atom = line[12:16].strip()
        if atom != "CA":
            continue
        resname = line[17:20].strip()
        chain = line[21:22].strip() or "-"
        try:
            resid = int(line[22:26])
        except ValueError:
            continue
        key = (chain, resid, resname)
        if key in seen:
            continue
        seen.add(key)
        residues.append(
            {
                "chain": chain,
                "resid": resid,
                "resname": resname,
                "label": f"{resname}:{chain}:{resid}",
            }
        )
    return residues


def gro_residues(path: Path) -> list[dict[str, object]]:
    residues: list[dict[str, object]] = []
    seen: set[int] = set()
    for line in path.read_text(errors="replace").splitlines()[2:-1]:
        if len(line) < 20:
            continue
        try:
            resid = int(line[:5])
        except ValueError:
            continue
        resname = line[5:10].strip()
        atom = line[10:15].strip()
        if resname in EXCLUDE_RESNAMES or atom != "CA" or resid in seen:
            continue
        seen.add(resid)
        residues.append(
            {
                "md_resid": resid,
                "md_resname": resname,
                "md_label": f"{resname}:SYSTEM:{resid}",
            }
        )
    return residues


def build_mapping(reference: list[dict[str, object]], md: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    if len(reference) != len(md):
        raise SystemExit(f"Residue-count mismatch: reference={len(reference)} md={len(md)}")
    mapping: dict[str, dict[str, object]] = {}
    mismatches: list[str] = []
    for ref, got in zip(reference, md):
        md_label = str(got["md_label"])
        ref_label = str(ref["label"])
        if str(ref["resname"]) != str(got["md_resname"]):
            # Histidine protonation names and terminal normalization can differ; keep an audit line.
            if not (str(ref["resname"]).startswith("HI") and str(got["md_resname"]).startswith("HI")):
                mismatches.append(f"{md_label} -> {ref_label}")
        mapping[md_label] = {
            **got,
            "pdb_chain": ref["chain"],
            "pdb_resid": ref["resid"],
            "pdb_resname": ref["resname"],
            "pdb_label": ref_label,
        }
    return mapping


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def map_rdc_csv(path: Path, mapping: dict[str, dict[str, object]]) -> None:
    rows = read_rows(path)
    if not rows:
        return
    mapped_rows: list[dict[str, object]] = []
    for row in rows:
        old_label = row.get("residue", "")
        info = mapping.get(old_label)
        new_row: dict[str, object] = dict(row)
        if info:
            new_row["md_residue"] = old_label
            new_row["md_resid"] = info["md_resid"]
            new_row["md_resname"] = info["md_resname"]
            new_row["residue"] = info["pdb_label"]
            new_row["resname"] = info["pdb_resname"]
            new_row["segid"] = info["pdb_chain"]
            new_row["resid"] = info["pdb_resid"]
        mapped_rows.append(new_row)
    write_rows(path, mapped_rows)


def map_residue_column_csv(path: Path, mapping: dict[str, dict[str, object]]) -> None:
    rows = read_rows(path)
    if not rows:
        return
    mapped_rows: list[dict[str, object]] = []
    for row in rows:
        old_label = row.get("residue", "")
        info = mapping.get(old_label)
        new_row: dict[str, object] = dict(row)
        if info:
            new_row["md_residue"] = old_label
            new_row["residue"] = info["pdb_label"]
            new_row["pdb_resid"] = info["pdb_resid"]
            new_row["pdb_chain"] = info["pdb_chain"]
        mapped_rows.append(new_row)
    write_rows(path, mapped_rows)


def map_timeseries(path: Path, mapping: dict[str, dict[str, object]]) -> None:
    if not path.exists():
        return
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        rows = list(reader)
    if not rows:
        return
    header = rows[0]
    new_header = [header[0]]
    used: set[str] = set()
    for label in header[1:]:
        mapped = str(mapping.get(label, {}).get("pdb_label", label))
        if mapped in used:
            raise SystemExit(f"Duplicate mapped residue label {mapped!r} from {path}")
        used.add(mapped)
        new_header.append(mapped)
    rows[0] = new_header
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    ref = reference_residues(args.reference_pdb)
    md = gro_residues(args.gro)
    mapping = build_mapping(ref, md)

    args.rdc_dir.mkdir(parents=True, exist_ok=True)
    map_rows = list(mapping.values())
    write_rows(args.rdc_dir / "residue_numbering_map.csv", map_rows)

    for name in ["residue_ligand_rdc_all.csv", "residue_ligand_rdc_within_cutoff.csv"]:
        map_rdc_csv(args.rdc_dir / name, mapping)
    for name in [
        "residue_ligand_distance_timeseries.csv",
        "residue_ligand_distance_timeseries_all.csv",
    ]:
        map_timeseries(args.rdc_dir / name, mapping)
    for name in [
        "residue_ligand_velocity_distribution.csv",
        "residue_ligand_velocity_distribution_summary.csv",
    ]:
        map_residue_column_csv(args.rdc_dir / name, mapping)

    metadata_path = args.rdc_dir / "rdc_metadata.json"
    metadata: dict[str, object] = {}
    if metadata_path.exists():
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            metadata = {}
    metadata["residue_numbering"] = {
        "display": "reference_pdb",
        "reference_pdb": str(args.reference_pdb),
        "gro": str(args.gro),
        "map_csv": str(args.rdc_dir / "residue_numbering_map.csv"),
        "note": "RDC labels were mapped from GROMACS/Amber sequential residue ids back to the reference PDB chain/residue numbering.",
    }
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print(f"mapped_residues={len(mapping)}")
    print(args.rdc_dir / "residue_numbering_map.csv")


if __name__ == "__main__":
    main()
