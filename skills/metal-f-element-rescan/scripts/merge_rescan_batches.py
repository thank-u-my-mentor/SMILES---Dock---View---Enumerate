#!/usr/bin/env python3
"""Merge Metal-F element-rescan batch XLSX/FASTA outputs."""

from __future__ import annotations

import argparse
import glob
import textwrap
from pathlib import Path

import pandas as pd


def concat_sheets(files: list[Path], sheet: str) -> pd.DataFrame:
    frames = []
    for path in files:
        xls = pd.ExcelFile(path)
        if sheet not in xls.sheet_names:
            continue
        frame = pd.read_excel(path, sheet_name=sheet)
        frame["source_batch_file"] = path.name
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def dedupe_frame(frame: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    if frame.empty:
        return frame
    usable = [key for key in keys if key in frame.columns]
    if not usable:
        return frame.drop_duplicates()
    return frame.drop_duplicates(subset=usable, keep="first")


def write_fasta_from_sheet(path: Path, frame: pd.DataFrame) -> int:
    if frame.empty or "header" not in frame.columns or "sequence" not in frame.columns:
        path.write_text("", encoding="utf-8")
        return 0
    count = 0
    seen = set()
    with open(path, "w", encoding="utf-8") as handle:
        for _, row in frame.iterrows():
            header = str(row.get("header", "")).strip()
            seq = str(row.get("sequence", "")).replace(" ", "").replace("\n", "").strip()
            if not header or not seq or seq.lower() == "nan":
                continue
            key = (header, seq)
            if key in seen:
                continue
            seen.add(key)
            handle.write(f">{header}\n")
            handle.write("\n".join(textwrap.wrap(seq, 80)) + "\n")
            count += 1
    return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-glob", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    files = [Path(item) for item in sorted(glob.glob(args.input_glob)) if Path(item).suffix.lower() == ".xlsx"]
    if not files:
        raise SystemExit(f"No XLSX files matched: {args.input_glob}")

    prefix = Path(args.prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    geometry = dedupe_frame(concat_sheets(files, "geometry_hits"), ["pdb_id", "nearest_F_atom", "nearest_metal_atom"])
    processed = dedupe_frame(concat_sheets(files, "processed_candidates"), ["pdb_id"])
    protein = dedupe_frame(concat_sheets(files, "protein_fasta"), ["header", "sequence"])
    dna = dedupe_frame(concat_sheets(files, "dna_fasta"), ["header", "sequence"])
    ccd_f = dedupe_frame(concat_sheets(files, "ccd_f_ligands"), ["comp_id"])
    ccd_metal = dedupe_frame(concat_sheets(files, "ccd_metal_ligands"), ["comp_id"])

    if "nearest_F_metal_A" in geometry.columns:
        geometry["_sort_dist"] = pd.to_numeric(geometry["nearest_F_metal_A"], errors="coerce")
        geometry = geometry.sort_values(["_sort_dist", "pdb_id"], na_position="last").drop(columns=["_sort_dist"])
    if "status" in processed.columns:
        processed = processed.sort_values(["status", "pdb_id"], na_position="last")

    protein_fasta = prefix.with_name(prefix.name + "_protein.fasta")
    dna_fasta = prefix.with_name(prefix.name + "_cds_dna.fasta")
    protein_count = write_fasta_from_sheet(protein_fasta, protein)
    dna_count = write_fasta_from_sheet(dna_fasta, dna)

    summary = pd.DataFrame(
        [
            {"metric": "batch_files", "value": ";".join(path.name for path in files)},
            {"metric": "batch_file_count", "value": len(files)},
            {"metric": "geometry_hit_count", "value": len(geometry)},
            {"metric": "processed_candidate_count", "value": len(processed)},
            {"metric": "protein_fasta_records", "value": protein_count},
            {"metric": "dna_fasta_records", "value": dna_count},
            {"metric": "protein_fasta", "value": str(protein_fasta)},
            {"metric": "dna_fasta", "value": str(dna_fasta)},
        ]
    )

    xlsx = prefix.with_suffix(".xlsx")
    with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
        geometry.to_excel(writer, sheet_name="geometry_hits", index=False)
        processed.to_excel(writer, sheet_name="processed_candidates", index=False)
        protein.to_excel(writer, sheet_name="protein_fasta", index=False)
        dna.to_excel(writer, sheet_name="dna_fasta", index=False)
        ccd_f.to_excel(writer, sheet_name="ccd_f_ligands", index=False)
        ccd_metal.to_excel(writer, sheet_name="ccd_metal_ligands", index=False)
        summary.to_excel(writer, sheet_name="run_summary", index=False)

    print(f"merged_xlsx: {xlsx}")
    print(f"protein_fasta: {protein_fasta} ({protein_count})")
    print(f"dna_fasta: {dna_fasta} ({dna_count})")
    print(f"geometry_hits: {len(geometry)}")
    print(f"processed_candidates: {len(processed)}")


if __name__ == "__main__":
    main()
