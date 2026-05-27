#!/usr/bin/env python3
"""
Generate iTOL annotations from the curated base.xlsx table.

Rules:
- rows whose notes contain "missed candidate" in Chinese are audited for merging
  into the existing Core Sequences highlight.
- rows whose notes contain "extra candidate" in Chinese are marked as yellow stars.
- core/curated rows get short enzyme-name text labels, using text inside
  half-width or full-width parentheses when present.
"""

from __future__ import annotations

import argparse
import csv
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, Iterable, List, Optional


NS = {
    "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}

HEADER_ALIASES = {
    "Family": "family",
    "\u9176\u540d\u79f0": "enzyme_name",
    "\u6765\u6e90\u751f\u7269": "source_organism",
    "\u5df2\u77e5\u6df7\u6742\u6027": "known_promiscuity",
    "Uniprot id": "uniprot_id",
    "UniProt id": "uniprot_id",
    "UniProt ID": "uniprot_id",
    "\u5907\u6ce8": "notes",
    "\u5173\u952eDOI": "key_doi",
    "Plasmid \uff08pET\uff09": "plasmid_pet",
    "Plasmid (pET)": "plasmid_pet",
}

MISSED_CANDIDATE_MARKER = "\u6f0f\u8865\u5019\u9009"
EXTRA_CANDIDATE_MARKER = "\u989d\u5916\u5019\u9009"


def col_index(cell_ref: str) -> int:
    letters = re.match(r"[A-Z]+", cell_ref).group(0)
    out = 0
    for ch in letters:
        out = out * 26 + (ord(ch) - ord("A") + 1)
    return out - 1


def load_shared_strings(zf: zipfile.ZipFile) -> List[str]:
    try:
        root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    except KeyError:
        return []
    return ["".join(t.text or "" for t in si.findall(".//m:t", NS)) for si in root.findall("m:si", NS)]


def cell_text(cell: ET.Element, shared: List[str]) -> str:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(t.text or "" for t in cell.findall(".//m:t", NS)).strip()

    v = cell.find("m:v", NS)
    val = "" if v is None else v.text or ""
    if cell_type == "s" and val:
        return shared[int(val)].strip()
    return val.strip()


def read_xlsx_rows(path: Path) -> List[Dict[str, str]]:
    with zipfile.ZipFile(path) as zf:
        shared = load_shared_strings(zf)
        root = ET.fromstring(zf.read("xl/worksheets/sheet1.xml"))
        raw_rows: List[List[str]] = []
        for row in root.findall(".//m:sheetData/m:row", NS):
            vals: List[str] = []
            for cell in row.findall("m:c", NS):
                idx = col_index(cell.attrib.get("r", "A1"))
                while len(vals) <= idx:
                    vals.append("")
                vals[idx] = cell_text(cell, shared)
            if any(vals):
                raw_rows.append(vals)
    if not raw_rows:
        return []
    headers = [HEADER_ALIASES.get(h.strip(), re.sub(r"[^A-Za-z0-9]+", "_", h.strip()).strip("_").lower()) for h in raw_rows[0]]
    while headers and not headers[-1]:
        headers.pop()
    rows = []
    for raw in raw_rows[1:]:
        row = {}
        for idx, header in enumerate(headers):
            row[header] = raw[idx].strip() if idx < len(raw) else ""
        rows.append(row)
    return rows


def read_nodes(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def clean_accession(raw: str) -> str:
    text = str(raw or "").strip()
    if "|" in text:
        for part in text.split("|"):
            if re.match(r"^[A-Z0-9]{6,10}(?:-\d+)?$", part.strip()):
                return part.strip()
    text = re.sub(r"^UniRef\d+_", "", text)
    text = re.sub(r"\.\d+$", "", text)
    match = re.search(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})", text)
    return match.group(1) if match else text.split()[0] if text else ""


def node_aliases(nodes: List[Dict[str, str]]) -> Dict[str, str]:
    aliases: Dict[str, str] = {}
    for node in nodes:
        node_id = node.get("id", "").strip()
        uid = node.get("uniprot_id", "").strip()
        for value in [node_id, uid, clean_accession(node_id), clean_accession(uid)]:
            if value:
                aliases[value] = node_id
    return aliases


def short_label(enzyme_name: str) -> str:
    name = (enzyme_name or "").strip()
    if not name:
        return ""
    matches = re.findall(r"[\(\uff08]([^()\uff08\uff09]+)[\)\uff09]", name)
    if matches:
        return matches[-1].strip()
    return re.sub(r"\s+", " ", name)


def row_kind(row: Dict[str, str]) -> str:
    notes = row.get("notes", "")
    if MISSED_CANDIDATE_MARKER in notes:
        return "missed_candidate"
    if EXTRA_CANDIDATE_MARKER in notes:
        return "extra_candidate"
    return "core"


def write_extra_star_dataset(rows: List[Dict[str, str]], aliases: Dict[str, str], outpath: Path) -> List[Dict[str, str]]:
    written: List[Dict[str, str]] = []
    with outpath.open("w", encoding="utf-8") as handle:
        handle.write("DATASET_BINARY\n")
        handle.write("SEPARATOR TAB\n")
        handle.write("DATASET_LABEL\tExtra unvalidated candidates\n")
        handle.write("COLOR\t#000000\n")
        handle.write("FIELD_SHAPES\t3\n")
        handle.write("FIELD_COLORS\t#FFD21F\n")
        handle.write("FIELD_LABELS\textra_unvalidated_candidate\n")
        handle.write("MARGIN\t18\n")
        handle.write("LEGEND_TITLE\tBase.xlsx extra candidates\n")
        handle.write("LEGEND_SHAPES\t3\n")
        handle.write("LEGEND_COLORS\t#FFD21F\n")
        handle.write("LEGEND_LABELS\tExtra unvalidated candidate\n")
        handle.write("DATA\n")
        for row in rows:
            kind = row_kind(row)
            if kind != "extra_candidate":
                continue
            uid = clean_accession(row.get("uniprot_id", ""))
            node_id = aliases.get(uid) or aliases.get(row.get("uniprot_id", "").strip())
            if not node_id:
                written.append({**row, "annotation_kind": kind, "node_id": "", "matched": "false"})
                continue
            handle.write(f"{node_id}\t1\n")
            written.append({**row, "annotation_kind": kind, "node_id": node_id, "matched": "true"})
    return written


def collect_missed_core_rows(rows: List[Dict[str, str]], aliases: Dict[str, str]) -> List[Dict[str, str]]:
    written: List[Dict[str, str]] = []
    for row in rows:
        kind = row_kind(row)
        if kind != "missed_candidate":
            continue
        uid = clean_accession(row.get("uniprot_id", ""))
        node_id = aliases.get(uid) or aliases.get(row.get("uniprot_id", "").strip())
        written.append({**row, "annotation_kind": kind, "node_id": node_id or "", "matched": "true" if node_id else "false"})
    return written


def write_text_dataset(rows: List[Dict[str, str]], aliases: Dict[str, str], outpath: Path) -> List[Dict[str, str]]:
    written: List[Dict[str, str]] = []
    with outpath.open("w", encoding="utf-8") as handle:
        handle.write("DATASET_TEXT\n")
        handle.write("SEPARATOR TAB\n")
        handle.write("DATASET_LABEL\tCore enzyme short names\n")
        handle.write("COLOR\t#111111\n")
        handle.write("MARGIN\t8\n")
        handle.write("DATA\n")
        for row in rows:
            if row_kind(row) == "extra_candidate":
                continue
            label = short_label(row.get("enzyme_name", ""))
            if not label:
                continue
            uid = clean_accession(row.get("uniprot_id", ""))
            node_id = aliases.get(uid) or aliases.get(row.get("uniprot_id", "").strip())
            if not node_id:
                written.append({**row, "annotation_kind": "core_label", "short_label": label, "node_id": "", "matched": "false"})
                continue
            handle.write(f"{node_id}\t{label}\t-1\t#111111\tbold\t1.1\t0\n")
            written.append({**row, "annotation_kind": "core_label", "short_label": label, "node_id": node_id, "matched": "true"})
    return written


def write_summary(path: Path, rows: List[Dict[str, str]]) -> None:
    headers = [
        "annotation_kind", "node_id", "matched", "family", "enzyme_name", "short_label",
        "source_organism", "known_promiscuity", "uniprot_id", "notes", "key_doi", "plasmid_pet",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate iTOL annotations from base.xlsx curation flags.")
    parser.add_argument("base_xlsx")
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--outdir", default="base_xlsx_itol_annotations")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = read_xlsx_rows(Path(args.base_xlsx))
    aliases = node_aliases(read_nodes(Path(args.nodes)))
    missed_rows = collect_missed_core_rows(rows, aliases)
    symbol_rows = write_extra_star_dataset(rows, aliases, outdir / "itol_extra_unvalidated_stars.txt")
    text_rows = write_text_dataset(rows, aliases, outdir / "itol_core_short_name_text.txt")
    write_summary(outdir / "base_itol_annotation_summary.csv", missed_rows + symbol_rows + text_rows)
    print(f"Wrote {outdir / 'itol_extra_unvalidated_stars.txt'}")
    print(f"Wrote {outdir / 'itol_core_short_name_text.txt'}")
    print(f"Wrote {outdir / 'base_itol_annotation_summary.csv'}")


if __name__ == "__main__":
    main()
