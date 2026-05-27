#!/usr/bin/env python3
"""
Extract a seed metadata/provenance template from a simple .xlsx workbook.

This intentionally uses only the Python standard library so it works before
pandas/openpyxl are installed in the user's md environment.
"""

from __future__ import annotations

import argparse
import csv
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple


NS = {
    "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}


EXTRA_COLUMNS = [
    "paper_title",
    "pubmed_id",
    "expression_host",
    "expression_strain",
    "plasmid_vector",
    "promoter",
    "antibiotic",
    "tag_type",
    "tag_position",
    "tag_cleaved",
    "expression_temperature_c",
    "inducer",
    "soluble_expression_reported",
    "purification_method",
    "evidence_quote",
    "evidence_source",
    "provenance_confidence",
    "notes_for_cloning",
]


HEADER_ALIASES = {
    "Family": "family",
    "酶名称": "enzyme_name",
    "來源生物": "source_organism",
    "来源生物": "source_organism",
    "已知混杂性": "known_promiscuity",
    "Uniprot id": "uniprot_id",
    "UniProt id": "uniprot_id",
    "UniProt ID": "uniprot_id",
    "备注": "seed_notes",
    "关键DOI": "key_doi",
}


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
    values: List[str] = []
    for si in root.findall("m:si", NS):
        values.append("".join(t.text or "" for t in si.findall(".//m:t", NS)))
    return values


def first_sheet_path(zf: zipfile.ZipFile) -> str:
    workbook = ET.fromstring(zf.read("xl/workbook.xml"))
    rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    sheet = workbook.find(".//m:sheet", NS)
    if sheet is None:
        raise ValueError("No worksheet found.")
    rel_id = sheet.attrib[f"{{{NS['r']}}}id"]
    rel_map: Dict[str, str] = {}
    for rel in rels:
        rel_map[rel.attrib["Id"]] = rel.attrib["Target"]
    target = rel_map[rel_id]
    return "xl/" + target.lstrip("/")


def read_xlsx_rows(path: Path) -> List[List[str]]:
    with zipfile.ZipFile(path) as zf:
        shared = load_shared_strings(zf)
        sheet_xml = first_sheet_path(zf)
        root = ET.fromstring(zf.read(sheet_xml))
        rows: List[List[str]] = []
        for row in root.findall(".//m:sheetData/m:row", NS):
            values: List[str] = []
            for cell in row.findall("m:c", NS):
                ref = cell.attrib.get("r", "A1")
                idx = col_index(ref)
                while len(values) <= idx:
                    values.append("")
                value_node = cell.find("m:v", NS)
                value = "" if value_node is None else value_node.text or ""
                if cell.attrib.get("t") == "s" and value:
                    value = shared[int(value)]
                values[idx] = value.strip()
            rows.append(values)
    return rows


def pad(row: List[str], width: int) -> List[str]:
    return row + [""] * max(0, width - len(row))


def normalize_header(header: str, style: str) -> str:
    header = header.strip()
    if style == "original":
        return header
    return HEADER_ALIASES.get(header, re.sub(r"[^A-Za-z0-9]+", "_", header).strip("_").lower() or "column")


def unique_headers(headers: List[str]) -> List[str]:
    counts: Dict[str, int] = {}
    out: List[str] = []
    for header in headers:
        base = header or "column"
        counts[base] = counts.get(base, 0) + 1
        out.append(base if counts[base] == 1 else f"{base}_{counts[base]}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Create seed provenance template CSV from base.xlsx.")
    parser.add_argument("xlsx", help="Input workbook, e.g. base.xlsx.")
    parser.add_argument("--output", default="seed_provenance_template.csv")
    parser.add_argument("--header-style", choices=["english", "original"], default="english",
                        help="Use stable English headers by default to avoid Chinese-header encoding issues.")
    args = parser.parse_args()

    rows = read_xlsx_rows(Path(args.xlsx))
    rows = [row for row in rows if any(cell.strip() for cell in row)]
    if not rows:
        raise ValueError("Workbook contains no non-empty rows.")

    raw_headers = [h.strip() for h in rows[0]]
    while raw_headers and raw_headers[-1] == "":
        raw_headers.pop()
    headers = [normalize_header(h, args.header_style) for h in raw_headers]
    headers = unique_headers(headers)
    out_headers = headers + [c for c in EXTRA_COLUMNS if c not in headers]

    with Path(args.output).open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_headers)
        writer.writeheader()
        for row in rows[1:]:
            row = pad(row, len(headers))
            record = {headers[i]: row[i].strip() for i in range(len(headers))}
            for col in EXTRA_COLUMNS:
                record.setdefault(col, "")
            writer.writerow(record)

    print(f"Wrote {args.output} with {len(rows) - 1} seed rows.")


if __name__ == "__main__":
    main()
