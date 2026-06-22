#!/usr/bin/env python3
"""Build a PyMOL PSE from merged Metal-F element-rescan XLSX."""

from __future__ import annotations

import argparse
import csv
import os
import re
import urllib.request
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd
from pymol import cmd


PROJECT = Path("/home/qin/Metal-F_project")
INTERACTIONS = PROJECT / "Metal-F_interactions"
CIF_DIR = PROJECT / ".cache" / "metal_f_element_ligand_rescan" / "cif"


def xlsx_rows(path: Path, sheet_name: str = "geometry_hits") -> list[dict[str, str]]:
    try:
        frame = pd.read_excel(path, sheet_name=sheet_name)
        frame = frame.fillna("")
        return [{key: str(value) for key, value in row.items()} for row in frame.to_dict("records")]
    except Exception:
        pass

    ns = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

    def col_index(ref: str) -> int:
        letters = "".join(ch for ch in ref if ch.isalpha())
        n = 0
        for ch in letters:
            n = n * 26 + ord(ch.upper()) - 64
        return n - 1

    with zipfile.ZipFile(path) as zf:
        shared = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall("a:si", ns):
                shared.append("".join(t.text or "" for t in si.findall(".//a:t", ns)))
        workbook = ET.fromstring(zf.read("xl/workbook.xml"))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rel_map = {rel.get("Id"): rel.get("Target") for rel in rels}
        sheet_path = None
        for sheet in workbook.findall(".//a:sheet", ns):
            if sheet.get("name") == sheet_name:
                rel_id = sheet.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
                target = rel_map[rel_id]
                target = target.lstrip("/")
                sheet_path = target if target.startswith("xl/") else "xl/" + target
                break
        if not sheet_path:
            raise ValueError(f"Sheet not found: {sheet_name}")
        sheet = ET.fromstring(zf.read(sheet_path))
        table = []
        for row in sheet.findall(".//a:row", ns):
            values = []
            for cell in row.findall("a:c", ns):
                idx = col_index(cell.get("r", "A1"))
                while len(values) <= idx:
                    values.append("")
                value_node = cell.find("a:v", ns)
                value = "" if value_node is None else value_node.text or ""
                if cell.get("t") == "s" and value:
                    value = shared[int(value)]
                values[idx] = value
            table.append(values)
    if not table:
        return []
    header = table[0]
    out = []
    for row in table[1:]:
        if not any(row):
            continue
        row = row + [""] * (len(header) - len(row))
        out.append(dict(zip(header, row)))
    return out


def ids_from_csv(path: Path, column: str, value: str, exclude_column: str = "", exclude_value: str = "") -> list[str]:
    if not path.exists():
        return []
    frame = pd.read_csv(path).fillna("")
    if "pdb_id" not in frame.columns:
        return []
    if column and column in frame.columns:
        frame = frame[frame[column].astype(str) == value]
    if exclude_column and exclude_column in frame.columns:
        excluded = {item.strip() for item in str(exclude_value).split(",") if item.strip()}
        if excluded:
            frame = frame[~frame[exclude_column].astype(str).isin(excluded)]
    return sorted(dict.fromkeys(frame["pdb_id"].astype(str).str.upper()))


def download_cif(pdb_id: str) -> Path:
    CIF_DIR.mkdir(parents=True, exist_ok=True)
    path = CIF_DIR / f"{pdb_id.upper()}.cif"
    if path.exists() and path.stat().st_size > 0:
        return path
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    request = urllib.request.Request(url, headers={"User-Agent": "metal-f-element-pse/1.0"})
    with urllib.request.urlopen(request, timeout=180) as response:
        path.write_bytes(response.read())
    return path


def atom_parts(label: str):
    parts = str(label or "").split(":")
    if len(parts) < 5:
        return None
    return parts[:5]


def atom_selection(obj: str, label: str) -> str:
    parts = atom_parts(label)
    if not parts:
        return "none"
    comp, chain, resi, atom, elem = parts
    clauses = [f"model {obj}", f"resn {comp}", f"resi {resi}", f"name {atom}", f"elem {elem}"]
    if chain:
        clauses.append(f"chain {chain}")
    return "(" + " and ".join(clauses) + ")"


def has_atoms(selection: str) -> bool:
    try:
        return cmd.count_atoms(selection) > 0
    except Exception:
        return False


def safe_name(text: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_]+", "_", str(text)).strip("_")
    if not value:
        value = "site"
    if value[0].isdigit():
        value = "pdb_" + value
    return value[:80]


def style_distance(name: str) -> None:
    cmd.set("dash_color", "yellow", name)
    cmd.set("dash_radius", 0.03, name)
    cmd.set("dash_gap", 0.35, name)
    cmd.set("label_color", "yellow", name)
    cmd.show("dashes", name)
    cmd.show("labels", name)


def draw_pairs(row: dict[str, str], obj: str, prefix: str, max_pairs: int = 12) -> int:
    count = 0
    for item in str(row.get("passing_pairs_le_cutoff_A", "")).split(";"):
        if count >= max_pairs:
            break
        if "--" not in item or ":" not in item:
            continue
        _, atoms = item.split(":", 1)
        f_label, metal_label = atoms.split("--", 1)
        f_sel = atom_selection(obj, f_label)
        metal_sel = atom_selection(obj, metal_label)
        if not (has_atoms(f_sel) and has_atoms(metal_sel)):
            continue
        name = f"{prefix}_{count}"
        cmd.distance(name, f_sel, metal_sel)
        style_distance(name)
        count += 1
    if count == 0:
        f_sel = atom_selection(obj, row.get("nearest_F_atom", ""))
        metal_sel = atom_selection(obj, row.get("nearest_metal_atom", ""))
        if has_atoms(f_sel) and has_atoms(metal_sel):
            name = f"{prefix}_nearest"
            cmd.distance(name, f_sel, metal_sel)
            style_distance(name)
            count = 1
    return count


def style_site(obj: str, f_sel: str, metal_sel: str) -> None:
    obj_sel = f"(model {obj})"
    core = f"(({f_sel}) or ({metal_sel}))"
    cmd.disable("all")
    cmd.enable(obj)
    cmd.hide("everything", obj_sel)
    cmd.show("cartoon", f"{obj_sel} and polymer.protein")
    cmd.color("gray80", f"{obj_sel} and polymer.protein and elem C")

    sidechains = f"byres ((polymer.protein within 4.0 of {core})) and not name N+C+O+OXT"
    if has_atoms(sidechains):
        cmd.show("sticks", sidechains)
        cmd.color("gray65", f"({sidechains}) and elem C")

    hetero = f"byres (((organic or inorganic) and not solvent) within 4.0 of {core})"
    if has_atoms(hetero):
        cmd.show("sticks", hetero)
        cmd.color("gray70", f"({hetero}) and elem C")

    cmd.show("spheres", f_sel)
    cmd.show("spheres", metal_sel)
    cmd.set("sphere_scale", 0.32, f_sel)
    cmd.set("sphere_scale", 0.32, metal_sel)


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "pdb_id",
        "scene",
        "nearest_F_metal_A",
        "nearest_F_atom",
        "nearest_metal_atom",
        "fluorine_comp_ids",
        "metal_comp_ids",
        "metal_elements",
        "protein_entity_count",
        "complex_warning",
        "fasta_seed_included",
        "distance_object_count",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build(
    xlsx: Path,
    output_pse: Path,
    summary_csv: Path,
    max_scenes: int = 0,
    include_pdb_ids: set[str] | None = None,
) -> int:
    rows = xlsx_rows(xlsx, "geometry_hits")
    if include_pdb_ids is not None:
        rows = [row for row in rows if str(row.get("pdb_id", "")).upper() in include_pdb_ids]
    if max_scenes:
        rows = rows[:max_scenes]
    cmd.delete("all")
    for scene in list(cmd.get_scene_list()):
        try:
            cmd.scene(scene, "clear")
        except Exception:
            pass

    made = []
    summary = []
    for index, row in enumerate(rows, 1):
        pdb_id = str(row.get("pdb_id", "")).upper()
        if not pdb_id:
            continue
        print(f"[Metal-F element PSE] {index}/{len(rows)} {pdb_id}", flush=True)
        cif = download_cif(pdb_id)
        obj = f"pdb_{pdb_id.lower()}"
        cmd.load(str(cif), obj)
        f_sel = atom_selection(obj, row.get("nearest_F_atom", ""))
        metal_sel = atom_selection(obj, row.get("nearest_metal_atom", ""))
        if not (has_atoms(f_sel) and has_atoms(metal_sel)):
            print(f"[Metal-F element PSE] atom lookup failed for {pdb_id}")
            continue
        style_site(obj, f_sel, metal_sel)
        scene_safe = safe_name(
            f"{pdb_id}_{row.get('fluorine_comp_ids', '')}_{row.get('metal_comp_ids', '')}_{row.get('nearest_F_metal_A', '')}"
        )
        distance_count = draw_pairs(row, obj, f"mfe_{scene_safe}")
        focus = f"(({f_sel}) or ({metal_sel})) expand 7"
        if has_atoms(focus):
            cmd.orient(focus)
            cmd.zoom(focus, 3.0)
        scene = f"MFE_{scene_safe}"
        cmd.scene(scene, "store")
        made.append(scene)
        out_row = dict(row)
        out_row["scene"] = scene
        out_row["distance_object_count"] = str(distance_count)
        summary.append(out_row)

    if made:
        cmd.scene(made[0], "recall")
        cmd.set("scene_buttons", 1)
        output_pse.parent.mkdir(parents=True, exist_ok=True)
        cmd.save(str(output_pse))
        write_summary(summary_csv, summary)
    print(f"[Metal-F element PSE] scenes={len(made)} output={output_pse}")
    return len(made)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xlsx", default=str(INTERACTIONS / "Metal-F_element_ligand_rescan_full_merged.xlsx"))
    parser.add_argument("--output-pse", default=str(INTERACTIONS / "Metal-F_element_ligand_rescan_full_merged.pse"))
    parser.add_argument("--summary-csv", default=str(INTERACTIONS / "Metal-F_element_ligand_rescan_full_merged_pse_summary.csv"))
    parser.add_argument("--max-scenes", type=int, default=0)
    parser.add_argument("--pdb-ids", default="", help="Comma/space separated PDB IDs to include.")
    parser.add_argument("--pdb-ids-from-csv", default="", help="Read include PDB IDs from CSV.")
    parser.add_argument("--csv-filter-column", default="")
    parser.add_argument("--csv-filter-value", default="")
    parser.add_argument("--csv-exclude-column", default="")
    parser.add_argument("--csv-exclude-value", default="")
    args = parser.parse_args()
    ids = set()
    if args.pdb_ids:
        ids |= {item.upper() for item in re.findall(r"\b[0-9][A-Za-z0-9]{3}\b", args.pdb_ids)}
    if args.pdb_ids_from_csv:
        ids |= set(
            ids_from_csv(
                Path(args.pdb_ids_from_csv),
                args.csv_filter_column,
                args.csv_filter_value,
                args.csv_exclude_column,
                args.csv_exclude_value,
            )
        )
    build(
        Path(args.xlsx),
        Path(args.output_pse),
        Path(args.summary_csv),
        args.max_scenes,
        ids if ids else None,
    )


if __name__ == "__main__":
    main()
