#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List


COLORS = {
    "Plant": "#A8D5BA",
    "Animal": "#8FB8E6",
    "Fungi": "#59A14F",
    "Bacteria": "#E8A0B0",
    "Archaea": "#F4C2A1",
    "Protist": "#F9E79F",
    "Metagenome": "#D5D5D5",
    "Eukaryota": "#EDC948",
    "Unknown": "#AAAAAA",
}

DEFAULT_KINGDOM_STYLE = {
    "strip_width": "25",
    "margin": "5",
    "border_width": "1",
    "border_color": "#000000",
}


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def short_label(row: Dict[str, str]) -> str:
    pdb = row.get("pdb_id", "")
    uid = row.get("uniprot_id", "") or row.get("seed_id", "")
    org = row.get("organism", "")
    set_name = row.get("set_name", "")
    if pdb and uid:
        return f"{pdb}/{uid}"
    if uid:
        return uid
    if org:
        return re.sub(r"\s+", "_", org)[:24]
    return set_name or row.get("id", "")


def seed_aliases(seed_metadata: Path) -> Dict[str, Dict[str, str]]:
    aliases: Dict[str, Dict[str, str]] = {}
    if not seed_metadata.exists():
        return aliases
    for row in read_csv(seed_metadata):
        keys = {
            row.get("seed_id", ""),
            row.get("uniprot_ids", "").split(";")[0].strip(),
            row.get("pdb_id", ""),
        }
        for key in keys:
            if key:
                aliases[key] = row
    return aliases


def write_kingdom(
    path: Path,
    nodes: List[Dict[str, str]],
    *,
    strip_width: str = DEFAULT_KINGDOM_STYLE["strip_width"],
    margin: str = DEFAULT_KINGDOM_STYLE["margin"],
    border_width: str = DEFAULT_KINGDOM_STYLE["border_width"],
    border_color: str = DEFAULT_KINGDOM_STYLE["border_color"],
) -> None:
    labels = ["Plant", "Animal", "Fungi", "Bacteria", "Archaea", "Protist", "Metagenome", "Eukaryota", "Unknown"]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("DATASET_COLORSTRIP\n")
        handle.write("SEPARATOR TAB\n")
        handle.write("DATASET_LABEL\tKingdom\n")
        handle.write("COLOR\t#000000\n")
        handle.write(f"STRIP_WIDTH\t{strip_width}\n")
        handle.write(f"MARGIN\t{margin}\n")
        handle.write(f"BORDER_WIDTH\t{border_width}\n")
        handle.write(f"BORDER_COLOR\t{border_color}\n")
        handle.write("SHOW_STRIP_LABELS\t0\n")
        handle.write("LEGEND_TITLE\tKingdom\n")
        handle.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(labels)) + "\n")
        handle.write("LEGEND_COLORS\t" + "\t".join(COLORS[k] for k in labels) + "\n")
        handle.write("LEGEND_LABELS\t" + "\t".join(labels) + "\n")
        handle.write("DATA\n")
        for row in nodes:
            node_id = row.get("id", "")
            kingdom = row.get("kingdom") or "Unknown"
            handle.write(f"{node_id}\t{COLORS.get(kingdom, COLORS['Unknown'])}\t{kingdom}\n")


def write_core_highlight(path: Path, nodes: List[Dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("DATASET_BINARY\n")
        handle.write("SEPARATOR TAB\n")
        handle.write("DATASET_LABEL\tCore Sequences\n")
        handle.write("COLOR\t#FF0000\n")
        handle.write("FIELD_SHAPES\t2\n")
        handle.write("FIELD_LABELS\tcore\n")
        handle.write("DATA\n")
        for row in nodes:
            if str(row.get("is_core", "")).upper() == "TRUE":
                handle.write(f"{row.get('id','')}\t1\n")


def write_core_text(path: Path, nodes: List[Dict[str, str]], aliases: Dict[str, Dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("DATASET_TEXT\n")
        handle.write("SEPARATOR TAB\n")
        handle.write("DATASET_LABEL\tCore enzyme short names\n")
        handle.write("COLOR\t#111111\n")
        handle.write("MARGIN\t8\n")
        handle.write("DATA\n")
        for row in nodes:
            if str(row.get("is_core", "")).upper() != "TRUE":
                continue
            node_id = row.get("id", "")
            uid = row.get("uniprot_id", "")
            seed = aliases.get(node_id) or aliases.get(uid) or row
            label = short_label(seed)
            handle.write(f"{node_id}\t{label}\t-1\t#111111\tbold\t1.1\t0\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build iTOL annotations matching the hydrolase tree style.")
    parser.add_argument("--nodes", default="/mnt/e/Tree-Metal-F/ssn/ssn_nodes.csv")
    parser.add_argument("--seed-metadata", default="/mnt/e/Tree-Metal-F/tree_core_metadata.csv")
    parser.add_argument("--outdir", default="/mnt/e/Tree-Metal-F/tree_analysis")
    parser.add_argument("--kingdom-strip-width", default=DEFAULT_KINGDOM_STYLE["strip_width"])
    parser.add_argument("--kingdom-margin", default=DEFAULT_KINGDOM_STYLE["margin"])
    parser.add_argument("--kingdom-border-width", default=DEFAULT_KINGDOM_STYLE["border_width"])
    parser.add_argument("--kingdom-border-color", default=DEFAULT_KINGDOM_STYLE["border_color"])
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    nodes = read_csv(Path(args.nodes))
    aliases = seed_aliases(Path(args.seed_metadata))
    write_kingdom(
        outdir / "itol_kingdom_color_strip.txt",
        nodes,
        strip_width=args.kingdom_strip_width,
        margin=args.kingdom_margin,
        border_width=args.kingdom_border_width,
        border_color=args.kingdom_border_color,
    )
    write_core_highlight(outdir / "itol_core_highlight.txt", nodes)
    ann_dir = outdir / "annotations" / "base_xlsx"
    ann_dir.mkdir(parents=True, exist_ok=True)
    write_core_text(ann_dir / "itol_core_short_name_text.txt", nodes, aliases)
    print(f"Wrote {outdir / 'itol_kingdom_color_strip.txt'}")
    print(f"Wrote {outdir / 'itol_core_highlight.txt'}")
    print(f"Wrote {ann_dir / 'itol_core_short_name_text.txt'}")


if __name__ == "__main__":
    main()
