#!/usr/bin/env python3
"""
Convert NetSolP/SoluProt-style CSV predictions into iTOL annotation files.

The input CSV can use flexible column names. IDs are matched against SSN/tree
IDs, UniProt IDs, or cleaned accessions when ssn_nodes.csv is provided.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional


ID_COLUMNS = ["id", "fa_id", "sequence_id", "seq_id", "name", "entry", "accession", "uniprot_id", "Uniprot id"]
SCORE_COLUMNS = [
    "soluprot_probability",
    "probability",
    "solubility_probability",
    "solubility",
    "soluble",
    "score",
    "NetSolP_solubility",
    "netsolp_solubility",
    "SoluProt",
    "soluprot",
]
CALL_COLUMNS = ["prediction", "class", "call", "soluble", "binary", "NetSolP_call", "SoluProt_call"]


def clean_accession(raw: str) -> str:
    text = str(raw).strip()
    if "|" in text:
        for part in text.split("|"):
            if re.match(r"^[A-Z0-9]{6,10}(?:-\d+)?$", part):
                return part
    text = re.sub(r"^UniRef\d+_", "", text)
    text = re.sub(r"\.\d+$", "", text)
    match = re.search(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})", text)
    return match.group(1) if match else text.split()[0]


def simplify_tree_id(seq_id: str) -> str:
    text = str(seq_id or "").strip()
    if "|" in text:
        left, right = text.split("|", 1)
        if right == left or right.startswith(left + "_"):
            return right
    return text


def first_existing(row: Dict[str, str], names: Iterable[str]) -> str:
    lower = {k.lower(): k for k in row}
    for name in names:
        if name in row and str(row[name]).strip():
            return str(row[name]).strip()
        key = lower.get(name.lower())
        if key and str(row[key]).strip():
            return str(row[key]).strip()
    return ""


def parse_float(value: str) -> Optional[float]:
    if not value:
        return None
    match = re.search(r"-?\d+(?:\.\d+)?", str(value))
    if not match:
        return None
    num = float(match.group(0))
    if num > 1.0 and num <= 100.0:
        num = num / 100.0
    return max(0.0, min(1.0, num))


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_nodes(path: Optional[Path]) -> List[Dict[str, str]]:
    if not path:
        return []
    return read_csv(path)


def read_id_mapping(path: Optional[Path]) -> Dict[str, str]:
    if not path:
        return {}
    mapping: Dict[str, str] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        if "\t" in sample:
            reader = csv.DictReader(handle, delimiter="\t")
        else:
            reader = csv.DictReader(handle)
        for row in reader:
            old_id = (row.get("old_id") or row.get("original_id") or row.get("tree_id") or "").strip()
            clean_id = (row.get("clean_id") or row.get("new_id") or row.get("fa_id") or "").strip()
            if old_id and clean_id:
                mapping[clean_id] = old_id
    return mapping


def build_node_aliases(nodes: List[Dict[str, str]]) -> Dict[str, str]:
    aliases: Dict[str, str] = {}
    for row in nodes:
        node_id = row.get("id", "").strip()
        if not node_id:
            continue
        for value in [
            node_id,
            simplify_tree_id(node_id),
            row.get("uniprot_id", ""),
            clean_accession(node_id),
            clean_accession(row.get("uniprot_id", "")),
        ]:
            if value:
                aliases[value] = node_id
    return aliases


def normalize_call(call: str, score: Optional[float], threshold: float) -> str:
    text = str(call or "").strip().lower()
    if text in {"soluble", "yes", "true", "1", "s"}:
        return "Soluble"
    if text in {"insoluble", "no", "false", "0", "i"}:
        return "Insoluble"
    if score is None:
        return "Unknown"
    return "Soluble" if score >= threshold else "Insoluble"


def call_color(call: str, score: Optional[float]) -> str:
    if call == "Unknown":
        return "#BDBDBD"
    if score is not None and score >= 0.75:
        return "#1B9E3C"
    if score is not None and score >= 0.50:
        return "#9CCC65"
    if score is not None and score < 0.25:
        return "#D7191C"
    return "#EF9A9A" if call == "Insoluble" else "#9CCC65"


def _lerp(a: int, b: int, t: float) -> int:
    return int(round(a + (b - a) * t))


def score_color(score: float) -> str:
    """Red-yellow-green color for a 0..1 SoluProt/solubility score."""
    score = max(0.0, min(1.0, score))
    red = (215, 25, 28)
    yellow = (247, 247, 166)
    green = (27, 158, 60)
    if score < 0.5:
        t = score / 0.5
        rgb = tuple(_lerp(red[i], yellow[i], t) for i in range(3))
    else:
        t = (score - 0.5) / 0.5
        rgb = tuple(_lerp(yellow[i], green[i], t) for i in range(3))
    return "#{:02X}{:02X}{:02X}".format(*rgb)


def score_color_rescaled(score: float, low_threshold: float) -> str:
    """Red-yellow-green color after treating low_threshold as the faint red floor."""
    if score <= low_threshold:
        return "#D7191C"
    scaled = (score - low_threshold) / max(1e-9, 1.0 - low_threshold)
    return score_color(scaled)


def symbol_size(score: float, low_threshold: float) -> float:
    """Draw low-solubility values small/faint and high-solubility values visibly."""
    if score <= low_threshold:
        return 4.0
    scaled = (score - low_threshold) / max(1e-9, 1.0 - low_threshold)
    return 5.0 + 13.0 * scaled


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate iTOL solubility annotations from prediction CSV.")
    parser.add_argument("solubility_csv", help="CSV from NetSolP/SoluProt or a manually normalized table.")
    parser.add_argument("--nodes", default=None, help="Optional ssn_nodes.csv for ID/UniProt mapping and merged output.")
    parser.add_argument("--id-mapping", default=None, help="Optional TSV/CSV with old_id and clean_id columns from prepare_soluprot_tree_inputs.py.")
    parser.add_argument("--outdir", default="solubility_itol")
    parser.add_argument("--tool-label", default="Solubility")
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--low-threshold", type=float, default=0.2, help="Scores at/below this are drawn faint red in gradient datasets.")
    parser.add_argument("--no-color-strip", action="store_true", help="Skip categorical color strip output.")
    parser.add_argument("--symbols-only", action="store_true", help="Only write gradient dot annotation plus normalized/merged CSV files.")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    predictions = read_csv(Path(args.solubility_csv))
    nodes = read_nodes(Path(args.nodes)) if args.nodes else []
    aliases = build_node_aliases(nodes)
    clean_to_old = read_id_mapping(Path(args.id_mapping)) if args.id_mapping else {}

    by_node: Dict[str, Dict[str, str]] = {}
    unmatched: List[Dict[str, str]] = []
    for row in predictions:
        raw_id = first_existing(row, ID_COLUMNS)
        mapped_raw_id = clean_to_old.get(raw_id, raw_id)
        key_candidates = [mapped_raw_id, raw_id, simplify_tree_id(raw_id), clean_accession(raw_id)]
        node_id = next((aliases[k] for k in key_candidates if k in aliases), raw_id)
        score = parse_float(first_existing(row, SCORE_COLUMNS))
        call = normalize_call(first_existing(row, CALL_COLUMNS), score, args.threshold)
        record = {
            "id": node_id,
            "source_id": raw_id,
            "solubility_score": "" if score is None else f"{score:.4f}",
            "solubility_call": call,
        }
        if args.nodes and node_id == raw_id and raw_id not in aliases:
            unmatched.append(record)
        by_node[node_id] = record

    strip_path = outdir / "itol_solubility_color_strip.txt"
    if not args.no_color_strip:
        with strip_path.open("w", encoding="utf-8") as handle:
            handle.write("DATASET_COLORSTRIP\n")
            handle.write("SEPARATOR TAB\n")
            handle.write(f"DATASET_LABEL\t{args.tool_label}\n")
            handle.write("COLOR\t#000000\n")
            handle.write("STRIP_WIDTH\t25\n")
            handle.write("MARGIN\t5\n")
            handle.write("BORDER_WIDTH\t1\n")
            handle.write("BORDER_COLOR\t#000000\n")
            handle.write("SHOW_STRIP_LABELS\t0\n")
            handle.write("LEGEND_TITLE\tSolubility prediction\n")
            handle.write("LEGEND_SHAPES\t1\t1\t1\t1\t1\n")
            handle.write("LEGEND_COLORS\t#1B9E3C\t#9CCC65\t#F4A3A3\t#D7191C\t#BDBDBD\n")
            handle.write("LEGEND_LABELS\tHigh soluble\tSoluble\tInsoluble\tLow soluble\tUnknown\n")
            handle.write("DATA\n")
            source_rows = nodes if nodes else [{"id": k} for k in by_node]
            for node in source_rows:
                node_id = node.get("id", "").strip()
                pred = by_node.get(node_id)
                if not pred:
                    continue
                score = parse_float(pred["solubility_score"])
                label = pred["solubility_call"] if not pred["solubility_score"] else f"{pred['solubility_call']} ({pred['solubility_score']})"
                handle.write(f"{node_id}\t{call_color(pred['solubility_call'], score)}\t{label}\n")

    bar_path = outdir / "itol_solubility_simplebar.txt"
    grad_path = outdir / "itol_solubility_gradient_strip.txt"
    if not args.symbols_only:
        with bar_path.open("w", encoding="utf-8") as handle:
            handle.write("DATASET_SIMPLEBAR\n")
            handle.write("SEPARATOR TAB\n")
            handle.write(f"DATASET_LABEL\t{args.tool_label} score\n")
            handle.write("COLOR\t#1B9E3C\n")
            handle.write("WIDTH\t80\n")
            handle.write("MARGIN\t5\n")
            handle.write("HEIGHT_FACTOR\t1\n")
            handle.write("BORDER_WIDTH\t1\n")
            handle.write("BORDER_COLOR\t#444444\n")
            handle.write("DATA\n")
            source_rows = nodes if nodes else [{"id": k} for k in by_node]
            for node in source_rows:
                node_id = node.get("id", "").strip()
                pred = by_node.get(node_id)
                if not pred:
                    continue
                score = parse_float(pred["solubility_score"])
                if score is None:
                    continue
                handle.write(f"{node_id}\t{score:.4f}\n")

        with grad_path.open("w", encoding="utf-8") as handle:
            handle.write("DATASET_COLORSTRIP\n")
            handle.write("SEPARATOR TAB\n")
            handle.write(f"DATASET_LABEL\t{args.tool_label} red-green\n")
            handle.write("COLOR\t#000000\n")
            handle.write("STRIP_WIDTH\t25\n")
            handle.write("MARGIN\t5\n")
            handle.write("BORDER_WIDTH\t1\n")
            handle.write("BORDER_COLOR\t#000000\n")
            handle.write("SHOW_STRIP_LABELS\t0\n")
            handle.write("LEGEND_TITLE\tSoluProt score\n")
            handle.write("LEGEND_SHAPES\t1\t1\t1\n")
            handle.write("LEGEND_COLORS\t#D7191C\t#F7F7A6\t#1B9E3C\n")
            handle.write("LEGEND_LABELS\tLow/Insoluble\tMiddle\tHigh/Soluble\n")
            handle.write("DATA\n")
            source_rows = nodes if nodes else [{"id": k} for k in by_node]
            for node in source_rows:
                node_id = node.get("id", "").strip()
                pred = by_node.get(node_id)
                if not pred:
                    continue
                score = parse_float(pred["solubility_score"])
                if score is None:
                    continue
                label = f"{score:.4f}"
                handle.write(f"{node_id}\t{score_color(score)}\t{label}\n")

    symbol_path = outdir / "itol_soluprot_gradient_symbols.txt"
    with symbol_path.open("w", encoding="utf-8") as handle:
        handle.write("DATASET_SYMBOL\n")
        handle.write("SEPARATOR TAB\n")
        handle.write(f"DATASET_LABEL\t{args.tool_label} gradient dots\n")
        handle.write("COLOR\t#000000\n")
        handle.write("MAXIMUM_SIZE\t18\n")
        handle.write("MARGIN\t8\n")
        handle.write("LEGEND_TITLE\tSoluProt score\n")
        handle.write("LEGEND_SHAPES\t2\t2\t2\n")
        handle.write("LEGEND_COLORS\t#D7191C\t#F7F7A6\t#1B9E3C\n")
        handle.write("LEGEND_LABELS\t<=0.2 low/faint\t0.5 middle\t1.0 high\n")
        handle.write("DATA\n")
        source_rows = nodes if nodes else [{"id": k} for k in by_node]
        for node in source_rows:
            node_id = node.get("id", "").strip()
            pred = by_node.get(node_id)
            if not pred:
                continue
            score = parse_float(pred["solubility_score"])
            if score is None:
                continue
            size = symbol_size(score, args.low_threshold)
            color = score_color_rescaled(score, args.low_threshold)
            label = f"{score:.4f}"
            handle.write(f"{node_id}\t2\t{size:.2f}\t{color}\t1\t1\t{label}\n")

    norm_path = outdir / "solubility_normalized.csv"
    with norm_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["id", "source_id", "solubility_score", "solubility_call"])
        writer.writeheader()
        for record in by_node.values():
            writer.writerow(record)

    if nodes:
        merged_path = outdir / "ssn_nodes_with_solubility.csv"
        headers = list(nodes[0].keys()) if nodes else []
        for col in ["solubility_score", "solubility_call"]:
            if col not in headers:
                headers.append(col)
        with merged_path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=headers)
            writer.writeheader()
            for node in nodes:
                pred = by_node.get(node.get("id", ""))
                out = dict(node)
                out["solubility_score"] = pred["solubility_score"] if pred else ""
                out["solubility_call"] = pred["solubility_call"] if pred else ""
                writer.writerow(out)

    if unmatched:
        unmatched_path = outdir / "unmatched_solubility_rows.csv"
        with unmatched_path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["id", "source_id", "solubility_score", "solubility_call"])
            writer.writeheader()
            writer.writerows(unmatched)

    if not args.no_color_strip:
        print(f"Wrote {strip_path}")
    if not args.symbols_only:
        print(f"Wrote {bar_path}")
        print(f"Wrote {grad_path}")
    print(f"Wrote {symbol_path}")
    print(f"Wrote {norm_path}")
    if nodes:
        print(f"Wrote {outdir / 'ssn_nodes_with_solubility.csv'}")
    if unmatched:
        print(f"Unmatched rows: {len(unmatched)}")


if __name__ == "__main__":
    main()
