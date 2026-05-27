#!/usr/bin/env python3
"""
Simplify duplicated UniProt-style tree labels and keep iTOL files in sync.

Example:
  A0ABV6D1V5|A0ABV6D1V5_9SPHN -> A0ABV6D1V5_9SPHN

If simplification would collide with an existing leaf name, a stable suffix is
added. The script writes a mapping CSV so every change is auditable.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def read_tree_ids(newick: str) -> List[str]:
    ids: List[str] = []
    token = []
    reading_length = False
    for ch in newick:
        if reading_length:
            if ch in ",);":
                reading_length = False
            continue
        if ch == ":":
            name = "".join(token).strip()
            if name and not re.fullmatch(r"\d+(?:\.\d+)?", name):
                ids.append(name)
            token = []
            reading_length = True
        elif ch in "(),;":
            token = []
        else:
            token.append(ch)
    return list(dict.fromkeys(ids))


def simplify_id(seq_id: str) -> str:
    text = seq_id.strip()
    if "|" in text:
        left, right = text.split("|", 1)
        if right == left:
            return left
        if right.startswith(left + "_"):
            return right
    return text


def unique_mapping(ids: Iterable[str]) -> Dict[str, str]:
    used = set()
    mapping: Dict[str, str] = {}
    for old in ids:
        base = simplify_id(old)
        new = base
        if new in used and new != old:
            suffix = 2
            while f"{base}_dup{suffix}" in used:
                suffix += 1
            new = f"{base}_dup{suffix}"
        mapping[old] = new
        used.add(new)
    return mapping


def replace_newick_ids(text: str, mapping: Dict[str, str]) -> str:
    for old in sorted(mapping, key=len, reverse=True):
        new = mapping[old]
        if old == new:
            continue
        text = text.replace(old + ":", new + ":")
        text = text.replace(old + ",", new + ",")
        text = text.replace(old + ")", new + ")")
    return text


def simplify_itol_file(src: Path, dst: Path, mapping: Dict[str, str]) -> None:
    with src.open("r", encoding="utf-8", errors="replace") as handle, dst.open("w", encoding="utf-8") as out:
        in_data = False
        for line in handle:
            if line.strip() == "DATA":
                in_data = True
                out.write(line)
                continue
            if not in_data or not line.strip() or line.startswith("#"):
                out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0] in mapping:
                parts[0] = mapping[parts[0]]
                out.write("\t".join(parts) + "\n")
            else:
                out.write(line)


def main() -> None:
    parser = argparse.ArgumentParser(description="Simplify duplicated UniProt tree labels and iTOL IDs.")
    parser.add_argument("tree", help="Input Newick tree.")
    parser.add_argument("--itol", nargs="*", default=[], help="iTOL annotation files to rewrite with simplified IDs.")
    parser.add_argument("--outdir", default=None, help="Output directory. Defaults to <tree_dir>/simplified_labels.")
    args = parser.parse_args()

    tree_path = Path(args.tree)
    outdir = Path(args.outdir) if args.outdir else tree_path.parent / "simplified_labels"
    outdir.mkdir(parents=True, exist_ok=True)

    text = tree_path.read_text(encoding="utf-8", errors="replace")
    ids = read_tree_ids(text)
    mapping = unique_mapping(ids)

    new_tree = replace_newick_ids(text, mapping)
    simplified_tree = outdir / f"simplified_{tree_path.name}"
    simplified_tree.write_text(new_tree, encoding="utf-8")

    with (outdir / "tree_label_mapping.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["old_id", "new_id", "changed"])
        for old in ids:
            writer.writerow([old, mapping[old], str(old != mapping[old]).lower()])

    for itol in args.itol:
        src = Path(itol)
        dst = outdir / f"simplified_{src.name}"
        simplify_itol_file(src, dst, mapping)

    changed = sum(1 for old in ids if old != mapping[old])
    print(f"Wrote {simplified_tree}")
    print(f"Wrote {outdir / 'tree_label_mapping.csv'}")
    print(f"Simplified labels: {changed}/{len(ids)}")


if __name__ == "__main__":
    main()
