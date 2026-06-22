#!/usr/bin/env python3
"""Audit CCD components that contain target metal elements.

This script classifies formula-level metal-containing CCD components into rough
human-readable buckets. It is intentionally heuristic: use it to explain why
formula-based metal component counts are much larger than true enzyme metal
center counts, not as a chemical ontology.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd


SKILL_SCRIPT = Path(__file__).resolve().with_name("metal_f_element_ligand_rescan.py")
sys.path.insert(0, str(SKILL_SCRIPT.parent))
import metal_f_element_ligand_rescan as rescan  # noqa: E402


BIOLOGICAL_CENTER_IDS = {
    "HEM",
    "HEA",
    "HEB",
    "HEC",
    "HEO",
    "HEV",
    "HAS",
    "SRM",
    "SF4",
    "F3S",
    "FES",
    "FS3",
    "FS4",
    "FEO",
    "F43",
}
SIMPLE_ION_IDS = {"FE", "FE2", "F3S", "CO", "CO3", "NI", "CU", "CU1", "CUA", "MN", "MN3", "CR"}
SALT_WORDS = {"CHLORIDE", "SULFATE", "PHOSPHATE", "NITRATE", "ACETATE", "FORMATE", "OXALATE", "CYANIDE"}
SYNTHETIC_WORDS = {
    "BIS",
    "TRIS",
    "COMPLEX",
    "PORPHYRIN",
    "PHTHALOCYANINE",
    "FERROCENE",
    "TERPYRIDINE",
    "BIPYRIDINE",
    "PHENANTHROLINE",
}


def classify(row: dict) -> str:
    comp_id = row["comp_id"].upper()
    name = row.get("name", "").upper()
    formula = row.get("formula", "")
    elements = set((row.get("target_metal_elements") or "").split(";"))
    heavy_nonmetals = [
        elem
        for elem in rescan.formula_elements(formula)
        if elem not in elements and elem not in {"H", "C", "N", "O", "S", "P"}
    ]
    if comp_id in BIOLOGICAL_CENTER_IDS:
        return "known_biological_metal_cofactor"
    if any(word in name for word in SYNTHETIC_WORDS):
        return "synthetic_coordination_complex_or_probe"
    if comp_id in SIMPLE_ION_IDS or (row.get("type") == "non-polymer" and len(formula.split()) <= 4):
        return "simple_metal_ion_or_small_salt"
    if any(word in name for word in SALT_WORDS):
        return "metal_salt_or_counterion"
    if len(heavy_nonmetals) >= 2:
        return "heteroatom_rich_synthetic_or_special_ligand"
    if len(formula.split()) > 20:
        return "large_metal_organic_component"
    return "ambiguous_metal_component"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="")
    parser.add_argument("--examples-per-class", type=int, default=12)
    args = parser.parse_args()

    ccd = rescan.parse_ccd(rescan.download_ccd())
    rows = []
    for row in ccd.values():
        if not row["target_metal_elements"]:
            continue
        item = dict(row)
        item["audit_class"] = classify(item)
        rows.append(item)
    frame = pd.DataFrame(rows).sort_values(["audit_class", "target_metal_elements", "comp_id"])

    print(f"target-metal CCD components: {len(frame)}")
    print("\nBy metal element:")
    counter = Counter()
    for value in frame["target_metal_elements"]:
        for elem in str(value).split(";"):
            if elem:
                counter[elem] += 1
    for elem, count in counter.most_common():
        print(f"  {elem}: {count}")

    print("\nBy audit class:")
    for klass, count in frame["audit_class"].value_counts().items():
        print(f"  {klass}: {count}")

    print("\nExamples:")
    for klass, sub in frame.groupby("audit_class", sort=True):
        print(f"\n### {klass}")
        for _, row in sub.head(args.examples_per_class).iterrows():
            name = re.sub(r"\s+", " ", str(row["name"]))
            print(
                f"{row['comp_id']:>4} | {row['target_metal_elements']:<8} | "
                f"{row['formula']:<32} | {name[:90]}"
            )

    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        if out.suffix.lower() == ".xlsx":
            with pd.ExcelWriter(out, engine="openpyxl") as writer:
                frame.to_excel(writer, sheet_name="metal_components", index=False)
                frame["audit_class"].value_counts().rename_axis("audit_class").reset_index(name="count").to_excel(
                    writer, sheet_name="class_counts", index=False
                )
        else:
            frame.to_csv(out, index=False)
        print(f"\nWrote: {out}")


if __name__ == "__main__":
    main()
