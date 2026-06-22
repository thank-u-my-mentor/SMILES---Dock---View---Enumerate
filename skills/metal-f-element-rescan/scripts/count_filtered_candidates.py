#!/usr/bin/env python3
"""Count Metal-F candidate entries after CCD metal audit and metadata filters."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import audit_metal_components as audit  # noqa: E402
import metal_f_element_ligand_rescan as rescan  # noqa: E402


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DEFAULT_EXCLUDED = {"ambiguous_metal_component"}


def post_search(query: dict, rows: int = 1000) -> set[str]:
    out: set[str] = set()
    start = 0
    while True:
        payload = {
            "query": query,
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": start, "rows": rows},
                "results_content_type": ["experimental"],
            },
        }
        request = urllib.request.Request(
            SEARCH_URL,
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json", "User-Agent": "metal-f-count/1.0"},
        )
        last_error = None
        for attempt in range(5):
            try:
                with urllib.request.urlopen(request, timeout=120) as response:
                    data = json.loads(response.read().decode("utf-8"))
                break
            except (urllib.error.URLError, TimeoutError, ConnectionResetError) as exc:
                last_error = exc
                time.sleep(1.5 * (attempt + 1))
        else:
            raise last_error
        result_set = data.get("result_set") or []
        out.update(row["identifier"].upper() for row in result_set if row.get("identifier"))
        if start + len(result_set) >= int(data.get("total_count") or 0) or not result_set:
            break
        start += len(result_set)
    return out


def terminal(attribute: str, operator: str, value, *, negation: bool = False) -> dict:
    parameters = {"attribute": attribute, "operator": operator, "value": value}
    if negation:
        parameters["negation"] = True
    return {"type": "terminal", "service": "text", "parameters": parameters}


def group(nodes: list[dict], logical_operator: str = "and") -> dict:
    return {"type": "group", "logical_operator": logical_operator, "nodes": nodes}


def comp_id_entries(comp_ids: list[str], chunk_size: int) -> set[str]:
    out: set[str] = set()
    comp_ids = sorted(set(comp_ids))
    total = (len(comp_ids) + chunk_size - 1) // chunk_size
    for i in range(0, len(comp_ids), chunk_size):
        chunk = comp_ids[i : i + chunk_size]
        print(f"comp_id chunk {i // chunk_size + 1}/{total}: {len(chunk)}", flush=True)
        out |= post_search(terminal("rcsb_chem_comp_container_identifiers.comp_id", "in", chunk))
    return out


def metadata_filtered_entries(base_ids: set[str], *, no_nucleic_acid: bool, protein_entity_one: bool) -> set[str]:
    nodes = [terminal("rcsb_entry_container_identifiers.entry_id", "in", sorted(base_ids))]
    if no_nucleic_acid:
        nodes.extend(
            [
                terminal("rcsb_entry_info.polymer_entity_count_DNA", "equals", 0),
                terminal("rcsb_entry_info.polymer_entity_count_RNA", "equals", 0),
                terminal("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
            ]
        )
    if protein_entity_one:
        nodes.append(terminal("rcsb_entry_info.polymer_entity_count_protein", "equals", 1))
    return post_search(group(nodes))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exclude-classes", default=",".join(sorted(DEFAULT_EXCLUDED)))
    parser.add_argument("--chunk-size", type=int, default=100)
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    excluded = {item.strip() for item in args.exclude_classes.split(",") if item.strip()}
    ccd = rescan.parse_ccd(rescan.download_ccd())
    f_components = [cid for cid, row in ccd.items() if row["has_f"]]
    metal_rows = []
    kept_metal_ids = []
    excluded_metal_ids = []
    for row in ccd.values():
        if not row["target_metal_elements"]:
            continue
        item = dict(row)
        item["audit_class"] = audit.classify(item)
        metal_rows.append(item)
        if item["audit_class"] in excluded:
            excluded_metal_ids.append(item["comp_id"])
        else:
            kept_metal_ids.append(item["comp_id"])

    print(f"CCD F components: {len(f_components)}")
    print(f"CCD metal components total: {len(metal_rows)}")
    print(f"CCD metal components kept: {len(kept_metal_ids)}")
    print(f"CCD metal components excluded: {len(excluded_metal_ids)} ({';'.join(sorted(excluded))})")

    f_entries = comp_id_entries(f_components, args.chunk_size)
    metal_entries_kept = comp_id_entries(kept_metal_ids, args.chunk_size)
    raw_intersection = f_entries & metal_entries_kept
    no_na = metadata_filtered_entries(raw_intersection, no_nucleic_acid=True, protein_entity_one=False)
    no_na_one_protein = metadata_filtered_entries(raw_intersection, no_nucleic_acid=True, protein_entity_one=True)

    summary = [
        {"metric": "ccd_f_component_count", "value": len(f_components)},
        {"metric": "ccd_metal_component_total", "value": len(metal_rows)},
        {"metric": "ccd_metal_component_kept", "value": len(kept_metal_ids)},
        {"metric": "ccd_metal_component_excluded", "value": len(excluded_metal_ids)},
        {"metric": "excluded_classes", "value": ";".join(sorted(excluded))},
        {"metric": "f_entry_count", "value": len(f_entries)},
        {"metric": "kept_metal_entry_count", "value": len(metal_entries_kept)},
        {"metric": "intersection_entry_count", "value": len(raw_intersection)},
        {"metric": "intersection_no_dna_rna_count", "value": len(no_na)},
        {"metric": "intersection_no_dna_rna_one_protein_entity_count", "value": len(no_na_one_protein)},
    ]
    for row in summary:
        print(f"{row['metric']}: {row['value']}")

    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        metal_frame = pd.DataFrame(metal_rows).sort_values(["audit_class", "comp_id"])
        frames = {
            "summary": pd.DataFrame(summary),
            "metal_components": metal_frame,
            "candidate_pdb_ids": pd.DataFrame({"pdb_id": sorted(raw_intersection)}),
            "candidate_no_dna_rna": pd.DataFrame({"pdb_id": sorted(no_na)}),
            "candidate_no_dna_rna_one_protein": pd.DataFrame({"pdb_id": sorted(no_na_one_protein)}),
        }
        if out.suffix.lower() == ".xlsx":
            with pd.ExcelWriter(out, engine="openpyxl") as writer:
                for sheet, frame in frames.items():
                    frame.to_excel(writer, sheet_name=sheet[:31], index=False)
        else:
            with open(out, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["metric", "value"])
                writer.writeheader()
                writer.writerows(summary)
        print(f"Wrote: {out}")


if __name__ == "__main__":
    main()
