#!/usr/bin/env python3
"""Search RCSB for fluorinated organic ligands co-occurring with transition metals.

This only builds candidate entry sets. Geometric claims still require mmCIF
coordinate filtering downstream.
"""

import argparse
import csv
import json
import os
import re
import urllib.request


ENDPOINT = "https://search.rcsb.org/rcsbsearch/v2/query"
DEFAULT_PROJECT = "/home/qin/Metal-F_project"
DEFAULT_METALS = ["FE", "CO", "NI", "CU", "MN"]

MOTIFS = {
    "CF3": {
        "description": "trifluoromethyl substructure",
        "smiles": "C(F)(F)F",
        "match_type": "sub-struct-graph-relaxed",
    },
    "C-F": {
        "description": "organic carbon-fluorine bond substructure",
        "smiles": "CF",
        "match_type": "sub-struct-graph-relaxed",
    },
    "fluoroaryl": {
        "description": "aryl fluoride / fluorobenzene-like substructure",
        "smiles": "c1ccccc1F",
        "match_type": "sub-struct-graph-relaxed",
    },
}


def safe_name(text):
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip())
    return value.strip("._-") or "fluoro_ligand"


def build_payload(smiles, metals, match_type="sub-struct-graph-relaxed", rows=1000):
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "chemical",
                    "parameters": {
                        "type": "descriptor",
                        "descriptor_type": "SMILES",
                        "match_type": match_type,
                        "value": smiles,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_instance_annotation.comp_id",
                        "operator": "in",
                        "value": metals,
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": rows},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
            "scoring_strategy": "combined",
        },
    }


def post(payload):
    data = json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(
        ENDPOINT,
        data=data,
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        result = json.load(response)
    return result


def result_rows(result):
    rows = []
    for row in result.get("result_set", []):
        rows.append({"pdb_id": row["identifier"], "score": row.get("score", "")})
    return rows


def write_csv(path, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pdb_id", "score"])
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", default="CF3", help="Output name: creates Metal-<name> directory.")
    parser.add_argument("--smiles", default="", help="SMILES query. Overrides --name motif if set.")
    parser.add_argument(
        "--match-type",
        default="",
        help="RCSB chemical match_type. Default uses motif setting or sub-struct-graph-relaxed.",
    )
    parser.add_argument("--project", default=DEFAULT_PROJECT)
    parser.add_argument("--metals", default=",".join(DEFAULT_METALS))
    parser.add_argument("--rows", type=int, default=1000)
    args = parser.parse_args()

    motif = MOTIFS.get(args.name)
    smiles = args.smiles or (motif["smiles"] if motif else args.name)
    match_type = args.match_type or (motif["match_type"] if motif else "sub-struct-graph-relaxed")
    metals = [x.strip().upper() for x in args.metals.split(",") if x.strip()]

    output_name = f"Metal-{safe_name(args.name)}"
    outdir = os.path.join(os.path.expanduser(args.project), output_name)
    os.makedirs(outdir, exist_ok=True)

    payload = build_payload(smiles, metals, match_type, args.rows)
    result = post(payload)
    rows = result_rows(result)

    payload_path = os.path.join(outdir, f"{output_name}_rcsb_query.json")
    result_path = os.path.join(outdir, f"{output_name}_rcsb_result.json")
    csv_path = os.path.join(outdir, f"{output_name}_pdb_ids.csv")
    ids_path = os.path.join(outdir, f"{output_name}_pdb_ids.txt")
    readme_path = os.path.join(outdir, "README.md")

    with open(payload_path, "w") as handle:
        json.dump(payload, handle, indent=2)
    with open(result_path, "w") as handle:
        json.dump(result, handle, indent=2)
    write_csv(csv_path, rows)
    with open(ids_path, "w") as handle:
        handle.write("\n".join(row["pdb_id"] for row in rows) + ("\n" if rows else ""))
    with open(readme_path, "w") as handle:
        handle.write(
            f"# {output_name}\n\n"
            f"RCSB candidate search for fluorinated organic ligand motif `{smiles}` "
            f"and transition metal components `{', '.join(metals)}`.\n\n"
            f"Result count: {result.get('total_count', len(rows))}\n\n"
            "This is a metadata/substructure candidate set. Download mmCIF files and "
            "calculate ligand-F-to-metal distances before making geometric claims.\n\n"
            "## Hard Geometry Rule\n\n"
            "Downstream filtering uses a fixed cutoff: keep only entries where at least one "
            "organic ligand fluorine atom is within <= 3.5 Angstrom of a FE/CO/NI/CU/MN atom "
            "in the mmCIF coordinates. The API query itself is not distance evidence.\n\n"
            "## FASTA Rule\n\n"
            "For HMMER/hydrolase-homolog-tree handoff, only protein polymer entities with "
            "sequence length >= 150 aa are written to the HMMER seed FASTA by default. "
            "Shorter protein fragments/domains are kept in an excluded CSV.\n"
        )

    print(f"name: {output_name}")
    print(f"motif: {smiles}")
    print(f"metals: {', '.join(metals)}")
    print(f"count: {result.get('total_count', len(rows))}")
    print(f"outdir: {outdir}")
    print("first_ids:", " ".join(row["pdb_id"] for row in rows[:25]))


if __name__ == "__main__":
    main()
