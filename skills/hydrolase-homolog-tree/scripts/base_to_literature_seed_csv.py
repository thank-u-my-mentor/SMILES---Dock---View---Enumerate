#!/usr/bin/env python3
"""
Convert seed provenance CSV into a Zotero-like CSV accepted by literature-agent.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def first_value(row: Dict[str, str], names: Iterable[str]) -> str:
    lower = {k.lower(): k for k in row}
    for name in names:
        key = lower.get(name.lower())
        if key and str(row.get(key, "")).strip():
            return str(row[key]).strip()
    return ""


def main() -> None:
    parser = argparse.ArgumentParser(description="Make a literature-agent seed CSV from hydrolase base/provenance CSV.")
    parser.add_argument("seed_csv", help="seed_provenance_template_v2.csv or seed_provenance_uniprot_enriched.csv")
    parser.add_argument("--output", default="hydrolase_literature_seed.csv")
    args = parser.parse_args()

    rows = read_csv(Path(args.seed_csv))
    headers = ["Title", "Author", "Publication Year", "DOI", "Abstract Note", "Publication Title", "Url", "Tags"]
    out_rows: List[Dict[str, str]] = []
    for row in rows:
        enzyme = first_value(row, ["enzyme_name", "enzyme", "title"])
        organism = first_value(row, ["source_organism", "organism"])
        doi = first_value(row, ["key_doi", "doi", "DOI", "uniprot_candidate_doi"])
        title = first_value(row, ["paper_title"])
        if not title:
            title = f"{enzyme} {organism} expression plasmid His-tag"
        note_parts = [
            f"enzyme={enzyme}",
            f"organism={organism}",
            f"uniprot={first_value(row, ['uniprot_id', 'Uniprot id'])}",
            f"promiscuity={first_value(row, ['known_promiscuity', 'promiscuity'])}",
            f"seed_notes={first_value(row, ['seed_notes', 'notes'])}",
            "review_goal=extract expression host, E. coli strain, pET vector, His-tag position, induction, purification, and methods/supplement evidence",
        ]
        out_rows.append({
            "Title": title,
            "Author": "",
            "Publication Year": "",
            "DOI": doi,
            "Abstract Note": "\n".join(p for p in note_parts if p and not p.endswith("=")),
            "Publication Title": "",
            "Url": f"https://doi.org/{doi}" if doi else "",
            "Tags": "hydrolase_seed;expression_provenance",
        })

    with Path(args.output).open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"Wrote {args.output} with {len(out_rows)} rows.")


if __name__ == "__main__":
    main()
