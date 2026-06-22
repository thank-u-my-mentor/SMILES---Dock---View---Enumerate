#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from pathlib import Path


def load_fetch_module():
    path = Path(__file__).resolve().parent / "fetch_hmmer_hmmsearch.py"
    spec = importlib.util.spec_from_file_location("fetch_hmmer_hmmsearch", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["fetch_hmmer_hmmsearch"] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def read_hits(path: Path, module):
    rows = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(module.Hit(
                accession=row.get("accession", ""),
                evalue=module.maybe_float(row.get("evalue")),
                score=module.maybe_float(row.get("score")),
                description=row.get("description", ""),
                organism=row.get("organism", ""),
                kingdom=row.get("kingdom", ""),
                phylum=row.get("phylum", ""),
                taxid=row.get("taxid", ""),
            ))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Refresh homolog_metadata.csv from existing FASTA and hmmer_hits.csv.")
    parser.add_argument("--fasta", default="/mnt/e/Tree-Metal-F/hmmer/homologs_plus_core.fasta")
    parser.add_argument("--hits", default="/mnt/e/Tree-Metal-F/hmmer/hmmer_hits.csv")
    parser.add_argument("--out", default="/mnt/e/Tree-Metal-F/hmmer/homolog_metadata.csv")
    args = parser.parse_args()

    module = load_fetch_module()
    records = module.read_fasta(Path(args.fasta))
    hits = read_hits(Path(args.hits), module)
    module.write_metadata_csv(module.fasta_metadata(records, hits), Path(args.out))
    print(f"Wrote {args.out} for {len(records)} FASTA records using {len(hits)} HMMER hits.")


if __name__ == "__main__":
    main()
