#!/usr/bin/env python
"""Collect a kinase-biased ChEMBL SMILES reference table.

Requires chembl_webresource_client. This script intentionally writes only SMILES
and light metadata so the refinement script can use it as a structural prior.
"""

from __future__ import annotations

import argparse
import csv
import sys
from time import monotonic
from pathlib import Path


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = ["molecule_chembl_id", "canonical_smiles", "target_chembl_id", "target_pref_name", "pchembl_value", "standard_type"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--max-targets", type=int, default=100)
    parser.add_argument("--max-activities", type=int, default=5000)
    parser.add_argument("--min-pchembl", type=float, default=6.0)
    args = parser.parse_args()

    try:
        from chembl_webresource_client.new_client import new_client
    except Exception as exc:
        raise SystemExit(
            "chembl_webresource_client is required. Install with: "
            "pip install chembl_webresource_client"
        ) from exc

    start = monotonic()
    print(
        f"[chembl] searching targets term=kinase max_targets={args.max_targets} "
        f"max_activities={args.max_activities} min_pchembl={args.min_pchembl}",
        flush=True,
    )
    targets = list(new_client.target.search("kinase"))[: args.max_targets]
    print(f"[chembl] targets_loaded={len(targets)}", flush=True)
    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    for index, target in enumerate(targets, start=1):
        target_id = target.get("target_chembl_id")
        if not target_id:
            print(f"[chembl] target {index}/{len(targets)} missing target_chembl_id; skipped", flush=True)
            continue
        target_name = target.get("pref_name", "")
        print(f"[chembl] target {index}/{len(targets)} {target_id} {target_name}", flush=True)
        acts = new_client.activity.filter(
            target_chembl_id=target_id,
            pchembl_value__gte=args.min_pchembl,
        ).only(
            "molecule_chembl_id",
            "canonical_smiles",
            "pchembl_value",
            "standard_type",
        )
        target_rows = 0
        for act in acts:
            smiles = act.get("canonical_smiles")
            mol_id = act.get("molecule_chembl_id")
            if not smiles or not mol_id or mol_id in seen:
                continue
            seen.add(mol_id)
            target_rows += 1
            rows.append(
                {
                    "molecule_chembl_id": mol_id,
                    "canonical_smiles": smiles,
                    "target_chembl_id": target_id,
                    "target_pref_name": target_name,
                    "pchembl_value": act.get("pchembl_value", ""),
                    "standard_type": act.get("standard_type", ""),
                }
            )
            if target_rows % 50 == 0:
                print(
                    f"[chembl] target {index}/{len(targets)} rows={target_rows} total={len(rows)}",
                    flush=True,
                )
            if len(rows) >= args.max_activities:
                write_csv(args.out.expanduser(), rows)
                elapsed = monotonic() - start
                print(f"[chembl] wrote {args.out.expanduser()} rows={len(rows)} elapsed_s={elapsed:.1f}", flush=True)
                return
        print(f"[chembl] target {index}/{len(targets)} done rows_added={target_rows} total={len(rows)}", flush=True)
    write_csv(args.out.expanduser(), rows)
    elapsed = monotonic() - start
    print(f"[chembl] wrote {args.out.expanduser()} rows={len(rows)} elapsed_s={elapsed:.1f}", flush=True)


if __name__ == "__main__":
    main()
