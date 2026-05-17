#!/usr/bin/env python
"""Report or remove chemically suspicious SMILES rows from dock_history.csv."""

from __future__ import annotations

import argparse
import csv
import shutil
from datetime import datetime
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase

RDLogger.DisableLog("rdApp.*")
rdBase.DisableLog("rdApp.*")


def read_csv(path: Path) -> list[dict[str, str]]:
    for encoding in ("utf-8-sig", "utf-8", "gb18030", "cp936", "latin1"):
        try:
            with path.open("r", encoding=encoding, newline="") as handle:
                return list(csv.DictReader(handle))
        except UnicodeDecodeError:
            continue
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles or "").strip())
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def sanity_reasons(smiles: str) -> list[str]:
    mol = Chem.MolFromSmiles(str(smiles or "").strip())
    if mol is None:
        return ["rdkit_parse_failed"]
    reasons: list[str] = []
    if len(Chem.GetMolFrags(mol)) != 1:
        reasons.append("multiple_fragments")
    if any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
        reasons.append("radical_atoms")
    if any(abs(atom.GetFormalCharge()) > 1 for atom in mol.GetAtoms()):
        reasons.append("large_atom_formal_charge")
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        if bond.IsInRing() and not bond.GetIsAromatic() and a.GetIsAromatic() and b.GetIsAromatic():
            reasons.append("nonaromatic_ring_bond_between_aromatic_atoms")
            break
    rings = mol.GetRingInfo().AtomRings()
    if rings and max(len(ring) for ring in rings) > 8:
        reasons.append("large_ring_gt8")
    return reasons


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--history-csv", type=Path, default=Path("~/vina_task2/dock_history/dock_history.csv"))
    parser.add_argument("--report-csv", type=Path, default=Path("~/vina_task2/dock_history/smiles_sanity_rejects.csv"))
    parser.add_argument("--apply", action="store_true", help="Rewrite dock_history.csv after creating a timestamped backup")
    parser.add_argument("--smiles", help="Optional single SMILES sanity check without reading history")
    args = parser.parse_args()

    if args.smiles:
        reasons = sanity_reasons(args.smiles)
        print(f"canonical={canonical(args.smiles)}")
        print(f"status={'reject' if reasons else 'ok'} reasons={';'.join(reasons)}")
        return

    history_csv = args.history_csv.expanduser().resolve()
    rows = read_csv(history_csv)
    if not rows:
        raise SystemExit(f"No rows found in {history_csv}")
    fields = list(rows[0].keys())
    rejects: list[dict[str, str]] = []
    keep: list[dict[str, str]] = []
    seen_can: dict[str, str] = {}
    for row in rows:
        smiles = row.get("canonical_smiles") or row.get("smiles") or row.get("input_smiles") or ""
        can = canonical(smiles)
        reasons = sanity_reasons(can or smiles)
        if can and can in seen_can and seen_can[can] != row.get("seq_id", ""):
            reasons.append(f"duplicate_canonical_of_{seen_can[can]}")
        elif can:
            seen_can[can] = row.get("seq_id", "")
        if reasons:
            rejected = dict(row)
            rejected["canonical_smiles_checked"] = can
            rejected["smiles_sanity_reasons"] = ";".join(reasons)
            rejects.append(rejected)
        else:
            keep.append(row)

    report_csv = args.report_csv.expanduser().resolve()
    report_fields = [*fields, "canonical_smiles_checked", "smiles_sanity_reasons"]
    report_csv.parent.mkdir(parents=True, exist_ok=True)
    write_csv(report_csv, rejects, report_fields)
    print(f"checked={len(rows)} reject={len(rejects)} keep={len(keep)} report={report_csv}")
    if not args.apply:
        print("dry_run=true; rerun with --apply to rewrite dock_history.csv after backup")
        return
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = history_csv.with_name(f"{history_csv.stem}.backup_before_sanity_clean_{stamp}{history_csv.suffix}")
    shutil.copyfile(history_csv, backup)
    write_csv(history_csv, keep, fields)
    print(f"applied=true backup={backup} wrote={history_csv}")


if __name__ == "__main__":
    main()
