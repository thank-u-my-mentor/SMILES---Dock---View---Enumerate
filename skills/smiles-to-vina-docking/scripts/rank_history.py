#!/usr/bin/env python
"""Rank dock_history.csv against a chosen reference SMILES."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path

import dock_utils as dl


RANK_FIELDS = [
    "rank",
    "seq_id",
    "nickname",
    "ancestor_smiles",
    "parent_smiles",
    "edit_label",
    "smiles",
    "affinity_kcal_mol",
    "official_binding_score",
    "reference_rmsd",
    "reference_near_mode_count",
    "mode_count",
    "inner_rmsd",
    "inner_cluster_fraction",
    "whole_rmsd",
    "hbond_count",
    "hydrophobic_count",
    "vdw_contact_count",
    "pi_contact_count",
    "ch_pi_count",
    "affinity_score",
    "reference_rmsd_score",
    "inner_score",
    "whole_score",
    "interaction_score",
    "final_score",
    "pose_path",
    "log_path",
]


def read_history(history_dir: Path) -> list[dict[str, str]]:
    return dl.read_csv(dl.ledger_path(history_dir))


def write_rank(path: Path, rows: list[dict[str, str | int | float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RANK_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in RANK_FIELDS})


def find_row(rows: list[dict[str, str]], smiles: str) -> dict[str, str] | None:
    target = dl.smiles_identity(smiles, smiles)
    matches = [row for row in rows if dl.same_molecule(dl.row_identity(row), target)]
    matches.sort(key=lambda row: float(row.get("affinity_kcal_mol") or 999.0))
    return matches[0] if matches else None


def ensure_reference(args: argparse.Namespace) -> dict[str, str]:
    rows = read_history(args.history_dir)
    found = find_row(rows, args.smiles)
    if found:
        return found
    print("reference_missing=true; docking it first", flush=True)
    cmd = [
        sys.executable,
        str(Path(__file__).with_name("dock_smiles.py")),
        "--smiles",
        args.smiles,
        "--history-dir",
        str(args.history_dir),
    ]
    for name in ("receptor", "config", "vina", "meeko", "obabel", "cpu", "exhaustiveness", "num_modes", "energy_range"):
        value = getattr(args, name, None)
        if value not in (None, ""):
            cmd.extend([f"--{name.replace('_', '-')}", str(value)])
    proc = subprocess.run(cmd, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"dock_smiles.py failed with exit code {proc.returncode}")
    found = find_row(read_history(args.history_dir), args.smiles)
    if not found:
        raise RuntimeError("reference was docked but not found in history")
    return found


def normal_scores(values: list[float | None], lower_is_better: bool) -> list[float | None]:
    finite = [value for value in values if value is not None and math.isfinite(value)]
    if not finite:
        return [None for _ in values]
    mean = sum(finite) / len(finite)
    stdev = math.sqrt(sum((value - mean) ** 2 for value in finite) / len(finite))
    if stdev < 1e-9:
        return [0.5 if value is not None else None for value in values]
    out: list[float | None] = []
    for value in values:
        if value is None or not math.isfinite(value):
            out.append(None)
            continue
        z = (value - mean) / stdev
        cdf = 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
        out.append(round((1.0 - cdf) if lower_is_better else cdf, 6))
    return out


def reference_coords(row: dict[str, str]) -> list[tuple[float, float, float]]:
    pose = Path(row.get("pose_path", "")).expanduser()
    modes = dl.parse_pdbqt_modes(pose)
    if not modes:
        return []
    best = min(modes, key=lambda mode: mode.affinity if mode.affinity is not None else float("inf"))
    return best.coords


def rmsd_to_reference(row: dict[str, str], ref: list[tuple[float, float, float]], cutoff: float) -> dict[str, str | int | float]:
    modes = dl.parse_pdbqt_modes(Path(row.get("pose_path", "")).expanduser())
    if not modes or not ref:
        return {"reference_rmsd": "", "reference_near_mode_count": "", "mode_count": row.get("mode_count", "")}
    values = [dl.coord_rmsd(ref, mode.coords) for mode in modes]
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return {"reference_rmsd": "", "reference_near_mode_count": "", "mode_count": len(modes)}
    return {
        "reference_rmsd": round(min(finite), 4),
        "reference_near_mode_count": sum(1 for value in finite if value <= cutoff),
        "mode_count": len(modes),
    }


def rank(args: argparse.Namespace) -> Path:
    args.history_dir = args.history_dir.expanduser().resolve()
    reference = ensure_reference(args)
    ref = reference_coords(reference)
    if not ref:
        raise RuntimeError(f"reference pose is missing/unreadable: {reference.get('pose_path')}")
    rows = read_history(args.history_dir)
    enriched: list[dict[str, str | int | float]] = []
    for row in rows:
        enriched.append({**row, **rmsd_to_reference(row, ref, args.near_rmsd_cutoff)})

    affinity = [float(row["affinity_kcal_mol"]) if row.get("affinity_kcal_mol") else None for row in enriched]
    reference_rmsd = [float(row["reference_rmsd"]) if row.get("reference_rmsd") not in ("", None) else None for row in enriched]
    inner = [float(row["inner_rmsd"]) if row.get("inner_rmsd") not in ("", None) else None for row in enriched]
    whole = [float(row["whole_rmsd"]) if row.get("whole_rmsd") not in ("", None) else None for row in enriched]
    interactions: list[float | None] = []
    for row in enriched:
        vals = [row.get("hbond_count"), row.get("hydrophobic_count"), row.get("pi_contact_count"), row.get("ch_pi_count")]
        nums = [float(v) for v in vals if v not in ("", None)]
        interactions.append(sum(nums) if nums else None)

    score_columns = {
        "affinity_score": normal_scores(affinity, True),
        "reference_rmsd_score": normal_scores(reference_rmsd, True),
        "inner_score": normal_scores(inner, True),
        "whole_score": normal_scores(whole, True),
        "interaction_score": normal_scores(interactions, False),
    }
    for i, row in enumerate(enriched):
        scores = []
        for key, values in score_columns.items():
            value = values[i]
            row[key] = "" if value is None else value
            if value is not None:
                scores.append(value)
        row["final_score"] = round(sum(scores) / len(scores), 6) if scores else ""

    enriched.sort(key=lambda row: (-(float(row["final_score"]) if row.get("final_score") else -1), float(row["affinity_kcal_mol"] or 999.0)))
    for i, row in enumerate(enriched, start=1):
        row["rank"] = i
    out = args.output.expanduser().resolve() if args.output else args.history_dir / "ranked_history.csv"
    write_rank(out, enriched)
    print(f"reference_seq_id={reference.get('seq_id')} wrote {out} rows={len(enriched)}", flush=True)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smiles", required=True)
    parser.add_argument("--history-dir", "--ledger-dir", dest="history_dir", metavar="HISTORY_DIR", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--receptor", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--vina", default="vina")
    parser.add_argument("--meeko", default="mk_prepare_ligand.py")
    parser.add_argument("--obabel", default="obabel")
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--exhaustiveness", type=int)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument("--energy-range", type=int)
    parser.add_argument("--near-rmsd-cutoff", type=float, default=4.0)
    args = parser.parse_args()
    rank(args)


if __name__ == "__main__":
    main()
