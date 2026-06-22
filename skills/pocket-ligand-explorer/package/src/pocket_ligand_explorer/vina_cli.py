#!/usr/bin/env python
"""Small Vina-compatible CLI backed by the PyPI vina Python API."""

from __future__ import annotations

import argparse
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receptor", required=True)
    parser.add_argument("--ligand", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--log")
    parser.add_argument("--center_x", type=float, required=True)
    parser.add_argument("--center_y", type=float, required=True)
    parser.add_argument("--center_z", type=float, required=True)
    parser.add_argument("--size_x", type=float, required=True)
    parser.add_argument("--size_y", type=float, required=True)
    parser.add_argument("--size_z", type=float, required=True)
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--cpu", type=int, default=0)
    parser.add_argument("--num_modes", "--num-modes", dest="num_modes", type=int, default=9)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--energy_range", "--energy-range", dest="energy_range", type=float, default=3.0)
    return parser


def format_log(energies: object) -> str:
    lines = [
        "-----+------------+----------+----------",
        "mode |   affinity | dist from best mode",
        "     | (kcal/mol) | rmsd l.b.| rmsd u.b.",
        "-----+------------+----------+----------",
    ]
    try:
        rows = list(energies)
    except Exception:
        rows = []
    for index, row in enumerate(rows, start=1):
        values = list(row)
        affinity = float(values[0]) if values else 0.0
        rmsd_lb = float(values[1]) if len(values) > 1 else 0.0
        rmsd_ub = float(values[2]) if len(values) > 2 else 0.0
        lines.append(f"{index:>4} {affinity:>11.3f} {rmsd_lb:>10.3f} {rmsd_ub:>10.3f}")
    return "\n".join(lines) + "\n"


def main() -> None:
    args = build_parser().parse_args()
    try:
        from vina import Vina
    except Exception as exc:
        raise SystemExit(f"PyPI vina package is not available: {exc}") from exc

    out = Path(args.out).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log).expanduser() if args.log else None
    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)

    v = Vina(sf_name="vina", cpu=int(args.cpu), seed=int(args.seed), verbosity=1)
    v.set_receptor(str(Path(args.receptor).expanduser()))
    v.set_ligand_from_file(str(Path(args.ligand).expanduser()))
    v.compute_vina_maps(
        center=[args.center_x, args.center_y, args.center_z],
        box_size=[args.size_x, args.size_y, args.size_z],
    )
    v.dock(exhaustiveness=int(args.exhaustiveness), n_poses=int(args.num_modes))
    v.write_poses(str(out), n_poses=int(args.num_modes), energy_range=float(args.energy_range), overwrite=True)
    log_text = format_log(v.energies(n_poses=int(args.num_modes), energy_range=float(args.energy_range)))
    if log_path:
        log_path.write_text(log_text, encoding="utf-8")
    print(log_text, end="")


if __name__ == "__main__":
    main()
