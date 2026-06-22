#!/usr/bin/env python3
"""Compare paired GAMESS spin-screen optimizations for Fe radical MCPB setup."""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

HARTREE_TO_KCAL = 627.509474
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from export_gamess_optimized_model import parse_coords_for_nserch  # noqa: E402


@dataclass
class Atom:
    element: str
    xyz: tuple[float, float, float]
    note: str = ""


@dataclass
class LogSummary:
    label: str
    path: Path
    normal: bool
    not_converged: bool
    energy: float | None
    nserch: int | None
    grad_max: float | None
    grad_rms: float | None
    s2: float | None
    wall_seconds: float | None
    coords: list[Atom]


def parse_xyz(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    lines = path.read_text(errors="replace").splitlines()
    for line in lines[2:]:
        if not line.strip():
            continue
        body, _, note = line.partition("#")
        parts = body.split()
        if len(parts) < 4:
            continue
        atoms.append(
            Atom(
                parts[0].upper(),
                (float(parts[1]), float(parts[2]), float(parts[3])),
                note.strip(),
            )
        )
    return atoms


def parse_log(path: Path, atom_count: int) -> LogSummary:
    lines = path.read_text(errors="replace").splitlines()
    text = "\n".join(lines)
    energies = [
        float(match.group(1))
        for match in re.finditer(r"TOTAL ENERGY\s*=\s*(-?\d+\.\d+)", text)
    ]
    searches = [
        (
            int(match.group(1)),
            float(match.group(2)),
            float(match.group(3)),
            float(match.group(4)),
        )
        for match in re.finditer(
            r"NSERCH:\s*(\d+)\s+E=\s*(-?\d+\.\d+)\s+GRAD\. MAX=\s*([0-9.Ee+-]+)\s+R\.M\.S\.=\s*([0-9.Ee+-]+)",
            text,
        )
    ]
    s2_vals = [float(match.group(1)) for match in re.finditer(r"S-SQUARED\s*=\s*([0-9.]+)", text)]
    wall_vals = [
        float(match.group(1))
        for match in re.finditer(r"TOTAL WALL CLOCK TIME=\s*([0-9.]+)\s+SECONDS", text)
    ]
    nserch = grad_max = grad_rms = None
    if searches:
        nserch, energy_from_search, grad_max, grad_rms = searches[-1]
    else:
        energy_from_search = None
    label_match = re.search(r"_M(\d+)_", path.name)
    label = f"M{label_match.group(1)}" if label_match else path.stem
    evaluated_coords = []
    if nserch is not None:
        _, coords, _ = parse_coords_for_nserch(path, atom_count, nserch)
        evaluated_coords = [Atom(atom.element, atom.xyz, "") for atom in coords]
    return LogSummary(
        label=label,
        path=path,
        normal="EXECUTION OF GAMESS TERMINATED NORMALLY" in text,
        not_converged="FAILURE TO LOCATE STATIONARY POINT" in text,
        energy=energies[-1] if energies else energy_from_search,
        nserch=nserch,
        grad_max=grad_max,
        grad_rms=grad_rms,
        s2=s2_vals[-1] if s2_vals else None,
        wall_seconds=wall_vals[-1] if wall_vals else None,
        coords=evaluated_coords,
    )


def dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def fmt(value: float | None, digits: int = 4) -> str:
    if value is None:
        return "NA"
    return f"{value:.{digits}f}"


def find_fe_and_donors(atoms: list[Atom]) -> tuple[int, list[int]]:
    fe_idx = next((i for i, atom in enumerate(atoms) if atom.element == "FE"), -1)
    preferred_notes = {
        "H1:187:NE2",
        "H2:270:NE2",
        "EAC:349:OE1",
        "ACT:501:O2",
        "UNL:1:N1",
        "HOH:875:O",
    }
    donors = [i for i, atom in enumerate(atoms) if atom.note.split()[0] in preferred_notes]
    if not donors and fe_idx >= 0:
        donors = [
            i
            for i, atom in enumerate(atoms)
            if i != fe_idx and atom.element in {"N", "O", "S"} and dist(atoms[fe_idx].xyz, atom.xyz) <= 3.2
        ]
    return fe_idx, donors


def expected_s2(mult_label: str) -> float | None:
    match = re.fullmatch(r"M(\d+)", mult_label)
    if not match:
        return None
    mult = int(match.group(1))
    spin = (mult - 1) / 2
    return spin * (spin + 1)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", required=True, help="Directory containing logs/ and feiii_radical_truncated_model.xyz")
    args = parser.parse_args()

    workdir = Path(args.workdir)
    ref_atoms = parse_xyz(workdir / "feiii_radical_truncated_model.xyz")
    logs = sorted((workdir / "logs").glob("FEIII_radical_M*_fastopt.log"))
    summaries = [parse_log(log, len(ref_atoms)) for log in logs]
    energies = [s.energy for s in summaries if s.energy is not None]
    min_energy = min(energies) if energies else None

    print("# GAMESS Fe-Radical Spin Screen Comparison\n")
    print("| model | normal | optimizer | energy Hartree | relative kcal/mol | wall h | S^2 | expected S^2 | NSERCH | grad max | grad RMS |")
    print("|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for summary in summaries:
        rel = None
        if summary.energy is not None and min_energy is not None:
            rel = (summary.energy - min_energy) * HARTREE_TO_KCAL
        wall_h = summary.wall_seconds / 3600 if summary.wall_seconds is not None else None
        optimizer = "not converged" if summary.not_converged else "converged"
        print(
            f"| {summary.label} | {str(summary.normal).lower()} | {optimizer} | "
            f"{fmt(summary.energy, 10)} | {fmt(rel, 3)} | {fmt(wall_h, 2)} | "
            f"{fmt(summary.s2, 3)} | {fmt(expected_s2(summary.label), 3)} | "
            f"{summary.nserch if summary.nserch is not None else 'NA'} | "
            f"{fmt(summary.grad_max, 6)} | {fmt(summary.grad_rms, 6)} |"
        )

    fe_idx, donor_indices = find_fe_and_donors(ref_atoms)
    if fe_idx < 0 or not donor_indices:
        return 0

    print("\n## Fe-Donor Distances (Angstrom)\n")
    labels = [ref_atoms[i].note.split()[0] if ref_atoms[i].note else f"{ref_atoms[i].element}{i + 1}" for i in donor_indices]
    print("| donor | initial | " + " | ".join(summary.label for summary in summaries) + " |")
    print("|---|---:|" + "---:|" * len(summaries))
    for donor_label, donor_idx in zip(labels, donor_indices):
        row = [donor_label, fmt(dist(ref_atoms[fe_idx].xyz, ref_atoms[donor_idx].xyz), 3)]
        for summary in summaries:
            if summary.coords and len(summary.coords) == len(ref_atoms):
                row.append(fmt(dist(summary.coords[fe_idx].xyz, summary.coords[donor_idx].xyz), 3))
            else:
                row.append("NA")
        print("| " + " | ".join(row) + " |")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
