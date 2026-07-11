#!/usr/bin/env python3
"""Summarize GAMESS SCF logs without fragile shell grep pipelines."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


ITER_RE = re.compile(
    r"^\s*(?P<iter>\d+)\s+(?P<ex>\d+)\s+"
    r"(?P<energy>[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+"
    r"(?P<echange>[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+"
    r"(?P<density>[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+"
    r"(?P<metric>[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
    re.MULTILINE,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--last", type=int, default=20)
    args = parser.parse_args()

    text = args.log.read_text(encoding="utf-8", errors="replace")
    rows = [m.groupdict() for m in ITER_RE.finditer(text)]
    markers = []
    for marker in [
        "DENSITY CONVERGED",
        "SCF IS UNCONVERGED",
        "SCF HAS NOT CONVERGED",
        "FINAL U-B3LYP ENERGY",
        "S-SQUARED",
        "EXECUTION OF GAMESS TERMINATED NORMALLY",
        "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-",
    ]:
        if marker in text:
            markers.append(marker)

    print(f"log={args.log}")
    print(f"iterations={len(rows)}")
    print("markers=" + ("; ".join(markers) if markers else "none"))
    if rows:
        print("last_iterations:")
        print("iter energy e_change density metric")
        for row in rows[-args.last :]:
            print(
                f"{row['iter']} {row['energy']} {row['echange']} "
                f"{row['density']} {row['metric']}"
            )
        last = rows[-1]
        density = float(last["density"])
        metric = float(last["metric"])
        energies = [float(row["energy"]) for row in rows]
        densities = [float(row["density"]) for row in rows]
        metrics = [float(row["metric"]) for row in rows]
        early = rows[: min(25, len(rows))]
        early_energies = [float(row["energy"]) for row in early]
        early_density_min = min(float(row["density"]) for row in early)
        early_metric_min = min(float(row["metric"]) for row in early)
        energy_rise_from_first = energies[-1] - energies[0]
        energy_rise_from_best = energies[-1] - min(energies)
        monotone_bad_window = 0
        for prev, cur in zip(energies, energies[1:]):
            if cur > prev:
                monotone_bad_window += 1
            else:
                monotone_bad_window = 0
        early_cliff = (
            len(rows) >= 12
            and energy_rise_from_first > 50.0
            and early_density_min > 0.02
            and early_metric_min > 0.5
        )
        electronic_collapse = (
            len(rows) >= 8
            and (max(densities[-min(10, len(densities)) :]) > 10.0 or max(metrics[-min(10, len(metrics)) :]) > 2.0)
        )
        if "DENSITY CONVERGED" in markers:
            verdict = "converged"
        elif "SCF IS UNCONVERGED" in markers or "SCF HAS NOT CONVERGED" in markers:
            verdict = "unconverged"
        elif early_cliff:
            verdict = "poor_trend_early_energy_cliff"
        elif electronic_collapse:
            verdict = "poor_trend_electronic_collapse"
        elif len(rows) >= 20 and (density > 0.02 or metric > 0.2):
            verdict = "poor_trend"
        else:
            verdict = "running_or_promising"
        print(f"energy_rise_from_first={energy_rise_from_first:.6f}")
        print(f"energy_rise_from_best={energy_rise_from_best:.6f}")
        print(f"early_density_min={early_density_min:.6g}")
        print(f"early_metric_min={early_metric_min:.6g}")
        print(f"consecutive_energy_rises_to_end={monotone_bad_window}")
        print(f"verdict={verdict}")
    else:
        print("last_iterations=none")
        print("verdict=no_iterations_yet")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
