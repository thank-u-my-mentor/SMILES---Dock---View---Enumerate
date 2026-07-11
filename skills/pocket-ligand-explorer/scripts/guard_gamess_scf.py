#!/usr/bin/env python3
"""Guard a GAMESS SCF log and stop electronic-state cliffs early.

This helper is intentionally conservative for Fe-radical MCPB work. It watches
the iteration table in a GAMESS log. If the first few SCF iterations show a
large energy rise with very large density changes, the wavefunction is treated
as a failed electronic-state branch and the matching job processes can be
terminated before hours are wasted.
"""

from __future__ import annotations

import argparse
import os
import re
import signal
import subprocess
import sys
import time
from pathlib import Path


ITER_RE = re.compile(
    r"^\s*(?P<iter>\d+)\s+(?P<ex>\d+)\s+"
    r"(?P<energy>[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+"
    r"(?P<echange>[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+"
    r"(?P<density>[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+"
    r"(?P<metric>[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
    re.MULTILINE,
)


def parse_rows(log: Path) -> tuple[str, list[dict[str, str]]]:
    if not log.exists():
        return "", []
    text = log.read_text(encoding="utf-8", errors="replace")
    return text, [m.groupdict() for m in ITER_RE.finditer(text)]


def classify(
    text: str,
    rows: list[dict[str, str]],
    *,
    min_rows: int,
    energy_rise_threshold: float,
    first_window_threshold: float,
    density_threshold: float,
    metric_threshold: float,
    catastrophic_energy_rise: float,
    catastrophic_density: float,
    catastrophic_metric: float,
) -> tuple[str, str]:
    if "DENSITY CONVERGED" in text:
        return "converged", "DENSITY CONVERGED"
    if "SCF IS UNCONVERGED" in text or "SCF HAS NOT CONVERGED" in text:
        return "unconverged", "GAMESS reported unconverged SCF"
    if "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-" in text:
        return "abnormal", "GAMESS terminated abnormally"
    if "EXECUTION OF GAMESS TERMINATED NORMALLY" in text and not rows:
        return "normal_no_iterations", "normal termination without parsed SCF rows"

    if len(rows) < 5:
        return "watching", f"only {len(rows)} SCF row(s)"

    energies = [float(row["energy"]) for row in rows]
    densities = [float(row["density"]) for row in rows]
    metrics = [float(row["metric"]) for row in rows]
    first = rows[: min(8, len(rows))]
    first_energies = [float(row["energy"]) for row in first]
    first_densities = [float(row["density"]) for row in first]
    first_metrics = [float(row["metric"]) for row in first]

    rise_from_first = energies[-1] - energies[0]
    rise_from_best = energies[-1] - min(energies)
    first_window_rise = max(first_energies) - min(first_energies)
    density_min = min(first_densities)
    metric_max = max(first_metrics)

    if len(rows) < min_rows:
        if (
            (rise_from_first > catastrophic_energy_rise or first_window_rise > catastrophic_energy_rise)
            and max(densities) > catastrophic_density
            and max(metrics) > catastrophic_metric
        ):
            return (
                "electronic_cliff",
                (
                    "catastrophic early SCF divergence before minimum rows: "
                    f"rows={len(rows)} rise_from_first={rise_from_first:.3f} "
                    f"first_window_rise={first_window_rise:.3f} "
                    f"density_max={max(densities):.3f} metric_max={max(metrics):.3f}"
                ),
            )
        return (
            "watching",
            (
                f"rows={len(rows)} below_min_rows={min_rows} "
                f"rise_from_first={rise_from_first:.3f} "
                f"first_window_rise={first_window_rise:.3f} "
                f"last_density={densities[-1]:.6f} last_metric={metrics[-1]:.6f}"
            ),
        )

    if (
        len(rows) >= min_rows
        and (
            rise_from_first > energy_rise_threshold
            or rise_from_best > energy_rise_threshold
            or first_window_rise > first_window_threshold
        )
        and min(densities[-min_rows:]) > density_threshold
        and max(metrics[-min_rows:]) > metric_threshold
    ):
        return (
            "electronic_cliff",
            (
                "early energy/density cliff: "
                f"rows={len(rows)} rise_from_first={rise_from_first:.3f} "
                f"rise_from_best={rise_from_best:.3f} first_window_rise={first_window_rise:.3f} "
                f"first_density_min={density_min:.3f} first_metric_max={metric_max:.3f}"
            ),
        )

    if len(rows) >= min_rows and (max(densities[-min_rows:]) > catastrophic_density or max(metrics[-min_rows:]) > catastrophic_metric):
        return (
            "electronic_collapse",
            (
                "large recent density/metric values: "
                f"recent_density_max={max(densities[-min_rows:]):.3f} "
                f"recent_metric_max={max(metrics[-min_rows:]):.3f}"
            ),
        )

    return (
        "watching",
        (
            f"rows={len(rows)} last_energy={energies[-1]:.6f} "
            f"last_density={densities[-1]:.6f} last_metric={metrics[-1]:.6f}"
        ),
    )


def matching_pids(pattern: str) -> list[int]:
    result = subprocess.run(
        ["pgrep", "-f", pattern],
        text=True,
        capture_output=True,
        check=False,
    )
    pids: list[int] = []
    exclude = {os.getpid(), os.getppid()}
    for line in result.stdout.splitlines():
        try:
            pid = int(line.strip())
        except ValueError:
            continue
        if pid not in exclude:
            pids.append(pid)
    return pids


def kill_pattern(pattern: str) -> list[int]:
    pids = matching_pids(pattern)
    for pid in pids:
        try:
            os.kill(pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    time.sleep(3)
    for pid in matching_pids(pattern):
        if pid in {os.getpid(), os.getppid()}:
            continue
        try:
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    return pids


def write_status(path: Path | None, message: str) -> None:
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(message.rstrip() + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--kill-pattern", required=True)
    parser.add_argument("--status", type=Path, default=None)
    parser.add_argument("--poll", type=float, default=30.0)
    parser.add_argument("--kill", action="store_true")
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--min-rows", type=int, default=25)
    parser.add_argument("--energy-rise-threshold", type=float, default=300.0)
    parser.add_argument("--first-window-threshold", type=float, default=1200.0)
    parser.add_argument("--density-threshold", type=float, default=1.0)
    parser.add_argument("--metric-threshold", type=float, default=0.5)
    parser.add_argument("--catastrophic-energy-rise", type=float, default=2500.0)
    parser.add_argument("--catastrophic-density", type=float, default=50.0)
    parser.add_argument("--catastrophic-metric", type=float, default=5.0)
    args = parser.parse_args()

    while True:
        text, rows = parse_rows(args.log)
        verdict, reason = classify(
            text,
            rows,
            min_rows=args.min_rows,
            energy_rise_threshold=args.energy_rise_threshold,
            first_window_threshold=args.first_window_threshold,
            density_threshold=args.density_threshold,
            metric_threshold=args.metric_threshold,
            catastrophic_energy_rise=args.catastrophic_energy_rise,
            catastrophic_density=args.catastrophic_density,
            catastrophic_metric=args.catastrophic_metric,
        )
        line = f"[guard] verdict={verdict} {reason}"
        print(line, flush=True)
        write_status(args.status, line)

        if verdict in {"converged", "unconverged", "abnormal", "normal_no_iterations"}:
            return 0 if verdict == "converged" else 2

        if verdict in {"electronic_cliff", "electronic_collapse"}:
            if args.kill:
                pids = kill_pattern(args.kill_pattern)
                kill_line = f"[guard] killed_pattern={args.kill_pattern} pids={','.join(map(str, pids))}"
                print(kill_line, flush=True)
                write_status(args.status, kill_line)
            return 3

        if args.once:
            return 0
        time.sleep(args.poll)


if __name__ == "__main__":
    raise SystemExit(main())
