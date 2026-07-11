#!/usr/bin/env python3
"""Wait for a GAMESS ENERGY job and optionally run a MOREAD ENERGY retry."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path


def has_text(path: Path, needle: str) -> bool:
    if not path.exists():
        return False
    return needle in path.read_text(encoding="utf-8", errors="replace")


def has_any(path: Path, needles: list[str]) -> bool:
    if not path.exists():
        return False
    text = path.read_text(encoding="utf-8", errors="replace")
    return any(needle in text for needle in needles)


def wait_for_first(log: Path, status: Path, poll_seconds: int, timeout_hours: float) -> str:
    deadline = time.time() + timeout_hours * 3600
    while time.time() < deadline:
        if has_any(status, ["exit_code=", "status=normal", "status=not_normal"]):
            return "status_finished"
        if has_any(log, ["EXECUTION OF GAMESS TERMINATED NORMALLY", "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-"]):
            return "log_finished"
        time.sleep(poll_seconds)
    return "timeout"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run a second fixed-geometry GAMESS MOREAD ENERGY job if the first one fails to converge."
    )
    parser.add_argument("--mcpb-dir", type=Path, required=True)
    parser.add_argument("--first-job", required=True)
    parser.add_argument("--source-job", required=True)
    parser.add_argument("--moread-job", required=True)
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--mult", type=int, required=True)
    parser.add_argument("--norb", type=int, required=True)
    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--poll-seconds", type=int, default=120)
    parser.add_argument("--timeout-hours", type=float, default=12.0)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument(
        "--helper",
        default="/mnt/e/Codex/skills/pocket-ligand-explorer/scripts/prepare_gamess_moread_energy.py",
    )
    args = parser.parse_args()

    root = args.mcpb_dir
    logs = root / "gamess_logs"
    logs.mkdir(exist_ok=True)
    chain_status = logs / f"{args.first_job}_to_{args.moread_job}.chain_status"
    first_log = logs / f"{args.first_job}.log"
    first_status = logs / f"{args.first_job}.status"
    dat = Path(f"/home/qin/softwares/gamess/restart/{args.first_job}.dat")

    def write(line: str) -> None:
        with chain_status.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {line}\n")

    write(f"watch_start first={args.first_job} moread={args.moread_job}")
    finished_by = wait_for_first(first_log, first_status, args.poll_seconds, args.timeout_hours)
    write(f"first_finished_by={finished_by}")
    if finished_by == "timeout":
        write("action=none reason=timeout")
        return 2

    if has_text(first_log, "DENSITY CONVERGED") and not has_any(
        first_log, ["SCF IS UNCONVERGED", "SCF HAS NOT CONVERGED"]
    ):
        write("action=none reason=first_density_converged")
        return 0

    if not dat.exists() or "$VEC" not in dat.read_text(encoding="utf-8", errors="replace"):
        write(f"action=none reason=no_vec dat={dat}")
        return 3

    write(f"action=prepare_moread dat={dat}")
    prep = [
        args.python,
        args.helper,
        "--mcpb-dir",
        str(root),
        "--source-job",
        args.source_job,
        "--out-job",
        args.moread_job,
        "--dat",
        str(dat),
        "--charge",
        str(args.charge),
        "--mult",
        str(args.mult),
        "--norb",
        str(args.norb),
        "--cores",
        str(args.cores),
        "--session",
        args.moread_job.lower(),
    ]
    subprocess.run(prep, cwd=root, check=True)

    runner = root / f"run_{args.moread_job}_foreground.sh"
    if not runner.exists():
        write(f"action=none reason=missing_runner runner={runner}")
        return 4

    env = os.environ.copy()
    env["NCORES"] = str(args.cores)
    write(f"action=run_moread runner={runner}")
    rc = subprocess.call(["bash", str(runner)], cwd=root, env=env)
    write(f"moread_exit_code={rc}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
