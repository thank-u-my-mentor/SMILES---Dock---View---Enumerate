#!/usr/bin/env python3
"""Fill a GAMESS Hessian input's $DATA section from an optimized GAMESS log."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


COORD_RE = re.compile(
    r"^\s*[A-Z][A-Za-z]?\s+\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*$"
)


def extract_final_data(log: Path) -> list[str]:
    lines = log.read_text(errors="replace").splitlines()
    marker = "COORDINATES OF ALL ATOMS ARE (ANGS)"
    indices = [idx for idx, line in enumerate(lines) if marker in line]
    if not indices:
        raise RuntimeError(f"no coordinate blocks found in {log}")
    start = indices[-1]
    coords: list[str] = []
    for line in lines[start + 1 :]:
        if not line.strip():
            if coords:
                break
            continue
        if "ATOM" in line and "ATOMIC" in line:
            continue
        if "CHARGE" in line or "COORDINATES" in line:
            continue
        if COORD_RE.match(line):
            parts = line.split()
            coords.append(f"{parts[0]:<2s} {float(parts[1]):10.1f} {float(parts[2]):12.6f} {float(parts[3]):12.6f} {float(parts[4]):12.6f}")
        elif coords:
            break
    if not coords:
        raise RuntimeError(f"failed to parse final coordinates from {log}")
    return coords


def replace_data_section(template: Path, out: Path, coords: list[str]) -> None:
    lines = template.read_text(errors="replace").splitlines()
    data_idx = next((idx for idx, line in enumerate(lines) if line.strip().upper() == "$DATA"), None)
    if data_idx is None:
        raise RuntimeError(f"no $DATA section in {template}")
    end_idx = next((idx for idx in range(data_idx + 1, len(lines)) if lines[idx].strip().upper() == "$END"), None)
    if end_idx is None:
        raise RuntimeError(f"no $END after $DATA in {template}")
    header = lines[data_idx : min(data_idx + 3, len(lines))]
    if len(header) < 3:
        raise RuntimeError(f"incomplete $DATA header in {template}")
    new_lines = lines[: data_idx + 3] + coords + [" $END"] + lines[end_idx + 1 :]
    out.write_text("\n".join(new_lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--opt-log", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    coords = extract_final_data(args.opt_log)
    if args.out.exists():
        backup = args.out.with_suffix(args.out.suffix + ".bak_empty_data")
        if not backup.exists():
            backup.write_text(args.out.read_text(errors="replace"), encoding="utf-8")
    replace_data_section(args.template, args.out, coords)
    print(f"atoms={len(coords)}")
    print(f"out={args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
