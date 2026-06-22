#!/usr/bin/env python
"""Rewrite selected GROMACS mdp key-value lines."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--set", action="append", default=[], help="key=value")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    replacements = {}
    for item in args.set:
        if "=" not in item:
            raise SystemExit(f"--set must be key=value, got {item!r}")
        key, value = item.split("=", 1)
        replacements[key.strip()] = value.strip()

    seen = set()
    out = []
    for line in args.input.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith(";") and "=" in line:
            key = line.split("=", 1)[0].strip()
            if key in replacements:
                out.append(f"{key:<24} = {replacements[key]}")
                seen.add(key)
                continue
        out.append(line)
    for key, value in replacements.items():
        if key not in seen:
            out.append(f"{key:<24} = {value}")
    args.output.write_text("\n".join(out) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
