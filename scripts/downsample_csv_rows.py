#!/usr/bin/env python
"""Write every Nth data row from a CSV while preserving the header."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--stride", required=True, type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.stride < 1:
        raise SystemExit("--stride must be >= 1")
    args.output.parent.mkdir(parents=True, exist_ok=True)

    input_rows = 0
    output_rows = 0
    with args.input.open(newline="", encoding="utf-8") as src, args.output.open(
        "w", newline="", encoding="utf-8"
    ) as dst:
        reader = csv.reader(src)
        writer = csv.writer(dst)
        header = next(reader)
        writer.writerow(header)
        output_rows += 1
        input_rows += 1
        for index, row in enumerate(reader):
            input_rows += 1
            if index % args.stride == 0:
                writer.writerow(row)
                output_rows += 1

    print(f"input_rows={input_rows}")
    print(f"output_rows={output_rows}")
    print(f"wrote={args.output}")


if __name__ == "__main__":
    main()
