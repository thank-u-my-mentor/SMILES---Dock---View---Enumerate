#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def atom_key(line: str) -> tuple[str, str]:
    atom = line[12:16].strip()
    resn = line[17:20].strip()
    return atom, resn


def main() -> int:
    parser = argparse.ArgumentParser(description="Remove problematic PDB CONECT records for preflight checks.")
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--remove-cross-atom-hetatm", action="store_true")
    parser.add_argument("--remove-fe-conect", action="store_true")
    args = parser.parse_args()

    lines = args.pdb.read_text(encoding="utf-8", errors="replace").splitlines()
    record_by_serial: dict[int, str] = {}
    atom_by_serial: dict[int, tuple[str, str]] = {}
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and line[6:11].strip().isdigit():
            serial = int(line[6:11])
            record_by_serial[serial] = line[:6].strip()
            atom_by_serial[serial] = atom_key(line)

    kept: list[str] = []
    removed: list[str] = []
    for line in lines:
        if not line.startswith("CONECT"):
            kept.append(line)
            continue
        nums: list[int] = []
        for i in range(6, len(line), 5):
            token = line[i : i + 5].strip()
            if token.isdigit():
                nums.append(int(token))
        records = {record_by_serial.get(n, "?") for n in nums}
        has_cross = "ATOM" in records and "HETATM" in records
        has_fe = any(
            atom_by_serial.get(n, ("", ""))[1] == "FE"
            or atom_by_serial.get(n, ("", ""))[0].upper().startswith("FE")
            for n in nums
        )
        if (args.remove_cross_atom_hetatm and has_cross) or (args.remove_fe_conect and has_fe):
            removed.append(line)
            continue
        kept.append(line)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(kept) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote={args.out}")
    print(f"removed={len(removed)}")
    for line in removed:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
