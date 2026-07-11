#!/usr/bin/env python3
"""Build tleap index probes and optionally rewrite simple bond selectors.

This is for MCPB.py handoffs where PDB residue IDs such as 187/431 are not the
same as the residue/unit indices that tleap uses in selectors like mol.673.NE2.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


STOP_PREFIXES = (
    "bond ",
    "savepdb ",
    "saveamberparm ",
    "solvatebox ",
    "addions ",
)


def read_leap_prefix(path: Path) -> list[str]:
    lines: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip().lower()
        if stripped == "quit" or any(stripped.startswith(prefix) for prefix in STOP_PREFIXES):
            break
        lines.append(line)
    return lines


def residue_order(pdb: Path) -> list[dict[str, str | int]]:
    order: list[dict[str, str | int]] = []
    seen: set[tuple[str, str, str, str]] = set()
    for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        res = line[17:20].strip()
        chain = line[21:22].strip()
        resseq = line[22:26].strip()
        icode = line[26:27].strip()
        key = (chain, resseq, icode, res)
        if key in seen:
            continue
        seen.add(key)
        order.append(
            {
                "index": len(order) + 1,
                "chain": chain,
                "resseq": resseq,
                "icode": icode,
                "res": res,
            }
        )
    return order


def atom_names_by_residue(pdb: Path) -> dict[tuple[str, str, str, str], set[str]]:
    atoms: dict[tuple[str, str, str, str], set[str]] = {}
    for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        key = (
            line[21:22].strip(),
            line[22:26].strip(),
            line[26:27].strip(),
            line[17:20].strip(),
        )
        atoms.setdefault(key, set()).add(line[12:16].strip())
    return atoms


def write_probe(args: argparse.Namespace) -> int:
    lines = read_leap_prefix(args.base_leap)
    probes = ["list", f"desc {args.unit}"]
    for target in args.targets:
        probes.append(f"desc {args.unit}.{target}")
    for target in args.targets:
        for atom in args.atoms:
            probes.append(f"desc {args.unit}.{target}.{atom}")
    probes.extend([f"check {args.unit}", "quit"])
    args.out.write_text("\n".join(lines + probes) + "\n", encoding="utf-8", newline="\n")
    print(f"out={args.out}")
    return 0


def write_scan(args: argparse.Namespace) -> int:
    lines = read_leap_prefix(args.base_leap)
    for i in range(args.start, args.end + 1):
        for atom in args.atoms:
            lines.append(f"desc {args.unit}.{i}.{atom}")
    lines.append("quit")
    args.out.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"out={args.out}")
    return 0


def rewrite_bonds(args: argparse.Namespace) -> int:
    order = residue_order(args.pdb)
    atoms = atom_names_by_residue(args.pdb)
    by_resseq: dict[str, list[dict[str, str | int]]] = {}
    for row in order:
        if args.chain and row["chain"] != args.chain:
            continue
        by_resseq.setdefault(str(row["resseq"]), []).append(row)

    bond_re = re.compile(
        rf"^(bond\s+{re.escape(args.unit)}\.)(\d+)(\.[A-Za-z0-9']+\s+"
        rf"{re.escape(args.unit)}\.)(\d+)(\.[A-Za-z0-9']+(?:\s+.*)?)$"
    )
    replacements: list[str] = []
    out_lines: list[str] = []
    for line in args.leap_in.read_text(encoding="utf-8", errors="replace").splitlines():
        match = bond_re.match(line.strip())
        if not match:
            out_lines.append(line)
            continue
        prefix, r1, mid, r2, suffix = match.groups()
        rows1 = by_resseq.get(r1, [])
        rows2 = by_resseq.get(r2, [])
        if len(rows1) != 1 or len(rows2) != 1:
            out_lines.append(line)
            replacements.append(f"UNRESOLVED {line}")
            continue
        n1 = str(rows1[0]["index"])
        n2 = str(rows2[0]["index"])
        new_line = f"{prefix}{n1}{mid}{n2}{suffix}"
        out_lines.append(new_line)
        replacements.append(
            f"{line}  ==>  {new_line}    "
            f"[{r1}:{rows1[0]['res']}->{n1}, {r2}:{rows2[0]['res']}->{n2}]"
        )

    args.out.write_text("\n".join(out_lines) + "\n", encoding="utf-8", newline="\n")
    map_out = args.map_out or args.out.with_suffix(args.out.suffix + ".map.txt")
    with map_out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("index\tchain\tresseq\ticode\tres\tatoms\n")
        for row in order:
            key = (str(row["chain"]), str(row["resseq"]), str(row["icode"]), str(row["res"]))
            handle.write(
                f"{row['index']}\t{row['chain']}\t{row['resseq']}\t"
                f"{row['icode']}\t{row['res']}\t"
                f"{','.join(sorted(atoms.get(key, [])))}\n"
            )
        handle.write("\nREPLACEMENTS\n")
        for item in replacements:
            handle.write(item + "\n")

    print(f"out={args.out}")
    print(f"map={map_out}")
    for item in replacements:
        print(item)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    probe = sub.add_parser("probe", help="write a compact tleap desc/check probe")
    probe.add_argument("--base-leap", type=Path, required=True)
    probe.add_argument("--out", type=Path, required=True)
    probe.add_argument("--unit", default="mol")
    probe.add_argument(
        "--targets",
        nargs="+",
        default=["1", "187", "270", "349", "431", "501", "875"],
    )
    probe.add_argument("--atoms", nargs="+", default=["FE", "N1", "O", "O2", "NE2", "OE1"])
    probe.set_defaults(func=write_probe)

    scan = sub.add_parser("scan", help="write a broad tleap desc scan")
    scan.add_argument("--base-leap", type=Path, required=True)
    scan.add_argument("--out", type=Path, required=True)
    scan.add_argument("--unit", default="mol")
    scan.add_argument("--start", type=int, default=1)
    scan.add_argument("--end", type=int, default=900)
    scan.add_argument("--atoms", nargs="+", default=["FE", "N1", "O", "O2", "NE2", "OE1"])
    scan.set_defaults(func=write_scan)

    rewrite = sub.add_parser(
        "rewrite-bonds",
        help="rewrite simple bond mol.<PDB-resseq> selectors to loadpdb order indices",
    )
    rewrite.add_argument("--pdb", type=Path, required=True)
    rewrite.add_argument("--leap-in", type=Path, required=True)
    rewrite.add_argument("--out", type=Path, required=True)
    rewrite.add_argument("--unit", default="mol")
    rewrite.add_argument("--chain", default="A")
    rewrite.add_argument("--map-out", type=Path, default=None)
    rewrite.set_defaults(func=rewrite_bonds)

    args = parser.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
