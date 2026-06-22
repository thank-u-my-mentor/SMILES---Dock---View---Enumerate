from __future__ import annotations

import sys
from pathlib import Path


ACT_NAME_BY_SERIAL_ORDER = {
    1: "C1",
    2: "C2",
    3: "O1",
    4: "O2",
    5: "H1",
    6: "H2",
    7: "H3",
}
ACT_CHAIN = "A"
ACT_RESI = 501


def rewrite_atom_line(line: str, atom_name: str, element: str) -> str:
    # Preserve coordinates/occupancy/tempfactor exactly; only change atom/residue names.
    padded = (line + " " * 80)[:80]
    return (
        f"{padded[:12]}{atom_name:>4s}{padded[16:17]}ACT"
        f"{padded[20:21]}{ACT_CHAIN}{ACT_RESI:4d}{padded[26:76]}{element:>2s}  "
    )


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: rename_acetate_ace_to_act.py input.pdb output.pdb")
        return 2

    inp = Path(sys.argv[1])
    out = Path(sys.argv[2])
    lines = inp.read_text(errors="replace").splitlines()

    ace_serials: list[int] = []
    rewritten: list[str] = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "ACE":
            serial = int(line[6:11])
            ace_serials.append(serial)
            order = len(ace_serials)
            atom_name = ACT_NAME_BY_SERIAL_ORDER.get(order, line[12:16].strip() or f"X{order}")
            element = "H" if atom_name.startswith("H") else atom_name[0]
            rewritten.append(rewrite_atom_line(line, atom_name, element))
        else:
            rewritten.append(line)

    if not ace_serials:
        raise SystemExit("No ACE residue found to rename")

    out.write_text("\n".join(rewritten) + "\n", encoding="utf-8")
    print(f"wrote={out}")
    print(f"renamed_ace_atoms={len(ace_serials)}")
    print("new_resname=ACT")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
