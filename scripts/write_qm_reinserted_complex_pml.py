from __future__ import annotations

import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: write_qm_reinserted_complex_pml.py complex_pdb out_pml")
        return 2
    complex_pdb = Path(sys.argv[1])
    out_pml = Path(sys.argv[2])
    out_pml.parent.mkdir(parents=True, exist_ok=True)
    out_pml.write_text(
        "\n".join(
            [
                "delete all",
                f"load {complex_pdb}, qm_ligand_complex",
                "hide everything",
                "show cartoon, polymer.protein",
                "color gray80, polymer.protein",
                "show sticks, resn ACT+UNL",
                "color palegreen, resn ACT+UNL and elem C",
                "color blue, resn ACT+UNL and elem N",
                "color red, resn ACT+UNL and elem O",
                "color white, resn ACT+UNL and elem H",
                "show spheres, resn FE+FE2+FE3 or elem Fe",
                "set sphere_scale, 0.45, resn FE+FE2+FE3 or elem Fe",
                "show_metal",
                "zoom resn ACT+UNL+FE+FE2+FE3 or elem Fe, 8",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(f"pml={out_pml}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
