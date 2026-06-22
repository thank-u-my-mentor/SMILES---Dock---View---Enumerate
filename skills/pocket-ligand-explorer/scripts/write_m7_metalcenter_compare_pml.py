#!/usr/bin/env python3
"""Write a compact PyMOL scene comparing initial and M7-optimized metal centers."""

from __future__ import annotations

import argparse
from pathlib import Path


def pml_path(path: Path) -> str:
    return path.as_posix()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a PyMOL PML that overlays initial and M7-optimized Fe radical metal-center geometry."
    )
    parser.add_argument("--initial", required=True, type=Path, help="Initial full complex PDB used for M7 preflight.")
    parser.add_argument("--m7", required=True, type=Path, help="Full complex PDB patched with M7 optimized coordinates.")
    parser.add_argument("--out", required=True, type=Path, help="Output PML path.")
    args = parser.parse_args()

    initial = args.initial.resolve()
    m7 = args.m7.resolve()
    out = args.out.resolve()

    for path in (initial, m7):
        if not path.exists():
            raise FileNotFoundError(path)

    out.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "reinitialize",
        "set retain_order, 1",
        "set connect_mode, 1",
        "set auto_zoom, off",
        "set bg_rgb, [1, 1, 1]",
        "set ray_opaque_background, off",
        f"load {pml_path(initial)}, initial",
        f"load {pml_path(m7)}, m7",
        "hide everything",
        "",
        "# Only show the chemically relevant local region. Full protein remains loaded but transparent.",
        "select active_initial, initial and ((chain A and resi 187+270+349+431+501) or (chain A and resn UNL) or (chain B and resi 875))",
        "select active_m7, m7 and ((chain A and resi 187+270+349+431+501) or (chain A and resn UNL) or (chain B and resi 875))",
        "show cartoon, initial and polymer",
        "color gray85, initial and polymer",
        "set cartoon_transparency, 0.82, initial and polymer",
        "show sticks, active_initial or active_m7",
        "set stick_radius, 0.12, active_initial or active_m7",
        "color gray65, active_initial and elem C",
        "color cyan, active_m7 and elem C",
        "color red, elem O",
        "color blue, elem N",
        "color white, elem H",
        "color orange, elem Fe",
        "show spheres, (active_initial or active_m7) and elem Fe",
        "set sphere_scale, 0.42, (active_initial or active_m7) and elem Fe",
        "",
        "# Fe coordination distances: yellow = initial, magenta = M7 optimized.",
    ]

    donors = [
        ("H187_NE2", "chain A and resi 187 and name NE2"),
        ("H270_NE2", "chain A and resi 270 and name NE2"),
        ("E349_OE1", "chain A and resi 349 and name OE1"),
        ("ACT501_O2", "chain A and resi 501 and name O2"),
        ("UNL_N1", "chain A and resn UNL and name N1"),
        ("HOH875_O", "chain B and resi 875 and name O"),
    ]
    for label, donor_sel in donors:
        lines.extend(
            [
                f"distance init_{label}, initial and chain A and resi 431 and name FE, initial and {donor_sel}",
                f"distance m7_{label}, m7 and chain A and resi 431 and name FE, m7 and {donor_sel}",
                f"set dash_color, yellow, init_{label}",
                f"set dash_color, magenta, m7_{label}",
                f"set dash_width, 2.0, init_{label}",
                f"set dash_width, 2.0, m7_{label}",
            ]
        )

    lines.extend(
        [
            "",
            "# Same-atom displacement guides for the largest local changes.",
            "distance move_GLU349_OE2, initial and chain A and resi 349 and name OE2, m7 and chain A and resi 349 and name OE2",
            "distance move_UNL_O1, initial and chain A and resn UNL and name O1, m7 and chain A and resn UNL and name O1",
            "distance move_UNL_N1, initial and chain A and resn UNL and name N1, m7 and chain A and resn UNL and name N1",
            "distance move_ACT_O1, initial and chain A and resi 501 and name O1, m7 and chain A and resi 501 and name O1",
            "set dash_color, hotpink, move_GLU349_OE2",
            "set dash_color, hotpink, move_UNL_O1",
            "set dash_color, hotpink, move_UNL_N1",
            "set dash_color, hotpink, move_ACT_O1",
            "set dash_width, 2.5, move_GLU349_OE2 move_UNL_O1 move_UNL_N1 move_ACT_O1",
            "",
            "# Labels are intentionally small; use hide labels if they become visually busy.",
            "set label_size, 14",
            "set label_color, black",
            "set dash_gap, 0.25",
            "set dash_radius, 0.05",
            "zoom active_initial or active_m7, 5",
            "orient active_initial or active_m7",
            "group Fe_distance_initial, init_H187_NE2 init_H270_NE2 init_E349_OE1 init_ACT501_O2 init_UNL_N1 init_HOH875_O",
            "group Fe_distance_M7, m7_H187_NE2 m7_H270_NE2 m7_E349_OE1 m7_ACT501_O2 m7_UNL_N1 m7_HOH875_O",
            "group atom_displacements, move_GLU349_OE2 move_UNL_O1 move_UNL_N1 move_ACT_O1",
            "",
            "print 'Initial carbon = gray; M7 carbon = cyan; Fe distances: yellow initial, magenta M7; hotpink = same-atom displacement guides.'",
        ]
    )

    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote={out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
