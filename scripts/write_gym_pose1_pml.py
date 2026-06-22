#!/usr/bin/env python3
"""Write a PyMOL scene for GYM 4X8B pose1 with palegreen surface."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor", type=Path, required=True)
    parser.add_argument("--ligand", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    text = f"""# Auto-generated GYM pose1 view.
reinitialize
load {args.receptor.as_posix()}, receptor
load {args.ligand.as_posix()}, ligand

bg_color white
set ray_opaque_background, off
set orthoscopic, on
set valence, 0
set stick_radius, 0.16
set transparency, 0.62
set surface_quality, 1
hide everything

show cartoon, receptor and polymer.protein
color gray70, receptor and polymer.protein
show surface, receptor and polymer.protein within 8 of ligand
color palegreen, receptor and polymer.protein within 8 of ligand

show sticks, ligand
util.cbag ligand

show spheres, receptor and (resn FE+FE2+FE3+MG+CA+CL or elem Fe+Mg+Ca+Cl)
util.cbag receptor and (resn FE+FE2+FE3+MG+CA+CL or elem Fe+Mg+Ca+Cl)
set sphere_scale, 0.45, receptor and (resn FE+FE2+FE3+MG+CA+CL or elem Fe+Mg+Ca+Cl)

select fe_core, receptor and resn FE+FE2+FE3
select metal_shell, byres (receptor and polymer.protein within 4 of fe_core)
show sticks, metal_shell and not name N+C+O+CA
color green, metal_shell and elem C
color blue, metal_shell and elem N
color red, metal_shell and elem O
color yellow, metal_shell and elem S

set dash_color, yellow
set dash_width, 2.2
set dash_radius, 0.045
distance Fe_to_HID51, fe_core, receptor and chain A and resi 51 and name NE2
distance Fe_to_HID134, fe_core, receptor and chain A and resi 134 and name NE2
distance Fe_to_HID138, fe_core, receptor and chain A and resi 138 and name NE2
hide labels, Fe_to_*

zoom ligand or fe_core or metal_shell, 8
deselect
"""
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text, encoding="utf-8")
    print(f"pml={args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
