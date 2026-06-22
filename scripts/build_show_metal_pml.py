#!/usr/bin/env python3
"""Build a PyMOL metal-site visualization script for a GROMACS structure."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path)
    parser.add_argument("--load-trajectory", action="store_true", help="Load the full trajectory into PyMOL; off by default to avoid memory blowups")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--metal-selection", default="resn FE+FE2+FE3 or elem Fe")
    parser.add_argument("--radius", type=float, default=4.0)
    parser.add_argument("--carbon-color", default="green")
    args = parser.parse_args()

    traj_line = f"load_traj {args.trajectory.as_posix()}, complex\n" if args.trajectory and args.load_trajectory else ""
    text = f"""# Auto-generated metal-site view.
# Load this with: pymol {args.out.as_posix()}
reinitialize
load {args.structure.as_posix()}, complex
{traj_line}

bg_color white
set ray_opaque_background, off
set orthoscopic, on
set valence, 0
set stick_radius, 0.16
set dash_color, yellow
set dash_width, 2.2
set dash_radius, 0.045
set auto_zoom, off
set defer_builds_mode, 3
hide everything
show cartoon, polymer.protein
color gray70, polymer.protein
show sticks, resn LIG
util.cbag resn LIG
show spheres, ({args.metal_selection})
set sphere_scale, 0.45, ({args.metal_selection})

python
from pymol import cmd

def ple_show_metal(selection="{args.metal_selection}", radius={args.radius:.3f}):
    cmd.delete("metal_dist_*")
    cmd.select("metal_core", selection)
    if cmd.count_atoms("metal_core") == 0:
        print(f"[ple_show_metal] no metal atoms matched: {selection}")
        return
    cmd.select("metal_scope", f"byres ((not solvent and not resn LIG and not ({{selection}})) within {{radius}} of ({{selection}}))")
    cmd.select("metal_shell", "metal_scope and not resn SOL+HOH+WAT+NA+CL+K+MG+CA")
    cmd.select("metal_sidechains", "metal_shell and not name N+C+O+CA+H+HA")
    print(f"[ple_show_metal] metals={{cmd.count_atoms('metal_core')}} shell_atoms={{cmd.count_atoms('metal_shell')}} sidechain_atoms={{cmd.count_atoms('metal_sidechains')}}")
    cmd.show("sticks", "metal_sidechains")
    cmd.color("{args.carbon_color}", "metal_sidechains and elem C")
    # Keep PyMOL/default hetero-atom convention for N/O/S and metal colors.
    cmd.show("spheres", "metal_core")
    cmd.set("sphere_scale", 0.45, "metal_core")
    donors = "metal_shell and (elem N+O+S) and not name N+O"
    cmd.select("metal_donors", donors)
    donor_count = cmd.count_atoms("metal_donors")
    print(f"[ple_show_metal] donors={{donor_count}}")
    for atom in cmd.get_model("metal_donors", state=1).atom:
        donor_sel = f"index {{atom.index}}"
        dist_name = f"metal_dist_{{atom.resn}}{{atom.resi}}_{{atom.name}}"
        cmd.distance(dist_name, "metal_core", donor_sel, cutoff=radius, mode=2, state=1)
        cmd.hide("labels", dist_name)
        cmd.color("yellow", dist_name)
    cmd.zoom("metal_shell or metal_core or resn LIG", 8)
    cmd.deselect()

cmd.extend("ple_show_metal", ple_show_metal)
cmd.extend("show_metal", ple_show_metal)
python end

ple_show_metal
"""
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text, encoding="utf-8")
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
