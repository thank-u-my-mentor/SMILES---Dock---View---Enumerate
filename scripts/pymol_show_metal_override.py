"""Persistent PyMOL show_metal override used by PLE.

Load from ~/.pymolrc with:
    run E:/Codex/scripts/pymol_show_metal_override.py

It intentionally replaces any existing show_metal command.
"""

from pymol import cmd


def show_metal(selection="resn FE+FE2+FE3 or elem Fe", radius=4.0, carbon_color="palegreen"):
    """Show metal, nearby side chains, and donor distances without loading trajectories."""
    cmd.delete("metal_dist_*")
    cmd.select("metal_core", selection)
    metal_count = cmd.count_atoms("metal_core")
    if metal_count == 0:
        print(f"[show_metal:PLE] no metal atoms matched: {selection}")
        return

    cmd.select(
        "metal_scope",
        f"byres ((not solvent and not resn LIG and not ({selection})) within {float(radius)} of ({selection}))",
    )
    cmd.select("metal_shell", "metal_scope and not resn SOL+HOH+WAT+NA+CL+K+MG+CA")
    cmd.select("metal_sidechains", "metal_shell and not name N+C+O+CA+H+HA")
    print(
        "[show_metal:PLE] "
        f"metals={metal_count} "
        f"shell_atoms={cmd.count_atoms('metal_shell')} "
        f"sidechain_atoms={cmd.count_atoms('metal_sidechains')}"
    )

    cmd.show("spheres", "metal_core")
    cmd.set("sphere_scale", 0.45, "metal_core")
    cmd.show("sticks", "metal_sidechains")
    cmd.color(carbon_color, "metal_sidechains and elem C")

    donors = "metal_shell and (elem N+O+S) and not name N+O"
    cmd.select("metal_donors", donors)
    donor_count = cmd.count_atoms("metal_donors")
    print(f"[show_metal:PLE] donors={donor_count}")
    for atom in cmd.get_model("metal_donors", state=1).atom:
        donor_sel = f"index {atom.index}"
        dist_name = f"metal_dist_{atom.resn}{atom.resi}_{atom.name}"
        cmd.distance(dist_name, "metal_core", donor_sel, cutoff=float(radius), mode=2, state=1)
        cmd.hide("labels", dist_name)
        cmd.color("yellow", dist_name)

    cmd.zoom("metal_shell or metal_core or resn LIG", 8)
    cmd.deselect()


cmd.extend("show_metal", show_metal)
cmd.extend("ple_show_metal", show_metal)
print("[show_metal:PLE] override loaded; commands: show_metal, ple_show_metal")
