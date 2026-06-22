param(
    [string]$Structure = "E:/TJ/md_handoff/gromacs/md_centered_compact.gro",
    [string]$Out = "E:/TJ/pymol/show_metal_centered_from_ple.pml",
    [string]$MetalSelection = "resn FE+FE2+FE3 or elem Fe",
    [double]$Radius = 4.0,
    [string]$CarbonColor = "palegreen"
)

$text = @"
# Auto-generated PLE metal-site view for Windows PyMOL.
reinitialize
load $Structure, complex

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
show spheres, ($MetalSelection)
set sphere_scale, 0.45, ($MetalSelection)

python
from pymol import cmd

def ple_show_metal(selection="$MetalSelection", radius=$($Radius.ToString("0.000", [Globalization.CultureInfo]::InvariantCulture))):
    cmd.delete("metal_dist_*")
    cmd.select("metal_core", selection)
    if cmd.count_atoms("metal_core") == 0:
        print(f"[ple_show_metal] no metal atoms matched: {selection}")
        return
    cmd.select("metal_scope", f"byres ((not solvent and not resn LIG and not ({selection})) within {radius} of ({selection}))")
    cmd.select("metal_shell", "metal_scope and not resn SOL+HOH+WAT+NA+CL+K+MG+CA")
    cmd.select("metal_sidechains", "metal_shell and not name N+C+O+CA+H+HA")
    print(f"[ple_show_metal] metals={cmd.count_atoms('metal_core')} shell_atoms={cmd.count_atoms('metal_shell')} sidechain_atoms={cmd.count_atoms('metal_sidechains')}")
    cmd.show("sticks", "metal_sidechains")
    cmd.color("$CarbonColor", "metal_sidechains and elem C")
    cmd.show("spheres", "metal_core")
    cmd.set("sphere_scale", 0.45, "metal_core")
    donors = "metal_shell and (elem N+O+S) and not name N+O"
    cmd.select("metal_donors", donors)
    donor_count = cmd.count_atoms("metal_donors")
    print(f"[ple_show_metal] donors={donor_count}")
    for atom in cmd.get_model("metal_donors", state=1).atom:
        donor_sel = f"index {atom.index}"
        dist_name = f"metal_dist_{atom.resn}{atom.resi}_{atom.name}"
        cmd.distance(dist_name, "metal_core", donor_sel, cutoff=radius, mode=2, state=1)
        cmd.hide("labels", dist_name)
        cmd.color("yellow", dist_name)
    cmd.zoom("metal_shell or metal_core or resn LIG", 8)
    cmd.deselect()

cmd.extend("ple_show_metal", ple_show_metal)
cmd.extend("show_metal", ple_show_metal)
python end

ple_show_metal
"@

$outPath = [System.IO.Path]::GetFullPath($Out)
New-Item -ItemType Directory -Force -Path ([System.IO.Path]::GetDirectoryName($outPath)) | Out-Null
Set-Content -LiteralPath $outPath -Value $text -Encoding UTF8
Write-Host $outPath
