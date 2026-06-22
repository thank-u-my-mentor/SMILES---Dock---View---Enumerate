"""MD post-processing helpers for PLE."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path


def run_command(cmd: list[str], *, stdin: str | None = None) -> None:
    print("$ " + " ".join(str(part) for part in cmd), flush=True)
    subprocess.run(cmd, input=stdin, text=True, check=True)


def gmx_env_prefix(args: argparse.Namespace) -> str:
    lines = []
    for path in getattr(args, "extra_ld_library_path", []) or []:
        lines.append(f"export LD_LIBRARY_PATH={path}:${{LD_LIBRARY_PATH:-}}")
    if getattr(args, "gmxrc", None):
        lines.append(f"source {args.gmxrc}")
    return "\n".join(lines)


def fix_pbc_command(args: argparse.Namespace) -> None:
    gmxdir = args.gmxdir.expanduser().resolve()
    gmx = str(args.gmx)
    tpr = gmxdir / args.tpr
    xtc = gmxdir / args.xtc
    nojump = gmxdir / args.nojump
    centered = gmxdir / args.centered
    centered_gro = gmxdir / args.centered_gro

    if args.gmxrc:
        run_command(["bash", "-lc", f"{gmx_env_prefix(args)}\ncd {gmxdir}\nprintf '0\\n' | {gmx} trjconv -s {tpr.name} -f {xtc.name} -o {nojump.name} -pbc nojump"])
        run_command(["bash", "-lc", f"{gmx_env_prefix(args)}\ncd {gmxdir}\nprintf '{args.center_group}\\n{args.output_group}\\n' | {gmx} trjconv -s {tpr.name} -f {nojump.name} -o {centered.name} -pbc mol -center -ur compact"])
        run_command(["bash", "-lc", f"{gmx_env_prefix(args)}\ncd {gmxdir}\nprintf '{args.output_group}\\n' | {gmx} trjconv -s {tpr.name} -f {centered.name} -o {centered_gro.name} -dump 0"])
    else:
        run_command([gmx, "trjconv", "-s", str(tpr), "-f", str(xtc), "-o", str(nojump), "-pbc", "nojump"], stdin="0\n")
        run_command([gmx, "trjconv", "-s", str(tpr), "-f", str(nojump), "-o", str(centered), "-pbc", "mol", "-center", "-ur", "compact"], stdin=f"{args.center_group}\n{args.output_group}\n")
        run_command([gmx, "trjconv", "-s", str(tpr), "-f", str(centered), "-o", str(centered_gro), "-dump", "0"], stdin=f"{args.output_group}\n")

    print(f"nojump={nojump}", flush=True)
    print(f"centered={centered}", flush=True)
    print(f"centered_gro={centered_gro}", flush=True)


def write_show_metal_pml(
    *,
    structure: Path,
    trajectory: Path | None,
    out: Path,
    metal_selection: str,
    radius: float,
    carbon_color: str,
    load_trajectory: bool = False,
) -> None:
    traj_line = f"load_traj {trajectory.as_posix()}, complex\n" if trajectory and load_trajectory else ""
    text = f"""# Auto-generated PLE metal-site view.
reinitialize
load {structure.as_posix()}, complex
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
show spheres, ({metal_selection})
set sphere_scale, 0.45, ({metal_selection})

python
from pymol import cmd

def ple_show_metal(selection="{metal_selection}", radius={radius:.3f}):
    cmd.delete("metal_dist_*")
    cmd.select("metal_core", selection)
    cmd.select("metal_scope", f"byres ((not solvent and not resn LIG and not ({{selection}})) within {{radius}} of ({{selection}}))")
    cmd.select("metal_shell", "metal_scope and not resn SOL+HOH+WAT+NA+CL+K+MG+CA")
    cmd.select("metal_sidechains", "metal_shell and not name N+C+O+CA+H+HA")
    print(f"[ple_show_metal] metals={{cmd.count_atoms('metal_core')}} shell_atoms={{cmd.count_atoms('metal_shell')}} sidechain_atoms={{cmd.count_atoms('metal_sidechains')}}")
    cmd.show("sticks", "metal_sidechains")
    cmd.color("{carbon_color}", "metal_sidechains and elem C")
    cmd.show("spheres", "metal_core")
    cmd.set("sphere_scale", 0.45, "metal_core")
    cmd.select("metal_donors", "metal_shell and (elem N+O+S) and not name N+O")
    print(f"[ple_show_metal] donors={{cmd.count_atoms('metal_donors')}}")
    for atom in cmd.get_model("metal_donors", state=1).atom:
        dist_name = f"metal_dist_{{atom.resn}}{{atom.resi}}_{{atom.name}}"
        cmd.distance(dist_name, "metal_core", f"index {{atom.index}}", cutoff=radius, mode=2, state=1)
        cmd.hide("labels", dist_name)
        cmd.color("yellow", dist_name)
    cmd.zoom("metal_shell or metal_core or resn LIG", 8)
    cmd.deselect()

cmd.extend("ple_show_metal", ple_show_metal)
python end

ple_show_metal
"""
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")


def pymol_metal_command(args: argparse.Namespace) -> None:
    out = args.out.expanduser().resolve()
    write_show_metal_pml(
        structure=args.structure.expanduser().resolve(),
        trajectory=args.trajectory.expanduser().resolve() if args.trajectory else None,
        out=out,
        metal_selection=args.metal_selection,
        radius=args.radius,
        carbon_color=args.carbon_color,
        load_trajectory=args.load_trajectory,
    )
    print(f"pml={out}", flush=True)


def _python() -> str:
    return sys.executable


def find_codex_scripts_dir() -> Path:
    """Find the shared Codex scripts directory while PLE migration is in progress."""
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / "scripts" / "ple_md_rdc.py"
        if candidate.exists():
            return parent / "scripts"
    fallback = Path("/mnt/e/Codex/scripts")
    if (fallback / "ple_md_rdc.py").exists():
        return fallback
    raise FileNotFoundError("could not find ple_md_rdc.py; expected a Codex scripts directory")


def rdc_command(args: argparse.Namespace) -> None:
    scripts = find_codex_scripts_dir()
    gmxdir = args.gmxdir.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    topology = args.topology or (gmxdir / "md.gro")
    trajectory = args.trajectory or (gmxdir / "md.xtc")
    run_command([
        _python(), str(scripts / "ple_md_rdc.py"),
        "--topology", str(topology),
        "--trajectory", str(trajectory),
        "--outdir", str(outdir),
        "--ligand-selection", args.ligand_selection,
        "--distance-cutoff", str(args.distance_cutoff),
        "--stride", str(args.stride),
        "--top-n", str(args.top_n),
        "--reference", args.reference,
        "--plot-traces",
    ])
    if args.reference_pdb and args.gro:
        run_command([
            _python(), str(scripts / "map_rdc_labels_to_reference_pdb.py"),
            "--rdc-dir", str(outdir),
            "--reference-pdb", str(args.reference_pdb.expanduser().resolve()),
            "--gro", str(args.gro.expanduser().resolve()),
        ])
    run_command([
        _python(), str(scripts / "build_rdc_dashboard.py"),
        "--rdc-dir", str(outdir),
        "--title", args.title,
        "--subtitle", args.subtitle,
        "--velocity-threshold", str(args.velocity_threshold),
    ])
    print(f"dashboard={outdir / 'index.html'}", flush=True)


def add_md_tool_parsers(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    pbc = subparsers.add_parser("md-fix-pbc", help="Write nojump/centered compact trajectory files for visualization")
    pbc.add_argument("--gmxdir", type=Path, default=Path("md_handoff/gromacs"))
    pbc.add_argument("--gmx", default="gmx_mpi")
    pbc.add_argument("--gmxrc", type=Path)
    pbc.add_argument("--extra-ld-library-path", action="append", default=[])
    pbc.add_argument("--tpr", default="md.tpr")
    pbc.add_argument("--xtc", default="md.xtc")
    pbc.add_argument("--nojump", default="md_nojump.xtc")
    pbc.add_argument("--centered", default="md_centered_compact.xtc")
    pbc.add_argument("--centered-gro", default="md_centered_compact.gro")
    pbc.add_argument("--center-group", default="1", help="GROMACS index group used for centering; default 1 Protein")
    pbc.add_argument("--output-group", default="0", help="GROMACS output group; default 0 System")

    rdc = subparsers.add_parser("md-rdc", help="Build residue-distance-change dashboard from a GROMACS MD trajectory")
    rdc.add_argument("--gmxdir", type=Path, default=Path("md_handoff/gromacs"))
    rdc.add_argument("--outdir", type=Path, default=Path("md_analysis/rdc"))
    rdc.add_argument("--topology", type=Path)
    rdc.add_argument("--trajectory", type=Path)
    rdc.add_argument("--gro", type=Path, help="GRO used for residue-number mapping")
    rdc.add_argument("--reference-pdb", type=Path)
    rdc.add_argument("--ligand-selection", default="resname LIG")
    rdc.add_argument("--distance-cutoff", type=float, default=10.0)
    rdc.add_argument("--stride", type=int, default=1)
    rdc.add_argument("--top-n", type=int, default=40)
    rdc.add_argument("--reference", choices=["first", "mean-first-10pct"], default="mean-first-10pct")
    rdc.add_argument("--velocity-threshold", type=float, default=2.0)
    rdc.add_argument("--title", default="RDC 动态分析 Dashboard")
    rdc.add_argument("--subtitle", default="Residue-ligand distance-change analysis from production MD.")

    metal = subparsers.add_parser("pymol-metal", help="Write a PyMOL show_metal script for a centered MD structure/trajectory")
    metal.add_argument("--structure", type=Path, required=True)
    metal.add_argument("--trajectory", type=Path)
    metal.add_argument("--out", type=Path, default=Path("pymol/show_metal_centered.pml"))
    metal.add_argument("--metal-selection", default="resn FE+FE2+FE3 or elem Fe")
    metal.add_argument("--radius", type=float, default=4.0)
    metal.add_argument("--carbon-color", default="green")
    metal.add_argument("--load-trajectory", action="store_true", help="Load full trajectory into PyMOL; disabled by default to avoid memory blowups")
