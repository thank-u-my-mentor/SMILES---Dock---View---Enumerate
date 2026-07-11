#!/usr/bin/env python3
"""Build a PyMOL-ready bundle for closest reactive-distance frames."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import MDAnalysis as mda
import pandas as pd


DONOR_BONDS = [
    ("FE1", "FE", "UL1", "N1"),
    ("FE1", "FE", "HD1", "NE2"),
    ("FE1", "FE", "HD2", "NE2"),
    ("FE1", "FE", "GU1", "OE1"),
    ("FE1", "FE", "AT1", "O2"),
    ("FE1", "FE", "HH1", "O"),
]

ACTIVE_RESN = "UL1+FE1+HD1+HD2+GU1+AT1+HH1"
CONTEXT_RESI = "255+256+322+333"
ACTIVE_RESN_SET = set(ACTIVE_RESN.split("+"))
CONTEXT_RESI_SET = set(CONTEXT_RESI.split("+"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--top", required=True)
    parser.add_argument("--traj", required=True)
    parser.add_argument("--n", type=int, default=20)
    return parser.parse_args()


def atom_records(path: Path) -> list[str]:
    rows = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            rows.append(line.rstrip())
    return rows


def parse_pdb_serials(first_pdb: Path) -> dict[tuple[str, str], int]:
    serials = {}
    for line in atom_records(first_pdb):
        serial = int(line[6:11])
        name = line[12:16].strip()
        resn = line[17:20].strip()
        serials[(resn, name)] = serial
    return serials


def connect_lines(serials: dict[tuple[str, str], int]) -> list[str]:
    lines = []
    for fe_resn, fe_name, resn, name in DONOR_BONDS:
        fe = serials.get((fe_resn, fe_name))
        donor = serials.get((resn, name))
        if fe is not None and donor is not None:
            lines.append(f"CONECT{fe:5d}{donor:5d}")
            lines.append(f"CONECT{donor:5d}{fe:5d}")
    return lines


def build_multistate_pdb(frames: pd.DataFrame, out_pdb: Path) -> None:
    serials = parse_pdb_serials(Path(frames.iloc[0]["pdb"]))
    lines = [
        "REMARK Closest UL1 C4-O1 MD frames merged as PyMOL states.",
        "REMARK MODEL number follows closest-distance rank, not chronological order.",
    ]
    for row in frames.itertuples(index=False):
        lines.append(f"MODEL     {int(row.rank):4d}")
        lines.append(
            f"REMARK rank={int(row.rank)} frame={int(row.frame)} time_ns={float(row.time_ns):.3f} C4_O1_A={float(row.C4_O1_distance_A):.3f}"
        )
        lines.extend(atom_records(Path(row.pdb)))
        lines.append("ENDMDL")
    lines.extend(connect_lines(serials))
    lines.append("END")
    out_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")


def is_active_site_record(line: str) -> bool:
    resn = line[17:20].strip()
    resid = line[22:26].strip()
    return resn in ACTIVE_RESN_SET or resid in CONTEXT_RESI_SET


def build_active_site_multistate_pdb(frames: pd.DataFrame, out_pdb: Path) -> None:
    serials = parse_pdb_serials(Path(frames.iloc[0]["pdb"]))
    lines = [
        "REMARK Active-site-only closest UL1 C4-O1 MD frames.",
        "REMARK Includes UL1, Fe/donors, and context residues 255/256/322/333.",
    ]
    for row in frames.itertuples(index=False):
        lines.append(f"MODEL     {int(row.rank):4d}")
        lines.append(
            f"REMARK rank={int(row.rank)} frame={int(row.frame)} time_ns={float(row.time_ns):.3f} C4_O1_A={float(row.C4_O1_distance_A):.3f}"
        )
        lines.extend(line for line in atom_records(Path(row.pdb)) if is_active_site_record(line))
        lines.append("ENDMDL")
    lines.extend(connect_lines(serials))
    lines.append("END")
    out_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_subset_xtc(frames: pd.DataFrame, top: str, traj: str, out_gro: Path, out_xtc: Path) -> None:
    u = mda.Universe(top, traj)
    u.trajectory[int(frames.iloc[0]["frame"])]
    u.atoms.write(str(out_gro))
    with mda.Writer(str(out_xtc), n_atoms=u.atoms.n_atoms) as writer:
        for frame in frames["frame"]:
            u.trajectory[int(frame)]
            writer.write(u.atoms)


def pml_text(multistate_pdb: Path, active_site_pdb: Path, subset_gro: Path, subset_xtc: Path, outdir: Path) -> str:
    pdb_posix = multistate_pdb.as_posix().replace("/mnt/e/", "E:/")
    active_posix = active_site_pdb.as_posix().replace("/mnt/e/", "E:/")
    gro_posix = subset_gro.as_posix().replace("/mnt/e/", "E:/")
    xtc_posix = subset_xtc.as_posix().replace("/mnt/e/", "E:/")
    return f"""reinitialize

load {pdb_posix}, closest_C4O1_states
load {active_posix}, closest_C4O1_active_site

dss closest_C4O1_states and polymer
rebuild

hide everything
show cartoon, closest_C4O1_states and polymer
show sticks, closest_C4O1_active_site
show spheres, closest_C4O1_active_site and resn FE1

color gray80, closest_C4O1_states
color gray65, closest_C4O1_states and polymer
color gray82, closest_C4O1_active_site
color orange, resn FE1
color gray90, resn UL1
color gray75, resn HD1+HD2+GU1+AT1+HH1
color gray55, resi {CONTEXT_RESI}

set stick_radius, 0.15
set sphere_scale, 0.35, resn FE1
set dash_width, 2.5
set dash_radius, 0.045
set valence, 0
set two_sided_lighting, on
set cartoon_transparency, 0.68, closest_C4O1_states

bond closest_C4O1_active_site and resn FE1 and name FE, closest_C4O1_active_site and resn UL1 and name N1
bond closest_C4O1_active_site and resn FE1 and name FE, closest_C4O1_active_site and resn HD1 and name NE2
bond closest_C4O1_active_site and resn FE1 and name FE, closest_C4O1_active_site and resn HD2 and name NE2
bond closest_C4O1_active_site and resn FE1 and name FE, closest_C4O1_active_site and resn GU1 and name OE1
bond closest_C4O1_active_site and resn FE1 and name FE, closest_C4O1_active_site and resn AT1 and name O2
bond closest_C4O1_active_site and resn FE1 and name FE, closest_C4O1_active_site and resn HH1 and name O

distance C4_O1, closest_C4O1_active_site and resn UL1 and name C4, closest_C4O1_active_site and resn UL1 and name O1
color red, C4_O1

show sticks, closest_C4O1_active_site and resn UL1 and name C4+O1
label closest_C4O1_active_site and resn UL1 and name C4, "C4"
label closest_C4O1_active_site and resn UL1 and name O1, "O1"
set label_size, 16
set label_color, black

orient closest_C4O1_active_site
zoom closest_C4O1_active_site, 7

set movie_fps, 4
mset 1 -20

# Optional: load the same frames as an XTC trajectory instead of multi-state PDB.
# load {gro_posix}, closest_C4O1_xtc
# load_traj {xtc_posix}, closest_C4O1_xtc

# To save an editable PyMOL session after running this file:
# save {outdir.as_posix().replace('/mnt/e/', 'E:/')}/closest_C4O1_states.pse
"""


def xtc_pml_text(subset_gro: Path, subset_xtc: Path, outdir: Path) -> str:
    gro_posix = subset_gro.as_posix().replace("/mnt/e/", "E:/")
    xtc_posix = subset_xtc.as_posix().replace("/mnt/e/", "E:/")
    return f"""reinitialize

load {gro_posix}, closest_C4O1_xtc
load_traj {xtc_posix}, closest_C4O1_xtc

dss closest_C4O1_xtc and polymer
rebuild

hide everything
show cartoon, closest_C4O1_xtc and polymer
show sticks, closest_C4O1_xtc and (resn {ACTIVE_RESN} or resi {CONTEXT_RESI})
show spheres, closest_C4O1_xtc and resn FE1

color gray80, closest_C4O1_xtc
color gray65, closest_C4O1_xtc and polymer
color orange, closest_C4O1_xtc and resn FE1
color gray90, closest_C4O1_xtc and resn UL1
color gray75, closest_C4O1_xtc and resn HD1+HD2+GU1+AT1+HH1
color gray55, closest_C4O1_xtc and resi {CONTEXT_RESI}

set stick_radius, 0.15
set sphere_scale, 0.35, closest_C4O1_xtc and resn FE1
set dash_width, 2.5
set dash_radius, 0.045
set valence, 0
set two_sided_lighting, on
set cartoon_transparency, 0.68, closest_C4O1_xtc

# Manual Fe coordination bonds. This is necessary because GRO/XTC carries no CONECT/topology bonds.
# PyMOL preset guesses metal bonds from state-1 distance heuristics; HD1 NE2 is often just outside that cutoff.
bond closest_C4O1_xtc and resn FE1 and name FE, closest_C4O1_xtc and resn UL1 and name N1
bond closest_C4O1_xtc and resn FE1 and name FE, closest_C4O1_xtc and resn HD1 and name NE2
bond closest_C4O1_xtc and resn FE1 and name FE, closest_C4O1_xtc and resn HD2 and name NE2
bond closest_C4O1_xtc and resn FE1 and name FE, closest_C4O1_xtc and resn GU1 and name OE1
bond closest_C4O1_xtc and resn FE1 and name FE, closest_C4O1_xtc and resn AT1 and name O2
bond closest_C4O1_xtc and resn FE1 and name FE, closest_C4O1_xtc and resn HH1 and name O

distance C4_O1_xtc, closest_C4O1_xtc and resn UL1 and name C4, closest_C4O1_xtc and resn UL1 and name O1
color red, C4_O1_xtc

label closest_C4O1_xtc and resn UL1 and name C4, "C4"
label closest_C4O1_xtc and resn UL1 and name O1, "O1"
set label_size, 16
set label_color, black

orient closest_C4O1_xtc and (resn {ACTIVE_RESN} or resi {CONTEXT_RESI})
zoom closest_C4O1_xtc and (resn {ACTIVE_RESN} or resi {CONTEXT_RESI}), 7

set movie_fps, 4
mset 1 -20

# To save an editable PyMOL session after running this file:
# save {outdir.as_posix().replace('/mnt/e/', 'E:/')}/closest_C4O1_xtc_forced_bonds.pse
"""


def safe_pml_text(multistate_pdb: Path, active_site_pdb: Path, outdir: Path) -> str:
    pdb_posix = multistate_pdb.as_posix().replace("/mnt/e/", "E:/")
    active_posix = active_site_pdb.as_posix().replace("/mnt/e/", "E:/")
    return f"""reinitialize
bg_color white

load {pdb_posix}, closest_C4O1_states
load {active_posix}, closest_C4O1_active_site

# Safe fallback display first: do not depend on PyMOL recognizing polymer/cartoon.
hide everything
show lines, closest_C4O1_states
show spheres, closest_C4O1_active_site
show sticks, closest_C4O1_active_site

# Then try secondary-structure cartoon. If dss does not classify this PDB, lines/spheres remain visible.
dss closest_C4O1_states
show cartoon, closest_C4O1_states

color gray88, closest_C4O1_states
color gray72, closest_C4O1_states and polymer
color gray82, closest_C4O1_active_site
color orange, closest_C4O1_active_site and resn FE1
color gray90, closest_C4O1_active_site and resn UL1
color gray75, closest_C4O1_active_site and resn HD1+HD2+GU1+AT1+HH1
color gray45, closest_C4O1_active_site and resi {CONTEXT_RESI}

set stick_radius, 0.15
set sphere_scale, 0.22, closest_C4O1_active_site
set sphere_scale, 0.40, closest_C4O1_active_site and resn FE1
set line_width, 1.0
set cartoon_transparency, 0.65, closest_C4O1_states
set two_sided_lighting, on
set valence, 0

bond closest_C4O1_active_site and resn FE1 and name FE, closest_C4O1_active_site and resn UL1 and name N1
bond closest_C4O1_active_site and resn FE1 and name FE, closest_C4O1_active_site and resn HD1 and name NE2
bond closest_C4O1_active_site and resn FE1 and name FE, closest_C4O1_active_site and resn HD2 and name NE2
bond closest_C4O1_active_site and resn FE1 and name FE, closest_C4O1_active_site and resn GU1 and name OE1
bond closest_C4O1_active_site and resn FE1 and name FE, closest_C4O1_active_site and resn AT1 and name O2
bond closest_C4O1_active_site and resn FE1 and name FE, closest_C4O1_active_site and resn HH1 and name O

distance C4_O1, closest_C4O1_active_site and resn UL1 and name C4, closest_C4O1_active_site and resn UL1 and name O1
color red, C4_O1
set dash_width, 2.5
set dash_radius, 0.045

label closest_C4O1_active_site and resn UL1 and name C4, "C4"
label closest_C4O1_active_site and resn UL1 and name O1, "O1"
set label_size, 16
set label_color, black

set all_states, off
set state, 1
set movie_fps, 4
mset 1 -20
frame 1

orient closest_C4O1_active_site
zoom closest_C4O1_active_site, 8
clip slab, 16
rebuild

# save {outdir.as_posix().replace('/mnt/e/', 'E:/')}/closest_C4O1_safe_view.pse
"""


def write_metadata(frames: pd.DataFrame, outdir: Path, multistate_pdb: Path, active_site_pdb: Path, subset_gro: Path, subset_xtc: Path) -> None:
    text = [
        "# Closest C4-O1 PyMOL Bundle",
        "",
        "Files:",
        f"- `{multistate_pdb.name}`: 20 closest C4-O1 frames as PyMOL states.",
        f"- `{active_site_pdb.name}`: lightweight active-site-only states for fast inspection.",
        f"- `{subset_gro.name}` + `{subset_xtc.name}`: same frames as a small trajectory subset.",
        "- `view_closest_C4O1_states.pml`: load and style the bundle.",
        "- `view_closest_C4O1_xtc_forced_bonds.pml`: load GRO/XTC and force Fe-donor bonds.",
        "- `view_closest_C4O1_safe.pml`: most robust fallback view; keeps lines/spheres visible even if cartoon recognition fails.",
        "",
        "PyMOL:",
        "",
        "```pymol",
        f"@{(outdir / 'view_closest_C4O1_states.pml').as_posix().replace('/mnt/e/', 'E:/')}",
        "```",
        "",
        "The PML manually bonds Fe to UL1 N1, HD1 NE2, HD2 NE2, GU1 OE1, AT1 O2, and HH1 O.",
        "This avoids relying on PyMOL preset heuristics for Fe-His bonds.",
        "",
        "Closest frames:",
        "",
        frames[["rank", "frame", "time_ns", "C4_O1_distance_A"]].to_string(index=False),
        "",
    ]
    (outdir / "PYMOL_BUNDLE_README.md").write_text("\n".join(text), encoding="utf-8")


def main() -> None:
    args = parse_args()
    analysis_dir = Path(args.analysis_dir)
    csv = analysis_dir / "closest_C4_O1_frames.csv"
    frames = pd.read_csv(csv).head(args.n).copy()

    outdir = analysis_dir / "pymol_closest_C4O1_bundle"
    outdir.mkdir(parents=True, exist_ok=True)
    multistate_pdb = outdir / "closest_C4O1_top20_multistate_with_conect.pdb"
    active_site_pdb = outdir / "closest_C4O1_top20_active_site_multistate.pdb"
    subset_gro = outdir / "closest_C4O1_top20.gro"
    subset_xtc = outdir / "closest_C4O1_top20.xtc"
    pml = outdir / "view_closest_C4O1_states.pml"
    xtc_pml = outdir / "view_closest_C4O1_xtc_forced_bonds.pml"
    safe_pml = outdir / "view_closest_C4O1_safe.pml"

    build_multistate_pdb(frames, multistate_pdb)
    build_active_site_multistate_pdb(frames, active_site_pdb)
    build_subset_xtc(frames, args.top, args.traj, subset_gro, subset_xtc)
    pml.write_text(pml_text(multistate_pdb, active_site_pdb, subset_gro, subset_xtc, outdir), encoding="utf-8")
    xtc_pml.write_text(xtc_pml_text(subset_gro, subset_xtc, outdir), encoding="utf-8")
    safe_pml.write_text(safe_pml_text(multistate_pdb, active_site_pdb, outdir), encoding="utf-8")
    write_metadata(frames, outdir, multistate_pdb, active_site_pdb, subset_gro, subset_xtc)
    print(outdir)


if __name__ == "__main__":
    main()
