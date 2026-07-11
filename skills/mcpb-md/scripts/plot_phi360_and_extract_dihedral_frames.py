#!/usr/bin/env python3
"""Plot 0-360 prochirality torsions and extract representative frames."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ACTIVE_SELECTION = "(resname UL1 FE1 HD1 HD2 GU1 AT1 HH1) or (resid 255 256 322 333)"
DONOR_BONDS = [
    ("FE1", "FE", "UL1", "N1"),
    ("FE1", "FE", "HD1", "NE2"),
    ("FE1", "FE", "HD2", "NE2"),
    ("FE1", "FE", "GU1", "OE1"),
    ("FE1", "FE", "AT1", "O2"),
    ("FE1", "FE", "HH1", "O"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--top", required=True)
    parser.add_argument("--traj", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--short-cutoff", type=float, default=3.0)
    return parser.parse_args()


def win_path(path: Path) -> str:
    s = path.as_posix()
    if s.startswith("/mnt/e/"):
        return "E:/" + s[len("/mnt/e/") :]
    return s


def wsl_path(path: Path) -> str:
    s = path.as_posix()
    if len(s) > 3 and s[1:3] == ":/":
        drive = s[0].lower()
        return f"/mnt/{drive}/" + s[3:]
    return s


def circular_delta_deg(a: pd.Series, target: float) -> pd.Series:
    return ((a - target + 180.0) % 360.0) - 180.0


def prepare_df(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df["phi_C1O1C4C5_360_deg"] = (df["dihedral_C1_O1_C4_C5_deg"] + 360.0) % 360.0
    df["abs_C1O1C4H02_deg"] = df["dihedral_C1_O1_C4_H02_deg"].abs()
    df["phi_side"] = np.where(df["phi_C1O1C4C5_360_deg"] < 180.0, "phi<180", "phi>180")
    return df


def choose_nearest(
    df: pd.DataFrame,
    label: str,
    target_phi: float,
    target_abs_h: float,
    used_frames: set[int],
    require_short: bool = False,
    short_cutoff: float = 3.0,
) -> pd.Series:
    pool = df.copy()
    if require_short:
        pool = pool[pool["C4_O1_A"] <= short_cutoff]
    pool = pool[~pool["frame"].astype(int).isin(used_frames)]
    if pool.empty:
        raise SystemExit(f"No candidate frames left for {label}")
    dphi = circular_delta_deg(pool["phi_C1O1C4C5_360_deg"], target_phi) / 30.0
    dh = (pool["abs_C1O1C4H02_deg"] - target_abs_h) / 25.0
    # Mildly prefer closer C4-O1 frames when several points are geometrically similar.
    dc4o1 = (pool["C4_O1_A"] - pool["C4_O1_A"].min()) / 2.5
    score = dphi * dphi + dh * dh + 0.04 * dc4o1 * dc4o1
    row = pool.loc[score.idxmin()].copy()
    row["selection_label"] = label
    row["target_phi360_deg"] = target_phi
    row["target_abs_H02_deg"] = target_abs_h
    row["selection_mode"] = "nearest_short" if require_short else "nearest_all"
    return row


def choose_min_c4o1(df: pd.DataFrame, label: str, phi_min: float, phi_max: float, used_frames: set[int]) -> pd.Series:
    pool = df[(df["C4_O1_A"] <= 3.0) & (df["phi_C1O1C4C5_360_deg"] >= phi_min) & (df["phi_C1O1C4C5_360_deg"] < phi_max)]
    pool = pool[~pool["frame"].astype(int).isin(used_frames)]
    if pool.empty:
        raise SystemExit(f"No short C4-O1 candidate for {label}")
    row = pool.loc[pool["C4_O1_A"].idxmin()].copy()
    row["selection_label"] = label
    row["target_phi360_deg"] = (phi_min + phi_max) / 2.0
    row["target_abs_H02_deg"] = np.nan
    row["selection_mode"] = "min_C4O1_in_phi_window"
    return row


def select_representative_frames(df: pd.DataFrame, short_cutoff: float) -> pd.DataFrame:
    targets = [
        ("all_phi150_H60", 150.0, 60.0, False),
        ("all_phi150_H100", 150.0, 100.0, False),
        ("all_phi150_H150", 150.0, 150.0, False),
        ("all_phi210_H60", 210.0, 60.0, False),
        ("all_phi210_H100", 210.0, 100.0, False),
        ("all_phi210_H150", 210.0, 150.0, False),
        ("short_phi150_H80", 150.0, 80.0, True),
        ("short_phi210_H80", 210.0, 80.0, True),
    ]
    rows = []
    used: set[int] = set()
    for label, phi, abs_h, require_short in targets:
        row = choose_nearest(df, label, phi, abs_h, used, require_short=require_short, short_cutoff=short_cutoff)
        used.add(int(row["frame"]))
        rows.append(row)
    for label, low, high in [
        ("short_min_C4O1_phi_0_180", 0.0, 180.0),
        ("short_min_C4O1_phi_180_360", 180.0, 360.0),
    ]:
        row = choose_min_c4o1(df, label, low, high, used)
        used.add(int(row["frame"]))
        rows.append(row)
    out = pd.DataFrame(rows)
    keep = [
        "selection_label",
        "selection_mode",
        "frame",
        "time_ns",
        "C4_O1_A",
        "dihedral_C1_O1_C4_C5_deg",
        "phi_C1O1C4C5_360_deg",
        "dihedral_C1_O1_C4_H02_deg",
        "abs_C1O1C4H02_deg",
        "pseudo_CIP_O1_C5_C3_H02_chiral_volume_A3",
        "lactone_plane_H02_C5_same_side_product_A2",
        "PHE322_ring_centroid_to_C4_A",
        "C4_radical_plane_normal_vs_PHE322_ring_normal_angle_deg",
        "target_phi360_deg",
        "target_abs_H02_deg",
    ]
    return out[keep].sort_values("frame").reset_index(drop=True)


def scatter_plots(df: pd.DataFrame, selected: pd.DataFrame, outdir: Path, short_cutoff: float) -> None:
    x = df["phi_C1O1C4C5_360_deg"]
    y = df["abs_C1O1C4H02_deg"]
    short = df[df["C4_O1_A"] <= short_cutoff]

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.scatter(x, y, s=3.0, c="black", alpha=0.035, linewidths=0, marker="o")
    ax.scatter(short["phi_C1O1C4C5_360_deg"], short["abs_C1O1C4H02_deg"], s=5.0, c="#1b9e77", alpha=0.23, linewidths=0)
    ax.scatter(selected["phi_C1O1C4C5_360_deg"], selected["abs_C1O1C4H02_deg"], s=28, c="#d95f02", edgecolors="white", linewidths=0.45, zorder=5)
    for i, row in selected.reset_index(drop=True).iterrows():
        ax.text(row["phi_C1O1C4C5_360_deg"] + 2.0, row["abs_C1O1C4H02_deg"] + 1.2, str(i + 1), fontsize=7, color="#8c2d04")
    ax.axvline(180.0, color="0.70", lw=0.8, ls="--")
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 180)
    ax.set_xlabel("Wrapped dihedral C1-O1-C4-C5(Ph), phi360 (deg)")
    ax.set_ylabel("|Dihedral C1-O1-C4-H02| (deg)")
    ax.set_xticks(np.arange(0, 361, 60))
    ax.set_yticks(np.arange(0, 181, 30))
    ax.set_title("All frames with C4-O1 <= 3 A overlay")
    fig.tight_layout()
    fig.savefig(outdir / "scatter_phi360_all_frames_with_short_and_selected.png", dpi=320)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    colors = short["C4_O1_A"]
    sc = ax.scatter(short["phi_C1O1C4C5_360_deg"], short["abs_C1O1C4H02_deg"], s=8.0, c=colors, cmap="viridis_r", alpha=0.62, linewidths=0, marker="o")
    ax.scatter(selected["phi_C1O1C4C5_360_deg"], selected["abs_C1O1C4H02_deg"], s=34, c="#d95f02", edgecolors="white", linewidths=0.45, zorder=5)
    for i, row in selected.reset_index(drop=True).iterrows():
        ax.text(row["phi_C1O1C4C5_360_deg"] + 2.0, row["abs_C1O1C4H02_deg"] + 1.2, str(i + 1), fontsize=7, color="#8c2d04")
    ax.axvline(180.0, color="0.70", lw=0.8, ls="--")
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 180)
    ax.set_xlabel("Wrapped dihedral C1-O1-C4-C5(Ph), phi360 (deg)")
    ax.set_ylabel("|Dihedral C1-O1-C4-H02| (deg)")
    ax.set_xticks(np.arange(0, 361, 60))
    ax.set_yticks(np.arange(0, 181, 30))
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("C4-O1 distance (A)")
    ax.set_title(f"C4-O1 <= {short_cutoff:.1f} A frames")
    fig.tight_layout()
    fig.savefig(outdir / "scatter_phi360_C4O1_le_3A_selected.png", dpi=320)
    plt.close(fig)


def atom_lines(path: Path) -> list[str]:
    lines = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            lines.append(line.rstrip())
    return lines


def build_multistate_pdb(frame_pdbs: list[Path], selected: pd.DataFrame, out_pdb: Path) -> None:
    lines = [
        "REMARK Representative C1-O1-C4-C5 phi360 dihedral states.",
        "REMARK MODEL order follows selected_frames_phi360.csv order.",
    ]
    for idx, (pdb, row) in enumerate(zip(frame_pdbs, selected.itertuples(index=False)), start=1):
        lines.append(f"MODEL     {idx:4d}")
        lines.append(
            "REMARK "
            f"idx={idx} label={row.selection_label} frame={int(row.frame)} time_ns={float(row.time_ns):.3f} "
            f"C4_O1_A={float(row.C4_O1_A):.3f} phi360={float(row.phi_C1O1C4C5_360_deg):.2f} "
            f"absH={float(row.abs_C1O1C4H02_deg):.2f}"
        )
        lines.extend(atom_lines(pdb))
        lines.append("ENDMDL")
    lines.append("END")
    out_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")


def extract_frames(top: str, traj: str, selected: pd.DataFrame, outdir: Path) -> tuple[Path, Path]:
    u = mda.Universe(top, traj)
    active = u.select_atoms(ACTIVE_SELECTION)
    if len(active) == 0:
        raise SystemExit(f"Active-site selection matched 0 atoms: {ACTIVE_SELECTION}")
    frame_dir = outdir / "representative_active_site_pdbs"
    full_frame_dir = outdir / "representative_full_system_pdbs"
    frame_dir.mkdir(parents=True, exist_ok=True)
    full_frame_dir.mkdir(parents=True, exist_ok=True)
    active_frame_pdbs: list[Path] = []
    full_frame_pdbs: list[Path] = []
    for idx, row in enumerate(selected.itertuples(index=False), start=1):
        frame = int(row.frame)
        u.trajectory[frame]
        pdb = frame_dir / f"{idx:02d}_{row.selection_label}_frame{frame}.pdb"
        active.write(str(pdb))
        active_frame_pdbs.append(pdb)
        full_pdb = full_frame_dir / f"{idx:02d}_{row.selection_label}_frame{frame}_full_system.pdb"
        u.atoms.write(str(full_pdb))
        full_frame_pdbs.append(full_pdb)
    active_multistate = outdir / "representative_phi360_active_site_multistate.pdb"
    full_multistate = outdir / "representative_phi360_full_system_multistate.pdb"
    build_multistate_pdb(active_frame_pdbs, selected, active_multistate)
    build_multistate_pdb(full_frame_pdbs, selected, full_multistate)
    return active_multistate, full_multistate


def pml_text(multistate_pdb: Path, selected_csv: Path, outdir: Path, path_style: str = "windows") -> str:
    obj = "phi360_representative_states"
    path_fn = wsl_path if path_style == "wsl" else win_path
    bond_lines = []
    for fe_resn, fe_name, resn, atom in DONOR_BONDS:
        bond_lines.append(f"bond {obj} and resn {fe_resn} and name {fe_name}, {obj} and resn {resn} and name {atom}")
    return f"""reinitialize
bg_color white

load {path_fn(multistate_pdb)}, {obj}

hide everything
show sticks, {obj}
show spheres, {obj} and resn FE1

color gray82, {obj}
color orange, {obj} and resn FE1
color gray90, {obj} and resn UL1
color gray70, {obj} and resn HD1+HD2+GU1+AT1+HH1
color slate, {obj} and resi 255+256+322+333

set stick_radius, 0.15
set sphere_scale, 0.35, {obj} and resn FE1
set valence, 0
set two_sided_lighting, on
set dash_width, 2.2
set dash_radius, 0.04

{chr(10).join(bond_lines)}

distance C4_O1, {obj} and resn UL1 and name C4, {obj} and resn UL1 and name O1
color red, C4_O1

label {obj} and resn UL1 and name C4, "C4"
label {obj} and resn UL1 and name O1, "O1"
label {obj} and resn UL1 and name C5, "C5/Ph"
label {obj} and resn UL1 and name H02, "H02"
set label_size, 14
set label_color, black

orient {obj} and resn UL1+FE1+HD1+HD2+GU1+AT1+HH1
zoom {obj} and (resn UL1+FE1+HD1+HD2+GU1+AT1+HH1 or resi 255+256+322+333), 7

set movie_fps, 3
mset 1 -10

# Frame metadata:
# {path_fn(selected_csv)}
# To save a PyMOL session:
# save {path_fn(outdir / 'representative_phi360_states.pse')}
"""


def full_system_pml_text(multistate_pdb: Path, selected_csv: Path, outdir: Path, path_style: str = "windows") -> str:
    obj = "phi360_full_system_states"
    path_fn = wsl_path if path_style == "wsl" else win_path
    active_sel = f"{obj} and (resn UL1+FE1+HD1+HD2+GU1+AT1+HH1 or resi 255+256+322+333)"
    bond_lines = []
    for fe_resn, fe_name, resn, atom in DONOR_BONDS:
        bond_lines.append(f"bond {obj} and resn {fe_resn} and name {fe_name}, {obj} and resn {resn} and name {atom}")
    return f"""reinitialize
bg_color white

load {path_fn(multistate_pdb)}, {obj}

dss {obj} and polymer
rebuild

hide everything
show cartoon, {obj} and polymer
show sticks, {active_sel}
show spheres, {obj} and resn FE1

color gray80, {obj} and polymer
color orange, {obj} and resn FE1
color gray90, {obj} and resn UL1
color gray70, {obj} and resn HD1+HD2+GU1+AT1+HH1
color cyan, {obj} and resi 322
color slate, {obj} and resi 255+256+333
color salmon, {obj} and resn GU1

set cartoon_transparency, 0.68, {obj} and polymer
set stick_radius, 0.15
set sphere_scale, 0.35, {obj} and resn FE1
set valence, 0
set two_sided_lighting, on
set dash_width, 2.2
set dash_radius, 0.04

{chr(10).join(bond_lines)}

distance C4_O1, {obj} and resn UL1 and name C4, {obj} and resn UL1 and name O1
distance PHE322_C4, {obj} and resi 322 and name CG+CD1+CE1+CZ+CE2+CD2, {obj} and resn UL1 and name C4
color red, C4_O1
color cyan, PHE322_C4

label {obj} and resn UL1 and name C4, "C4"
label {obj} and resn UL1 and name O1, "O1"
label {obj} and resn UL1 and name C5, "C5/Ph"
label {obj} and resn UL1 and name H02, "H02"
label {obj} and resi 322 and name CZ, "PHE322/F336"
set label_size, 13
set label_color, black

orient {active_sel}
zoom {active_sel}, 8

set movie_fps, 3
mset 1 -10

# Full protein/protein+MCPB representative states.
# The original active-site-only PDBs are still kept for quick loading.
# Frame metadata:
# {path_fn(selected_csv)}
# To save a PyMOL session:
# save {path_fn(outdir / 'representative_phi360_full_system_states.pse')}
"""


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = prepare_df(args.csv)
    selected = select_representative_frames(df, args.short_cutoff)
    selected_csv = outdir / "selected_frames_phi360.csv"
    selected.to_csv(selected_csv, index=False)
    scatter_plots(df, selected, outdir, args.short_cutoff)
    active_multistate, full_multistate = extract_frames(args.top, args.traj, selected, outdir)
    pml = outdir / "view_representative_phi360_states.pml"
    pml.write_text(pml_text(active_multistate, selected_csv, outdir), encoding="utf-8")
    pml_wsl = outdir / "view_representative_phi360_states_wsl.pml"
    pml_wsl.write_text(pml_text(active_multistate, selected_csv, outdir, path_style="wsl"), encoding="utf-8")
    full_pml = outdir / "view_representative_phi360_full_system_states.pml"
    full_pml.write_text(full_system_pml_text(full_multistate, selected_csv, outdir), encoding="utf-8")
    full_pml_wsl = outdir / "view_representative_phi360_full_system_states_wsl.pml"
    full_pml_wsl.write_text(full_system_pml_text(full_multistate, selected_csv, outdir, path_style="wsl"), encoding="utf-8")
    readme = [
        "# Phi360 Prochirality Representative Frames",
        "",
        "This directory wraps `dihedral_C1_O1_C4_C5_deg` into 0-360 degrees:",
        "",
        "`phi360 = (dihedral_C1_O1_C4_C5_deg + 360) % 360`",
        "",
        "This avoids splitting geometries near -180/+180 across opposite edges of the plot.",
        "",
        "Important interpretation:",
        "- The dihedral is a product-like prochirality descriptor around the putative O1-C4 forming bond.",
        "- It does not by itself assign R/S. Map phi360 basins to pro-R/pro-S only after calibration against an explicit product model.",
        "- The selected frames intentionally cover both sides of the 180 degree boundary and several H02 orientations.",
        "",
        "Open in PyMOL:",
        "",
        "Recommended full-system view:",
        f"Windows PyMOL: `@{win_path(full_pml)}`",
        f"WSL PyMOL: `@{wsl_path(full_pml_wsl)}`",
        "",
        "Active-site-only quick view:",
        f"Windows PyMOL: `@{win_path(pml)}`",
        f"WSL PyMOL: `@{wsl_path(pml_wsl)}`",
        "",
    ]
    (outdir / "README_phi360_representatives.md").write_text("\n".join(readme), encoding="utf-8")
    print(outdir)


if __name__ == "__main__":
    main()
