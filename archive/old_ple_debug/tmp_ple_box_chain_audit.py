from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path("/mnt/e/Pose_create/ple_project_degree5_focused")
PDB = Path("/mnt/e/Pose_create/pdb2r5v.pdb")
HPP = Path("/mnt/e/Pose_create/pdb2r5v_chainA_GLU108_for_Hpp.pdb")
CONFIG = ROOT / "configs" / "focused_box_round1.txt"
OUTDIR = ROOT / "pose_analysis" / "box_chain_audit"


def dist(a, b):
    return math.dist(a, b)


def read_config(path: Path):
    vals = {}
    for line in path.read_text().splitlines():
        if "=" not in line:
            continue
        key, value = [part.strip() for part in line.split("=", 1)]
        try:
            vals[key] = float(value)
        except ValueError:
            pass
    return vals


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    vals = read_config(CONFIG)
    center = (vals["center_x"], vals["center_y"], vals["center_z"])
    size = (vals["size_x"], vals["size_y"], vals["size_z"])

    protein = []
    het = []
    for line in PDB.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = line[21].strip() or "_"
        resn = line[17:20].strip()
        atom = line[12:16].strip()
        elem = (line[76:78].strip() or atom[:1]).upper()
        resi = line[22:26].strip()
        xyz = tuple(float(line[i : i + 8]) for i in (30, 38, 46))
        if line.startswith("ATOM"):
            protein.append((chain, resn, resi, atom, elem, xyz))
        else:
            het.append((chain, resn, resi, atom, elem, xyz))

    rows = []
    for chain in sorted({item[0] for item in protein}):
        pts = [item[5] for item in protein if item[0] == chain]
        centroid = tuple(sum(point[i] for point in pts) / len(pts) for i in range(3))
        rows.append(
            {
                "object": f"chain_{chain}_protein",
                "chain": chain,
                "resname": "",
                "resseq": "",
                "atom": "",
                "distance_to_box_center": dist(center, centroid),
                "nearest_atom_to_box_center": min(dist(center, point) for point in pts),
                "x": centroid[0],
                "y": centroid[1],
                "z": centroid[2],
            }
        )

    for chain, resn, resi, atom, elem, xyz in het:
        if resn in {"FE", "HHH", "PO4"} or elem == "FE":
            rows.append(
                {
                    "object": f"{chain}_{resn}_{resi}_{atom}",
                    "chain": chain,
                    "resname": resn,
                    "resseq": resi,
                    "atom": atom,
                    "distance_to_box_center": dist(center, xyz),
                    "nearest_atom_to_box_center": "",
                    "x": xyz[0],
                    "y": xyz[1],
                    "z": xyz[2],
                }
            )
    rows.sort(key=lambda row: float(row["distance_to_box_center"]))

    with (OUTDIR / "box_chain_distance_audit.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "object",
                "chain",
                "resname",
                "resseq",
                "atom",
                "distance_to_box_center",
                "nearest_atom_to_box_center",
                "x",
                "y",
                "z",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    pose_rows = []
    with (ROOT / "dock_history" / "dock_history.csv").open(errors="replace") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("reason"):
                continue
            try:
                affinity = float(row.get("affinity_kcal_mol") or "nan")
            except ValueError:
                continue
            pose_path = Path(row.get("pose_path", ""))
            if pose_path.exists():
                pose_rows.append((affinity, row.get("seq_id"), pose_path))
    pose_rows.sort(key=lambda item: item[0])

    cx, cy, cz = center
    sx, sy, sz = size
    hx, hy, hz = sx / 2, sy / 2, sz / 2
    corners = [
        (cx + dx * hx, cy + dy * hy, cz + dz * hz)
        for dx in (-1, 1)
        for dy in (-1, 1)
        for dz in (-1, 1)
    ]
    edges = []
    for index_a, a in enumerate(corners):
        for index_b, b in enumerate(corners):
            if index_b <= index_a:
                continue
            if sum(1 for k in range(3) if abs(a[k] - b[k]) > 1e-6) == 1:
                edges.append((a, b))

    lines = [
        "delete all",
        "set mouse_selection_mode, 1",
        "set retain_order, 1",
        "set pdb_retain_ids, on",
        "bg_color white",
        f"load {PDB}, receptor_full",
        f"load {HPP}, hpp_chainA_repaired",
        "hide everything, all",
        "show cartoon, receptor_full",
        "color gray80, receptor_full",
        "color marine, receptor_full and chain A",
        "color wheat, receptor_full and chain B",
        "set cartoon_transparency, 0.35, receptor_full",
        "show cartoon, hpp_chainA_repaired",
        "color cyan, hpp_chainA_repaired",
        "set cartoon_transparency, 0.75, hpp_chainA_repaired",
        "show spheres, receptor_full and (resn FE or elem FE)",
        "color orange, receptor_full and (resn FE or elem FE)",
        "show sticks, receptor_full and resn HHH+PO4",
        "color green, receptor_full and resn HHH+PO4",
        f"pseudoatom focused_box_center, pos=[{cx:.3f},{cy:.3f},{cz:.3f}]",
        "show spheres, focused_box_center",
        "color magenta, focused_box_center",
        "set sphere_scale, 0.45, focused_box_center",
        "from pymol.cgo import CYLINDER",
        "from pymol import cmd",
        "box=[]",
    ]
    for a, b in edges:
        lines.append(
            "box.extend([CYLINDER,"
            f"{a[0]:.3f},{a[1]:.3f},{a[2]:.3f},"
            f"{b[0]:.3f},{b[1]:.3f},{b[2]:.3f},"
            "0.12,1.0,0.0,1.0,1.0,0.0,1.0])"
        )
    lines.extend(["cmd.load_cgo(box, 'focused_box_21A')", "set cgo_transparency, 0.15, focused_box_21A"])

    for rank, (affinity, seq_id, pose_path) in enumerate(pose_rows[:12], start=1):
        obj = f"pose_{rank:02d}_{seq_id}"
        lines.extend(
            [
                f"load {pose_path}, {obj}",
                f"show sticks, {obj}",
                f"color yellow, {obj} and elem C",
                f"color nitrogen, {obj} and elem N",
                f"color oxygen, {obj} and elem O",
                f"color sulfur, {obj} and elem S",
            ]
        )
    lines.extend(
        [
            "select box_contacts_A, receptor_full and chain A and byres (receptor_full within 6 of focused_box_center)",
            "select box_contacts_B, receptor_full and chain B and byres (receptor_full within 6 of focused_box_center)",
            "show sticks, box_contacts_A or box_contacts_B",
            "zoom focused_box_center or pose_01_S000006 or receptor_full and (resn FE or elem FE), 18",
            "center focused_box_center",
            "deselect",
            f"print('focused box center {center}, size {size}; see box_chain_distance_audit.csv')",
        ]
    )
    (OUTDIR / "focused_box_chain_pose_overview.pml").write_text("\n".join(lines) + "\n")

    print("wrote", OUTDIR / "box_chain_distance_audit.csv")
    print("wrote", OUTDIR / "focused_box_chain_pose_overview.pml")
    for row in rows[:12]:
        print(row["object"], round(float(row["distance_to_box_center"]), 3), "nearest", row["nearest_atom_to_box_center"])


if __name__ == "__main__":
    main()
