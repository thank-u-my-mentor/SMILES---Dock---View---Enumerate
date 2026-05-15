#!/usr/bin/env python
"""Infer binding-score-friendly Vina pose families and export pocket datasets."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

try:
    import pymol
    from pymol import cmd as pymol_cmd
except Exception:  # pragma: no cover - optional runtime dependency
    pymol = None
    pymol_cmd = None


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
HBOND_ELEMENTS = {"N", "O", "S"}
HYDROPHOBIC_ELEMENTS = {"C", "Cl", "Br", "I", "F"}
AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def safe_float(value: object) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except Exception:
        return None


def parse_vina_config(path: Path | None) -> dict[str, float]:
    if path is None:
        return {}
    values: dict[str, float] = {}
    for line in path.expanduser().read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = [part.strip() for part in line.split("=", 1)]
        try:
            values[key] = float(value)
        except ValueError:
            continue
    return values


VINA_BOX_KEYS = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")


def has_vina_box(box: dict[str, float]) -> bool:
    return all(key in box for key in VINA_BOX_KEYS)


def clamp_vina_box(box: dict[str, float], min_size: float) -> dict[str, float]:
    clamped = {key: float(box[key]) for key in VINA_BOX_KEYS if key in box}
    for key in ("size_x", "size_y", "size_z"):
        if key in clamped:
            clamped[key] = max(float(min_size), clamped[key])
    return clamped


def vina_box_signature(box: dict[str, float]) -> tuple[float, ...]:
    return tuple(round(float(box[key]), 6) for key in VINA_BOX_KEYS)


def grid_box_objective(correlation: float, target: float) -> tuple[float, float]:
    # The first optimization target is positive correlation:
    # official binding score should rise as predicted binding-surface coverage rises.
    objective = correlation
    loss = max(0.0, target - correlation)
    return objective, loss


def pearson_correlation(pairs: list[tuple[float, float]]) -> float:
    if len(pairs) < 3:
        return 0.0
    xs = [item[0] for item in pairs]
    ys = [item[1] for item in pairs]
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in pairs)
    x_var = sum((x - x_mean) ** 2 for x in xs)
    y_var = sum((y - y_mean) ** 2 for y in ys)
    if x_var <= 0.0 or y_var <= 0.0:
        return 0.0
    return numerator / math.sqrt(x_var * y_var)


def rank_values(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        average_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[order[k]] = average_rank
        i = j
    return ranks


def spearman_correlation(pairs: list[tuple[float, float]]) -> float:
    if len(pairs) < 3:
        return 0.0
    xs = rank_values([item[0] for item in pairs])
    ys = rank_values([item[1] for item in pairs])
    return pearson_correlation(list(zip(xs, ys)))


def vina_box_neighbors(
    box: dict[str, float],
    center_step: float,
    size_step: float,
    min_size: float,
    base_box: dict[str, float] | None = None,
    max_center_shift: float | None = None,
    max_size: float | None = None,
) -> list[tuple[str, dict[str, float]]]:
    neighbors: list[tuple[str, dict[str, float]]] = [("current", dict(box))]
    for axis in ("x", "y", "z"):
        for direction in (-1.0, 1.0):
            candidate = dict(box)
            candidate[f"center_{axis}"] += direction * center_step
            neighbors.append((f"center_{axis}{direction:+g}", clamp_vina_box_bounds(candidate, base_box, min_size, max_center_shift, max_size)))
    for axis in ("x", "y", "z"):
        for direction in (-1.0, 1.0):
            candidate = dict(box)
            key = f"size_{axis}"
            candidate[key] = max(min_size, candidate[key] + direction * size_step)
            neighbors.append((f"{key}{direction:+g}", clamp_vina_box_bounds(candidate, base_box, min_size, max_center_shift, max_size)))
    unique: list[tuple[str, dict[str, float]]] = []
    seen: set[tuple[float, ...]] = set()
    for label, candidate in neighbors:
        signature = tuple(round(candidate[key], 6) for key in VINA_BOX_KEYS)
        if signature in seen:
            continue
        seen.add(signature)
        unique.append((label, candidate))
    return unique


def clamp_vina_box_bounds(
    box: dict[str, float],
    base_box: dict[str, float] | None,
    min_size: float,
    max_center_shift: float | None,
    max_size: float | None,
) -> dict[str, float]:
    candidate = clamp_vina_box(box, min_size)
    if max_size is not None:
        for key in ("size_x", "size_y", "size_z"):
            candidate[key] = min(max_size, candidate[key])
    if base_box is not None and max_center_shift is not None:
        for axis in ("x", "y", "z"):
            key = f"center_{axis}"
            low = base_box[key] - max_center_shift
            high = base_box[key] + max_center_shift
            candidate[key] = min(max(candidate[key], low), high)
    return candidate


def vina_box_starting_points(
    base_box: dict[str, float],
    restarts: int,
    center_step: float,
    size_step: float,
    min_size: float,
) -> list[tuple[str, dict[str, float]]]:
    starts: list[tuple[str, dict[str, float]]] = [("base", clamp_vina_box(base_box, min_size))]
    for axis in ("x", "y", "z"):
        for direction in (-1.0, 1.0):
            candidate = dict(starts[0][1])
            candidate[f"center_{axis}"] += direction * center_step
            starts.append((f"restart_center_{axis}{direction:+g}", candidate))
    for axis in ("x", "y", "z"):
        for direction in (-1.0, 1.0):
            candidate = dict(starts[0][1])
            key = f"size_{axis}"
            candidate[key] = max(min_size, candidate[key] + direction * size_step)
            starts.append((f"restart_{key}{direction:+g}", candidate))

    unique: list[tuple[str, dict[str, float]]] = []
    seen: set[tuple[float, ...]] = set()
    for label, candidate in starts:
        signature = vina_box_signature(candidate)
        if signature in seen:
            continue
        seen.add(signature)
        unique.append((label, candidate))
    return unique[: max(1, restarts + 1)]


def receptor_contacts_for_ligand(
    ligand_atoms: list[dict[str, object]],
    receptor_atoms: list[dict[str, object]],
    cutoff: float = 4.0,
) -> dict[str, set]:
    residues: set[str] = set()
    atom_indices: set[int] = set()
    for ligand_atom in ligand_atoms:
        for receptor_index, receptor_atom in enumerate(receptor_atoms):
            if distance(ligand_atom, receptor_atom) <= cutoff:
                residues.add(residue_label(receptor_atom))
                atom_indices.add(receptor_index)
    return {"residues": residues, "atom_indices": atom_indices}


def grid_surface_atom_indices(
    receptor_atoms: list[dict[str, object]],
    surface_labels: set[str],
    box: dict[str, float],
    padding: float,
    include_buried: bool,
) -> set[int]:
    if not box:
        return {
            index
            for index, atom in enumerate(receptor_atoms)
            if include_buried or residue_label(atom) in surface_labels
        }
    selected: set[int] = set()
    for index, atom in enumerate(receptor_atoms):
        label = residue_label(atom)
        if not include_buried and label not in surface_labels:
            continue
        if atom_in_box(atom, box, padding):
            selected.add(index)
    return selected


def evaluate_vina_box_surface_correlation(
    *,
    box: dict[str, float],
    receptor_atoms: list[dict[str, object]],
    receptor_residues: dict[str, list[dict[str, object]]],
    surface_labels: set[str],
    scored_contact_rows: list[tuple[dict[str, str], list[dict[str, set]]]],
    padding: float,
    include_buried: bool,
    pose_strategy: str,
    metric: str,
) -> tuple[float, float, int, int]:
    binding_surface_labels = grid_surface_residues(receptor_residues, surface_labels, box, padding, include_buried)
    binding_surface_atoms = grid_surface_atom_indices(receptor_atoms, surface_labels, box, padding, include_buried)
    binding_surface_total = len(binding_surface_atoms) if metric == "atom" else len(binding_surface_labels)
    pairs: list[tuple[float, float]] = []
    for row, mode_contacts in scored_contact_rows:
        score = safe_float(row.get("official_binding_score"))
        if score is None or not mode_contacts or not binding_surface_total:
            continue
        coverages: list[float] = []
        selected_contacts = mode_contacts if pose_strategy == "best-coverage" else mode_contacts[:1]
        for contacts in selected_contacts:
            if metric == "atom":
                coverages.append(len(contacts["atom_indices"] & binding_surface_atoms) / binding_surface_total)
            else:
                coverages.append(len(contacts["residues"] & binding_surface_labels) / binding_surface_total)
        if coverages:
            pairs.append((score, max(coverages)))
    return pearson_correlation(pairs), spearman_correlation(pairs), len(pairs), len(binding_surface_labels)


def tune_vina_grid_box(
    *,
    base_box: dict[str, float],
    receptor_atoms: list[dict[str, object]],
    receptor_residues: dict[str, list[dict[str, object]]],
    surface_labels: set[str],
    scored_contact_rows: list[tuple[dict[str, str], list[dict[str, set]]]],
    padding: float,
    include_buried: bool,
    iterations: int,
    center_step: float,
    size_step: float,
    min_size: float,
    max_size: float,
    max_center_shift: float,
    restarts: int,
    max_steps: int,
    target_pearson: float,
    pose_strategy: str,
    metric: str,
) -> tuple[dict[str, float], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    eval_cache: dict[tuple[float, ...], tuple[float, float, int, int]] = {}
    base_box = clamp_vina_box(base_box, min_size)

    def evaluate_cached(candidate: dict[str, float]) -> tuple[float, float, int, int]:
        candidate = clamp_vina_box_bounds(candidate, base_box, min_size, max_center_shift, max_size)
        signature = vina_box_signature(candidate)
        if signature not in eval_cache:
            eval_cache[signature] = evaluate_vina_box_surface_correlation(
                box=candidate,
                receptor_atoms=receptor_atoms,
                receptor_residues=receptor_residues,
                surface_labels=surface_labels,
                scored_contact_rows=scored_contact_rows,
                padding=padding,
                include_buried=include_buried,
                pose_strategy=pose_strategy,
                metric=metric,
            )
        return eval_cache[signature]

    global_best_box = dict(base_box)
    global_best_correlation, global_best_spearman, global_best_pair_count, global_best_surface_count = evaluate_cached(global_best_box)
    global_best_start_id = 0
    starts = vina_box_starting_points(base_box, restarts, center_step, size_step, min_size)
    frontier: list[tuple[float, int, str, int, dict[str, float]]] = []
    expanded: set[tuple[float, ...]] = set()
    step = 0
    for start_id, (start_label, start_box) in enumerate(starts):
        candidate = clamp_vina_box_bounds(start_box, base_box, min_size, max_center_shift, max_size)
        correlation, spearman, pair_count, surface_count = evaluate_cached(candidate)
        objective, loss = grid_box_objective(correlation, target_pearson)
        rows.append(
            {
                "step": step,
                "start_id": start_id,
                "start_label": start_label,
                "iteration": 0,
                "move": "start",
                "metric": metric,
                "pearson_r": round(correlation, 6),
                "spearman_r": round(spearman, 6),
                "objective_score": round(objective, 6),
                "loss": round(loss, 6),
                "pair_count": pair_count,
                "binding_surface_residue_count": surface_count,
                **{key: round(candidate[key], 4) for key in VINA_BOX_KEYS},
            }
        )
        step += 1
        frontier.append((correlation, start_id, start_label, 0, candidate))
        if correlation > global_best_correlation + 1e-9 or (
            abs(correlation - global_best_correlation) <= 1e-9 and surface_count < global_best_surface_count
        ):
            global_best_box = dict(candidate)
            global_best_correlation = correlation
            global_best_spearman = spearman
            global_best_pair_count = pair_count
            global_best_surface_count = surface_count
            global_best_start_id = start_id

    while frontier and step < max_steps and global_best_correlation < target_pearson:
        frontier.sort(key=lambda item: item[0], reverse=True)
        _, start_id, start_label, depth, current = frontier.pop(0)
        signature = vina_box_signature(current)
        if signature in expanded or depth >= iterations:
            continue
        expanded.add(signature)
        for move, candidate in vina_box_neighbors(
            current,
            center_step,
            size_step,
            min_size,
            base_box=base_box,
            max_center_shift=max_center_shift,
            max_size=max_size,
        ):
            if step >= max_steps:
                break
            correlation, spearman, pair_count, surface_count = evaluate_cached(candidate)
            objective, loss = grid_box_objective(correlation, target_pearson)
            rows.append(
                {
                    "step": step,
                    "start_id": start_id,
                    "start_label": start_label,
                    "iteration": depth + 1,
                    "move": move,
                    "metric": metric,
                    "pearson_r": round(correlation, 6),
                    "spearman_r": round(spearman, 6),
                    "objective_score": round(objective, 6),
                    "loss": round(loss, 6),
                    "pair_count": pair_count,
                    "binding_surface_residue_count": surface_count,
                    **{key: round(candidate[key], 4) for key in VINA_BOX_KEYS},
                }
            )
            step += 1
            candidate_signature = vina_box_signature(candidate)
            if candidate_signature not in expanded:
                frontier.append((correlation, start_id, start_label, depth + 1, dict(candidate)))
            if correlation > global_best_correlation + 1e-9 or (
                abs(correlation - global_best_correlation) <= 1e-9 and surface_count < global_best_surface_count
            ):
                global_best_box = dict(candidate)
                global_best_correlation = correlation
                global_best_spearman = spearman
                global_best_pair_count = pair_count
                global_best_surface_count = surface_count
                global_best_start_id = start_id
    objective, loss = grid_box_objective(global_best_correlation, target_pearson)
    rows.append(
        {
            "step": step,
            "start_id": global_best_start_id,
            "start_label": "selected",
            "iteration": "selected",
            "move": "selected",
            "metric": metric,
            "pearson_r": round(global_best_correlation, 6),
            "spearman_r": round(global_best_spearman, 6),
            "objective_score": round(objective, 6),
            "loss": round(loss, 6),
            "pair_count": global_best_pair_count,
            "binding_surface_residue_count": global_best_surface_count,
            **{key: round(global_best_box[key], 4) for key in VINA_BOX_KEYS},
        }
    )
    return global_best_box, rows


def atom_element(atom_name: str, atom_type: str = "") -> str:
    source = atom_type or atom_name
    text = re.sub(r"[^A-Za-z]", "", source)
    if not text:
        return ""
    element = text[:2].title()
    if element in {"Cl", "Br"}:
        return element
    return element[0].upper()


def parse_atom_line(line: str, source: str) -> dict[str, object] | None:
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None
    atom_name = line[12:16].strip()
    residue_name = line[17:20].strip()
    chain = line[21].strip() or "A"
    residue_id = line[22:26].strip()
    atom_type = line[77:].strip().split()[0] if len(line) >= 78 and line[77:].strip() else ""
    return {
        "source": source,
        "serial": line[6:11].strip(),
        "atom_name": atom_name,
        "residue_name": residue_name,
        "chain": chain,
        "residue_id": residue_id,
        "x": x,
        "y": y,
        "z": z,
        "atom_type": atom_type,
        "element": atom_element(atom_name, atom_type),
        "line": line,
    }


def residue_label(atom: dict[str, object]) -> str:
    return f"{atom.get('residue_name', '')}:{atom.get('chain', '')}:{atom.get('residue_id', '')}"


def residue_key(atom: dict[str, object]) -> tuple[str, str, str]:
    return (
        str(atom.get("residue_name", "")),
        str(atom.get("chain", "")),
        str(atom.get("residue_id", "")),
    )


def distance(a: dict[str, object], b: dict[str, object]) -> float:
    return math.sqrt(
        (float(a["x"]) - float(b["x"])) ** 2
        + (float(a["y"]) - float(b["y"])) ** 2
        + (float(a["z"]) - float(b["z"])) ** 2
    )


def read_receptor_atoms(path: Path) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atom = parse_atom_line(line, "receptor")
            if atom and atom.get("element") != "H":
                atoms.append(atom)
    return atoms


def receptor_residue_index(receptor_atoms: list[dict[str, object]]) -> dict[str, list[dict[str, object]]]:
    residues: dict[str, list[dict[str, object]]] = defaultdict(list)
    for atom in receptor_atoms:
        label = residue_label(atom)
        if label:
            residues[label].append(atom)
    return residues


def centroid(atoms: list[dict[str, object]]) -> tuple[float, float, float]:
    n = max(len(atoms), 1)
    return (
        sum(float(atom["x"]) for atom in atoms) / n,
        sum(float(atom["y"]) for atom in atoms) / n,
        sum(float(atom["z"]) for atom in atoms) / n,
    )


def point_distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def surface_residue_table(
    receptor_atoms: list[dict[str, object]],
    surface_neighbor_cutoff: float = 10.0,
) -> tuple[list[dict[str, object]], set[str], dict[str, float]]:
    residues = receptor_residue_index(receptor_atoms)
    centers = {label: centroid(atoms) for label, atoms in residues.items()}
    rows: list[dict[str, object]] = []
    raw_scores: dict[str, float] = {}
    labels = sorted(residues)
    for label in labels:
        center = centers[label]
        neighbor_count = sum(
            1
            for other in labels
            if other != label and point_distance(center, centers[other]) <= surface_neighbor_cutoff
        )
        atom_count = len(residues[label])
        # Lower local residue density implies more exposed/surface-like.
        surface_score = 1.0 / (1.0 + neighbor_count)
        raw_scores[label] = surface_score
        resn, chain, resi = label.split(":", 2)
        rows.append(
            {
                "residue": label,
                "residue_name": resn,
                "chain": chain,
                "residue_id": resi,
                "atom_count": atom_count,
                "neighbor_count_10a": neighbor_count,
                "surface_score": round(surface_score, 6),
            }
        )
    values = sorted(raw_scores.values(), reverse=True)
    if not values:
        return rows, set(), raw_scores
    threshold_index = max(0, int(len(values) * 0.30) - 1)
    threshold = values[threshold_index]
    surface_labels = {label for label, value in raw_scores.items() if value >= threshold}
    for row in rows:
        row["surface_like"] = 1 if row["residue"] in surface_labels else 0
    return rows, surface_labels, raw_scores


def atom_in_box(atom: dict[str, object], box: dict[str, float], padding: float = 0.0) -> bool:
    required = {"center_x", "center_y", "center_z", "size_x", "size_y", "size_z"}
    if not required <= set(box):
        return False
    return (
        abs(float(atom["x"]) - box["center_x"]) <= box["size_x"] / 2.0 + padding
        and abs(float(atom["y"]) - box["center_y"]) <= box["size_y"] / 2.0 + padding
        and abs(float(atom["z"]) - box["center_z"]) <= box["size_z"] / 2.0 + padding
    )


def grid_surface_residues(
    receptor_residues: dict[str, list[dict[str, object]]],
    surface_labels: set[str],
    box: dict[str, float],
    padding: float,
    include_buried: bool,
) -> set[str]:
    if not box:
        return set(surface_labels)
    selected: set[str] = set()
    for label, atoms in receptor_residues.items():
        if not include_buried and label not in surface_labels:
            continue
        if any(atom_in_box(atom, box, padding) for atom in atoms):
            selected.add(label)
    return selected


def pymol_surface_residue_table(
    receptor_path: Path,
    fallback_atoms: list[dict[str, object]],
    surface_area_quantile: float,
) -> tuple[list[dict[str, object]], set[str], dict[str, float]]:
    if pymol_cmd is None:
        raise RuntimeError("PyMOL is not available")
    pymol.finish_launching(["pymol", "-cq"])
    pymol_cmd.reinitialize()
    pymol_cmd.load(str(receptor_path), "rec")
    pymol_cmd.remove("hydrogens")
    pymol_cmd.set("dot_solvent", 1)
    pymol_cmd.set("dot_density", 3)
    pymol_cmd.get_area("rec", load_b=1)

    atom_area_rows: list[tuple[str, float]] = []
    space = {"atom_area_rows": atom_area_rows}
    pymol_cmd.iterate(
        "rec and polymer",
        "atom_area_rows.append((resn + ':' + chain + ':' + resi, float(b)))",
        space=space,
    )
    pymol_cmd.reinitialize()
    residue_areas: dict[str, float] = defaultdict(float)
    residue_atom_counts: dict[str, int] = defaultdict(int)
    for label, area in atom_area_rows:
        residue_areas[label] += area
        residue_atom_counts[label] += 1
    if not residue_areas:
        raise RuntimeError("PyMOL returned no residue surface areas")
    values = sorted(residue_areas.values(), reverse=True)
    q = min(max(surface_area_quantile, 0.01), 1.0)
    threshold_index = max(0, min(len(values) - 1, int(len(values) * q) - 1))
    threshold = values[threshold_index]
    surface_labels = {label for label, area in residue_areas.items() if area >= threshold}
    rows: list[dict[str, object]] = []
    for label, area in sorted(residue_areas.items()):
        parts = label.split(":", 2)
        resn, chain, resi = parts if len(parts) == 3 else ("", "", "")
        rows.append(
            {
                "residue": label,
                "residue_name": resn,
                "chain": chain,
                "residue_id": resi,
                "atom_count": residue_atom_counts[label],
                "neighbor_count_10a": "",
                "surface_score": round(area, 6),
                "surface_like": 1 if label in surface_labels else 0,
                "surface_method": "pymol_sasa",
            }
        )
    return rows, surface_labels, dict(residue_areas)


def parse_pdbqt_modes(path: Path) -> list[dict[str, object]]:
    modes: list[dict[str, object]] = []
    if not path.exists():
        return modes
    atoms: list[dict[str, object]] = []
    mode_index = 1
    affinity: float | None = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            if atoms:
                modes.append({"mode": mode_index, "affinity": affinity, "atoms": atoms})
                atoms = []
                affinity = None
            try:
                mode_index = int(line.split()[1])
            except Exception:
                mode_index = len(modes) + 1
        elif line.startswith("REMARK VINA RESULT:"):
            parts = line.split()
            if len(parts) >= 4:
                affinity = safe_float(parts[3])
        elif line.startswith(("ATOM", "HETATM")):
            atom = parse_atom_line(line, "ligand")
            if atom and atom.get("element") != "H":
                atoms.append(atom)
        elif line.startswith("ENDMDL") and atoms:
            modes.append({"mode": mode_index, "affinity": affinity, "atoms": atoms})
            atoms = []
            affinity = None
    if atoms:
        modes.append({"mode": mode_index, "affinity": affinity, "atoms": atoms})
    return modes


def contact_features(
    ligand_atoms: list[dict[str, object]],
    receptor_atoms: list[dict[str, object]],
    cutoff: float,
    surface_labels: set[str] | None = None,
    binding_surface_labels: set[str] | None = None,
    receptor_residue_count: int | None = None,
) -> dict[str, object]:
    surface_labels = surface_labels or set()
    binding_surface_labels = binding_surface_labels or surface_labels
    residues4: set[str] = set()
    residues5: set[str] = set()
    residues6: set[str] = set()
    hbonds: set[str] = set()
    hydrophobic: set[str] = set()
    aromatic: set[str] = set()
    ch_pi: set[str] = set()
    contacted_ligand_atoms: set[int] = set()
    min_distance = float("inf")
    residue_min: dict[str, float] = {}
    for i, ligand_atom in enumerate(ligand_atoms):
        le = str(ligand_atom.get("element", ""))
        atom_type = str(ligand_atom.get("atom_type", ""))
        for receptor_atom in receptor_atoms:
            d = distance(ligand_atom, receptor_atom)
            label = residue_label(receptor_atom)
            relem = str(receptor_atom.get("element", ""))
            resn = str(receptor_atom.get("residue_name", ""))
            min_distance = min(min_distance, d)
            if label:
                previous = residue_min.get(label)
                if previous is None or d < previous:
                    residue_min[label] = d
            if d <= 4.0:
                residues4.add(label)
                contacted_ligand_atoms.add(i)
            if d <= 5.0:
                residues5.add(label)
            if d <= 6.0:
                residues6.add(label)
            if le in HBOND_ELEMENTS and relem in HBOND_ELEMENTS and d <= 3.5:
                hbonds.add(label)
            if le in HYDROPHOBIC_ELEMENTS and relem in HYDROPHOBIC_ELEMENTS and d <= 4.5:
                hydrophobic.add(label)
            if le == "C" and resn in AROMATIC_RESIDUES and d <= 5.0:
                aromatic.add(label)
            if le == "C" and resn in AROMATIC_RESIDUES and d <= 4.2 and not atom_type.startswith("A"):
                ch_pi.add(label)
    closest = sorted(residue_min.items(), key=lambda item: item[1])[:12]
    surface_contact4 = residues4 & surface_labels
    surface_contact5 = residues5 & surface_labels
    binding_surface_contact4 = residues4 & binding_surface_labels
    binding_surface_contact5 = residues5 & binding_surface_labels
    receptor_total = receptor_residue_count or 0
    binding_surface_total = len(binding_surface_labels)
    return {
        "residues4": residues4,
        "residues5": residues5,
        "residues6": residues6,
        "hbonds": hbonds,
        "hydrophobic": hydrophobic,
        "aromatic": aromatic,
        "ch_pi": ch_pi,
        "contact_residue_count_4a": len(residues4),
        "contact_residue_count_5a": len(residues5),
        "contact_residue_count_6a": len(residues6),
        "surface_contact_residue_count_4a": len(surface_contact4),
        "surface_contact_residue_count_5a": len(surface_contact5),
        "surface_contact_fraction_4a": round(len(surface_contact4) / len(residues4), 4) if residues4 else "",
        "binding_surface_contact_count_4a": len(binding_surface_contact4),
        "binding_surface_contact_count_5a": len(binding_surface_contact5),
        "binding_surface_coverage_4a": round(len(binding_surface_contact4) / binding_surface_total, 4) if binding_surface_total else "",
        "binding_surface_coverage_5a": round(len(binding_surface_contact5) / binding_surface_total, 4) if binding_surface_total else "",
        "binding_surface_total_residues": binding_surface_total,
        "receptor_residue_coverage_4a": round(len(residues4) / receptor_total, 6) if receptor_total else "",
        "receptor_residue_coverage_5a": round(len(residues5) / receptor_total, 6) if receptor_total else "",
        "ligand_atom_contact_fraction_4a": round(len(contacted_ligand_atoms) / len(ligand_atoms), 4) if ligand_atoms else "",
        "min_ligand_receptor_distance": round(min_distance, 4) if math.isfinite(min_distance) else "",
        "closest_residues": ";".join(f"{label}:{d:.2f}" for label, d in closest),
    }


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def mode_similarity(a: dict[str, object], b: dict[str, object]) -> float:
    weights = {
        "residues4": 0.35,
        "residues5": 0.2,
        "hbonds": 0.15,
        "hydrophobic": 0.1,
        "aromatic": 0.1,
        "ch_pi": 0.1,
    }
    return round(sum(weight * jaccard(a[key], b[key]) for key, weight in weights.items()), 6)


def best_score_row(rows: list[dict[str, str]]) -> dict[str, str]:
    scored = [row for row in rows if safe_float(row.get("official_binding_score")) is not None]
    if not scored:
        raise ValueError("analysis CSV has no official_binding_score values")
    # In current Task2 score table, larger official score is treated as better.
    return max(scored, key=lambda row: safe_float(row.get("official_binding_score")) or float("-inf"))


def best_score_row_with_pose(rows: list[dict[str, str]]) -> tuple[dict[str, str], list[dict[str, object]]]:
    scored = [row for row in rows if safe_float(row.get("official_binding_score")) is not None]
    scored.sort(key=lambda row: safe_float(row.get("official_binding_score")) or float("-inf"), reverse=True)
    for row in scored:
        modes = parse_pdbqt_modes(Path(row.get("pose_path", "")).expanduser())
        if modes:
            return row, modes
    raise ValueError("no scored row has a readable pose_path")


def pdb_atom_line(atom: dict[str, object], serial: int, record: str = "ATOM") -> str:
    element = str(atom.get("element", ""))[:2].rjust(2)
    atom_name = str(atom.get("atom_name", "X"))[:4].rjust(4)
    resn = str(atom.get("residue_name", "MOL"))[:3].rjust(3)
    chain = str(atom.get("chain", "A") or "A")[:1]
    try:
        resi = int(str(atom.get("residue_id", "1")).strip())
    except Exception:
        resi = 1
    return (
        f"{record:<6}{serial:5d} {atom_name} {resn} {chain}{resi:4d}    "
        f"{float(atom['x']):8.3f}{float(atom['y']):8.3f}{float(atom['z']):8.3f}"
        f"  1.00  0.00          {element}\n"
    )


def export_pocket(
    outdir: Path,
    standard_mode: dict[str, object],
    features: dict[str, object],
    receptor_atoms: list[dict[str, object]],
    cutoff: float,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    ligand_atoms = standard_mode["atoms"]
    pocket_atoms = [
        atom for atom in receptor_atoms
        if any(distance(atom, ligand_atom) <= cutoff for ligand_atom in ligand_atoms)
    ]
    pocket_residues = sorted({residue_label(atom) for atom in pocket_atoms})
    with (outdir / "ligand_mode.pdb").open("w", encoding="utf-8") as handle:
        for serial, atom in enumerate(ligand_atoms, start=1):
            handle.write(pdb_atom_line(atom, serial, "HETATM"))
        handle.write("END\n")
    with (outdir / "pocket_atoms.pdb").open("w", encoding="utf-8") as handle:
        for serial, atom in enumerate(pocket_atoms, start=1):
            handle.write(pdb_atom_line(atom, serial, "ATOM"))
        handle.write("END\n")
    graph = {
        "standard_mode": standard_mode["mode"],
        "standard_affinity": standard_mode.get("affinity"),
        "features": {
            key: sorted(value) if isinstance(value, set) else value
            for key, value in features.items()
        },
        "ligand_atoms": [
            {key: atom[key] for key in ("atom_name", "element", "x", "y", "z")}
            for atom in ligand_atoms
        ],
        "pocket_atoms": [
            {key: atom[key] for key in ("atom_name", "residue_name", "chain", "residue_id", "element", "x", "y", "z")}
            for atom in pocket_atoms
        ],
        "pocket_residues": pocket_residues,
    }
    (outdir / "pocket_graph.json").write_text(json.dumps(graph, indent=2), encoding="utf-8")


def write_binding_surface_pdb(
    path: Path,
    receptor_residues: dict[str, list[dict[str, object]]],
    binding_surface_labels: set[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        serial = 1
        for label in sorted(binding_surface_labels):
            for atom in receptor_residues.get(label, []):
                handle.write(pdb_atom_line(atom, serial, "ATOM"))
                serial += 1
        handle.write("END\n")


def residue_selection(residues: set[str] | list[str]) -> str:
    grouped: dict[tuple[str, str], list[str]] = defaultdict(list)
    for residue in residues:
        parts = str(residue).split(":", 2)
        if len(parts) != 3:
            continue
        resn, chain, resi = parts
        grouped[(chain, resn)].append(resi)
    pieces = []
    for (chain, resn), resis in sorted(grouped.items()):
        resi_text = "+".join(sorted(set(resis), key=lambda item: int(item) if item.isdigit() else item))
        pieces.append(f"(chain {chain} and resn {resn} and resi {resi_text})")
    return " or ".join(pieces) if pieces else "none"


def write_pymol_view(
    *,
    pml_path: Path,
    pse_path: Path,
    receptor_path: Path,
    ligand_pdb_path: Path,
    pocket_pdb_path: Path,
    binding_surface_pdb_path: Path,
    features: dict[str, object],
    title: str,
    save_pse: bool,
) -> None:
    lines = [
        "reinitialize",
        f"load {receptor_path}, receptor",
        f"load {ligand_pdb_path}, ligand_mode",
        f"load {pocket_pdb_path}, pocket_atoms",
        f"load {binding_surface_pdb_path}, predicted_binding_surface",
        "hide everything, all",
        "show cartoon, receptor",
        "color gray80, receptor",
        "set cartoon_transparency, 0.15, receptor",
        "show sticks, ligand_mode",
        "show spheres, ligand_mode",
        "color yellow, ligand_mode and element C",
        "color nitrogen, ligand_mode and element N",
        "color oxygen, ligand_mode and element O",
        "color sulfur, ligand_mode and element S",
        "set sphere_scale, 0.18, ligand_mode",
        "show surface, predicted_binding_surface",
        "color lime, predicted_binding_surface",
        "set transparency, 0.25, predicted_binding_surface",
        "show sticks, predicted_binding_surface and (sidechain or name CA)",
        "set stick_radius, 0.15, predicted_binding_surface",
        "color nitrogen, predicted_binding_surface and element N",
        "color oxygen, predicted_binding_surface and element O",
        "color sulfur, predicted_binding_surface and element S",
        "select contact_4a, byres (receptor within 4.0 of ligand_mode)",
        "show sticks, contact_4a and (sidechain or name CA)",
        "color tv_orange, contact_4a",
        "select predicted_surface_contact_4a, byres (predicted_binding_surface within 4.0 of ligand_mode)",
        "show sticks, predicted_surface_contact_4a and (sidechain or name CA)",
        "color tv_orange, predicted_surface_contact_4a",
        "set field_of_view, 18",
        "zoom ligand_mode or predicted_binding_surface, 10",
        "center ligand_mode",
        f"print('{title}')",
    ]
    if save_pse:
        lines.append(f"save {pse_path}")
    pml_path.parent.mkdir(parents=True, exist_ok=True)
    pml_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_pymol_session(pml_path: Path) -> bool:
    if pymol_cmd is None:
        return False
    try:
        pymol.finish_launching(["pymol", "-cq"])
        pymol_cmd.do(f"run {pml_path}")
        return True
    except Exception as exc:
        print(f"pse_save_failed={pml_path} reason={exc}", flush=True)
        return False


def infer(args: argparse.Namespace) -> None:
    analysis_rows = read_csv(args.analysis_csv.expanduser().resolve())
    receptor_path = args.receptor.expanduser().resolve()
    receptor_atoms = read_receptor_atoms(receptor_path)
    receptor_residues = receptor_residue_index(receptor_atoms)
    try:
        surface_rows, surface_labels, _ = pymol_surface_residue_table(
            receptor_path,
            receptor_atoms,
            args.surface_area_quantile,
        )
        surface_method = "pymol_sasa"
    except Exception as exc:
        print(f"surface_method=fallback_neighbor_density reason={exc}", flush=True)
        surface_rows, surface_labels, _ = surface_residue_table(receptor_atoms, args.surface_neighbor_cutoff)
        for row in surface_rows:
            row["surface_method"] = "neighbor_density"
        surface_method = "neighbor_density"
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    vina_box = parse_vina_config(args.config.expanduser().resolve()) if args.config else {}
    scored_rows = [row for row in analysis_rows if safe_float(row.get("official_binding_score")) is not None and row.get("pose_path")]
    scored_mode_rows: list[tuple[dict[str, str], list[dict[str, object]]]] = [
        (row, parse_pdbqt_modes(Path(row["pose_path"]).expanduser()))
        for row in scored_rows
    ]
    scored_contact_rows: list[tuple[dict[str, str], list[dict[str, set]]]] = []
    for row, modes in scored_mode_rows:
        mode_contacts: list[dict[str, set]] = []
        for mode in modes:
            mode_contacts.append(receptor_contacts_for_ligand(mode["atoms"], receptor_atoms, cutoff=4.0))
        scored_contact_rows.append((row, mode_contacts))
    if args.tune_grid_box and has_vina_box(vina_box):
        vina_box, tuning_rows = tune_vina_grid_box(
            base_box=vina_box,
            receptor_atoms=receptor_atoms,
            receptor_residues=receptor_residues,
            surface_labels=surface_labels,
            scored_contact_rows=scored_contact_rows,
            padding=args.grid_surface_padding,
            include_buried=args.include_buried_grid_residues,
            iterations=args.grid_tune_iterations,
            center_step=args.grid_center_step,
            size_step=args.grid_size_step,
            min_size=args.grid_min_size,
            max_size=args.grid_max_size,
            max_center_shift=args.grid_max_center_shift,
            restarts=args.grid_tune_restarts,
            max_steps=args.grid_tune_max_steps,
            target_pearson=args.grid_target_pearson,
            pose_strategy=args.grid_tune_pose_strategy,
            metric=args.grid_tune_metric,
        )
        write_csv(
            outdir / "grid_box_tuning.csv",
            tuning_rows,
            ["step", "start_id", "start_label", "iteration", "move", "metric", "pearson_r", "spearman_r", "objective_score", "loss", "pair_count", "binding_surface_residue_count", *VINA_BOX_KEYS],
        )
        write_csv(
            outdir / "selected_grid_box.csv",
            [tuning_rows[-1]],
            ["step", "start_id", "start_label", "iteration", "move", "metric", "pearson_r", "spearman_r", "objective_score", "loss", "pair_count", "binding_surface_residue_count", *VINA_BOX_KEYS],
        )
        selected = tuning_rows[-1]
        print(
            "grid_tuning_result="
            f"metric={selected.get('metric')} "
            f"pearson_r={selected.get('pearson_r')} "
            f"spearman_r={selected.get('spearman_r')} "
            f"loss={selected.get('loss')} "
            f"steps={selected.get('step')} "
            f"target_pearson={args.grid_target_pearson} "
            f"surface_residues={selected.get('binding_surface_residue_count')} "
            f"center=({selected.get('center_x')},{selected.get('center_y')},{selected.get('center_z')}) "
            f"size=({selected.get('size_x')},{selected.get('size_y')},{selected.get('size_z')})",
            flush=True,
        )
    elif args.tune_grid_box:
        print("warning: --tune-grid-box requested but no valid --config Vina box was provided; using untuned receptor surface", flush=True)
    binding_surface_labels = grid_surface_residues(
        receptor_residues,
        surface_labels,
        vina_box,
        args.grid_surface_padding,
        include_buried=args.include_buried_grid_residues,
    )
    binding_surface_source = "vina_grid_box" if vina_box else f"full_receptor_{surface_method}"
    binding_surface_pdb_path = outdir / "predicted_binding_surface.pdb"
    write_binding_surface_pdb(binding_surface_pdb_path, receptor_residues, binding_surface_labels)
    write_csv(
        outdir / "receptor_surface_residues.csv",
        surface_rows,
        ["residue", "residue_name", "chain", "residue_id", "atom_count", "neighbor_count_10a", "surface_score", "surface_like", "surface_method"],
    )
    binding_surface_rows: list[dict[str, object]] = []
    for residue in sorted(binding_surface_labels):
        resn, chain, resi = residue.split(":", 2)
        binding_surface_rows.append(
            {
                "residue": residue,
                "residue_name": resn,
                "chain": chain,
                "residue_id": resi,
                "surface_like": 1 if residue in surface_labels else 0,
                "atom_count": len(receptor_residues.get(residue, [])),
                "source": binding_surface_source,
                **{key: round(vina_box[key], 4) for key in VINA_BOX_KEYS if key in vina_box},
                "grid_padding": args.grid_surface_padding if vina_box else "",
            }
        )
    write_csv(
        outdir / "predicted_binding_surface_residues.csv",
        binding_surface_rows,
        ["residue", "residue_name", "chain", "residue_id", "surface_like", "atom_count", "source", *VINA_BOX_KEYS, "grid_padding"],
    )

    standard_row, standard_modes = best_score_row_with_pose(analysis_rows)
    standard_features = {
        int(mode["mode"]): contact_features(
            mode["atoms"],
            receptor_atoms,
            args.pocket_cutoff,
            surface_labels=surface_labels,
            binding_surface_labels=binding_surface_labels,
            receptor_residue_count=len(receptor_residues),
        )
        for mode in standard_modes
    }

    match_rows: list[dict[str, object]] = []
    pocket_residue_rows: list[dict[str, object]] = []
    family_support: Counter[int] = Counter()
    family_weighted_support: defaultdict[int, float] = defaultdict(float)
    residue_support: defaultdict[str, float] = defaultdict(float)

    for row, modes in scored_mode_rows:
        score = safe_float(row.get("official_binding_score")) or 0.0
        for mode in modes:
            features = contact_features(
                mode["atoms"],
                receptor_atoms,
                args.pocket_cutoff,
                surface_labels=surface_labels,
                binding_surface_labels=binding_surface_labels,
                receptor_residue_count=len(receptor_residues),
            )
            best_mode_index = None
            best_similarity = -1.0
            for standard_index, standard_feature in standard_features.items():
                similarity = mode_similarity(features, standard_feature)
                if similarity > best_similarity:
                    best_similarity = similarity
                    best_mode_index = standard_index
            if best_mode_index is None:
                continue
            family_support[best_mode_index] += 1
            family_weighted_support[best_mode_index] += score * max(best_similarity, 0.0)
            for residue in features["residues4"]:
                residue_support[str(residue)] += score
            match_rows.append(
                {
                    "seq_id": row.get("seq_id", ""),
                    "nickname": row.get("nickname", ""),
                    "official_binding_score": row.get("official_binding_score", ""),
                    "ligand_pose_path": row.get("pose_path", ""),
                    "ligand_mode": mode["mode"],
                    "ligand_mode_affinity": mode.get("affinity", ""),
                    "standard_seq_id": standard_row.get("seq_id", ""),
                    "standard_mode": best_mode_index,
                    "mode_similarity": best_similarity,
                    "contact_residue_count_4a": features["contact_residue_count_4a"],
                    "surface_contact_fraction_4a": features["surface_contact_fraction_4a"],
                    "binding_surface_contact_count_4a": features["binding_surface_contact_count_4a"],
                    "binding_surface_coverage_4a": features["binding_surface_coverage_4a"],
                    "binding_surface_total_residues": features["binding_surface_total_residues"],
                    "receptor_residue_coverage_4a": features["receptor_residue_coverage_4a"],
                    "ligand_atom_contact_fraction_4a": features["ligand_atom_contact_fraction_4a"],
                    "min_ligand_receptor_distance": features["min_ligand_receptor_distance"],
                    "closest_residues": features["closest_residues"],
                }
            )

    summary_rows: list[dict[str, object]] = []
    for mode in standard_modes:
        index = int(mode["mode"])
        features = standard_features[index]
        pocket_dir = outdir / "pockets" / f"standard_mode_{index}"
        export_pocket(pocket_dir, mode, features, receptor_atoms, args.pocket_cutoff)
        for residue in sorted(features["residues5"]):
            atoms = receptor_residues.get(residue, [])
            resn, chain, resi = residue.split(":", 2)
            pocket_residue_rows.append(
                {
                    "standard_mode": index,
                    "residue": residue,
                    "residue_name": resn,
                    "chain": chain,
                    "residue_id": resi,
                    "within_4a": 1 if residue in features["residues4"] else 0,
                    "within_5a": 1,
                    "surface_like": 1 if residue in surface_labels else 0,
                    "predicted_binding_surface": 1 if residue in binding_surface_labels else 0,
                    "atom_count": len(atoms),
                }
            )
        summary_rows.append(
            {
                "standard_seq_id": standard_row.get("seq_id", ""),
                "standard_score": standard_row.get("official_binding_score", ""),
                "standard_input_smiles": standard_row.get("input_smiles", ""),
                "standard_canonical_smiles": standard_row.get("canonical_smiles", ""),
                "standard_pose_path": standard_row.get("pose_path", ""),
                "standard_mode": index,
                "standard_mode_affinity": mode.get("affinity", ""),
                "support_count": family_support[index],
                "weighted_support": round(family_weighted_support[index], 6),
                "standard_contact_residue_count_4a": features["contact_residue_count_4a"],
                "standard_contact_residue_count_5a": features["contact_residue_count_5a"],
                "standard_surface_contact_fraction_4a": features["surface_contact_fraction_4a"],
                "standard_binding_surface_contact_count_4a": features["binding_surface_contact_count_4a"],
                "standard_binding_surface_coverage_4a": features["binding_surface_coverage_4a"],
                "standard_binding_surface_total_residues": features["binding_surface_total_residues"],
                "standard_receptor_residue_coverage_4a": features["receptor_residue_coverage_4a"],
                "standard_receptor_residue_coverage_5a": features["receptor_residue_coverage_5a"],
                "standard_ligand_atom_contact_fraction_4a": features["ligand_atom_contact_fraction_4a"],
                "standard_min_ligand_receptor_distance": features["min_ligand_receptor_distance"],
                "standard_closest_residues": features["closest_residues"],
                "pocket_dir": str(pocket_dir),
                "_features": features,
            }
        )
    summary_rows.sort(
        key=lambda row: (
            -float(row.get("standard_binding_surface_coverage_4a") or 0.0),
            -float(row.get("standard_surface_contact_fraction_4a") or 0.0),
            -int(row["support_count"]),
            -float(row["weighted_support"]),
        )
    )
    visualization_dir = outdir / "visualizations"
    visualization_dir.mkdir(parents=True, exist_ok=True)
    for old_view in list(visualization_dir.glob("*.pml")) + list(visualization_dir.glob("*.pse")):
        old_view.unlink()
    for rank, row in enumerate(summary_rows, start=1):
        pocket_dir = Path(str(row["pocket_dir"]))
        features = row.pop("_features")
        surface_fraction = str(row.get("standard_surface_contact_fraction_4a", "")).replace(".", "p")
        binding_coverage = str(row.get("standard_binding_surface_coverage_4a", "")).replace(".", "p")
        view_base = f"rank{rank:02d}_bind{binding_coverage}_surface{surface_fraction}_seq{row['standard_seq_id']}_mode{row['standard_mode']}"
        local_pml = pocket_dir / "view.pml"
        local_pse = pocket_dir / "view.pse"
        ranked_pml = visualization_dir / f"{view_base}.pml"
        ranked_pse = visualization_dir / f"{view_base}.pse"
        title = (
            f"{view_base} score={row['standard_score']} support={row['support_count']} "
            f"pose={row['standard_pose_path']}"
        )
        write_pymol_view(
            pml_path=local_pml,
            pse_path=local_pse,
            receptor_path=receptor_path,
            ligand_pdb_path=pocket_dir / "ligand_mode.pdb",
            pocket_pdb_path=pocket_dir / "pocket_atoms.pdb",
            binding_surface_pdb_path=binding_surface_pdb_path,
            features=features,
            title=title,
            save_pse=args.save_pse,
        )
        write_pymol_view(
            pml_path=ranked_pml,
            pse_path=ranked_pse,
            receptor_path=receptor_path,
            ligand_pdb_path=pocket_dir / "ligand_mode.pdb",
            pocket_pdb_path=pocket_dir / "pocket_atoms.pdb",
            binding_surface_pdb_path=binding_surface_pdb_path,
            features=features,
            title=title,
            save_pse=args.save_pse,
        )
        if args.save_pse:
            save_pymol_session(ranked_pml)
        row["view_pml"] = str(ranked_pml)
        row["view_pse"] = str(ranked_pse) if args.save_pse else ""
    residue_rows = [
        {"residue": residue, "score_weighted_support": round(value, 6)}
        for residue, value in sorted(residue_support.items(), key=lambda item: -item[1])
    ]

    write_csv(
        outdir / "mode_family_summary.csv",
        summary_rows,
        [
            "standard_seq_id",
            "standard_score",
            "standard_input_smiles",
            "standard_canonical_smiles",
            "standard_pose_path",
            "standard_mode",
            "standard_mode_affinity",
            "support_count",
            "weighted_support",
            "standard_contact_residue_count_4a",
            "standard_contact_residue_count_5a",
            "standard_surface_contact_fraction_4a",
            "standard_binding_surface_contact_count_4a",
            "standard_binding_surface_coverage_4a",
            "standard_binding_surface_total_residues",
            "standard_receptor_residue_coverage_4a",
            "standard_receptor_residue_coverage_5a",
            "standard_ligand_atom_contact_fraction_4a",
            "standard_min_ligand_receptor_distance",
            "standard_closest_residues",
            "pocket_dir",
            "view_pml",
            "view_pse",
        ],
    )
    write_csv(
        outdir / "scored_mode_matches.csv",
        match_rows,
        [
            "seq_id",
            "nickname",
            "official_binding_score",
            "ligand_pose_path",
            "ligand_mode",
            "ligand_mode_affinity",
            "standard_seq_id",
            "standard_mode",
            "mode_similarity",
            "contact_residue_count_4a",
            "surface_contact_fraction_4a",
            "binding_surface_contact_count_4a",
            "binding_surface_coverage_4a",
            "binding_surface_total_residues",
            "receptor_residue_coverage_4a",
            "ligand_atom_contact_fraction_4a",
            "min_ligand_receptor_distance",
            "closest_residues",
        ],
    )
    write_csv(
        outdir / "pocket_residues.csv",
        pocket_residue_rows,
        ["standard_mode", "residue", "residue_name", "chain", "residue_id", "within_4a", "within_5a", "surface_like", "predicted_binding_surface", "atom_count"],
    )
    write_csv(outdir / "consensus_residue_contacts.csv", residue_rows, ["residue", "score_weighted_support"])
    print(
        f"standard_seq_id={standard_row.get('seq_id')} standard_score={standard_row.get('official_binding_score')} "
        f"standard_modes={len(standard_modes)} scored_rows={len(scored_rows)} surface_method={surface_method} "
        f"binding_surface_residues={len(binding_surface_labels)} binding_surface_source={binding_surface_source} wrote={outdir}",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-csv", type=Path, required=True)
    parser.add_argument("--history-csv", type=Path)
    parser.add_argument("--receptor", type=Path, required=True)
    parser.add_argument("--config", type=Path, help="Optional Vina config; center/size define the predicted binding surface region")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--pocket-cutoff", type=float, default=5.0)
    parser.add_argument("--surface-neighbor-cutoff", type=float, default=10.0)
    parser.add_argument("--surface-area-quantile", type=float, default=0.35, help="Top fraction of receptor residues by PyMOL SASA treated as surface-like")
    parser.add_argument("--grid-surface-padding", type=float, default=4.0, help="Extra Angstrom padding around the Vina grid when selecting predicted binding-surface residues")
    parser.add_argument("--surface-only-grid-residues", dest="include_buried_grid_residues", action="store_false", help="Restrict predicted binding surface to SASA surface-like residues inside the Vina grid")
    parser.add_argument("--tune-grid-box", action="store_true", help="Locally tune the Vina grid box to maximize correlation between binding-surface coverage and official binding score")
    parser.add_argument("--grid-tune-iterations", type=int, default=12, help="Maximum local-search iterations for --tune-grid-box")
    parser.add_argument("--grid-tune-restarts", type=int, default=6, help="Number of extra deterministic starting points for --tune-grid-box")
    parser.add_argument("--grid-tune-max-steps", type=int, default=1000, help="Maximum candidate boxes to evaluate for --tune-grid-box")
    parser.add_argument("--grid-target-pearson", type=float, default=0.7, help="Stop --tune-grid-box early once this Pearson correlation is reached")
    parser.add_argument("--grid-tune-metric", choices=["atom", "residue"], default="atom", help="Optimize binding-surface coverage by non-hydrogen receptor atoms or by residues")
    parser.add_argument("--grid-center-step", type=float, default=0.5, help="Center step in Angstrom for --tune-grid-box")
    parser.add_argument("--grid-size-step", type=float, default=1.0, help="Size step in Angstrom for --tune-grid-box")
    parser.add_argument("--grid-min-size", type=float, default=18.0, help="Minimum size_x/size_y/size_z allowed during --tune-grid-box")
    parser.add_argument("--grid-max-size", type=float, default=30.0, help="Maximum size_x/size_y/size_z allowed during --tune-grid-box")
    parser.add_argument("--grid-max-center-shift", type=float, default=6.0, help="Maximum center shift in Angstrom from the input config during --tune-grid-box")
    parser.add_argument("--grid-tune-pose-strategy", choices=["mode1", "best-coverage"], default="best-coverage", help="Use only Vina mode 1 or each ligand's best coverage among modes when tuning the box")
    parser.add_argument("--save-pse", action="store_true", help="Also execute generated PML files in PyMOL and save PSE sessions")
    parser.set_defaults(include_buried_grid_residues=True)
    args = parser.parse_args()
    infer(args)


if __name__ == "__main__":
    main()
