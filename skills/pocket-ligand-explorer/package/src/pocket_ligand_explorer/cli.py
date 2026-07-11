#!/usr/bin/env python
"""CLI for project-directory ligand docking and pose-guided box focusing."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from pocket_ligand_explorer.hpp import write_protonation_map
from pocket_ligand_explorer.md_tools import add_md_tool_parsers, fix_pbc_command, pymol_metal_command, rdc_command
from pocket_ligand_explorer.rosetta_loop_repair import add_repair_loops_parser, repair_loops_command


RECEPTOR_SUFFIXES = {".pdb", ".pdbqt"}
DEFAULT_STATE_DIRNAME = "ple_project"
DEFAULT_INITIAL_PADDING = 8.0
DEFAULT_INITIAL_MAX_SIZE = 80.0
DEFAULT_FOCUS_PADDING = 6.0
DEFAULT_FOCUS_MIN_SIZE = 18.0
DEFAULT_FOCUS_TARGET_SIZE = 21.0
DEFAULT_FOCUS_MAX_SIZE = 21.0
DEFAULT_FOCUS_CLUSTER_RADIUS = 18.0
DEFAULT_FOCUS_MIN_CLUSTER_SIZE = 3
DEFAULT_FOCUS_MAX_POSES = 24
DEFAULT_WIDE_SCOUT_BATCH_SIZE = 18
ENGINE_NAMES = ("vina", "vina.exe", "gnina", "gnina.exe", "glina", "glina.exe", "ple-vina")
GNINA_NAMES = ("gnina", "gnina.exe", "glina", "glina.exe")
DEFAULT_AMBER_LOAD_PARAMS = ("frcmod.ions234lm_126_tip3p",)
STANDARD_MD_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH",
    "GLY", "HID", "HIE", "HIP", "HIS", "ILE", "LEU", "LYS", "LYN", "MET", "PHE",
    "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
DEFAULT_MD_KEEP_RESNAMES = {"FE"}
DEFAULT_MD_KEEP_ELEMENTS = {"FE"}


@dataclass(frozen=True)
class LigandInput:
    smiles: str
    nickname: str = ""


@dataclass(frozen=True)
class Box:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float


@dataclass(frozen=True)
class DegreePreset:
    initial_rounds: int
    initial_batch_size: int
    focus_rounds: int
    focus_batch_size: int
    edit_mode: str
    history_parent_count: int = 0
    history_rounds: int = 1
    history_batch_size: int = 0
    history_nickname_suffix: str = "history"


@dataclass(frozen=True)
class PoseRecord:
    seq_id: str
    nickname: str
    smiles: str
    pose_path: Path
    coords: list[tuple[float, float, float]]
    centroid: tuple[float, float, float]
    focus_score: float
    affinity: float | None
    inner_rmsd: float | None
    cnn_pose_score: float | None


@dataclass(frozen=True)
class PdbAtom:
    line: str
    serial: int
    name: str
    resname: str
    chain: str
    resseq: str
    icode: str
    element: str
    x: float
    y: float
    z: float


DEGREE_PRESETS: dict[int, DegreePreset] = {
    1: DegreePreset(1, 24, 1, 48, "add,delete,shrink,replace"),
    2: DegreePreset(1, 40, 1, 60, "add,delete,shrink,replace"),
    3: DegreePreset(1, 40, 1, 60, "add,delete,shrink,replace,drastic", history_parent_count=4, history_batch_size=16, history_nickname_suffix="degree3"),
    4: DegreePreset(2, 50, 2, 75, "add,delete,shrink,replace,drastic", history_parent_count=8, history_rounds=2, history_batch_size=25, history_nickname_suffix="degree4"),
    5: DegreePreset(2, 50, 2, 75, "add,delete,shrink,replace,drastic", history_parent_count=12, history_rounds=2, history_batch_size=35, history_nickname_suffix="degree5"),
}


def clean_name(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip())
    return text.strip("_") or "ligand"


def safe_float(text: object) -> float | None:
    try:
        if text in ("", None):
            return None
        value = float(text)
    except Exception:
        return None
    return value if math.isfinite(value) else None


def degree_preset(degree: int) -> DegreePreset:
    if degree not in DEGREE_PRESETS:
        raise ValueError("--degree must be an integer from 1 to 5")
    return DEGREE_PRESETS[degree]


def effective_value(args: argparse.Namespace, name: str, default: object) -> object:
    value = getattr(args, name, None)
    return default if value is None else value


def run(cmd: list[str], *, dry_run: bool = False) -> None:
    try:
        print("$ " + shlex.join(cmd), flush=True)
    except BrokenPipeError:
        return
    if dry_run:
        return
    subprocess.run(cmd, check=True)


def emit(message: str) -> None:
    try:
        print(message, flush=True)
    except BrokenPipeError:
        raise SystemExit(0)


def project_dir_from_args(value: Path | None) -> Path:
    return (value or Path.cwd()).expanduser().resolve()


def detect_receptor(project_dir: Path, explicit: Path | None) -> Path:
    if explicit is not None:
        path = explicit.expanduser()
        if not path.is_absolute():
            path = project_dir / path
        if not path.exists():
            raise FileNotFoundError(f"--pdb does not exist: {path}")
        return path.resolve()
    candidates = sorted(
        path for path in project_dir.iterdir()
        if path.is_file() and path.suffix.lower() in RECEPTOR_SUFFIXES
    )
    if not candidates:
        raise FileNotFoundError(f"no receptor *.pdb or *.pdbqt found in {project_dir}; pass --pdb")
    if len(candidates) > 1:
        names = ", ".join(path.name for path in candidates)
        raise ValueError(f"multiple receptor files found: {names}; pass --pdb to choose one")
    return candidates[0].resolve()


def parse_atom_line(line: str) -> tuple[float, float, float] | None:
    if not line.startswith(("ATOM", "HETATM")):
        return None
    try:
        return float(line[30:38]), float(line[38:46]), float(line[46:54])
    except Exception:
        return None


def read_atom_coords(path: Path) -> list[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parsed = parse_atom_line(line)
            if parsed is not None:
                coords.append(parsed)
    return coords


def split_option_values(values: list[str] | None) -> list[str]:
    items: list[str] = []
    for value in values or []:
        for part in re.split(r"[,\s]+", value.strip()):
            if part:
                items.append(part)
    return items


def read_pdb_lines(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    if "\n" not in text and "\\n" in text:
        text = text.replace("\\n", "\n")
    return text.splitlines()


def pdb_resname(line: str) -> str:
    return line[17:20].strip().upper() if len(line) >= 20 else ""


def pdb_chain(line: str) -> str:
    return line[21:22].strip() if len(line) >= 22 else ""


def pdb_resseq(line: str) -> str:
    return line[22:26].strip() if len(line) >= 26 else ""


def pdb_atom_name(line: str) -> str:
    return line[12:16].strip().upper() if len(line) >= 16 else ""


def pdb_residue_label(line: str) -> str:
    chain = pdb_chain(line) or "-"
    return f"{chain}:{pdb_resseq(line)}:{pdb_resname(line)}:{pdb_atom_name(line)}"


def pdb_atom_from_line(line: str) -> PdbAtom | None:
    if not line.startswith(("ATOM  ", "HETATM")):
        return None
    try:
        serial = int(line[6:11])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None
    return PdbAtom(
        line=line.rstrip(),
        serial=serial,
        name=pdb_atom_name(line),
        resname=pdb_resname(line),
        chain=pdb_chain(line),
        resseq=pdb_resseq(line),
        icode=line[26:27].strip() if len(line) >= 27 else "",
        element=pdb_element(line),
        x=x,
        y=y,
        z=z,
    )


def pdb_distance(a: PdbAtom, b: PdbAtom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def atom_residue_id(atom: PdbAtom) -> str:
    chain = atom.chain or "-"
    return f"{atom.resname}:{chain}:{atom.resseq}{atom.icode}"


def prepare_receptor_site_copy(args: argparse.Namespace, receptor: Path, state_dir: Path) -> Path:
    drop_resnames = {item.upper() for item in split_option_values(getattr(args, "drop_receptor_resname", []))}
    if not drop_resnames:
        return receptor
    if receptor.suffix.lower() != ".pdb":
        raise ValueError("--drop-receptor-resname requires a receptor PDB before PDBQT conversion")

    suffix = "_".join(sorted(drop_resnames))
    output = state_dir / "receptor" / f"{receptor.stem}_drop_{suffix}.pdb"
    output.parent.mkdir(parents=True, exist_ok=True)
    kept: list[str] = []
    removed_coord = 0
    removed_meta = 0

    for line in read_pdb_lines(receptor):
        record = line[:6].strip()
        line_resname = pdb_resname(line)
        upper_line = line.upper()
        if record in {"ATOM", "HETATM", "ANISOU"} and line_resname in drop_resnames:
            removed_coord += 1
            continue
        if any(name in upper_line for name in drop_resnames) and (
            record in {"HET", "HETNAM", "FORMUL", "LINK", "SITE"} or line.startswith("REMARK")
        ):
            removed_meta += 1
            continue
        kept.append(line.rstrip())

    if not kept:
        raise RuntimeError(f"receptor filtering removed every record from {receptor}")
    output.write_text("\n".join(kept) + "\n", encoding="utf-8")
    remaining = sum(
        1
        for line in kept
        if line.startswith(("ATOM", "HETATM", "ANISOU")) and pdb_resname(line) in drop_resnames
    )
    print(
        f"receptor_drop_resname={','.join(sorted(drop_resnames))} "
        f"removed_coord_records={removed_coord} removed_meta_records={removed_meta} "
        f"remaining_coord_records={remaining} output={output}",
        flush=True,
    )
    return output.resolve()


def site_box_from_args(args: argparse.Namespace, receptor: Path) -> Box | None:
    site_resname = getattr(args, "site_resname", None)
    if not site_resname:
        return None
    target_resname = re.sub(r"[^A-Za-z0-9]", "", site_resname).upper()
    if not target_resname:
        raise ValueError("--site-resname must contain at least one letter or number")
    target_chain = (getattr(args, "site_chain", "") or "").strip()
    target_resseq = (getattr(args, "site_resseq", "") or "").strip()
    target_atom = (getattr(args, "site_atom_name", "") or "").strip().upper()

    matches: list[tuple[str, tuple[float, float, float]]] = []
    residue_keys: set[tuple[str, str, str]] = set()
    for line in read_pdb_lines(receptor):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if pdb_resname(line) != target_resname:
            continue
        if target_chain and pdb_chain(line) != target_chain:
            continue
        if target_resseq and pdb_resseq(line) != target_resseq:
            continue
        if target_atom and pdb_atom_name(line) != target_atom:
            continue
        parsed = parse_atom_line(line)
        if parsed is None:
            continue
        matches.append((line, parsed))
        residue_keys.add((pdb_chain(line), pdb_resseq(line), pdb_resname(line)))

    selector = f"resname={target_resname}"
    if target_chain:
        selector += f" chain={target_chain}"
    if target_resseq:
        selector += f" resseq={target_resseq}"
    if target_atom:
        selector += f" atom={target_atom}"
    if not matches:
        raise ValueError(f"no receptor atom matched --site-{selector}")
    if len(residue_keys) > 1 and not target_resseq:
        examples = ", ".join(sorted(f"{chain or '-'}:{resseq}:{resname}" for chain, resseq, resname in residue_keys)[:8])
        raise ValueError(f"multiple receptor residues matched --site-{selector}: {examples}; pass --site-chain and --site-resseq")

    center = centroid([coord for _, coord in matches])
    side = float(getattr(args, "site_box_size", DEFAULT_FOCUS_TARGET_SIZE) or DEFAULT_FOCUS_TARGET_SIZE)
    if side <= 0:
        raise ValueError("--site-box-size must be positive")
    box = Box(center[0], center[1], center[2], side, side, side)
    labels = sorted({pdb_residue_label(line) for line, _ in matches})
    print(
        f"site_box={selector} atoms={len(matches)} residues={len(residue_keys)} "
        f"center=({box.center_x:.3f},{box.center_y:.3f},{box.center_z:.3f}) "
        f"size=({box.size_x:.1f},{box.size_y:.1f},{box.size_z:.1f}) "
        f"matches={','.join(labels[:6])}",
        flush=True,
    )
    return box


def bounding_box(coords: list[tuple[float, float, float]], padding: float, min_size: float = 0.0, max_size: float | None = None) -> Box:
    if not coords:
        raise ValueError("no atom coordinates found for box generation")
    xs = [item[0] for item in coords]
    ys = [item[1] for item in coords]
    zs = [item[2] for item in coords]

    def side(values: list[float]) -> float:
        size = max(values) - min(values) + 2.0 * padding
        if min_size:
            size = max(size, min_size)
        if max_size is not None:
            size = min(size, max_size)
        return size

    return Box(
        center_x=(min(xs) + max(xs)) / 2.0,
        center_y=(min(ys) + max(ys)) / 2.0,
        center_z=(min(zs) + max(zs)) / 2.0,
        size_x=side(xs),
        size_y=side(ys),
        size_z=side(zs),
    )


def write_vina_config(path: Path, receptor: Path, box: Box, exhaustiveness: int, num_modes: int, cpu: int | None = None) -> None:
    lines = [
        f"receptor = {receptor}",
        f"center_x = {box.center_x:.4f}",
        f"center_y = {box.center_y:.4f}",
        f"center_z = {box.center_z:.4f}",
        f"size_x = {box.size_x:.4f}",
        f"size_y = {box.size_y:.4f}",
        f"size_z = {box.size_z:.4f}",
        f"exhaustiveness = {exhaustiveness}",
        f"num_modes = {num_modes}",
    ]
    if cpu is not None:
        lines.append(f"cpu = {cpu}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_smiles_file(path: Path) -> list[LigandInput]:
    ligands: list[LigandInput] = []
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            parts = text.split()
            smiles = parts[0].strip()
            nickname = parts[1].strip() if len(parts) > 1 else f"ligand_{line_number}"
            ligands.append(LigandInput(smiles=smiles, nickname=nickname))
    return ligands


def parse_csv_file(path: Path, smiles_column: str, nickname_column: str | None) -> list[LigandInput]:
    ligands: list[LigandInput] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if smiles_column not in (reader.fieldnames or []):
            raise ValueError(f"CSV has no --smiles-column {smiles_column!r}")
        for index, row in enumerate(reader, start=1):
            smiles = str(row.get(smiles_column, "")).strip()
            if not smiles:
                continue
            nickname = str(row.get(nickname_column, "")).strip() if nickname_column else ""
            ligands.append(LigandInput(smiles=smiles, nickname=nickname or f"ligand_{index}"))
    return ligands


def parse_rdkit_structure(path: Path) -> list[LigandInput]:
    try:
        from rdkit import Chem
    except Exception as exc:
        raise RuntimeError("RDKit is required for --sdf/--mol input; use --smiles or --smiles-file instead") from exc
    ligands: list[LigandInput] = []
    if path.suffix.lower() == ".sdf":
        supplier = Chem.SDMolSupplier(str(path), removeHs=False)
        for index, mol in enumerate(supplier, start=1):
            if mol is None:
                continue
            smiles = Chem.MolToSmiles(mol, canonical=True)
            nickname = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else f"ligand_{index}"
            ligands.append(LigandInput(smiles=smiles, nickname=nickname))
    else:
        mol = Chem.MolFromMolFile(str(path), removeHs=False)
        if mol is None:
            raise ValueError(f"RDKit could not parse {path}")
        ligands.append(LigandInput(smiles=Chem.MolToSmiles(mol, canonical=True), nickname=path.stem))
    return ligands


def collect_ligands(args: argparse.Namespace, project_dir: Path) -> list[LigandInput]:
    ligands: list[LigandInput] = []
    for smiles in args.smiles or []:
        ligands.append(LigandInput(smiles=smiles.strip(), nickname=args.nickname or "lead"))
    if args.smiles_file:
        path = resolve_project_path(project_dir, args.smiles_file)
        ligands.extend(parse_smiles_file(path))
    if args.csv:
        path = resolve_project_path(project_dir, args.csv)
        ligands.extend(parse_csv_file(path, args.smiles_column, args.nickname_column))
    if args.sdf:
        ligands.extend(parse_rdkit_structure(resolve_project_path(project_dir, args.sdf)))
    if args.mol:
        ligands.extend(parse_rdkit_structure(resolve_project_path(project_dir, args.mol)))
    cleaned = [ligand for ligand in ligands if ligand.smiles]
    if not cleaned:
        raise ValueError("provide --smiles, --smiles-file, --csv, --sdf, or --mol")
    return cleaned


def resolve_project_path(project_dir: Path, path: Path) -> Path:
    expanded = path.expanduser()
    if not expanded.is_absolute():
        expanded = project_dir / expanded
    if not expanded.exists():
        raise FileNotFoundError(expanded)
    return expanded.resolve()


def find_skills_dir(explicit: Path | None) -> Path:
    if explicit is not None:
        return explicit.expanduser().resolve()
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "smiles-to-vina-docking").exists() and (parent / "pose-analyzer").exists():
            return parent
    env = os.environ.get("CODEX_SKILLS_DIR")
    if env:
        return Path(env).expanduser().resolve()
    # Editable install from skills/pocket-ligand-explorer/package usually reaches
    # the skill directory through the source path parents.
    candidate = here.parents[5] if len(here.parents) > 5 else here.parent
    return candidate


def sibling_scripts(skills_dir: Path) -> Path:
    dock_script = skills_dir / "smiles-to-vina-docking" / "scripts" / "dock_smiles.py"
    infer_script = skills_dir / "pose-analyzer" / "scripts" / "infer_binding_mode_families.py"
    report_script = skills_dir / "pose-analyzer" / "scripts" / "build_pose_report.py"
    for path in (dock_script, infer_script, report_script):
        if not path.exists():
            raise FileNotFoundError(f"required sibling script not found: {path}")
    return dock_script


def report_script_path(skills_dir: Path) -> Path:
    return skills_dir / "pose-analyzer" / "scripts" / "build_pose_report.py"


def score_model_script_path(skills_dir: Path) -> Path:
    return skills_dir / "pose-analyzer" / "scripts" / "score_space_model.py"


def score_dashboard_script_path(skills_dir: Path) -> Path:
    return skills_dir / "pose-analyzer" / "scripts" / "build_score_space_dashboard.py"


def dock_ligand_batch(
    *,
    ligands: list[LigandInput],
    dock_script: Path,
    history_dir: Path,
    receptor: Path,
    config: Path,
    engine_exe: str,
    max_rounds: int,
    batch_size: int,
    exhaustiveness: int,
    num_modes: int,
    cpu: int | None,
    edit_mode: str,
    engine_args: list[str],
    dry_run: bool,
    nickname_suffix: str,
) -> None:
    for index, ligand in enumerate(ligands, start=1):
        nickname = clean_name(ligand.nickname or f"ligand_{index}")
        if nickname_suffix:
            nickname = f"{nickname}_{nickname_suffix}"
        cmd = [
            sys.executable,
            str(dock_script),
            "--smiles",
            ligand.smiles,
            "--nickname",
            nickname,
            "--history-dir",
            str(history_dir),
            "--receptor",
            str(receptor),
            "--config",
            str(config),
            "--vina",
            engine_exe,
            "--max-rounds",
            str(max_rounds),
            "--batch-size",
            str(batch_size),
            "--exhaustiveness",
            str(exhaustiveness),
            "--num-modes",
            str(num_modes),
            "--edit-mode",
            edit_mode,
        ]
        if cpu is not None:
            cmd.extend(["--cpu", str(cpu)])
        for engine_arg in engine_args:
            cmd.append(f"--engine-arg={engine_arg}")
        run(cmd, dry_run=dry_run)


def history_pose_paths(history_csv: Path) -> list[Path]:
    if not history_csv.exists():
        return []
    paths: list[Path] = []
    with history_csv.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("reason"):
                continue
            affinity = safe_float(row.get("affinity_kcal_mol"))
            if affinity is None:
                continue
            text = str(row.get("pose_path", "")).strip()
            if not text:
                continue
            path = Path(text).expanduser()
            if path.exists():
                paths.append(path.resolve())
    return paths


def clean_smiles_text(text: object) -> str:
    return re.sub(r"\s+", "", str(text or "").strip())


def row_smiles(row: dict[str, str]) -> str:
    return clean_smiles_text(row.get("smiles") or row.get("canonical_smiles") or row.get("input_smiles") or row.get("mol_smiles") or "")


def history_score(row: dict[str, str]) -> float | None:
    affinity = safe_float(row.get("affinity_kcal_mol"))
    if affinity is None:
        return None
    inner = safe_float(row.get("inner_rmsd"))
    cnn_pose = safe_float(row.get("cnn_pose_score"))
    penalty = 0.08 * (inner if inner is not None else 6.0)
    cnn_reward = 0.8 * (cnn_pose if cnn_pose is not None else 0.0)
    if cnn_pose is None:
        penalty += 0.25
    elif cnn_pose < 0.2:
        penalty += 0.25
    elif cnn_pose < 0.4:
        penalty += 0.10
    return affinity + penalty - cnn_reward


def is_gnina_engine(engine_exe: str) -> bool:
    name = Path(engine_exe).name.lower()
    return name in GNINA_NAMES or "gnina" in name or "glina" in name


def engine_args_with_cnn(engine_args: list[str], cnn_scoring: str | None) -> list[str]:
    args = [item for item in engine_args if item]
    if cnn_scoring and cnn_scoring != "none":
        args.extend(["--cnn_scoring", cnn_scoring])
    return args


def select_history_parents(history_csv: Path, count: int, exclude_smiles: set[str]) -> list[LigandInput]:
    if count <= 0 or not history_csv.exists():
        return []
    rows: list[tuple[float, dict[str, str]]] = []
    seen = {clean_smiles_text(item) for item in exclude_smiles if item}
    with history_csv.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("reason"):
                continue
            smiles = row_smiles(row)
            if not smiles or smiles in seen:
                continue
            score = history_score(row)
            if score is None:
                continue
            rows.append((score, row))
            seen.add(smiles)
    rows.sort(key=lambda item: item[0])
    parents: list[LigandInput] = []
    for _, row in rows[:count]:
        nickname = row.get("nickname") or row.get("seq_id") or "history_parent"
        parents.append(LigandInput(smiles=row_smiles(row), nickname=nickname))
    return parents


def centroid(coords: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    return (
        sum(item[0] for item in coords) / len(coords),
        sum(item[1] for item in coords) / len(coords),
        sum(item[2] for item in coords) / len(coords),
    )


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def first_pose_coords(path: Path) -> list[tuple[float, float, float]]:
    file_coords: list[tuple[float, float, float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("MODEL") and file_coords:
                break
            parsed = parse_atom_line(line)
            if parsed is not None:
                file_coords.append(parsed)
    return file_coords


def focus_score(row: dict[str, str]) -> float | None:
    score = history_score(row)
    if score is None:
        return None
    cnn_pose = safe_float(row.get("cnn_pose_score"))
    inner = safe_float(row.get("inner_rmsd"))
    if cnn_pose is not None and cnn_pose < 0.5:
        score += 0.4
    if inner is not None and inner > 25.0:
        score += 0.4
    return score


def history_pose_records(history_csv: Path, max_poses: int = DEFAULT_FOCUS_MAX_POSES) -> list[PoseRecord]:
    if not history_csv.exists():
        return []
    rows: list[tuple[float, dict[str, str]]] = []
    with history_csv.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("reason"):
                continue
            text = str(row.get("pose_path", "")).strip()
            if not text:
                continue
            score = focus_score(row)
            if score is None:
                continue
            path = Path(text).expanduser()
            if not path.exists():
                continue
            rows.append((score, row))
    rows.sort(key=lambda item: item[0])
    records: list[PoseRecord] = []
    for score, row in rows[:max_poses]:
        path = Path(str(row.get("pose_path", "")).strip()).expanduser()
        coords = first_pose_coords(path)
        if not coords:
            continue
        records.append(
            PoseRecord(
                seq_id=str(row.get("seq_id", "")),
                nickname=str(row.get("nickname", "")),
                smiles=row_smiles(row),
                pose_path=path.resolve(),
                coords=coords,
                centroid=centroid(coords),
                focus_score=score,
                affinity=safe_float(row.get("affinity_kcal_mol")),
                inner_rmsd=safe_float(row.get("inner_rmsd")),
                cnn_pose_score=safe_float(row.get("cnn_pose_score")),
            )
        )
    return records


def select_focus_records(
    records: list[PoseRecord],
    *,
    radius: float,
    min_cluster_size: int,
) -> list[PoseRecord]:
    if not records:
        return []
    best: tuple[float, int, float, list[PoseRecord]] | None = None
    for anchor in records:
        members = [record for record in records if distance(anchor.centroid, record.centroid) <= radius]
        if len(members) < min_cluster_size:
            members = [anchor]
        avg_score = sum(record.focus_score for record in members) / len(members)
        key = (-len(members), 0, avg_score)
        candidate = (key[0], key[1], key[2], members)
        if best is None or candidate[:3] < best[:3]:
            best = candidate
    selected = best[3] if best else records[:1]
    selected.sort(key=lambda record: record.focus_score)
    return selected


def focus_box_from_records(
    records: list[PoseRecord],
    *,
    padding: float,
    min_size: float,
    target_size: float,
    max_size: float,
) -> Box:
    coords = [coord for record in records for coord in record.coords]
    raw = bounding_box(coords, padding=padding, min_size=min_size)
    center = centroid([record.centroid for record in records])
    sizes = [raw.size_x, raw.size_y, raw.size_z]
    if max(sizes) > max_size:
        side = min(max_size, max(min_size, target_size))
        return Box(
            center_x=center[0],
            center_y=center[1],
            center_z=center[2],
            size_x=side,
            size_y=side,
            size_z=side,
        )

    def side(value: float) -> float:
        return min(max_size, max(min_size, target_size, value))

    return Box(
        center_x=center[0],
        center_y=center[1],
        center_z=center[2],
        size_x=side(raw.size_x),
        size_y=side(raw.size_y),
        size_z=side(raw.size_z),
    )


def pose_coords(paths: list[Path], max_poses: int = 200) -> list[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    for path in paths[:max_poses]:
        file_coords: list[tuple[float, float, float]] = []
        file_coords = first_pose_coords(path)
        coords.extend(file_coords)
    return coords


def box_is_scattered(box: Box, initial: Box, scatter_fraction: float) -> bool:
    focused_volume = max(box.size_x * box.size_y * box.size_z, 1.0)
    initial_volume = max(initial.size_x * initial.size_y * initial.size_z, 1.0)
    return focused_volume / initial_volume > scatter_fraction


def choose_focused_box(history_csv: Path, args: argparse.Namespace) -> tuple[Box | None, int, int]:
    focus_records = history_pose_records(history_csv, max_poses=args.focus_max_poses)
    selected_records = select_focus_records(
        focus_records,
        radius=args.focus_cluster_radius,
        min_cluster_size=args.focus_min_cluster_size,
    )
    if not selected_records:
        return None, 0, len(focus_records)
    focused_box = focus_box_from_records(
        selected_records,
        padding=args.focus_padding,
        min_size=args.focus_min_size,
        target_size=args.focus_target_size,
        max_size=args.focus_max_size,
    )
    return focused_box, len(selected_records), len(focus_records)


def run_pose_analysis(
    *,
    skills_dir: Path,
    history_dir: Path,
    receptor: Path,
    config: Path,
    state_dir: Path,
    save_pse: bool,
    dry_run: bool,
) -> None:
    infer_script = skills_dir / "pose-analyzer" / "scripts" / "infer_binding_mode_families.py"
    report_script = report_script_path(skills_dir)
    pose_analysis_dir = state_dir / "pose_analysis"
    report_dir = state_dir / "pose_report"
    history_csv = history_dir / "dock_history.csv"
    infer_cmd = [
        sys.executable,
        str(infer_script),
        "--analysis-csv",
        str(history_csv),
        "--history-csv",
        str(history_csv),
        "--receptor",
        str(receptor),
        "--config",
        str(config),
        "--outdir",
        str(pose_analysis_dir),
        "--score-column",
        "affinity_kcal_mol",
        "--score-direction",
        "lower-is-better",
    ]
    if save_pse:
        infer_cmd.append("--save-pse")
    run(infer_cmd, dry_run=dry_run)
    report_cmd = [
        sys.executable,
        str(report_script),
        "--history-csv",
        str(history_csv),
        "--pose-analysis-dir",
        str(pose_analysis_dir),
        "--outdir",
        str(report_dir),
    ]
    run(report_cmd, dry_run=dry_run)
    print(f"report={report_dir / 'index.html'}", flush=True)


def run_score_space_dashboard(
    *,
    skills_dir: Path,
    history_dir: Path,
    receptor: Path,
    state_dir: Path,
    dry_run: bool,
    train_size: int,
    test_size: int,
) -> None:
    model_script = score_model_script_path(skills_dir)
    dashboard_script = score_dashboard_script_path(skills_dir)
    for path in (model_script, dashboard_script):
        if not path.exists():
            raise FileNotFoundError(f"required score-space script not found: {path}")
    history_csv = history_dir / "dock_history.csv"
    model_dir = state_dir / "score_space_model"
    dashboard_dir = state_dir / "score_space_dashboard"
    model_cmd = [
        sys.executable,
        str(model_script),
        "--history-csv",
        str(history_csv),
        "--receptor",
        str(receptor),
        "--outdir",
        str(model_dir),
        "--analysis-mode",
        "unsupervised",
        "--score-column",
        "affinity_kcal_mol",
        "--score-direction",
        "lower-is-better",
        "--score-label",
        "PLE docking analysis score",
        "--train-size",
        str(train_size),
        "--test-size",
        str(test_size),
    ]
    run(model_cmd, dry_run=dry_run)
    dashboard_cmd = [
        sys.executable,
        str(dashboard_script),
        "--feature-matrix",
        str(model_dir / "binding_score_feature_matrix.csv"),
        "--metrics-json",
        str(model_dir / "binding_score_model_metrics.json"),
        "--outdir",
        str(dashboard_dir),
    ]
    run(dashboard_cmd, dry_run=dry_run)
    print(f"score_space_dashboard={dashboard_dir / 'index.html'}", flush=True)


def check_engine(engine_exe: str) -> None:
    path = Path(engine_exe).expanduser()
    if path.exists() and not path.is_dir():
        return
    if shutil.which(engine_exe):
        return
    raise FileNotFoundError(f"docking engine not found: {engine_exe}; pass --engine-exe /path/to/vina")


def maybe_existing_executable(path: Path) -> str | None:
    expanded = path.expanduser()
    if expanded.exists() and not expanded.is_dir():
        return str(expanded.resolve())
    return None


def engine_search_paths(project_dir: Path, state_dir: Path) -> list[Path]:
    paths: list[Path] = []
    env_paths = [
        os.environ.get("PLE_VINA"),
        os.environ.get("VINA_EXE"),
        os.environ.get("VINA"),
        os.environ.get("GNINA"),
    ]
    paths.extend(Path(value).expanduser() for value in env_paths if value)
    home = Path.home()
    prefixes = [
        project_dir,
        state_dir,
        project_dir / "vina_bin",
        project_dir / "bin",
        state_dir / "vina_bin",
        home / "vina_task2" / "vina_bin",
        home / "bin",
    ]
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        conda = Path(conda_prefix)
        prefixes.extend([conda / "bin", conda / "Scripts", conda])
    for prefix in prefixes:
        for name in ENGINE_NAMES:
            paths.append(prefix / name)
    return paths


def prefer_gnina_engine(project_dir: Path, state_dir: Path) -> str | None:
    conda_env_name = Path(os.environ.get("CONDA_PREFIX", "")).name.lower()
    should_prefer = conda_env_name in {"gnina", "glina"} or os.environ.get("PLE_ENGINE", "").lower() in {"gnina", "glina"}
    if not should_prefer:
        return None
    for key in ("GNINA", "GLINA"):
        value = os.environ.get(key)
        if not value:
            continue
        existing = maybe_existing_executable(Path(value))
        if existing:
            return existing
        found = shutil.which(value)
        if found:
            return found
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        conda = Path(conda_prefix)
        for prefix in (conda / "bin", conda / "Scripts", conda):
            for name in GNINA_NAMES:
                existing = maybe_existing_executable(prefix / name)
                if existing:
                    return existing
    for name in GNINA_NAMES:
        found = shutil.which(name)
        if found:
            return found
    for path in engine_search_paths(project_dir, state_dir):
        if path.name.lower() in GNINA_NAMES:
            existing = maybe_existing_executable(path)
            if existing:
                return existing
    return None


def resolve_engine_exe(value: str, project_dir: Path, state_dir: Path) -> str:
    if value and value != "auto":
        explicit = maybe_existing_executable(Path(value))
        if explicit:
            return explicit
        found = shutil.which(value)
        if found:
            return found
        return value
    preferred_gnina = prefer_gnina_engine(project_dir, state_dir)
    if preferred_gnina:
        return preferred_gnina
    env_values = [
        os.environ.get("PLE_VINA"),
        os.environ.get("VINA_EXE"),
        os.environ.get("VINA"),
        os.environ.get("GNINA"),
        os.environ.get("GLINA"),
    ]
    for env_value in env_values:
        if not env_value:
            continue
        existing = maybe_existing_executable(Path(env_value))
        if existing:
            return existing
        found = shutil.which(env_value)
        if found:
            return found
    preferred = [
        project_dir / "vina_bin" / "vina",
        project_dir / "vina_bin" / "vina.exe",
        project_dir / "bin" / "vina",
        project_dir / "bin" / "vina.exe",
        state_dir / "vina_bin" / "vina",
        state_dir / "vina_bin" / "vina.exe",
        Path.home() / "vina_task2" / "vina_bin" / "vina",
        Path.home() / "vina_task2" / "vina_bin" / "vina.exe",
    ]
    for path in preferred:
        existing = maybe_existing_executable(path)
        if existing:
            return existing
    for name in ("vina", "vina.exe"):
        found = shutil.which(name)
        if found:
            return found
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        conda = Path(conda_prefix)
        for prefix in (conda / "bin", conda / "Scripts", conda):
            for name in ("vina", "vina.exe"):
                existing = maybe_existing_executable(prefix / name)
                if existing:
                    return existing
    for name in ENGINE_NAMES:
        found = shutil.which(name)
        if found:
            return found
    for path in engine_search_paths(project_dir, state_dir):
        existing = maybe_existing_executable(path)
        if existing:
            return existing
    wrapper = shutil.which("ple-vina")
    if wrapper:
        return wrapper
    return "ple-vina"


def command_exists(command: str) -> bool:
    path = Path(command).expanduser()
    return (path.exists() and not path.is_dir()) or shutil.which(command) is not None


def receptor_pdbqt_path(receptor: Path, state_dir: Path) -> Path:
    return state_dir / "receptor" / f"{receptor.stem}.pdbqt"


def prepare_receptor_with_obabel(
    *,
    receptor: Path,
    output: Path,
    obabel: str,
    add_hydrogens: bool,
    dry_run: bool,
) -> Path:
    if not command_exists(obabel):
        raise FileNotFoundError(f"OpenBabel not found: {obabel}; install obabel or pass --obabel /path/to/obabel")
    output.parent.mkdir(parents=True, exist_ok=True)
    cmd = [obabel, str(receptor)]
    if add_hydrogens:
        cmd.append("-h")
    cmd.extend(["-O", str(output), "-xr", "-xc", "-xn"])
    run(cmd, dry_run=dry_run)
    if not dry_run and (not output.exists() or output.stat().st_size == 0):
        raise RuntimeError(f"OpenBabel did not create receptor PDBQT: {output}")
    if not dry_run:
        clean_receptor_pdbqt(output)
    return output


def clean_receptor_pdbqt(path: Path) -> None:
    allowed_prefixes = ("ATOM", "HETATM", "TER", "END")
    raw_lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    cleaned = [line for line in raw_lines if line.startswith(allowed_prefixes)]
    if not cleaned:
        raise RuntimeError(f"prepared receptor PDBQT has no atom records after cleaning: {path}")
    path.write_text("\n".join(cleaned) + "\n", encoding="utf-8")


def ensure_receptor_for_docking(args: argparse.Namespace, receptor: Path, state_dir: Path) -> Path:
    if receptor.suffix.lower() == ".pdbqt":
        return receptor
    if receptor.suffix.lower() != ".pdb":
        raise ValueError(f"unsupported receptor suffix for docking: {receptor.suffix}")
    if args.receptor_prep == "none":
        return receptor
    output = receptor_pdbqt_path(receptor, state_dir)
    if output.exists() and output.stat().st_size > 0 and not args.force_receptor_prep:
        print(f"receptor_pdbqt={output}", flush=True)
        return output.resolve()
    print(f"preparing_receptor_pdbqt={output}", flush=True)
    if args.receptor_prep not in {"auto", "obabel"}:
        raise ValueError(f"unsupported --receptor-prep {args.receptor_prep!r}")
    return prepare_receptor_with_obabel(
        receptor=receptor,
        output=output,
        obabel=args.obabel,
        add_hydrogens=args.receptor_add_hydrogens,
        dry_run=args.dry_run,
    ).resolve()


def run_workflow(args: argparse.Namespace) -> None:
    preset = degree_preset(args.degree)
    project_dir = project_dir_from_args(args.project_dir)
    source_receptor = detect_receptor(project_dir, args.pdb)
    ligands = [] if args.analysis_only else collect_ligands(args, project_dir)
    skills_dir = find_skills_dir(args.skills_dir)
    dock_script = sibling_scripts(skills_dir)
    state_dir = project_dir / args.state_dir
    config_dir = state_dir / "configs"
    history_dir = state_dir / "dock_history"
    state_dir.mkdir(parents=True, exist_ok=True)
    site_receptor = prepare_receptor_site_copy(args, source_receptor, state_dir)
    site_box = site_box_from_args(args, site_receptor)
    receptor = ensure_receptor_for_docking(args, site_receptor, state_dir)

    initial_rounds = int(effective_value(args, "initial_rounds", preset.initial_rounds))
    initial_batch_size = int(effective_value(args, "initial_batch_size", preset.initial_batch_size))
    focus_rounds = int(effective_value(args, "focus_rounds", preset.focus_rounds))
    focus_batch_size = int(effective_value(args, "focus_batch_size", preset.focus_batch_size))
    edit_mode = args.edit_mode or preset.edit_mode
    history_parent_count = int(effective_value(args, "history_parent_count", preset.history_parent_count))
    history_rounds = int(effective_value(args, "history_rounds", preset.history_rounds))
    history_batch_size = int(effective_value(args, "history_batch_size", preset.history_batch_size))
    base_engine_args = list(args.engine_arg or [])

    if args.analysis_only:
        if args.config:
            analysis_config = resolve_project_path(project_dir, args.config)
        else:
            existing_configs = [
                config_dir / "focused_box_round1.txt",
                config_dir / "initial_large_box.txt",
            ]
            analysis_config = next((path for path in existing_configs if path.exists()), config_dir / "initial_large_box.txt")
            if not analysis_config.exists():
                initial_box = site_box or bounding_box(
                    read_atom_coords(site_receptor),
                    padding=args.initial_padding,
                    max_size=args.initial_max_size,
                )
                write_vina_config(analysis_config, receptor, initial_box, args.exhaustiveness, args.num_modes, args.cpu)
                print(f"initial_box_config={analysis_config}", flush=True)
        if not args.skip_pose_analysis:
            run_pose_analysis(
                skills_dir=skills_dir,
                history_dir=history_dir,
                receptor=receptor,
                config=analysis_config,
                state_dir=state_dir,
                save_pse=args.save_pse,
                dry_run=args.dry_run,
            )
        if args.build_space_dashboard:
            run_score_space_dashboard(
                skills_dir=skills_dir,
                history_dir=history_dir,
                receptor=receptor,
                state_dir=state_dir,
                dry_run=args.dry_run,
                train_size=args.space_train_size,
                test_size=args.space_test_size,
            )
        if args.skip_pose_analysis and not args.build_space_dashboard:
            print("analysis_only_nothing_to_run=use_pose_analysis_or_build_space_dashboard", flush=True)
        return

    engine_exe = resolve_engine_exe(args.engine_exe, project_dir, state_dir)
    print(f"docking_engine={engine_exe}", flush=True)
    if not args.skip_engine_check and not args.dry_run:
        check_engine(engine_exe)
    gnina_engine = is_gnina_engine(engine_exe)
    requested_cnn_scoring = args.gnina_cnn_scoring
    wide_cnn_scoring = requested_cnn_scoring
    focused_cnn_scoring = requested_cnn_scoring
    if gnina_engine and requested_cnn_scoring == "refinement" and args.gnina_wide_scout:
        wide_cnn_scoring = "rescore"
        focused_cnn_scoring = "refinement"
        if args.initial_rounds is None:
            initial_rounds = 1
        if args.initial_batch_size is None:
            initial_batch_size = min(initial_batch_size, args.wide_scout_batch_size)
        print(
            f"gnina_wide_scout=enabled wide_cnn_scoring={wide_cnn_scoring} "
            f"focused_cnn_scoring={focused_cnn_scoring}",
            flush=True,
        )
    wide_engine_args = engine_args_with_cnn(base_engine_args, wide_cnn_scoring)
    focused_engine_args = engine_args_with_cnn(base_engine_args, focused_cnn_scoring)

    print(
        f"degree={args.degree} initial_batch={initial_batch_size} focus_batch={focus_batch_size} "
        f"edit_mode={edit_mode} history_parents={history_parent_count}",
        flush=True,
    )

    initial_config = resolve_project_path(project_dir, args.config) if args.config else config_dir / "initial_large_box.txt"
    initial_box: Box | None = None
    if args.config is None:
        initial_box = site_box or bounding_box(
            read_atom_coords(site_receptor),
            padding=args.initial_padding,
            max_size=args.initial_max_size,
        )
        write_vina_config(initial_config, receptor, initial_box, args.exhaustiveness, args.num_modes, args.cpu)
        print(f"initial_box_config={initial_config}", flush=True)

    dock_ligand_batch(
        ligands=ligands,
        dock_script=dock_script,
        history_dir=history_dir,
        receptor=receptor,
        config=initial_config,
        engine_exe=engine_exe,
        max_rounds=initial_rounds,
        batch_size=initial_batch_size,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        cpu=args.cpu,
        edit_mode=edit_mode,
        engine_args=wide_engine_args,
        dry_run=args.dry_run,
        nickname_suffix="wide",
    )

    focused_config = initial_config
    if not args.no_focus:
        focused_box, selected_count, record_count = choose_focused_box(history_dir / "dock_history.csv", args)
        if focused_box is None:
            print("focus_box_skipped=no_successful_pose_coords", flush=True)
        else:
            if initial_box is not None and box_is_scattered(focused_box, initial_box, args.scatter_volume_fraction):
                print("focus_box_skipped=poses_scattered_across_initial_box", flush=True)
            else:
                focused_config = config_dir / "focused_box_round1.txt"
                write_vina_config(focused_config, receptor, focused_box, args.exhaustiveness, args.num_modes, args.cpu)
                print(
                    f"focused_box_config={focused_config} "
                    f"cluster_poses={selected_count}/{record_count} "
                    f"center=({focused_box.center_x:.3f},{focused_box.center_y:.3f},{focused_box.center_z:.3f}) "
                    f"size=({focused_box.size_x:.1f},{focused_box.size_y:.1f},{focused_box.size_z:.1f})",
                    flush=True,
                )
                dock_ligand_batch(
                    ligands=ligands,
                    dock_script=dock_script,
                    history_dir=history_dir,
                    receptor=receptor,
                    config=focused_config,
                    engine_exe=engine_exe,
                    max_rounds=focus_rounds,
                    batch_size=focus_batch_size,
                    exhaustiveness=args.exhaustiveness,
                    num_modes=args.num_modes,
                    cpu=args.cpu,
                    edit_mode=edit_mode,
                    engine_args=focused_engine_args,
                    dry_run=args.dry_run,
                    nickname_suffix="focused",
                )

    history_parents = select_history_parents(
        history_dir / "dock_history.csv",
        count=history_parent_count,
        exclude_smiles={ligand.smiles for ligand in ligands},
    )
    if history_parents and history_batch_size > 0:
        print(
            "history_parent_smiles="
            + ",".join(parent.nickname for parent in history_parents),
            flush=True,
        )
        dock_ligand_batch(
            ligands=history_parents,
            dock_script=dock_script,
            history_dir=history_dir,
            receptor=receptor,
            config=focused_config,
            engine_exe=engine_exe,
            max_rounds=history_rounds,
            batch_size=history_batch_size,
            exhaustiveness=args.exhaustiveness,
            num_modes=args.num_modes,
            cpu=args.cpu,
            edit_mode=edit_mode,
            engine_args=focused_engine_args,
            dry_run=args.dry_run,
            nickname_suffix=preset.history_nickname_suffix,
        )

    if not args.skip_pose_analysis:
        run_pose_analysis(
            skills_dir=skills_dir,
            history_dir=history_dir,
            receptor=receptor,
            config=focused_config,
            state_dir=state_dir,
            save_pse=args.save_pse,
            dry_run=args.dry_run,
        )
    if args.build_space_dashboard:
        run_score_space_dashboard(
            skills_dir=skills_dir,
            history_dir=history_dir,
            receptor=receptor,
            state_dir=state_dir,
            dry_run=args.dry_run,
            train_size=args.space_train_size,
            test_size=args.space_test_size,
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ple", description="Pocket ligand explorer")
    subparsers = parser.add_subparsers(dest="command")

    def add_focus_box_args(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument("--focus-padding", type=float, default=DEFAULT_FOCUS_PADDING)
        subparser.add_argument("--focus-min-size", type=float, default=DEFAULT_FOCUS_MIN_SIZE)
        subparser.add_argument("--focus-target-size", type=float, default=DEFAULT_FOCUS_TARGET_SIZE, help="Preferred focused box side length in Angstroms when poses are scattered")
        subparser.add_argument("--focus-max-size", type=float, default=DEFAULT_FOCUS_MAX_SIZE, help="Largest automatic focused box side length in Angstroms")
        subparser.add_argument("--focus-cluster-radius", type=float, default=DEFAULT_FOCUS_CLUSTER_RADIUS, help="Centroid radius used to find the common pose landing region")
        subparser.add_argument("--focus-min-cluster-size", type=int, default=DEFAULT_FOCUS_MIN_CLUSTER_SIZE)
        subparser.add_argument("--focus-max-poses", type=int, default=DEFAULT_FOCUS_MAX_POSES, help="Top pose records considered when choosing the focused box")

    def add_site_args(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument("--drop-receptor-resname", action="append", default=[], help="Remove this receptor residue name before docking, e.g. HHH; comma-separated or repeatable")
        subparser.add_argument("--site-resname", help="Use atoms from this receptor residue name as the initial box center, e.g. FE")
        subparser.add_argument("--site-chain", help="Restrict --site-resname to one chain, e.g. A")
        subparser.add_argument("--site-resseq", help="Restrict --site-resname to one residue number")
        subparser.add_argument("--site-atom-name", help="Restrict --site-resname to one atom name")
        subparser.add_argument("--site-box-size", type=float, default=DEFAULT_FOCUS_TARGET_SIZE, help="Cubic initial box side length around --site-resname, in Angstroms")

    init_parser = subparsers.add_parser("init", help="Create a simple project skeleton")
    init_parser.add_argument("--project-dir", type=Path, help="Project root; defaults to current directory")
    init_parser.add_argument("--force", action="store_true", help="Overwrite existing template files")
    prep_parser = subparsers.add_parser("prepare-receptor", help="Prepare the project receptor PDB as PDBQT")
    prep_parser.add_argument("--project-dir", type=Path, help="Project root; defaults to current directory")
    prep_parser.add_argument("--state-dir", default=DEFAULT_STATE_DIRNAME, help="Workflow output folder inside project-dir")
    prep_parser.add_argument("--pdb", type=Path, help="Receptor PDB/PDBQT; required when project has multiple receptor files")
    prep_parser.add_argument("--obabel", default="obabel", help="OpenBabel executable name or path")
    prep_parser.add_argument("--receptor-add-hydrogens", action="store_true", help="Ask OpenBabel to add hydrogens before PDBQT conversion")
    prep_parser.add_argument("--force", action="store_true", help="Overwrite an existing prepared receptor PDBQT")
    prep_parser.add_argument("--dry-run", action="store_true")

    repair_parser = subparsers.add_parser("repair-loops", help="Audit SEQRES/ATOM gaps and prepare isolated one-residue Rosetta loop repairs")
    add_repair_loops_parser(repair_parser)

    hpp_parser = subparsers.add_parser("hpp-map", help="Convert H++ pkout records into an Amber protonation-map CSV")
    hpp_parser.add_argument("--project-dir", type=Path, help="Project root; defaults to current directory")
    hpp_parser.add_argument("--pkout", type=Path, required=True, help="H++ .pkout or copied pkout.txt file")
    hpp_parser.add_argument("--pdb", type=Path, required=True, help="Reference receptor PDB used for chain/residue numbering and metal-contact audit")
    hpp_parser.add_argument("--out", type=Path, help="Output CSV; defaults to <project-dir>/protonation_map.csv")
    hpp_parser.add_argument("--chain", default="A", help="Protein chain used in H++; default A")
    hpp_parser.add_argument("--ph", type=float, default=7.4, help="pH used for heuristic protonation choices")
    hpp_parser.add_argument("--metal-cutoff", type=float, default=2.8, help="Fe coordination contact cutoff in Angstroms")

    metal_parser = subparsers.add_parser("metal-model", help="Prepare a first-pass MCPB.py metal-center model draft")
    metal_parser.add_argument("--project-dir", type=Path, help="Project root; defaults to current directory")
    metal_parser.add_argument("--pdb", type=Path, help="Receptor/MD receptor PDB containing the metal atom")
    metal_parser.add_argument("--ligand-pdb", type=Path, help="Optional docked ligand PDB to include as a possible metal ligand")
    metal_parser.add_argument("--ligand-mol2", type=Path, help="Optional antechamber mol2 for the docked ligand; copied into MCPB.py naa_mol2files")
    metal_parser.add_argument("--outdir", type=Path, help="Output directory for metal-site report and MCPB.py draft inputs")
    metal_parser.add_argument("--metal-element", default="FE", help="Metal element to model, default FE")
    metal_parser.add_argument("--metal-serial", type=int, help="Restrict to a specific PDB atom serial")
    metal_parser.add_argument("--metal-chain", help="Restrict metal by chain ID")
    metal_parser.add_argument("--metal-resseq", help="Restrict metal by residue number")
    metal_parser.add_argument("--metal-charge", type=float, default=3.0, help="Formal metal charge for the generated metal mol2; default 3.0 for Fe(III)")
    metal_parser.add_argument("--cutoff", type=float, default=2.8, help="Coordination search cutoff in Angstroms")
    metal_parser.add_argument("--coord-elements", default="N,O,S", help="Elements treated as likely coordination atoms, comma-separated")
    metal_parser.add_argument("--ligand-resname", default="LIG", help="Residue name to assign to the docked ligand in MCPB draft PDB")
    metal_parser.add_argument("--ligand-chain", default="Z", help="Chain ID to assign to the docked ligand in MCPB draft PDB")
    metal_parser.add_argument(
        "--ligand-role",
        choices=["nonbonded", "coordinating"],
        default="nonbonded",
        help="Use nonbonded to keep the docked ligand out of MCPB bonded metal-site fitting; coordinating includes it in the metal model",
    )
    metal_parser.add_argument("--group-name", help="MCPB.py group_name; defaults to a metal/chain/residue label")
    metal_parser.add_argument("--water-model", default="TIP3P", help="MCPB.py water_model line")
    metal_parser.add_argument("--force-field", default="ff19SB", help="MCPB.py force_field value; default ff19SB")
    metal_parser.add_argument("--gaff", type=int, choices=[0, 1, 2], default=2, help="MCPB.py gaff setting: 0 off, 1 GAFF, 2 GAFF2")
    metal_parser.add_argument("--software-version", choices=["g03", "g09", "g16", "gau", "gms"], default="g16", help="QM input flavor for MCPB.py; default g16")
    metal_parser.add_argument("--smmodel-chg", type=int, help="Charge for the MCPB small QM model; omit for MCPB.py auto charge")
    metal_parser.add_argument("--smmodel-spin", type=int, default=6, help="Spin multiplicity for the MCPB small QM model; default 6 for high-spin Fe(III)")
    metal_parser.add_argument("--lgmodel-chg", type=int, help="Charge for the MCPB large RESP model; omit for MCPB.py auto charge")
    metal_parser.add_argument("--lgmodel-spin", type=int, default=6, help="Spin multiplicity for the MCPB large RESP model; default 6 for high-spin Fe(III)")
    metal_parser.add_argument("--large-opt", type=int, choices=[0, 1, 2], default=1, help="MCPB.py large_opt: 0 none, 1 H-only, 2 full large-model optimization")
    metal_parser.add_argument("--add-redcrd", type=int, choices=[0, 1, 2], default=1, help="MCPB.py add_redcrd setting for small-model redundant coordinates")
    metal_parser.add_argument("--scale-factor", type=float, default=1.0, help="MCPB.py force-constant frequency scale factor")
    metal_parser.add_argument("--ambertools-bin", type=Path, help="AmberTools bin directory; used to infer AMBERHOME for run_mcpb_step1.sh")
    metal_parser.add_argument("--amberhome", type=Path, help="Explicit AMBERHOME for run_mcpb_step1.sh")

    focus_parser = subparsers.add_parser("focus-box", help="Infer a focused docking box from existing dock_history poses")
    focus_parser.add_argument("--project-dir", type=Path, help="Project root; defaults to current directory")
    focus_parser.add_argument("--state-dir", default=DEFAULT_STATE_DIRNAME, help="Workflow output folder inside project-dir")
    focus_parser.add_argument("--pdb", type=Path, help="Receptor PDB/PDBQT; required when project has multiple receptor files")
    focus_parser.add_argument("--out-config", type=Path, help="Output Vina/GNINA config path")
    focus_parser.add_argument("--receptor-prep", choices=["auto", "obabel", "none"], default="auto", help="How to handle a PDB receptor before docking")
    focus_parser.add_argument("--obabel", default="obabel", help="OpenBabel executable name or path for receptor preparation")
    focus_parser.add_argument("--receptor-add-hydrogens", action="store_true")
    focus_parser.add_argument("--force-receptor-prep", action="store_true")
    focus_parser.add_argument("--exhaustiveness", type=int, default=8)
    focus_parser.add_argument("--num-modes", type=int, default=9)
    focus_parser.add_argument("--cpu", type=int)
    focus_parser.add_argument("--dry-run", action="store_true")
    add_site_args(focus_parser)
    add_focus_box_args(focus_parser)

    md_parser = subparsers.add_parser("md-handoff", help="Prepare AmberTools/ACPYPE/GROMACS files from a selected docked pose")
    md_parser.add_argument("--project-dir", type=Path, help="Project root; defaults to current directory")
    md_parser.add_argument("--state-dir", default=DEFAULT_STATE_DIRNAME, help="Workflow output folder inside project-dir")
    md_parser.add_argument("--pdb", type=Path, help="Receptor PDB; required when project has multiple receptor files")
    md_parser.add_argument("--seq-id", help="Dock-history seq_id to send to MD; default selects the best triage row")
    md_parser.add_argument("--pose-mode", type=int, default=1, help="Docked pose MODEL number to export for the ligand coordinates")
    md_parser.add_argument("--rank-by", choices=["triage", "affinity"], default="triage", help="Selection rule when --seq-id is not supplied")
    md_parser.add_argument("--outdir", type=Path, help="Output MD handoff directory")
    md_parser.add_argument("--ligand-charge", type=int, default=0, help="Formal ligand charge passed to antechamber -nc")
    md_parser.add_argument("--ligand-resname", default="LIG", help="Ligand residue name passed to antechamber")
    md_parser.add_argument("--protein-forcefield", default="leaprc.protein.ff19SB", help="AmberTools leaprc protein force field")
    md_parser.add_argument("--ligand-forcefield", default="leaprc.gaff2", help="AmberTools leaprc ligand force field")
    md_parser.add_argument("--ambertools-bin", type=Path, help="Directory containing AmberTools/ACPYPE commands such as antechamber, parmchk2, tleap, acpype, and parmed")
    md_parser.add_argument("--receptor-md-mode", choices=["protein-only", "full"], default="protein-only", help="Use a protein-only receptor copy for first MD smoke tests, or keep the full PDB when all cofactors are parameterized")
    md_parser.add_argument("--protonation-map", type=Path, help="Manual H++/curation residue rename map, e.g. CSV chain,resseq,resname or lines like A:123=HIE")
    md_parser.add_argument("--protonate", action="append", default=[], help="One manual Amber residue rename, e.g. A:123=HID, A:45=ASH, or 128=GLH; repeat as needed")
    md_parser.add_argument("--keep-hetatm-resname", action="append", default=[], help="Keep this HETATM residue name in protein-only MD mode; comma-separated or repeatable")
    md_parser.add_argument("--keep-hetatm-element", action="append", default=[], help="Keep this HETATM element in protein-only MD mode; comma-separated or repeatable")
    md_parser.add_argument("--drop-default-md-metals", action="store_true", help="Do not automatically keep Fe/FE records in protein-only MD mode")
    md_parser.add_argument("--amber-frcmod", action="append", default=[], help="Extra Amber frcmod file copied into handoff and loaded by tleap; repeat as needed")
    md_parser.add_argument("--amber-lib", action="append", default=[], help="Extra Amber .lib/.off file copied into handoff and loaded by tleap")
    md_parser.add_argument("--amber-prep", action="append", default=[], help="Extra Amber prepi/prep file copied into handoff and loaded by tleap")
    md_parser.add_argument("--amber-mol2", action="append", default=[], help="Extra Amber mol2 residue/cofactor file copied into handoff and loaded by tleap")
    md_parser.add_argument("--water-box", default="TIP3PBOX", help="tleap solvent box name")
    md_parser.add_argument("--box-padding", type=float, default=10.0, help="tleap solvateBox padding in Angstroms")
    md_parser.add_argument("--gmx", default="gmx_mpi", help="GROMACS executable name or path")
    md_parser.add_argument("--gmxrc", type=Path, help="Optional GROMACS GMXRC file to source before running gmx")
    md_parser.add_argument("--extra-ld-library-path", action="append", default=[], help="Extra library path exported before GROMACS commands; repeat as needed")
    md_parser.add_argument("--time", "--production-ns", dest="production_ns", type=float, default=20.0, help="Production MD length in ns; for example --time 10 writes a 10 ns production run")
    md_parser.add_argument("--dt-ps", type=float, default=0.002, help="MD timestep in ps")
    md_parser.add_argument("--maxwarn", type=int, default=1, help="GROMACS grompp -maxwarn")
    md_parser.add_argument("--mdrun-args", default="", help="Extra tokens appended to every gmx mdrun command, e.g. '-nb gpu -pme gpu'")
    md_parser.add_argument("--mdrun-arg", action="append", default=[], help="One extra token appended to every gmx mdrun command; repeat for tokens that start with '-'")

    md_run_parser = subparsers.add_parser("md-run", help="Prepare, run, fix-PBC, and analyze a selected docked pose with a flat MD project layout")
    md_run_parser.add_argument("--project-dir", type=Path, help="Project root; defaults to current directory")
    md_run_parser.add_argument("--state-dir", default=DEFAULT_STATE_DIRNAME, help="Docking state folder containing dock_history; output stays flat in project-dir/md_handoff")
    md_run_parser.add_argument("--pdb", type=Path, help="H++/curated receptor PDB for MD; required when project has multiple receptor files")
    md_run_parser.add_argument("--seq-id", help="Dock-history seq_id to send to MD; default selects the best triage row")
    md_run_parser.add_argument("--pose-mode", type=int, default=1, help="Docked pose MODEL number to export for the ligand coordinates")
    md_run_parser.add_argument("--rank-by", choices=["triage", "affinity"], default="triage")
    md_run_parser.add_argument("--outdir", type=Path, default=Path("md_handoff"), help="Flat handoff directory; default <project-dir>/md_handoff")
    md_run_parser.add_argument("--ligand-charge", type=int, default=0)
    md_run_parser.add_argument("--ligand-resname", default="LIG")
    md_run_parser.add_argument("--protein-forcefield", default="leaprc.protein.ff19SB")
    md_run_parser.add_argument("--ligand-forcefield", default="leaprc.gaff2")
    md_run_parser.add_argument("--ambertools-bin", type=Path)
    md_run_parser.add_argument("--receptor-md-mode", choices=["protein-only", "full"], default="protein-only")
    md_run_parser.add_argument("--protonation-map", type=Path, help="CSV from ple hpp-map or manual Amber residue rename map")
    md_run_parser.add_argument("--protonate", action="append", default=[], help="Manual Amber residue rename, e.g. A:123=HIE")
    md_run_parser.add_argument("--keep-hetatm-resname", action="append", default=[])
    md_run_parser.add_argument("--keep-hetatm-element", action="append", default=[])
    md_run_parser.add_argument("--drop-default-md-metals", action="store_true")
    md_run_parser.add_argument("--amber-frcmod", action="append", default=[])
    md_run_parser.add_argument("--amber-lib", action="append", default=[])
    md_run_parser.add_argument("--amber-prep", action="append", default=[])
    md_run_parser.add_argument("--amber-mol2", action="append", default=[])
    md_run_parser.add_argument("--water-box", default="TIP3PBOX")
    md_run_parser.add_argument("--box-padding", type=float, default=10.0)
    md_run_parser.add_argument("--gmx", default="gmx_mpi")
    md_run_parser.add_argument("--gmxrc", type=Path)
    md_run_parser.add_argument("--extra-ld-library-path", action="append", default=[])
    md_run_parser.add_argument("--time", "--production-ns", dest="production_ns", type=float, default=20.0)
    md_run_parser.add_argument("--dt-ps", type=float, default=0.002)
    md_run_parser.add_argument("--maxwarn", type=int, default=1)
    md_run_parser.add_argument("--mdrun-args", default="")
    md_run_parser.add_argument("--mdrun-arg", action="append", default=[])
    md_run_parser.add_argument("--prepare-only", action="store_true", help="Only write the coherent MD scripts; do not execute them")
    md_run_parser.add_argument("--background", action="store_true", help="Run 02_run_gromacs.sh in the background after preparation; postprocess is chained unless --skip-postprocess is set")
    md_run_parser.add_argument("--skip-postprocess", action="store_true", help="Skip automatic PBC/RDC/PyMOL postprocessing after a foreground run")
    md_run_parser.add_argument("--dry-run", action="store_true")
    md_run_parser.add_argument("--center-group", default="Protein", help="GROMACS index group or name for centering during postprocess; default Protein")
    md_run_parser.add_argument("--output-group", default="System", help="GROMACS output group or name for postprocess; default System")
    md_run_parser.add_argument("--metal-radius", type=float, default=4.0)
    md_run_parser.add_argument("--carbon-color", default="green")
    md_run_parser.add_argument("--ligand-selection", default="resname LIG")
    md_run_parser.add_argument("--distance-cutoff", type=float, default=10.0)
    md_run_parser.add_argument("--rdc-stride", type=int, default=1)
    md_run_parser.add_argument("--top-n", type=int, default=80)
    md_run_parser.add_argument("--velocity-threshold", type=float, default=1.0)

    run_parser = subparsers.add_parser("run", help="Run wide-box docking, pose-guided focusing, and report generation")
    run_parser.add_argument("--project-dir", type=Path, help="Project root; defaults to current directory")
    run_parser.add_argument("--state-dir", default=DEFAULT_STATE_DIRNAME, help="Workflow output folder inside project-dir")
    run_parser.add_argument("--pdb", type=Path, help="Receptor PDB/PDBQT; required when project has multiple receptor files")
    run_parser.add_argument("--smiles", action="append", help="Starting ligand SMILES; can be repeated")
    run_parser.add_argument("--nickname", default="", help="Nickname for manually supplied --smiles")
    run_parser.add_argument("--smiles-file", type=Path, help="SMI/TXT file: SMILES [nickname] per line")
    run_parser.add_argument("--csv", type=Path, help="CSV file containing SMILES")
    run_parser.add_argument("--smiles-column", default="smiles")
    run_parser.add_argument("--nickname-column")
    run_parser.add_argument("--sdf", type=Path, help="SDF ligand file converted to SMILES with RDKit")
    run_parser.add_argument("--mol", type=Path, help="MOL ligand file converted to SMILES with RDKit")
    run_parser.add_argument("--config", type=Path, help="Existing Vina config; skips automatic large-box creation")
    run_parser.add_argument("--engine-exe", default="auto", help="Docking executable name/path; auto searches PATH, env vars, conda, project, and ~/vina_task2/vina_bin")
    run_parser.add_argument("--engine-arg", action="append", default=[], help="Extra docking-engine token passed through to dock_smiles.py; repeat for multi-token options")
    run_parser.add_argument("--gnina-cnn-scoring", choices=["none", "rescore", "refinement", "all"], help="Pass GNINA --cnn_scoring value, for example refinement")
    run_parser.add_argument("--gnina-wide-scout", action=argparse.BooleanOptionalAction, default=True, help="For GNINA refinement runs, use a cheaper wide-box scout before focused refinement")
    run_parser.add_argument("--wide-scout-batch-size", type=int, default=DEFAULT_WIDE_SCOUT_BATCH_SIZE)
    run_parser.add_argument("--receptor-prep", choices=["auto", "obabel", "none"], default="auto", help="How to handle a PDB receptor before docking")
    run_parser.add_argument("--obabel", default="obabel", help="OpenBabel executable name or path for receptor preparation")
    run_parser.add_argument("--receptor-add-hydrogens", action="store_true", help="Ask OpenBabel to add hydrogens during receptor PDBQT conversion")
    run_parser.add_argument("--force-receptor-prep", action="store_true", help="Regenerate receptor PDBQT even if it already exists")
    run_parser.add_argument("--skills-dir", type=Path, help="Directory containing smiles-to-vina-docking and pose-analyzer")
    run_parser.add_argument("--initial-padding", type=float, default=DEFAULT_INITIAL_PADDING)
    run_parser.add_argument("--initial-max-size", type=float, default=DEFAULT_INITIAL_MAX_SIZE)
    add_site_args(run_parser)
    add_focus_box_args(run_parser)
    run_parser.add_argument("--scatter-volume-fraction", type=float, default=0.55)
    run_parser.add_argument("--degree", type=int, choices=sorted(DEGREE_PRESETS), default=1, help="Exploration strength: 1 about 50 rows, 2 about 100, 3-5 add more drastic/history-guided exploration")
    run_parser.add_argument("--edit-mode", help="Override degree edit families, e.g. add,delete,shrink,replace,drastic")
    run_parser.add_argument("--initial-rounds", type=int)
    run_parser.add_argument("--initial-batch-size", type=int)
    run_parser.add_argument("--focus-rounds", type=int)
    run_parser.add_argument("--focus-batch-size", type=int)
    run_parser.add_argument("--history-parent-count", type=int, help="Override number of successful dock_history rows selected as extra parents")
    run_parser.add_argument("--history-rounds", type=int, help="Override history-guided analog rounds")
    run_parser.add_argument("--history-batch-size", type=int, help="Override analog batch per selected history parent")
    run_parser.add_argument("--exhaustiveness", type=int, default=8)
    run_parser.add_argument("--num-modes", type=int, default=9)
    run_parser.add_argument("--cpu", type=int)
    run_parser.add_argument("--no-focus", action="store_true")
    run_parser.add_argument("--analysis-only", action="store_true", help="Skip docking and rebuild pose-analysis/report/dashboard outputs from the existing dock_history.csv")
    run_parser.add_argument("--skip-pose-analysis", action="store_true")
    run_parser.add_argument("--save-pse", action="store_true")
    run_parser.add_argument("--build-space-dashboard", action="store_true", help="Build the pose-analyzer molecule-space dashboard from dock_history for clustering/dimensionality triage")
    run_parser.add_argument("--space-train-size", type=int, default=60, help="Reserved for supervised score-space runs; ignored by PLE unsupervised dashboards")
    run_parser.add_argument("--space-test-size", type=int, default=20, help="Reserved for supervised score-space runs; ignored by PLE unsupervised dashboards")
    run_parser.add_argument("--skip-engine-check", action="store_true")
    run_parser.add_argument("--dry-run", action="store_true")
    add_md_tool_parsers(subparsers)
    return parser


def init_project(args: argparse.Namespace) -> None:
    project_dir = project_dir_from_args(args.project_dir)
    project_dir.mkdir(parents=True, exist_ok=True)
    state_dir = project_dir / DEFAULT_STATE_DIRNAME
    (state_dir / "configs").mkdir(parents=True, exist_ok=True)
    ligand_file = project_dir / "ligands.smi"
    note_file = project_dir / "PLE_PROJECT.txt"
    if args.force or not ligand_file.exists():
        ligand_file.write_text("# SMILES nickname\n", encoding="utf-8")
    if args.force or not note_file.exists():
        note_file.write_text(
            "\n".join(
                [
                    "Pocket Ligand Explorer project",
                    "",
                    "Put exactly one receptor .pdbqt or .pdb file in this directory,",
                    "or run ple run --pdb path/to/receptor.pdbqt.",
                    "",
                    "Add ligands to ligands.smi as: SMILES nickname",
                    "Run: ple run --smiles-file ligands.smi",
                    "",
                ]
            ),
            encoding="utf-8",
        )
    print(f"initialized={project_dir}", flush=True)


def prepare_receptor_command(args: argparse.Namespace) -> None:
    project_dir = project_dir_from_args(args.project_dir)
    receptor = detect_receptor(project_dir, args.pdb)
    state_dir = project_dir / args.state_dir
    state_dir.mkdir(parents=True, exist_ok=True)
    if receptor.suffix.lower() == ".pdbqt":
        print(f"receptor_already_pdbqt={receptor}", flush=True)
        return
    output = receptor_pdbqt_path(receptor, state_dir)
    if output.exists() and output.stat().st_size > 0 and not args.force:
        print(f"receptor_pdbqt={output}", flush=True)
        return
    prepared = prepare_receptor_with_obabel(
        receptor=receptor,
        output=output,
        obabel=args.obabel,
        add_hydrogens=args.receptor_add_hydrogens,
        dry_run=args.dry_run,
    )
    print(f"receptor_pdbqt={prepared}", flush=True)


def focus_box_command(args: argparse.Namespace) -> None:
    project_dir = project_dir_from_args(args.project_dir)
    source_receptor = detect_receptor(project_dir, args.pdb)
    state_dir = project_dir / args.state_dir
    config_dir = state_dir / "configs"
    history_dir = state_dir / "dock_history"
    state_dir.mkdir(parents=True, exist_ok=True)
    site_receptor = prepare_receptor_site_copy(args, source_receptor, state_dir)
    receptor = ensure_receptor_for_docking(args, site_receptor, state_dir)
    history_csv = history_dir / "dock_history.csv"
    if not history_csv.exists():
        raise FileNotFoundError(f"no dock history found: {history_csv}")
    focused_box, selected_count, record_count = choose_focused_box(history_csv, args)
    if focused_box is None:
        raise RuntimeError("could not infer focused box: no successful pose coordinates found")
    output = args.out_config
    if output is None:
        output = config_dir / "focused_box_round1.txt"
    else:
        output = output.expanduser()
        if not output.is_absolute():
            output = project_dir / output
    write_vina_config(output, receptor, focused_box, args.exhaustiveness, args.num_modes, args.cpu)
    print(
        f"focused_box_config={output} "
        f"cluster_poses={selected_count}/{record_count} "
        f"center=({focused_box.center_x:.3f},{focused_box.center_y:.3f},{focused_box.center_z:.3f}) "
        f"size=({focused_box.size_x:.1f},{focused_box.size_y:.1f},{focused_box.size_z:.1f})",
        flush=True,
    )


def hpp_map_command(args: argparse.Namespace) -> None:
    project_dir = project_dir_from_args(args.project_dir)
    pkout = args.pkout.expanduser()
    if not pkout.is_absolute():
        pkout = project_dir / pkout
    pdb = args.pdb.expanduser()
    if not pdb.is_absolute():
        pdb = project_dir / pdb
    out = args.out.expanduser() if args.out else project_dir / "protonation_map.csv"
    if not out.is_absolute():
        out = project_dir / out
    mapped_count, metal_contact_count = write_protonation_map(
        pkout=pkout,
        pdb=pdb,
        out_csv=out,
        chain=args.chain,
        ph=args.ph,
        metal_cutoff_a=args.metal_cutoff,
    )
    print(f"protonation_map={out}", flush=True)
    print(f"mapped_residues={mapped_count}", flush=True)
    print(f"metal_contacts={metal_contact_count}", flush=True)
    print(f"metal_contacts_csv={out.with_name(out.stem + '_metal_contacts.csv')}", flush=True)


def successful_history_rows(history_csv: Path) -> list[dict[str, str]]:
    if not history_csv.exists():
        return []
    rows: list[dict[str, str]] = []
    with history_csv.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("reason"):
                continue
            if safe_float(row.get("affinity_kcal_mol")) is None:
                continue
            rows.append(row)
    return rows


def select_md_row(history_csv: Path, seq_id: str | None, rank_by: str) -> dict[str, str]:
    rows = successful_history_rows(history_csv)
    if not rows:
        raise RuntimeError(f"no successful docked rows found in {history_csv}")
    if seq_id:
        for row in rows:
            if row.get("seq_id") == seq_id:
                return row
        raise ValueError(f"--seq-id {seq_id!r} not found in successful rows")
    if rank_by == "affinity":
        rows.sort(key=lambda row: safe_float(row.get("affinity_kcal_mol")) or float("inf"))
    else:
        rows.sort(key=lambda row: history_score(row) if history_score(row) is not None else float("inf"))
    return rows[0]


def copy_existing(src_text: str, dst: Path) -> Path | None:
    if not src_text:
        return None
    src = Path(src_text).expanduser()
    if not src.exists():
        return None
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return dst


def pdb_count_key(line: str) -> str:
    record = line[:6].strip() or "OTHER"
    if line.startswith(("ATOM  ", "HETATM")):
        residue = line[17:20].strip() or "UNK"
        return f"{record}:{residue}"
    return record


def pdb_element(line: str) -> str:
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element:
        return element
    letters = "".join(ch for ch in line[12:16].strip() if ch.isalpha()).upper()
    if len(letters) >= 2 and letters[:2] in {"FE", "ZN", "MG", "CA", "MN", "CU", "CO", "NI", "NA", "CL"}:
        return letters[:2]
    return letters[:1]


def pdb_residue_key(line: str) -> tuple[str, str, str]:
    return (line[21:22].strip(), line[22:26].strip(), line[26:27].strip())


def normalize_resname(value: str) -> str:
    name = re.sub(r"[^A-Za-z0-9]", "", value).upper()
    if not name or len(name) > 3:
        raise ValueError(f"Amber residue names must be 1-3 letters/numbers, got {value!r}")
    return name


def parse_residue_assignment(text: str) -> tuple[tuple[str, str, str], str]:
    raw = text.strip()
    if not raw or raw.startswith("#"):
        raise ValueError("empty residue assignment")
    if "=" in raw:
        selector, new_name = raw.split("=", 1)
        selector = selector.strip()
        new_name = normalize_resname(new_name)
        parts = selector.replace("/", ":").split(":")
        if len(parts) == 1:
            chain = ""
            resid = parts[0].strip()
        elif len(parts) == 2:
            chain = parts[0].strip()
            resid = parts[1].strip()
        elif len(parts) == 3:
            chain = parts[0].strip()
            resid = parts[1].strip()
            icode = parts[2].strip()
            return (chain, resid, icode), new_name
        else:
            raise ValueError(f"could not parse residue selector {selector!r}")
        match = re.fullmatch(r"(-?\d+)([A-Za-z]?)", resid)
        if match:
            return (chain, match.group(1), match.group(2)), new_name
        return (chain, resid, ""), new_name

    fields = [field for field in re.split(r"[\s,]+", raw) if field]
    if len(fields) == 2:
        return ("", fields[0], ""), normalize_resname(fields[1])
    if len(fields) == 3:
        return (fields[0], fields[1], ""), normalize_resname(fields[2])
    if len(fields) == 4:
        # Accept either "chain resid icode new_name" or the common manual
        # curation form "chain resid old_name new_name".
        icode = fields[2] if len(fields[2]) == 1 else ""
        return (fields[0], fields[1], icode), normalize_resname(fields[3])
    raise ValueError(f"could not parse residue assignment {text!r}")


def load_protonation_map(path: Path | None, direct_assignments: list[str] | None, project_dir: Path) -> dict[tuple[str, str, str], str]:
    assignments: dict[tuple[str, str, str], str] = {}
    if path is not None:
        resolved = path.expanduser()
        if not resolved.is_absolute():
            resolved = project_dir / resolved
        if not resolved.exists():
            raise FileNotFoundError(resolved)
        lines = [line for line in resolved.read_text(encoding="utf-8", errors="replace").splitlines() if line.strip() and not line.lstrip().startswith("#")]
        if lines:
            header = [item.strip().lower() for item in next(csv.reader([lines[0]]))]
            has_resname = any(name in header for name in ("resname", "amber_resname", "name"))
            has_resid = any(name in header for name in ("resseq", "resid", "residue_id"))
            if "chain" in header and has_resname and has_resid:
                for row in csv.DictReader(lines):
                    lowered = {str(key).strip().lower(): str(value or "") for key, value in row.items()}
                    chain = lowered.get("chain", "").strip()
                    resid = (lowered.get("resseq") or lowered.get("resid") or lowered.get("residue_id") or "").strip()
                    icode = (lowered.get("icode") or lowered.get("insertion_code") or "").strip()
                    resname = normalize_resname(lowered.get("resname") or lowered.get("amber_resname") or lowered.get("name") or "")
                    assignments[(chain, resid, icode)] = resname
            else:
                for line in lines:
                    key, resname = parse_residue_assignment(line)
                    assignments[key] = resname
    for item in direct_assignments or []:
        key, resname = parse_residue_assignment(item)
        assignments[key] = resname
    return assignments


def md_keep_set(values: list[str] | None, defaults: set[str], drop_defaults: bool) -> set[str]:
    kept: set[str] = set() if drop_defaults else set(defaults)
    for value in values or []:
        for part in re.split(r"[,\s]+", value):
            if part:
                kept.add(part.upper())
    return kept


def apply_residue_rename(
    line: str,
    protonation_map: dict[tuple[str, str, str], str],
    renamed: dict[tuple[str, str, str, str, str], int],
) -> str:
    if not line.startswith(("ATOM  ", "HETATM")) or not protonation_map:
        return line
    chain, resid, icode = pdb_residue_key(line)
    new_name = protonation_map.get((chain, resid, icode)) or protonation_map.get(("", resid, icode))
    if not new_name:
        return line
    old_name = line[17:20].strip()
    if old_name == new_name:
        return line
    renamed[(chain, resid, icode, old_name, new_name)] = renamed.get((chain, resid, icode, old_name, new_name), 0) + 1
    return f"{line[:17]}{new_name:>3}{line[20:]}"


def write_md_receptor_pdb(
    source: Path,
    output: Path,
    mode: str,
    *,
    protonation_map: dict[tuple[str, str, str], str] | None = None,
    keep_hetatm_resnames: set[str] | None = None,
    keep_hetatm_elements: set[str] | None = None,
) -> dict[str, object]:
    counts: dict[str, int] = {}
    dropped: dict[str, int] = {}
    renamed: dict[tuple[str, str, str, str, str], int] = {}
    keep_resnames = {item.upper() for item in (keep_hetatm_resnames or set())}
    keep_elements = {item.upper() for item in (keep_hetatm_elements or set())}
    kept_lines: list[str] = [
        f"REMARK PLE receptor_md_mode={mode}",
        "REMARK Original receptor is preserved as input/receptor.pdb",
    ]
    atom_count = 0
    last_kept_was_atom = False
    last_kept_record = ""

    for raw in source.read_text(encoding="utf-8", errors="replace").splitlines():
        key = pdb_count_key(raw)
        counts[key] = counts.get(key, 0) + 1
        raw = apply_residue_rename(raw, protonation_map or {}, renamed)
        keep = False
        if mode == "full":
            keep = raw.startswith(("ATOM  ", "HETATM", "TER", "END"))
        elif mode == "protein-only":
            if raw.startswith("ATOM  "):
                keep = raw[17:20].strip() in STANDARD_MD_RESNAMES
            elif raw.startswith("HETATM"):
                resname = raw[17:20].strip().upper()
                element = pdb_element(raw)
                keep = resname in keep_resnames or element in keep_elements
            elif raw.startswith("TER"):
                keep = last_kept_was_atom
            elif raw.startswith("END"):
                keep = True
        else:
            raise ValueError(f"unsupported --receptor-md-mode {mode!r}")

        if keep:
            if raw.startswith("HETATM") and last_kept_record == "ATOM":
                kept_lines.append("TER")
                last_kept_was_atom = False
            kept_lines.append(raw)
            if raw.startswith(("ATOM  ", "HETATM")):
                atom_count += 1
                last_kept_was_atom = True
                last_kept_record = raw[:6].strip()
            elif raw.startswith("TER"):
                last_kept_was_atom = False
                last_kept_record = "TER"
            elif raw.startswith("END"):
                last_kept_record = "END"
        elif raw.startswith(("ATOM  ", "HETATM")):
            dropped[key] = dropped.get(key, 0) + 1

    if atom_count == 0:
        raise RuntimeError(f"MD receptor preparation kept no atoms from {source}")
    if not any(line.startswith("END") for line in kept_lines):
        kept_lines.append("END")
    output.write_text("\n".join(kept_lines) + "\n", encoding="utf-8")
    return {
        "mode": mode,
        "input": str(source),
        "output": str(output),
        "kept_atoms": atom_count,
        "kept_hetatm_resnames": sorted(keep_resnames),
        "kept_hetatm_elements": sorted(keep_elements),
        "manual_residue_renames": [
            {
                "chain": chain,
                "resseq": resid,
                "icode": icode,
                "from": old_name,
                "to": new_name,
                "atom_records": count,
            }
            for (chain, resid, icode, old_name, new_name), count in sorted(renamed.items())
        ],
        "input_counts": counts,
        "dropped_counts": dropped,
    }


def mdp_common() -> str:
    return """cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.0
coulombtype              = PME
rcoulomb                 = 1.0
rvdw                     = 1.0
pbc                      = xyz
constraints              = h-bonds
constraint-algorithm     = lincs
lincs-iter               = 1
lincs-order              = 4
"""


def md_shell_env(args: argparse.Namespace) -> str:
    lines: list[str] = []
    ambertools_bin = getattr(args, "ambertools_bin", None)
    if ambertools_bin:
        lines.append(f"export PATH=$PATH:{shlex.quote(str(Path(ambertools_bin).expanduser()))}")
    gmxrc = getattr(args, "gmxrc", None)
    if gmxrc:
        path = Path(gmxrc).expanduser()
        lines.append(f"if [[ -f {shlex.quote(str(path))} ]]; then set +u; source {shlex.quote(str(path))}; set -u; fi")
    else:
        gmx_path = Path(str(args.gmx)).expanduser()
        candidate = gmx_path.parent / "GMXRC"
        if gmx_path.parent != Path(".") and candidate.exists():
            lines.append(f"if [[ -f {shlex.quote(str(candidate))} ]]; then set +u; source {shlex.quote(str(candidate))}; set -u; fi")
    for path in getattr(args, "extra_ld_library_path", []) or []:
        lines.append(f"export LD_LIBRARY_PATH={shlex.quote(str(path))}:${{LD_LIBRARY_PATH:-}}")
    return "\n".join(lines)


def md_mdrun_array(args: argparse.Namespace) -> str:
    tokens: list[str] = []
    text = str(getattr(args, "mdrun_args", "") or "").strip()
    if text:
        try:
            tokens.extend(shlex.split(text))
        except ValueError:
            tokens.extend(text.split())
    tokens.extend(str(item) for item in (getattr(args, "mdrun_arg", []) or []) if str(item))
    return " ".join(shlex.quote(item) for item in tokens)


def ns_label(value: float) -> str:
    if abs(value - round(value)) < 1e-6:
        return str(int(round(value)))
    return f"{value:g}".replace(".", "p")


def copy_custom_amber_params(args: argparse.Namespace, amber_dir: Path, project_dir: Path) -> tuple[list[str], dict[str, list[str]]]:
    param_dir = amber_dir / "custom_params"
    leap_lines: list[str] = []
    copied: dict[str, list[str]] = {"frcmod": [], "lib": [], "prep": [], "mol2": []}

    for name in DEFAULT_AMBER_LOAD_PARAMS:
        leap_lines.append(f"loadamberparams {name}")

    def resolve(path: Path) -> Path:
        expanded = path.expanduser()
        if not expanded.is_absolute():
            expanded = project_dir / expanded
        if not expanded.exists():
            raise FileNotFoundError(expanded)
        return expanded.resolve()

    for arg_name, leap_command, key in (
        ("amber_frcmod", "loadamberparams", "frcmod"),
        ("amber_lib", "loadoff", "lib"),
        ("amber_prep", "loadamberprep", "prep"),
        ("amber_mol2", "loadmol2", "mol2"),
    ):
        for item in getattr(args, arg_name, []) or []:
            src = resolve(Path(item))
            param_dir.mkdir(parents=True, exist_ok=True)
            dst = param_dir / src.name
            shutil.copy2(src, dst)
            copied[key].append(str(dst))
            leap_lines.append(f"{leap_command} custom_params/{dst.name}")
    return leap_lines, copied


def read_pdb_atoms(path: Path) -> list[PdbAtom]:
    atoms: list[PdbAtom] = []
    for line in read_pdb_lines(path):
        atom = pdb_atom_from_line(line)
        if atom is not None:
            atoms.append(atom)
    return atoms


def pdb_line_with_serial(line: str, serial: int, *, resname: str | None = None, chain: str | None = None) -> str:
    text = line.rstrip()
    if len(text) < 80:
        text = text.ljust(80)
    text = f"{text[:6]}{serial:5d}{text[11:]}"
    if resname:
        text = f"{text[:17]}{resname[:3]:>3}{text[20:]}"
    if chain is not None:
        text = f"{text[:21]}{chain[:1] or ' '}{text[22:]}"
    return text.rstrip()


def write_metal_mol2(path: Path, metal: PdbAtom, *, charge: float) -> None:
    name = metal.element.upper() or metal.name.upper()
    path.write_text(
        f"""@<TRIPOS>MOLECULE
{name}
    1     0     1     0     0
SMALL
No Charge or Current Charge


@<TRIPOS>ATOM
      1 {name:<4s} {metal.x:10.4f} {metal.y:10.4f} {metal.z:10.4f} {name:<6s}    1 {name:<4s} {charge:12.6f}
@<TRIPOS>BOND
@<TRIPOS>SUBSTRUCTURE
     1 {name:<4s}       1 TEMP              0 ****  ****    0 ROOT
""",
        encoding="utf-8",
    )


def read_mol2_atom_names(path: Path) -> list[str]:
    names: list[str] = []
    in_atoms = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if line.startswith("@<TRIPOS>") and in_atoms:
            break
        if not in_atoms or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 2:
            names.append(parts[1])
    return names


def pdb_line_with_atom_name(line: str, atom_name: str) -> str:
    text = line.rstrip()
    if len(text) < 80:
        text = text.ljust(80)
    name = atom_name[:4]
    return f"{text[:12]}{name:>4}{text[16:]}".rstrip()


def metal_model_command(args: argparse.Namespace) -> None:
    project_dir = project_dir_from_args(args.project_dir)
    receptor = detect_receptor(project_dir, args.pdb)
    if receptor.suffix.lower() != ".pdb":
        raise ValueError("metal-model needs a PDB receptor/MD receptor; pass --pdb path/to/receptor_md.pdb")

    outdir = args.outdir.expanduser() if args.outdir else project_dir / "metal_model"
    if not outdir.is_absolute():
        outdir = project_dir / outdir
    outdir.mkdir(parents=True, exist_ok=True)

    receptor_atoms = read_pdb_atoms(receptor)
    ligand_atoms: list[PdbAtom] = []
    ligand_path = args.ligand_pdb.expanduser() if args.ligand_pdb else None
    if ligand_path:
        if not ligand_path.is_absolute():
            ligand_path = project_dir / ligand_path
        ligand_atoms = read_pdb_atoms(ligand_path)
    ligand_mol2_path = args.ligand_mol2.expanduser() if args.ligand_mol2 else None
    ligand_mol2_name = ""
    if ligand_mol2_path:
        if not ligand_mol2_path.is_absolute():
            ligand_mol2_path = project_dir / ligand_mol2_path
        if not ligand_mol2_path.exists():
            raise FileNotFoundError(ligand_mol2_path)
        ligand_mol2_name = f"{clean_name(args.ligand_resname.upper())}.mol2"
        shutil.copy2(ligand_mol2_path, outdir / ligand_mol2_name)
    ligand_mol2_atom_names = read_mol2_atom_names(ligand_mol2_path) if ligand_mol2_path else []
    if ligand_mol2_atom_names and ligand_atoms and len(ligand_mol2_atom_names) != len(ligand_atoms):
        raise RuntimeError(
            f"ligand mol2 atom count {len(ligand_mol2_atom_names)} does not match ligand PDB atom count {len(ligand_atoms)}"
        )
    include_ligand_in_mcpb = args.ligand_role == "coordinating"

    metal_element = args.metal_element.upper()
    metals = [atom for atom in receptor_atoms if atom.element.upper() == metal_element or atom.name.upper() == metal_element]
    if args.metal_serial is not None:
        metals = [atom for atom in metals if atom.serial == args.metal_serial]
    if args.metal_chain:
        metals = [atom for atom in metals if atom.chain == args.metal_chain]
    if args.metal_resseq:
        metals = [atom for atom in metals if atom.resseq == str(args.metal_resseq)]
    if not metals:
        raise RuntimeError(f"no {metal_element} atom matched the metal selector")
    if len(metals) > 1:
        labels = ", ".join(f"serial={atom.serial} {atom_residue_id(atom)} {atom.name}" for atom in metals[:12])
        raise RuntimeError(f"multiple {metal_element} atoms matched; add --metal-serial/--metal-chain/--metal-resseq. Matches: {labels}")
    metal = metals[0]

    coord_elements = {item.upper() for item in split_option_values([args.coord_elements])}
    all_atoms = [(atom, "receptor") for atom in receptor_atoms] + [(atom, "ligand") for atom in ligand_atoms]
    candidates: list[dict[str, object]] = []
    residue_keys: set[tuple[str, str, str, str, str]] = set()
    for atom, source in all_atoms:
        if source == "receptor" and atom.serial == metal.serial:
            continue
        if atom.element.upper().startswith("H"):
            continue
        distance = pdb_distance(metal, atom)
        if distance > args.cutoff:
            continue
        likely = atom.element.upper() in coord_elements
        candidates.append(
            {
                "source": source,
                "serial": atom.serial,
                "atom_name": atom.name,
                "element": atom.element,
                "resname": atom.resname,
                "chain": atom.chain,
                "resseq": atom.resseq,
                "icode": atom.icode,
                "residue_id": atom_residue_id(atom),
                "distance_A": round(distance, 4),
                "likely_coordination_atom": likely,
            }
        )
        if likely and (source == "receptor" or include_ligand_in_mcpb):
            residue_keys.add((source, atom.resname, atom.chain, atom.resseq, atom.icode))

    candidates.sort(key=lambda row: (float(row["distance_A"]), str(row["source"]), int(row["serial"])))

    combined_lines: list[str] = []
    serial_map: dict[tuple[str, int], int] = {}
    next_serial = 1
    for atom in receptor_atoms:
        serial_map[("receptor", atom.serial)] = next_serial
        combined_lines.append(pdb_line_with_serial(atom.line, next_serial))
        next_serial += 1
    if include_ligand_in_mcpb and ligand_atoms:
        combined_lines.append("TER")
        for ligand_idx, atom in enumerate(ligand_atoms):
            serial_map[("ligand", atom.serial)] = next_serial
            line = atom.line
            if ligand_mol2_atom_names:
                line = pdb_line_with_atom_name(line, ligand_mol2_atom_names[ligand_idx])
            combined_lines.append(pdb_line_with_serial(line, next_serial, resname=args.ligand_resname, chain=args.ligand_chain))
            next_serial += 1
    combined_lines.extend(["TER", "END"])
    original_pdb = outdir / "mcpb_original.pdb"
    original_pdb.write_text("\n".join(combined_lines) + "\n", encoding="utf-8")

    site_lines: list[str] = []
    site_next_serial = 1
    selected_receptor_residues = {(resname, chain, resseq, icode) for source, resname, chain, resseq, icode in residue_keys if source == "receptor"}
    selected_ligand_residues = {(resname, chain, resseq, icode) for source, resname, chain, resseq, icode in residue_keys if source == "ligand"}
    for atom in receptor_atoms:
        if atom.serial == metal.serial or (atom.resname, atom.chain, atom.resseq, atom.icode) in selected_receptor_residues:
            site_lines.append(pdb_line_with_serial(atom.line, site_next_serial))
            site_next_serial += 1
    if selected_ligand_residues:
        site_lines.append("TER")
        for ligand_idx, atom in enumerate(ligand_atoms):
            if (atom.resname, atom.chain, atom.resseq, atom.icode) in selected_ligand_residues:
                line = atom.line
                if ligand_mol2_atom_names:
                    line = pdb_line_with_atom_name(line, ligand_mol2_atom_names[ligand_idx])
                site_lines.append(pdb_line_with_serial(line, site_next_serial, resname=args.ligand_resname, chain=args.ligand_chain))
                site_next_serial += 1
    site_lines.extend(["TER", "END"])
    (outdir / "metal_site.pdb").write_text("\n".join(site_lines) + "\n", encoding="utf-8")

    mol2_name = f"{metal_element}.mol2"
    write_metal_mol2(outdir / mol2_name, metal, charge=float(args.metal_charge))
    mcpb_group = clean_name(args.group_name or f"{metal_element}_{metal.chain or 'chain'}{metal.resseq}_{metal.serial}")
    mcpb_input = f"""original_pdb {original_pdb.name}
group_name {mcpb_group}
cut_off {args.cutoff:g}
ion_ids {serial_map[("receptor", metal.serial)]}
ion_mol2files {mol2_name}
water_model {args.water_model}
force_field {args.force_field}
gaff {args.gaff}
software_version {args.software_version}
large_opt {args.large_opt}
add_redcrd {args.add_redcrd}
scale_factor {args.scale_factor:g}
"""
    if args.smmodel_chg is not None:
        mcpb_input += f"smmodel_chg {args.smmodel_chg}\n"
    if args.smmodel_spin is not None:
        mcpb_input += f"smmodel_spin {args.smmodel_spin}\n"
    if args.lgmodel_chg is not None:
        mcpb_input += f"lgmodel_chg {args.lgmodel_chg}\n"
    if args.lgmodel_spin is not None:
        mcpb_input += f"lgmodel_spin {args.lgmodel_spin}\n"
    if include_ligand_in_mcpb and ligand_mol2_name:
        mcpb_input += f"naa_mol2files {ligand_mol2_name}\n"
    (outdir / "mcpb.in").write_text(mcpb_input, encoding="utf-8")

    audit_model_args = ""
    if args.smmodel_chg is not None and args.smmodel_spin is not None:
        audit_model_args = f" \\\n  --model-charge {args.smmodel_chg} \\\n  --mult {args.smmodel_spin}"
    visual_model_args = ""
    if args.smmodel_chg is not None and args.smmodel_spin is not None:
        visual_model_args = f" \\\n  --model-charge {args.smmodel_chg} \\\n  --mult {args.smmodel_spin}"

    with (outdir / "metal_coordination_candidates.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "source",
            "serial",
            "atom_name",
            "element",
            "resname",
            "chain",
            "resseq",
            "icode",
            "residue_id",
            "distance_A",
            "likely_coordination_atom",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(candidates)

    run_script = f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
export AMBERHOME={shlex.quote(str(args.amberhome or (args.ambertools_bin.parent if args.ambertools_bin else Path('/mnt/l/WSL/conda_envs/AmberTools25'))))}
export PATH="$AMBERHOME/bin:$PATH"
MCPB.py -i mcpb.in -s 1
/mnt/l/WSL/softwares/conda_envs/md/bin/python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/build_mcpb_visual_check_pdb.py \\
  --small-pdb {mcpb_group}_small.pdb \\
  --mol2-dir "$PWD" \\
  --out {mcpb_group}_small_visual_check.pdb \\
  --report {mcpb_group}_small_visual_check_report.tsv{visual_model_args}
/mnt/l/WSL/softwares/conda_envs/md/bin/python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/audit_mcpb_small_model_hydrogens.py \\
  --small-pdb {mcpb_group}_small.pdb \\
  --mol2-dir "$PWD" \\
  --out mcpb_small_model_hydrogen_audit.tsv \\
  --max-warnings "${{PLE_MCPB_H_WARNING_MAX:-2}}"{audit_model_args}
"""
    (outdir / "run_mcpb_step1.sh").write_text(run_script, encoding="utf-8")
    (outdir / "run_mcpb_step1.sh").chmod(0o755)

    pipeline_script = f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
export AMBERHOME={shlex.quote(str(args.amberhome or (args.ambertools_bin.parent if args.ambertools_bin else Path('/mnt/l/WSL/conda_envs/AmberTools25'))))}
export PATH="$AMBERHOME/bin:$PATH"

# Step 1 writes small/standard/large model PDBs and QM input files.
MCPB.py -i mcpb.in -s 1
/mnt/l/WSL/softwares/conda_envs/md/bin/python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/build_mcpb_visual_check_pdb.py \\
  --small-pdb {mcpb_group}_small.pdb \\
  --mol2-dir "$PWD" \\
  --out {mcpb_group}_small_visual_check.pdb \\
  --report {mcpb_group}_small_visual_check_report.tsv{visual_model_args}
/mnt/l/WSL/softwares/conda_envs/md/bin/python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/audit_mcpb_small_model_hydrogens.py \\
  --small-pdb {mcpb_group}_small.pdb \\
  --mol2-dir "$PWD" \\
  --out mcpb_small_model_hydrogen_audit.tsv \\
  --max-warnings "${{PLE_MCPB_H_WARNING_MAX:-2}}"{audit_model_args}

# Step 2 needs completed QM output for the small-model optimization/frequency job.
# For Gaussian, run {mcpb_group}_small_opt.com and {mcpb_group}_small_fc.com first,
# then pass the frequency log/fchk back to MCPB.py, for example:
# MCPB.py -i mcpb.in -s 2 --logf {mcpb_group}_small_fc.log --fchk {mcpb_group}_small_opt.fchk

# Step 3 needs completed large-model electrostatic-potential output.
# For Gaussian, run {mcpb_group}_large_mk.com first, then:
# MCPB.py -i mcpb.in -s 3 --logf {mcpb_group}_large_mk.log

# Step 4 builds the bonded-model leap input after steps 2 and 3 succeed:
# MCPB.py -i mcpb.in -s 4
"""
    (outdir / "run_mcpb_workflow_skeleton.sh").write_text(pipeline_script, encoding="utf-8")
    (outdir / "run_mcpb_workflow_skeleton.sh").chmod(0o755)

    summary = {
        "receptor": str(receptor),
        "ligand_pdb": str(ligand_path) if ligand_path else "",
        "ligand_mol2": str(ligand_mol2_path) if ligand_mol2_path else "",
        "ligand_role": args.ligand_role,
        "metal": {
            "serial": metal.serial,
            "element": metal.element,
            "atom_name": metal.name,
            "residue_id": atom_residue_id(metal),
            "x": metal.x,
            "y": metal.y,
            "z": metal.z,
            "mcpb_original_serial": serial_map[("receptor", metal.serial)],
        },
        "cutoff_A": args.cutoff,
        "coord_elements": sorted(coord_elements),
        "candidate_count": len(candidates),
        "likely_candidate_count": sum(1 for row in candidates if row["likely_coordination_atom"]),
        "mcpb": {
            "force_field": args.force_field,
            "gaff": args.gaff,
            "software_version": args.software_version,
            "smmodel_chg": args.smmodel_chg if args.smmodel_chg is not None else "auto",
            "smmodel_spin": args.smmodel_spin,
            "lgmodel_chg": args.lgmodel_chg if args.lgmodel_chg is not None else "auto",
            "lgmodel_spin": args.lgmodel_spin,
            "large_opt": args.large_opt,
            "add_redcrd": args.add_redcrd,
            "scale_factor": args.scale_factor,
        },
        "important_caveat": "MCPB files are preparation drafts; confirm Fe oxidation state, spin state, ligand protonation, and coordination chemistry before production MD.",
    }
    (outdir / "metal_model_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    emit(f"metal_model={outdir.resolve()}")
    emit(f"metal={atom_residue_id(metal)} serial={metal.serial} mcpb_serial={serial_map[('receptor', metal.serial)]}")
    for row in candidates[:12]:
        mark = "*" if row["likely_coordination_atom"] else "-"
        emit(
            f"{mark} {row['distance_A']:.3f}A {row['source']} {row['residue_id']} "
            f"{row['atom_name']} {row['element']} serial={row['serial']}"
        )


def write_md_handoff_files(
    *,
    outdir: Path,
    row: dict[str, str],
    receptor: Path,
    args: argparse.Namespace,
) -> None:
    project_dir = project_dir_from_args(getattr(args, "project_dir", None))
    input_dir = outdir / "input"
    amber_dir = outdir / "amber"
    gmx_dir = outdir / "gromacs"
    logs_dir = outdir / "logs"
    for path in (input_dir, amber_dir, gmx_dir, logs_dir):
        path.mkdir(parents=True, exist_ok=True)

    shell_env = md_shell_env(args)
    mdrun_array = md_mdrun_array(args)
    python_exe = sys.executable
    pose_mode = max(1, int(getattr(args, "pose_mode", 1) or 1))

    receptor_copy = input_dir / "receptor.pdb"
    shutil.copy2(receptor, receptor_copy)
    receptor_md_copy = input_dir / "receptor_md.pdb"
    protonation_map = load_protonation_map(args.protonation_map, args.protonate, project_dir)
    keep_resnames = md_keep_set(args.keep_hetatm_resname, DEFAULT_MD_KEEP_RESNAMES, args.drop_default_md_metals)
    keep_elements = md_keep_set(args.keep_hetatm_element, DEFAULT_MD_KEEP_ELEMENTS, args.drop_default_md_metals)
    receptor_report = write_md_receptor_pdb(
        receptor_copy,
        receptor_md_copy,
        args.receptor_md_mode,
        protonation_map=protonation_map,
        keep_hetatm_resnames=keep_resnames,
        keep_hetatm_elements=keep_elements,
    )
    (input_dir / "receptor_md_report.json").write_text(json.dumps(receptor_report, indent=2), encoding="utf-8")
    pose_copy = copy_existing(row.get("pose_path", ""), input_dir / "ligand_pose.pdbqt")
    pdbqt_copy = copy_existing(row.get("pdbqt_path", ""), input_dir / "ligand_prepared.pdbqt")
    sdf_copy = copy_existing(row.get("sdf_path", ""), input_dir / "ligand_input.sdf")
    (input_dir / "ligand.smi").write_text(row.get("smiles", "") + "\n", encoding="utf-8")
    (outdir / "selected_candidate.json").write_text(json.dumps(row, indent=2), encoding="utf-8")

    nsteps = max(1, int(round(float(args.production_ns) * 1000.0 / float(args.dt_ps))))
    md_ns_label = ns_label(float(args.production_ns))
    production_mdp = f"md_{md_ns_label}ns.mdp"
    custom_leap_lines, custom_param_report = copy_custom_amber_params(args, amber_dir, project_dir)
    (amber_dir / "custom_param_report.json").write_text(json.dumps(custom_param_report, indent=2), encoding="utf-8")
    custom_leap = "\n".join(custom_leap_lines)
    if custom_leap:
        custom_leap += "\n"
    leap = f"""source {args.protein_forcefield}
source {args.ligand_forcefield}
source leaprc.water.tip3p
{custom_leap}

LIG = loadmol2 ligand.mol2
loadamberparams ligand.frcmod
REC = loadpdb receptor_md.pdb
COMPLEX = combine {{ REC LIG }}
solvateBox COMPLEX {args.water_box} {args.box_padding:.3f}
addIonsRand COMPLEX Na+ 0
addIonsRand COMPLEX Cl- 0
saveamberparm COMPLEX complex.prmtop complex.inpcrd
savepdb COMPLEX complex_amber.pdb
quit
"""
    (amber_dir / "leap.in").write_text(leap, encoding="utf-8")

    (outdir / "build_ligand_pose_sdf.py").write_text(
        """from pathlib import Path
import math
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D


def pdbqt_element(line):
    atom_type = (line[70:].split() or [""])[-1]
    if atom_type in {"A", "C"}:
        return "C"
    if atom_type.startswith("O"):
        return "O"
    if atom_type.startswith("N"):
        return "N"
    if atom_type.startswith("H"):
        return "H"
    name = line[12:16].strip()
    letters = "".join(ch for ch in name if ch.isalpha())
    if letters:
        two = letters[:2].capitalize()
        if two in {"Cl", "Br", "Fe", "Zn", "Mg", "Ca", "Na", "Mn", "Cu", "Co", "Ni"}:
            return two
        return letters[0].upper()
    return atom_type[:1].upper()


def read_pdbqt_atoms(path, model_index=1):
    atoms = []
    serial_to_mol_index = {}
    model_count = 0
    in_target_model = False
    saw_model = False
    for line in Path(path).read_text(errors="replace").splitlines():
        if line.startswith("MODEL"):
            saw_model = True
            try:
                model_count = int(line.split()[1])
            except Exception:
                model_count += 1
            in_target_model = model_count == model_index
            continue
        if line.startswith("ENDMDL"):
            if in_target_model:
                break
            in_target_model = False
            continue
        if saw_model and not in_target_model:
            continue
        if line.startswith("REMARK SMILES IDX"):
            values = [int(item) for item in line.split()[3:]]
            for idx in range(0, len(values) - 1, 2):
                mol_index = values[idx] - 1
                serial = values[idx + 1]
                serial_to_mol_index[serial] = mol_index
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        serial = int(line[6:11])
        element = pdbqt_element(line)
        if element.upper().startswith("H"):
            continue
        atoms.append(
            {
                "serial": serial,
                "element": element,
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "name": line[12:16].strip(),
            }
        )
    if saw_model and model_index > 1 and not atoms:
        raise RuntimeError(f"MODEL {model_index} was not found in {path}")
    return atoms, serial_to_mol_index


def dist3(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def heavy_positions(mol):
    conf = mol.GetConformer()
    rows = []
    heavy_idx = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol().upper().startswith("H"):
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        rows.append((heavy_idx, atom.GetSymbol(), (pos.x, pos.y, pos.z)))
        heavy_idx += 1
    return rows


def map_prepared_serials_to_sdf(prepared_atoms, mol):
    sdf_atoms = heavy_positions(mol)
    unused = set(range(len(sdf_atoms)))
    mapping = {}
    report = []
    for pdbqt_atom in prepared_atoms:
        candidates = []
        for row_index in unused:
            mol_idx, symbol, coords = sdf_atoms[row_index]
            if symbol.upper() != pdbqt_atom["element"].upper():
                continue
            d = dist3(coords, (pdbqt_atom["x"], pdbqt_atom["y"], pdbqt_atom["z"]))
            candidates.append((d, row_index, mol_idx, symbol))
        if not candidates:
            raise RuntimeError(
                f"could not map prepared PDBQT serial {pdbqt_atom['serial']} element {pdbqt_atom['element']} to SDF"
            )
        d, row_index, mol_idx, symbol = min(candidates)
        unused.remove(row_index)
        mapping[pdbqt_atom["serial"]] = mol_idx
        report.append((pdbqt_atom["serial"], mol_idx, symbol, pdbqt_atom["element"], d))
    if unused:
        missing = ",".join(str(sdf_atoms[row_index][0] + 1) for row_index in sorted(unused))
        raise RuntimeError(f"prepared PDBQT did not map all SDF heavy atoms; missing {missing}")
    return mapping, report


def map_serials_with_remark(serial_to_mol_index):
    return {serial: mol_index for serial, mol_index in serial_to_mol_index.items()}


def heavy_bond_lengths(mol):
    conf = mol.GetConformer()
    rows = []
    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        pa = conf.GetAtomPosition(a)
        pb = conf.GetAtomPosition(b)
        d = math.sqrt((pa.x - pb.x) ** 2 + (pa.y - pb.y) ** 2 + (pa.z - pb.z) ** 2)
        rows.append((a + 1, b + 1, d))
    return rows


def main():
    pose_mode = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    input_sdf = Path("input/ligand_input.sdf")
    pose_pdbqt = Path("input/ligand_pose.pdbqt")
    prepared_pdbqt = Path("input/ligand_prepared.pdbqt")
    output_sdf = Path("amber/ligand_pose.sdf")
    mapping_tsv = Path("amber/ligand_pose_mapping.tsv")
    if not input_sdf.exists() or not pose_pdbqt.exists():
        return 2
    supplier = Chem.SDMolSupplier(str(input_sdf), removeHs=False)
    mol = next((item for item in supplier if item is not None), None)
    if mol is None:
        raise RuntimeError(f"could not read {input_sdf}")
    heavy = Chem.RemoveHs(mol)
    pose_atoms, serial_to_mol_index = read_pdbqt_atoms(pose_pdbqt, pose_mode)
    pose_by_serial = {item["serial"]: item for item in pose_atoms}
    if len(pose_by_serial) != heavy.GetNumAtoms():
        raise RuntimeError(
            f"docked heavy atom count {len(pose_by_serial)} does not match SDF heavy atom count {heavy.GetNumAtoms()}"
        )
    mapping_report = []
    if prepared_pdbqt.exists():
        prepared_atoms, _ = read_pdbqt_atoms(prepared_pdbqt, 1)
        if len(prepared_atoms) != heavy.GetNumAtoms():
            raise RuntimeError(
                f"prepared heavy atom count {len(prepared_atoms)} does not match SDF heavy atom count {heavy.GetNumAtoms()}"
            )
        serial_to_mol_index, mapping_report = map_prepared_serials_to_sdf(prepared_atoms, mol)
    elif serial_to_mol_index:
        serial_to_mol_index = map_serials_with_remark(serial_to_mol_index)
    else:
        print("pose_sdf_mapping_warning=no prepared PDBQT or REMARK SMILES IDX; falling back to PDBQT serial order", file=sys.stderr)
        serial_to_mol_index = {
            serial: idx for idx, serial in enumerate(sorted(pose_by_serial))
        }
    if heavy.GetNumConformers() == 0:
        conf = Chem.Conformer(heavy.GetNumAtoms())
        heavy.AddConformer(conf)
    conf = heavy.GetConformer()
    mismatches = []
    assignments = [None] * heavy.GetNumAtoms()
    for serial, mol_index in serial_to_mol_index.items():
        if serial not in pose_by_serial:
            continue
        if 0 <= mol_index < heavy.GetNumAtoms():
            assignments[mol_index] = pose_by_serial[serial]
    missing = [str(idx + 1) for idx, item in enumerate(assignments) if item is None]
    if missing:
        raise RuntimeError("missing docked coordinates for heavy atom indices: " + ",".join(missing))
    for idx, item in enumerate(assignments):
        x, y, z, element = item["x"], item["y"], item["z"], item["element"]
        atom = heavy.GetAtomWithIdx(idx)
        if element and atom.GetSymbol().upper() != element.upper():
            mismatches.append(f"{idx + 1}:{atom.GetSymbol()}!={element}")
        conf.SetAtomPosition(idx, Point3D(x, y, z))
    if mismatches:
        print("pose_sdf_element_warnings=" + ",".join(mismatches[:12]), file=sys.stderr)
    bond_lengths = heavy_bond_lengths(heavy)
    max_bond = max((row[2] for row in bond_lengths), default=0.0)
    print(f"pose_sdf_mode={pose_mode}")
    print(f"pose_sdf_heavy_bond_max_A={max_bond:.3f}")
    long_bonds = [row for row in bond_lengths if row[2] > 2.1]
    if long_bonds:
        text = ",".join(f"{a}-{b}:{d:.3f}" for a, b, d in long_bonds[:12])
        print(f"pose_sdf_long_bond_warning={text}", file=sys.stderr)
    mapping_lines = ["pdbqt_serial\\tsdf_atom_index\\tsdf_element\\tpdbqt_element\\tprepared_distance_A"]
    if mapping_report:
        for serial, mol_idx, sdf_element, pdbqt_element_value, d in mapping_report:
            mapping_lines.append(f"{serial}\\t{mol_idx + 1}\\t{sdf_element}\\t{pdbqt_element_value}\\t{d:.4f}")
    else:
        for serial, mol_index in sorted(serial_to_mol_index.items()):
            element = pose_by_serial.get(serial, {}).get("element", "")
            symbol = heavy.GetAtomWithIdx(mol_index).GetSymbol() if 0 <= mol_index < heavy.GetNumAtoms() else ""
            mapping_lines.append(f"{serial}\\t{mol_index + 1}\\t{symbol}\\t{element}\\t")
    mapping_tsv.write_text("\\n".join(mapping_lines) + "\\n", encoding="utf-8")
    posed = Chem.AddHs(heavy, addCoords=True)
    try:
        ff = AllChem.UFFGetMoleculeForceField(posed, confId=0)
        for idx in range(heavy.GetNumAtoms()):
            ff.AddFixedPoint(idx)
        ff.Minimize(maxIts=200)
    except Exception as exc:
        print(f"pose_sdf_hydrogen_minimize_warning={type(exc).__name__}:{exc}", file=sys.stderr)
    writer = Chem.SDWriter(str(output_sdf))
    writer.write(posed)
    writer.close()
    print(f"pose_sdf_from_docked_heavy_atoms={output_sdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
""",
        encoding="utf-8",
    )

    (outdir / "build_posres_topology.py").write_text(
        """from pathlib import Path
import sys


def is_section(line, name=None):
    text = line.strip()
    if not (text.startswith("[") and text.endswith("]")):
        return False
    section = text.strip("[]").strip().lower()
    return section == name if name else True


def parse_first_molecule_atoms(lines):
    molecule_indices = [idx for idx, line in enumerate(lines) if is_section(line, "moleculetype")]
    if not molecule_indices:
        raise RuntimeError("topology has no [ moleculetype ] section")
    start = molecule_indices[0]
    end = molecule_indices[1] if len(molecule_indices) > 1 else len(lines)
    atoms_section = None
    for idx in range(start, end):
        if is_section(lines[idx], "atoms"):
            atoms_section = idx
            break
    if atoms_section is None:
        raise RuntimeError("first molecule type has no [ atoms ] section")

    atom_numbers = []
    for line in lines[atoms_section + 1:end]:
        if is_section(line):
            break
        fields = line.split(";", 1)[0].split()
        if len(fields) < 8:
            continue
        try:
            atom_number = int(fields[0])
            atom_mass = float(fields[7])
        except ValueError:
            continue
        atom_name = fields[4].upper()
        if atom_mass > 1.5 and not atom_name.startswith("H"):
            atom_numbers.append(atom_number)
    if not atom_numbers:
        raise RuntimeError("no heavy atoms found for position restraints")
    return atom_numbers, end


def main():
    top = Path(sys.argv[1])
    force = float(sys.argv[2]) if len(sys.argv) > 2 else 1000.0
    lines = top.read_text(errors="replace").splitlines()
    atom_numbers, insert_at = parse_first_molecule_atoms(lines)
    posre = top.with_name("posre_complex.itp")
    posre_lines = [
        "[ position_restraints ]",
        "; atom  type      fx      fy      fz",
    ]
    for atom_number in atom_numbers:
        posre_lines.append(f"{atom_number:6d}     1  {force:7.1f} {force:7.1f} {force:7.1f}")
    posre.write_text("\\n".join(posre_lines) + "\\n", encoding="utf-8")

    include_line = '#include "posre_complex.itp"'
    if include_line not in lines:
        insert = [
            "",
            "#ifdef POSRES",
            include_line,
            "#endif",
            "",
        ]
        lines[insert_at:insert_at] = insert
        top.write_text("\\n".join(lines) + "\\n", encoding="utf-8")
    print(f"position_restraints={posre} heavy_atoms={len(atom_numbers)} force={force:g}")


if __name__ == "__main__":
    raise SystemExit(main())
""",
        encoding="utf-8",
    )

    (gmx_dir / "minim.mdp").write_text(
        """integrator               = steep
emtol                    = 1000.0
emstep                   = 0.01
nsteps                   = 50000
""" + mdp_common(),
        encoding="utf-8",
    )
    (gmx_dir / "nvt.mdp").write_text(
        f"""define                   = -DPOSRES
integrator               = md
dt                       = {args.dt_ps}
nsteps                   = 50000
nstxout-compressed       = 5000
nstenergy                = 1000
nstlog                   = 1000
continuation             = no
gen-vel                  = yes
gen-temp                 = 300
gen-seed                 = -1
tcoupl                   = V-rescale
tc-grps                  = System
tau-t                    = 0.1
ref-t                    = 300
pcoupl                   = no
""" + mdp_common(),
        encoding="utf-8",
    )
    (gmx_dir / "npt.mdp").write_text(
        f"""define                   = -DPOSRES
integrator               = md
dt                       = {args.dt_ps}
nsteps                   = 250000
nstxout-compressed       = 5000
nstenergy                = 1000
nstlog                   = 1000
continuation             = yes
gen-vel                  = no
tcoupl                   = V-rescale
tc-grps                  = System
tau-t                    = 0.1
ref-t                    = 300
pcoupl                   = C-rescale
pcoupltype               = isotropic
tau-p                    = 2.0
ref-p                    = 1.0
compressibility          = 4.5e-5
refcoord_scaling         = com
""" + mdp_common(),
        encoding="utf-8",
    )
    (gmx_dir / production_mdp).write_text(
        f"""integrator               = md
dt                       = {args.dt_ps}
nsteps                   = {nsteps}
nstxout-compressed       = 5000
nstenergy                = 1000
nstlog                   = 1000
continuation             = yes
gen-vel                  = no
tcoupl                   = V-rescale
tc-grps                  = System
tau-t                    = 0.1
ref-t                    = 300
pcoupl                   = Parrinello-Rahman
pcoupltype               = isotropic
tau-p                    = 2.0
ref-p                    = 1.0
compressibility          = 4.5e-5
""" + mdp_common(),
        encoding="utf-8",
    )

    check_script = f"""#!/usr/bin/env bash
set -u
{shell_env}
echo "PLE MD handoff environment check"
GMX={shlex.quote(str(args.gmx))}
PYTHON={shlex.quote(str(python_exe))}
commands=("$GMX" obabel antechamber parmchk2 tleap acpype)
missing=0
for cmd in "${{commands[@]}}"; do
  printf '%s\\t' "$cmd"
  if command -v "$cmd" >/dev/null 2>&1; then
    command -v "$cmd"
  else
    echo "missing"
    missing=1
  fi
done
echo "gmx_version_probe"
"$GMX" --version 2>&1 | sed -n '1,8p' || true
"$PYTHON" - <<'PY'
mods = ["MDAnalysis", "mdtraj", "matplotlib", "parmed", "rdkit"]
for mod in mods:
    try:
        __import__(mod)
        print(f"{{mod}} ok")
    except Exception as exc:
        print(f"{{mod}} missing: {{type(exc).__name__}}: {{exc}}")
PY
if [[ "$missing" -ne 0 ]]; then
  echo "install_hint=conda install -c conda-forge ambertools acpype parmed matplotlib"
fi
"""
    (outdir / "00_check_env.sh").write_text(check_script, encoding="utf-8")

    prep_script = f"""#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
cd "$ROOT"
mkdir -p logs amber gromacs
exec > >(tee -a logs/prepare_amber_to_gmx.log) 2>&1

{shell_env}
GMX="{args.gmx}"
PYTHON={shlex.quote(str(python_exe))}
LIG_CHARGE="{args.ligand_charge}"
LIG_RESNAME="{args.ligand_resname}"

echo "[1/5] Preparing ligand coordinates"
if [[ -s input/ligand_input.sdf && -s input/ligand_pose.pdbqt ]]; then
  if "$PYTHON" build_ligand_pose_sdf.py "{pose_mode}"; then
    true
  else
    echo "pose_sdf_fallback=original_sdf"
    cp input/ligand_input.sdf amber/ligand_pose.sdf
  fi
  obabel amber/ligand_pose.sdf -O amber/ligand_pose.pdb --title "$LIG_RESNAME" || true
  obabel input/ligand_pose.pdbqt -O amber/docked_pose_reference.pdb -f {pose_mode} -l {pose_mode} --title "$LIG_RESNAME" || true
elif [[ -s input/ligand_input.sdf ]]; then
  cp input/ligand_input.sdf amber/ligand_pose.sdf
  obabel amber/ligand_pose.sdf -O amber/ligand_pose.pdb --title "$LIG_RESNAME" || true
elif [[ -s input/ligand_pose.pdbqt ]]; then
  obabel input/ligand_pose.pdbqt -O amber/ligand_pose.sdf -f {pose_mode} -l {pose_mode} -h --title "$LIG_RESNAME"
  obabel input/ligand_pose.pdbqt -O amber/ligand_pose.pdb -f {pose_mode} -l {pose_mode} --title "$LIG_RESNAME"
else
  echo "No ligand_pose.pdbqt or ligand_input.sdf found" >&2
  exit 2
fi

echo "[2/5] AmberTools ligand parameters"
cd amber
antechamber -i ligand_pose.sdf -fi sdf -o ligand.mol2 -fo mol2 -c bcc -nc "$LIG_CHARGE" -rn "$LIG_RESNAME" -s 2
parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod

echo "[3/5] Receptor PDB copy"
cp ../input/receptor_md.pdb receptor_md.pdb

echo "[4/5] tleap complex build"
tleap -f leap.in

echo "[5/5] ACPYPE Amber to GROMACS conversion"
acpype -p complex.prmtop -x complex.inpcrd -b complex
cd ..
GMX_GRO="$(find amber -name '*GMX.gro' | head -n 1)"
GMX_TOP="$(find amber -name '*GMX.top' | head -n 1)"
if [[ -z "$GMX_GRO" || -z "$GMX_TOP" ]]; then
  echo "ACPYPE did not produce *GMX.gro/*GMX.top" >&2
  exit 3
fi
cp "$GMX_GRO" gromacs/complex_GMX.gro
cp "$GMX_TOP" gromacs/complex_GMX.top
find amber -name '*.itp' -exec cp {{}} gromacs/ \\;
"$PYTHON" build_posres_topology.py gromacs/complex_GMX.top 1000
echo "Prepared gromacs/complex_GMX.gro and gromacs/complex_GMX.top"
"""
    (outdir / "01_prepare_amber_to_gmx.sh").write_text(prep_script, encoding="utf-8")

    run_script = f"""#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
cd "$ROOT/gromacs"
mkdir -p ../logs
exec > >(tee -a ../logs/run_gromacs.log) 2>&1
{shell_env}
GMX="{args.gmx}"
MDRUN_ARGS=({mdrun_array})

"$GMX" grompp -f minim.mdp -c complex_GMX.gro -p complex_GMX.top -o minim.tpr -maxwarn {args.maxwarn} -v
"$GMX" mdrun -deffnm minim -v
"$GMX" grompp -f nvt.mdp -c minim.gro -r minim.gro -p complex_GMX.top -o nvt.tpr -maxwarn {args.maxwarn}
"$GMX" mdrun -deffnm nvt -v "${{MDRUN_ARGS[@]}}"
"$GMX" grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p complex_GMX.top -o npt.tpr -maxwarn {args.maxwarn}
"$GMX" mdrun -deffnm npt -v "${{MDRUN_ARGS[@]}}"
"$GMX" grompp -f {production_mdp} -c npt.gro -t npt.cpt -p complex_GMX.top -o md.tpr -maxwarn {args.maxwarn}
"$GMX" mdrun -deffnm md -v "${{MDRUN_ARGS[@]}}"
"""
    (outdir / "02_run_gromacs.sh").write_text(run_script, encoding="utf-8")

    analyze_script = f"""#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
cd "$ROOT/gromacs"
mkdir -p ../analysis ../logs
exec > >(tee -a ../logs/analyze_basic.log) 2>&1
{shell_env}
GMX="{args.gmx}"
PYTHON={shlex.quote(str(python_exe))}
printf "System\\n" | "$GMX" trjconv -s md.tpr -f md.xtc -o ../analysis/md_nojump.xtc -pbc nojump
printf "Backbone\\nSystem\\n" | "$GMX" trjconv -s md.tpr -f ../analysis/md_nojump.xtc -o ../analysis/md_fit.xtc -fit rot+trans
printf "Backbone\\nBackbone\\n" | "$GMX" rms -s md.tpr -f ../analysis/md_fit.xtc -o ../analysis/rmsd.xvg -tu ns
printf "Protein\\n" | "$GMX" gyrate -s md.tpr -f ../analysis/md_fit.xtc -o ../analysis/rg.xvg || true
cd ../analysis
"$PYTHON" plot_fel_rmsd_rg.py
"""
    (outdir / "03_analyze_basic.sh").write_text(analyze_script, encoding="utf-8")
    (outdir / "analysis").mkdir(exist_ok=True)
    (outdir / "analysis" / "plot_fel_rmsd_rg.py").write_text(
        """from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def load_xvg(path):
    rows = []
    for line in Path(path).read_text(errors="replace").splitlines():
        if not line or line[0] in "#@":
            continue
        parts = line.split()
        if len(parts) >= 2:
            rows.append((float(parts[0]), float(parts[1])))
    return np.asarray(rows, float)

rmsd = load_xvg("rmsd.xvg")
rg = load_xvg("rg.xvg")
n = min(len(rmsd), len(rg))
if n == 0:
    raise SystemExit("missing rmsd.xvg or rg.xvg data")
x = rmsd[:n, 1]
y = rg[:n, 1]
hist, xedges, yedges = np.histogram2d(x, y, bins=60)
prob = hist / max(hist.sum(), 1.0)
free = -np.log(prob + 1e-12)
free -= np.nanmin(free)
plt.figure(figsize=(6, 5))
plt.imshow(
    free.T,
    origin="lower",
    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
    aspect="auto",
    cmap="viridis",
)
plt.colorbar(label="Relative free energy (kT)")
plt.xlabel("Backbone RMSD (nm)")
plt.ylabel("Radius of gyration (nm)")
plt.tight_layout()
plt.savefig("fel_rmsd_rg.png", dpi=180)
print("wrote fel_rmsd_rg.png")
""",
        encoding="utf-8",
    )

    for script in ("00_check_env.sh", "01_prepare_amber_to_gmx.sh", "02_run_gromacs.sh", "03_analyze_basic.sh"):
        (outdir / script).chmod(0o755)

    ambertools_text = str(getattr(args, "ambertools_bin", "") or "")
    gmxrc_text = str(getattr(args, "gmxrc", "") or "")
    ld_text = ", ".join(str(path) for path in (getattr(args, "extra_ld_library_path", []) or []))
    mdrun_text = md_mdrun_array(args)
    dropped_counts = receptor_report.get("dropped_counts", {})
    dropped_text = ", ".join(f"{key}={value}" for key, value in sorted(dropped_counts.items())) if dropped_counts else "none"
    keep_residue_text = ", ".join(receptor_report.get("kept_hetatm_resnames", [])) or "none"
    keep_element_text = ", ".join(receptor_report.get("kept_hetatm_elements", [])) or "none"
    rename_report = receptor_report.get("manual_residue_renames", [])
    rename_text = ", ".join(
        f"{item.get('chain') or '-'}:{item.get('resseq')}{item.get('icode') or ''} {item.get('from')}->{item.get('to')}"
        for item in rename_report
    ) or "none"
    custom_params = json.loads((amber_dir / "custom_param_report.json").read_text(encoding="utf-8"))
    custom_param_text = ", ".join(str(path) for values in custom_params.values() for path in values) or "none"
    readme = f"""# PLE MD Handoff

Selected candidate: `{row.get('seq_id')}`

- affinity_kcal_mol: `{row.get('affinity_kcal_mol')}`
- cnn_pose_score: `{row.get('cnn_pose_score')}`
- inner_rmsd: `{row.get('inner_rmsd')}`
- smiles: `{row.get('smiles')}`
- selected docked pose mode: `{pose_mode}`

Environment written into the generated scripts:

- AmberTools bin: `{ambertools_text or 'PATH'}`
- GROMACS: `{args.gmx}`
- GMXRC: `{gmxrc_text or 'auto/none'}`
- extra LD library path: `{ld_text or 'none'}`
- mdrun args: `{mdrun_text or 'none'}`
- receptor MD mode: `{args.receptor_md_mode}`
- production MD time: `{args.production_ns:g} ns`
- production MDP: `{production_mdp}`
- receptor atoms kept: `{receptor_report.get('kept_atoms')}`
- kept HETATM residue names: `{keep_residue_text}`
- kept HETATM elements: `{keep_element_text}`
- manual protonation/residue renames: `{rename_text}`
- custom Amber parameter files: `{custom_param_text}`
- receptor records dropped from MD copy: `{dropped_text}`

Run from WSL inside the conda MD environment:

```bash
cd {outdir}
./00_check_env.sh
./01_prepare_amber_to_gmx.sh
./02_run_gromacs.sh
./03_analyze_basic.sh
```

If the check script reports missing AmberTools/ACPYPE commands, install them in
the MD environment with:

```bash
conda install -c conda-forge ambertools acpype parmed matplotlib
```

The critical conversion step copied from the working notebook is:

```bash
acpype -p complex.prmtop -x complex.inpcrd -b complex
```

Expected hard failure points are AmberTools/tleap residue naming, ligand charge,
metal/heme/nonstandard residue parameters, or GROMACS topology warnings. Keep the
first failing log under `logs/`; it is the next thing to fix.

By default, `ple md-handoff` uses `--receptor-md-mode protein-only`: the original
PDB is kept as `input/receptor.pdb`, while `input/receptor_md.pdb` removes waters,
heme-like cofactors, phosphate, and other HETATM records so a first MD smoke test
can reach GROMACS. Fe/FE records are kept by default because catalytic iron is
not a disposable solvent component. Use `--drop-default-md-metals` only for a
temporary protein-only smoke test.

Manual H++/curation residue names can be applied with either repeated one-line
rules:

```bash
ple md-handoff --protonate A:123=HIE --protonate A:45=ASH
```

or a CSV/line file via `--protonation-map`. Amber names should be used directly:
`ASH`, `GLH`, `HID`, `HIE`, `HIP`, `CYM`, `CYX`, and related residue names.

The generated `leap.in` loads `frcmod.ions234lm_126_tip3p` by default, matching
the earlier working notebook style. That helps ordinary ions, but it is not a
complete model for a catalytic iron center with defined oxidation state, spin
state, and coordinating residues. For a real metalloenzyme MD system, add
literature/MCPB.py-derived files with `--amber-frcmod`, `--amber-lib`,
`--amber-prep`, or `--amber-mol2`, then rerun `./01_prepare_amber_to_gmx.sh`.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")
    print(f"md_handoff={outdir}", flush=True)
    print(f"selected_seq_id={row.get('seq_id')} affinity={row.get('affinity_kcal_mol')} cnn_pose={row.get('cnn_pose_score')}", flush=True)
    if pose_copy is None:
        print("warning=no_pose_pdbqt_copied", flush=True)
    if pdbqt_copy is None:
        print("warning=no_prepared_ligand_pdbqt_copied", flush=True)
    if sdf_copy is None:
        print("warning=no_input_sdf_copied", flush=True)


def md_handoff_command(args: argparse.Namespace) -> None:
    project_dir = project_dir_from_args(args.project_dir)
    receptor = detect_receptor(project_dir, args.pdb)
    if receptor.suffix.lower() != ".pdb":
        raise ValueError("MD handoff needs a receptor PDB; pass --pdb path/to/receptor.pdb")
    state_dir = project_dir / args.state_dir
    history_csv = state_dir / "dock_history" / "dock_history.csv"
    row = select_md_row(history_csv, args.seq_id, args.rank_by)
    outdir = args.outdir.expanduser() if args.outdir else state_dir / "md_handoff" / str(row.get("seq_id") or "selected")
    if not outdir.is_absolute():
        outdir = project_dir / outdir
    outdir.mkdir(parents=True, exist_ok=True)
    write_md_handoff_files(outdir=outdir.resolve(), row=row, receptor=receptor, args=args)


def write_md_postprocess_script(args: argparse.Namespace, handoff_dir: Path, project_dir: Path) -> Path:
    gmx_dir = handoff_dir / "gromacs"
    analysis_dir = project_dir / "md_analysis" / "rdc"
    pymol_dir = project_dir / "pymol"
    script = handoff_dir / "03_postprocess.sh"
    env_lines = md_shell_env(args)
    gmxrc_arg = f" --gmxrc {shlex.quote(str(args.gmxrc))}" if getattr(args, "gmxrc", None) else ""
    ld_args = "".join(
        f" --extra-ld-library-path {shlex.quote(str(path))}"
        for path in (getattr(args, "extra_ld_library_path", []) or [])
    )
    ref_pdb = handoff_dir / "input" / "receptor_md.pdb"
    script.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
cd "$ROOT"
{env_lines}

echo "[PLE] fixing trajectory PBC: raw md.xtc stays untouched; writing md_nojump.xtc and md_centered_compact.xtc/gro"
ple md-fix-pbc --gmxdir {shlex.quote(str(gmx_dir))} --gmx {shlex.quote(str(args.gmx))}{gmxrc_arg}{ld_args} --center-group {shlex.quote(str(args.center_group))} --output-group {shlex.quote(str(args.output_group))}
ple pymol-metal --structure {shlex.quote(str(gmx_dir / "md_centered_compact.gro"))} --trajectory {shlex.quote(str(gmx_dir / "md_centered_compact.xtc"))} --out {shlex.quote(str(pymol_dir / "show_metal_centered.pml"))} --radius {args.metal_radius:g} --carbon-color {shlex.quote(str(args.carbon_color))}
ple md-rdc --gmxdir {shlex.quote(str(gmx_dir))} --outdir {shlex.quote(str(analysis_dir))} --topology {shlex.quote(str(gmx_dir / "md_centered_compact.gro"))} --trajectory {shlex.quote(str(gmx_dir / "md_centered_compact.xtc"))} --gro {shlex.quote(str(gmx_dir / "md_centered_compact.gro"))} --reference-pdb {shlex.quote(str(ref_pdb))} --ligand-selection {shlex.quote(str(args.ligand_selection))} --distance-cutoff {args.distance_cutoff:g} --stride {args.rdc_stride:d} --top-n {args.top_n:d} --velocity-threshold {args.velocity_threshold:g}
""",
        encoding="utf-8",
    )
    script.chmod(0o755)
    return script


def md_run_command(args: argparse.Namespace) -> None:
    """Prepare and optionally run the flat PLE MD workflow."""
    project_dir = project_dir_from_args(args.project_dir)
    receptor = detect_receptor(project_dir, args.pdb)
    if receptor.suffix.lower() != ".pdb":
        raise ValueError("md-run needs a receptor PDB; pass --pdb path/to/receptor.pdb")
    state_dir = project_dir / args.state_dir
    history_csv = state_dir / "dock_history" / "dock_history.csv"
    row = select_md_row(history_csv, args.seq_id, args.rank_by)
    outdir = args.outdir.expanduser() if args.outdir else project_dir / "md_handoff"
    if not outdir.is_absolute():
        outdir = project_dir / outdir
    outdir.mkdir(parents=True, exist_ok=True)

    write_md_handoff_files(outdir=outdir.resolve(), row=row, receptor=receptor, args=args)
    postprocess_script = write_md_postprocess_script(args, outdir.resolve(), project_dir)

    print(f"md_handoff={outdir.resolve()}", flush=True)
    print(f"postprocess_script={postprocess_script}", flush=True)
    if args.prepare_only:
        print("status=prepared_only", flush=True)
        return

    for script_name in ("00_check_env.sh", "01_prepare_amber_to_gmx.sh"):
        run(["bash", str(outdir / script_name)], dry_run=args.dry_run)

    run_script = outdir / "02_run_gromacs.sh"
    background_script = run_script
    if not args.skip_postprocess:
        background_script = outdir / "04_run_gromacs_then_postprocess.sh"
        background_script.write_text(
            f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
bash ./02_run_gromacs.sh
bash {shlex.quote(postprocess_script.name)}
""",
            encoding="utf-8",
        )
        background_script.chmod(0o755)
    if args.background:
        logs_dir = outdir / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        if args.dry_run:
            print("$ nohup bash " + shlex.quote(str(background_script)) + " > " + shlex.quote(str(logs_dir / "run_gromacs_background.log")) + " 2>&1 &", flush=True)
            return
        handle = (logs_dir / "run_gromacs_background.log").open("ab")
        proc = subprocess.Popen(["bash", str(background_script)], stdout=handle, stderr=subprocess.STDOUT)
        (logs_dir / "md.pid").write_text(str(proc.pid) + "\n", encoding="utf-8")
        print(f"started_pid={proc.pid}", flush=True)
        print(f"monitor=tail -f {logs_dir / 'run_gromacs_background.log'}", flush=True)
        return

    run(["bash", str(run_script)], dry_run=args.dry_run)
    if not args.skip_postprocess:
        run(["bash", str(postprocess_script)], dry_run=args.dry_run)


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    raw_args = list(argv if argv is not None else sys.argv[1:])
    if raw_args and raw_args[0].startswith("-"):
        raw_args = ["run", *raw_args]
    args = parser.parse_args(raw_args)
    if args.command == "init":
        init_project(args)
        return
    if args.command == "prepare-receptor":
        prepare_receptor_command(args)
        return
    if args.command == "repair-loops":
        repair_loops_command(args)
        return
    if args.command == "hpp-map":
        hpp_map_command(args)
        return
    if args.command == "metal-model":
        metal_model_command(args)
        return
    if args.command == "focus-box":
        focus_box_command(args)
        return
    if args.command == "md-handoff":
        md_handoff_command(args)
        return
    if args.command == "md-run":
        md_run_command(args)
        return
    if args.command == "md-fix-pbc":
        fix_pbc_command(args)
        return
    if args.command == "md-rdc":
        rdc_command(args)
        return
    if args.command == "pymol-metal":
        pymol_metal_command(args)
        return
    if args.command == "run":
        run_workflow(args)
        return
    if args.command is None:
        parser.print_help()
        return
    parser.error(f"unknown command: {args.command}")


if __name__ == "__main__":
    main()
