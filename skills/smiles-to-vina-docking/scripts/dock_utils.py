#!/usr/bin/env python
"""Simple chronological SMILES -> Vina docking history.

This script collects facts. It does not rank molecules against a global lead.
Each molecule is recorded by time sequence, ancestor_smiles, parent_smiles,
SMILES identities, Vina affinity, pose stability, simple interaction counts,
and file paths.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, QED, rdMolDescriptors

try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError:  # pragma: no cover
    rdMolStandardize = None  # type: ignore[assignment]


RDLogger.DisableLog("rdApp.warning")
RDLogger.DisableLog("rdApp.error")

LEDGER_FILENAME = "dock_history.csv"
CONFIG_FILENAME = "history_config.json"
DEFAULT_TASK_ROOT = Path("~/vina_task2")
DEFAULT_LEDGER_DIR = DEFAULT_TASK_ROOT / "dock_history"
AFFINITY_RE = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+", re.MULTILINE)
PDBQT_AFFINITY_RE = re.compile(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)")
HBOND_ELEMENTS = {"N", "O", "S"}
HYDROPHOBIC_ELEMENTS = {"C", "Cl", "Br", "I", "F"}
AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS"}

LEDGER_FIELDS = [
    "seq_id",
    "timestamp",
    "nickname",
    "ancestor_smiles",
    "parent_smiles",
    "edit_label",
    "refinement_source",
    "druglike_refinement_score",
    "refinement_qed",
    "synthetic_score_proxy",
    "property_window_score",
    "qed_component_score",
    "synthetic_component_score",
    "property_component_score",
    "parent_similarity_score",
    "expert_similarity_score",
    "reference_similarity_score",
    "reference_partial_similarity_score",
    "qve_component_score",
    "inherited_structure_score",
    "structural_alert_penalty_score",
    "refinement_parent_similarity",
    "refinement_expert_similarity",
    "kinase_reference_similarity",
    "kinase_reference_partial_similarity",
    "qve_delta",
    "structural_alert_penalty",
    "refinement_generation",
    "input_smiles",
    "smiles",
    "canonical_smiles",
    "affinity_kcal_mol",
    "official_binding_score",
    "binding_score_source",
    "mode_count",
    "best_affinity_mode",
    "inner_rmsd",
    "inner_cluster_fraction",
    "whole_rmsd",
    "hbond_count",
    "hydrophobic_count",
    "vdw_contact_count",
    "pi_contact_count",
    "ch_pi_count",
    "mw",
    "logp",
    "hbd",
    "hba",
    "tpsa",
    "rot_bonds",
    "heavy_atoms",
    "formal_charge",
    "sdf_path",
    "pdbqt_path",
    "pose_path",
    "log_path",
    "prep_log_path",
    "reason",
]

REFINEMENT_LEDGER_FIELDS = [
    "refinement_source",
    "druglike_refinement_score",
    "refinement_qed",
    "synthetic_score_proxy",
    "property_window_score",
    "qed_component_score",
    "synthetic_component_score",
    "property_component_score",
    "parent_similarity_score",
    "expert_similarity_score",
    "reference_similarity_score",
    "reference_partial_similarity_score",
    "qve_component_score",
    "inherited_structure_score",
    "structural_alert_penalty_score",
    "refinement_parent_similarity",
    "refinement_expert_similarity",
    "kinase_reference_similarity",
    "kinase_reference_partial_similarity",
    "qve_delta",
    "structural_alert_penalty",
    "refinement_generation",
]


@dataclass(frozen=True)
class VinaBox:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    exhaustiveness: int = 8
    cpu: int = 1
    num_modes: int = 9
    energy_range: int | None = None


@dataclass(frozen=True)
class PoseMode:
    mode_index: int
    affinity: float | None
    coords: list[tuple[float, float, float]]
    atoms: list[dict[str, str | float]]


def clean_smiles(smiles: str | None) -> str:
    return re.sub(r"\s+", "", str(smiles or "").strip())


def mol_from_smiles(smiles: str | None) -> Chem.Mol | None:
    text = clean_smiles(smiles)
    if not text:
        return None
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return mol


def canonical(smiles: str | None) -> str | None:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def standardized_mol(mol: Chem.Mol) -> Chem.Mol:
    out = Chem.Mol(mol)
    if rdMolStandardize is None:
        return out
    try:
        out = rdMolStandardize.Cleanup(out)
        out = rdMolStandardize.FragmentParent(out)
        out = rdMolStandardize.Uncharger().uncharge(out)
        Chem.SanitizeMol(out)
    except Exception:
        out = Chem.Mol(mol)
    return out


def smiles_identity(smiles: str | None, input_smiles: str | None = None) -> dict[str, str]:
    original = clean_smiles(input_smiles) or clean_smiles(smiles)
    mol = mol_from_smiles(smiles)
    if mol is None:
        return {
            "input_smiles": original,
            "smiles": clean_smiles(smiles),
            "canonical_smiles": "",
        }
    can = Chem.MolToSmiles(mol, canonical=True)
    return {
        "input_smiles": original,
        "smiles": can,
        "canonical_smiles": can,
    }


def same_molecule(a: dict[str, str], b: dict[str, str]) -> bool:
    return bool(a.get("canonical_smiles")) and a.get("canonical_smiles") == b.get("canonical_smiles")


def properties(smiles: str) -> dict[str, float | int]:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return {}
    return {
        "qed": round(float(QED.qed(mol)), 4),
        "mw": round(float(Descriptors.MolWt(mol)), 3),
        "logp": round(float(Crippen.MolLogP(mol)), 3),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(mol)), 3),
        "rot_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "formal_charge": int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms())),
    }


def analog_property_rejection(smiles: str, args: argparse.Namespace) -> str:
    props = properties(smiles)
    if not props:
        return "property_parse_failed"
    if float(props.get("mw", 0)) > float(args.max_analog_mw):
        return f"mw>{args.max_analog_mw}"
    if int(props.get("rot_bonds", 99)) > int(args.max_analog_rot_bonds):
        return f"rot_bonds>{args.max_analog_rot_bonds}"
    if float(props.get("tpsa", 0)) > float(args.max_analog_tpsa):
        return f"tpsa>{args.max_analog_tpsa}"
    if float(props.get("qed", 0)) < float(args.min_analog_qed):
        return f"qed<{args.min_analog_qed}"
    return ""


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    raw = path.read_bytes()
    for encoding in ("utf-8-sig", "utf-8", "gb18030", "cp936", "latin1"):
        try:
            text = raw.decode(encoding)
            return list(csv.DictReader(io.StringIO(text)))
        except UnicodeDecodeError:
            continue
    text = raw.decode("utf-8", errors="replace")
    return list(csv.DictReader(io.StringIO(text)))


def clean_float(value: object, default: float = 0.0) -> float:
    try:
        if value in ("", None):
            return default
        return float(value)
    except Exception:
        return default


def append_csv(path: Path, rows: list[dict[str, str | int | float]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LEDGER_FIELDS)
        if write_header:
            writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in LEDGER_FIELDS})


def write_csv(path: Path, rows: list[dict[str, str | int | float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LEDGER_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in LEDGER_FIELDS})


def refinement_meta_from_args(args: argparse.Namespace) -> dict[str, str]:
    meta: dict[str, str] = {}
    for field in REFINEMENT_LEDGER_FIELDS:
        value = getattr(args, field, "")
        if value not in ("", None):
            meta[field] = str(value)
    return meta


def merge_refinement_meta(row: dict[str, str | int | float], meta: dict[str, str]) -> None:
    if not meta:
        return
    current_score = clean_float(row.get("druglike_refinement_score"), default=float("-inf"))
    incoming_score = clean_float(meta.get("druglike_refinement_score"), default=float("-inf"))
    should_replace_scores = incoming_score >= current_score
    for field, value in meta.items():
        if value in ("", None):
            continue
        if field == "refinement_source":
            existing = str(row.get(field, ""))
            row[field] = value if not existing else existing
            continue
        if should_replace_scores or row.get(field, "") in ("", None):
            row[field] = value


def ledger_path(ledger_dir: Path) -> Path:
    return ledger_dir / LEDGER_FILENAME


def normalize_ledger_dir(path: Path | None) -> Path:
    raw = path or DEFAULT_LEDGER_DIR
    expanded = raw.expanduser()
    if expanded.name.lower() == "dock_hisotry" and not expanded.exists():
        corrected = expanded.with_name("dock_history")
        if corrected.exists():
            return corrected.resolve()
    return expanded.resolve()


def next_seq(rows: list[dict[str, str]]) -> int:
    found: list[int] = []
    for row in rows:
        match = re.search(r"(\d+)$", row.get("seq_id", ""))
        if match:
            found.append(int(match.group(1)))
    return max(found, default=0) + 1


def seq_id(number: int) -> str:
    return f"S{number:06d}"


def path_text(path: str | None) -> str:
    return str(Path(path).expanduser()) if path else ""


def default_vina_cpu() -> int:
    logical = os.cpu_count() or 1
    if logical <= 2:
        return logical
    # Vina often scales better around physical-core count than full SMT count.
    return max(1, min(16, logical // 2))


def load_saved_config(ledger_dir: Path) -> dict[str, str | int | float]:
    path = ledger_dir / CONFIG_FILENAME
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def save_config(ledger_dir: Path, args: argparse.Namespace) -> None:
    fields = {
        "receptor": path_text(str(args.receptor) if args.receptor else ""),
        "config": path_text(str(args.config) if args.config else ""),
        "vina": args.vina,
        "meeko": args.meeko,
        "obabel": args.obabel,
        "cpu": args.cpu,
        "exhaustiveness": args.exhaustiveness,
        "num_modes": args.num_modes,
        "energy_range": args.energy_range,
    }
    existing = load_saved_config(ledger_dir)
    existing.update({key: value for key, value in fields.items() if value not in ("", None)})
    (ledger_dir / CONFIG_FILENAME).write_text(json.dumps(existing, indent=2), encoding="utf-8")


def arg_or_config(args: argparse.Namespace, name: str, default: str | int | None = None):
    value = getattr(args, name, None)
    if value not in (None, ""):
        return value
    return load_saved_config(args.ledger_dir.expanduser().resolve()).get(name, default)


def resolve_executable(value: str, ledger_dir: Path, filename: str) -> str:
    if command_exists(value):
        return value
    local = ledger_dir.parent / "vina_bin" / filename
    if local.exists():
        return str(local)
    return value


def parse_vina_config(path: Path, args: argparse.Namespace) -> VinaBox:
    values: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = [part.strip() for part in line.split("=", 1)]
        try:
            values[key] = float(value)
        except ValueError:
            continue
    required = ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]
    missing = [key for key in required if key not in values]
    if missing:
        raise ValueError(f"Vina config missing: {', '.join(missing)}")
    saved = load_saved_config(args.ledger_dir.expanduser().resolve())

    def choose(name: str, default: int | float | None) -> int | float | None:
        cli_value = getattr(args, name, None)
        if cli_value not in (None, ""):
            return cli_value
        if name in values:
            return values[name]
        saved_value = saved.get(name, default)
        return saved_value

    return VinaBox(
        center_x=values["center_x"],
        center_y=values["center_y"],
        center_z=values["center_z"],
        size_x=values["size_x"],
        size_y=values["size_y"],
        size_z=values["size_z"],
        exhaustiveness=int(choose("exhaustiveness", 8) or 8),
        cpu=int(choose("cpu", default_vina_cpu()) or default_vina_cpu()),
        num_modes=int(choose("num_modes", 9) or 9),
        energy_range=choose("energy_range", None),
    )


def run(cmd: list[str]) -> tuple[int, str]:
    try:
        proc = subprocess.run(cmd, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        return proc.returncode, proc.stdout
    except FileNotFoundError as exc:
        return 127, str(exc)


def command_exists(command: str) -> bool:
    return Path(command).exists() or shutil.which(command) is not None


def embed_sdf(smiles: str, name: str, sdf_path: Path, seed: int) -> None:
    mol = mol_from_smiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise RuntimeError("RDKit 3D embedding failed")
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=250)
    except Exception:
        pass
    mol.SetProp("_Name", name)
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol)
    writer.close()


def prepare_pdbqt(sdf_path: Path, pdbqt_path: Path, meeko: str, obabel: str) -> tuple[bool, str]:
    commands: list[list[str]] = []
    if command_exists(meeko):
        commands.append([meeko, "-i", str(sdf_path), "-o", str(pdbqt_path)])
    commands.append([sys.executable, "-m", "meeko.cli.mk_prepare_ligand", "-i", str(sdf_path), "-o", str(pdbqt_path)])
    commands.append([obabel, str(sdf_path), "-O", str(pdbqt_path)])
    outputs: list[str] = []
    for cmd in commands:
        code, output = run(cmd)
        outputs.append("$ " + " ".join(cmd) + "\n" + output)
        if code == 0 and pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
            return True, "\n\n".join(outputs)
    return False, "\n\n".join(outputs)


def dock(vina: str, receptor: Path, ligand: Path, pose_path: Path, log_path: Path, box: VinaBox, seed: int) -> tuple[bool, float | None, str]:
    cmd = [
        vina,
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
        "--out",
        str(pose_path),
        "--log",
        str(log_path),
        "--center_x",
        str(box.center_x),
        "--center_y",
        str(box.center_y),
        "--center_z",
        str(box.center_z),
        "--size_x",
        str(box.size_x),
        "--size_y",
        str(box.size_y),
        "--size_z",
        str(box.size_z),
        "--exhaustiveness",
        str(box.exhaustiveness),
        "--cpu",
        str(box.cpu),
        "--num_modes",
        str(box.num_modes),
        "--seed",
        str(seed),
    ]
    if box.energy_range is not None:
        cmd.extend(["--energy_range", str(box.energy_range)])
    code, output = run(cmd)
    if not log_path.exists():
        log_path.write_text(output, encoding="utf-8", errors="replace")
    text = output
    if log_path.exists():
        text += "\n" + log_path.read_text(encoding="utf-8", errors="replace")
    if pose_path.exists():
        text += "\n" + pose_path.read_text(encoding="utf-8", errors="replace")
    match = AFFINITY_RE.search(text) or PDBQT_AFFINITY_RE.search(text)
    affinity = float(match.group(1)) if match else None
    return code == 0 and affinity is not None and pose_path.exists(), affinity, output


def parse_atom_line(line: str) -> dict[str, str | float] | None:
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None
    atom_name = line[12:16].strip()
    residue_name = line[17:20].strip()
    residue_id = line[22:26].strip()
    atom_type = line[77:].strip().split()[0] if len(line) >= 78 and line[77:].strip() else ""
    element = re.sub(r"[^A-Za-z]", "", atom_type or atom_name)[:2].title()
    if len(element) == 2 and element not in {"Cl", "Br"}:
        element = element[0]
    return {
        "x": x,
        "y": y,
        "z": z,
        "atom_name": atom_name,
        "atom_type": atom_type,
        "element": element,
        "residue_name": residue_name,
        "residue_id": residue_id,
    }


def parse_pdbqt_modes(path: Path) -> list[PoseMode]:
    modes: list[PoseMode] = []
    if not path.exists():
        return modes
    current_atoms: list[dict[str, str | float]] = []
    current_affinity: float | None = None
    current_index = 1
    saw_model = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            saw_model = True
            if current_atoms:
                modes.append(PoseMode(current_index, current_affinity, [(float(a["x"]), float(a["y"]), float(a["z"])) for a in current_atoms], current_atoms))
                current_atoms = []
                current_affinity = None
            try:
                current_index = int(line.split()[1])
            except Exception:
                current_index = len(modes) + 1
        elif line.startswith("REMARK VINA RESULT:"):
            match = PDBQT_AFFINITY_RE.search(line)
            if match:
                current_affinity = float(match.group(1))
        elif line.startswith(("ATOM", "HETATM")):
            atom = parse_atom_line(line)
            if atom and str(atom.get("element")) != "H":
                current_atoms.append(atom)
        elif line.startswith("ENDMDL") and current_atoms:
            modes.append(PoseMode(current_index, current_affinity, [(float(a["x"]), float(a["y"]), float(a["z"])) for a in current_atoms], current_atoms))
            current_atoms = []
            current_affinity = None
    if current_atoms:
        modes.append(PoseMode(current_index if saw_model else 1, current_affinity, [(float(a["x"]), float(a["y"]), float(a["z"])) for a in current_atoms], current_atoms))
    return modes


def parse_receptor_atoms(path: Path) -> list[dict[str, str | float]]:
    atoms: list[dict[str, str | float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atom = parse_atom_line(line)
            if atom:
                atoms.append(atom)
    return atoms


def coord_rmsd(a: list[tuple[float, float, float]], b: list[tuple[float, float, float]]) -> float:
    n = min(len(a), len(b))
    if n == 0:
        return float("nan")
    total = 0.0
    for i in range(n):
        total += (a[i][0] - b[i][0]) ** 2 + (a[i][1] - b[i][1]) ** 2 + (a[i][2] - b[i][2]) ** 2
    return math.sqrt(total / n)


def squared_distance(a: dict[str, str | float], b: dict[str, str | float]) -> float:
    return (float(a["x"]) - float(b["x"])) ** 2 + (float(a["y"]) - float(b["y"])) ** 2 + (float(a["z"]) - float(b["z"])) ** 2


def interaction_counts(ligand_atoms: list[dict[str, str | float]], receptor_atoms: list[dict[str, str | float]]) -> dict[str, int]:
    hbonds: set[tuple[int, int]] = set()
    hydrophobic: set[tuple[int, int]] = set()
    vdw: set[tuple[int, int]] = set()
    pi_contacts: set[tuple[int, str]] = set()
    ch_pi: set[tuple[int, str]] = set()
    for i, la in enumerate(ligand_atoms):
        le = str(la.get("element", ""))
        atom_type = str(la.get("atom_type", ""))
        for j, ra in enumerate(receptor_atoms):
            relem = str(ra.get("element", ""))
            d2 = squared_distance(la, ra)
            if le in HBOND_ELEMENTS and relem in HBOND_ELEMENTS and d2 <= 3.5**2:
                hbonds.add((i, j))
            if le in HYDROPHOBIC_ELEMENTS and relem in HYDROPHOBIC_ELEMENTS and d2 <= 4.5**2:
                hydrophobic.add((i, j))
            if d2 <= 4.0**2:
                vdw.add((i, j))
            residue = str(ra.get("residue_name", ""))
            if le == "C" and residue in AROMATIC_RESIDUES and d2 <= 5.0**2:
                pi_contacts.add((i, f"{residue}{ra.get('residue_id', '')}"))
            if le == "C" and residue in AROMATIC_RESIDUES and d2 <= 4.2**2 and not atom_type.startswith("A"):
                ch_pi.add((i, f"{residue}{ra.get('residue_id', '')}"))
    return {
        "hbond_count": len(hbonds),
        "hydrophobic_count": len(hydrophobic),
        "vdw_contact_count": len(vdw),
        "pi_contact_count": len(pi_contacts),
        "ch_pi_count": len(ch_pi),
    }


def pose_metrics(pose_path: Path, receptor_atoms: list[dict[str, str | float]], cluster_cutoff: float) -> dict[str, str | float | int]:
    modes = parse_pdbqt_modes(pose_path)
    if not modes:
        return {"mode_count": 0}
    best = min(modes, key=lambda mode: mode.affinity if mode.affinity is not None else float("inf"))
    rmsds = [coord_rmsd(best.coords, mode.coords) for mode in modes]
    finite = [r for r in rmsds if math.isfinite(r)]
    near = [r for r in finite if r <= cluster_cutoff]
    pairwise: list[float] = []
    for i, mode in enumerate(modes):
        for other in modes[i + 1 :]:
            r = coord_rmsd(mode.coords, other.coords)
            if math.isfinite(r):
                pairwise.append(r)
    return {
        "mode_count": len(modes),
        "best_affinity_mode": best.mode_index,
        "inner_rmsd": round(sum(finite) / len(finite), 4) if finite else "",
        "inner_cluster_fraction": round(len(near) / len(modes), 4) if modes else "",
        "whole_rmsd": round(sum(pairwise) / len(pairwise), 4) if pairwise else (0.0 if len(modes) == 1 else ""),
        **interaction_counts(best.atoms, receptor_atoms),
    }


def sanitize_product(mol: Chem.Mol) -> str | None:
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def add_substituent(mol: Chem.Mol, atom_idx: int, label: str) -> str | None:
    rw = Chem.RWMol(mol)
    def add_atom(symbol: str) -> int:
        return rw.AddAtom(Chem.Atom(symbol))
    try:
        if label in {"F", "Cl", "Br", "C", "N", "O"}:
            new = add_atom({"C": "C", "N": "N", "O": "O"}.get(label, label))
            rw.AddBond(atom_idx, new, Chem.BondType.SINGLE)
        elif label == "Me":
            new = add_atom("C")
            rw.AddBond(atom_idx, new, Chem.BondType.SINGLE)
        elif label == "OH":
            new = add_atom("O")
            rw.AddBond(atom_idx, new, Chem.BondType.SINGLE)
        elif label == "NH2":
            new = add_atom("N")
            rw.AddBond(atom_idx, new, Chem.BondType.SINGLE)
        elif label == "CN":
            c = add_atom("C")
            n = add_atom("N")
            rw.AddBond(atom_idx, c, Chem.BondType.SINGLE)
            rw.AddBond(c, n, Chem.BondType.TRIPLE)
        elif label == "CHO":
            c = add_atom("C")
            o = add_atom("O")
            rw.AddBond(atom_idx, c, Chem.BondType.SINGLE)
            rw.AddBond(c, o, Chem.BondType.DOUBLE)
        elif label == "OMe":
            o = add_atom("O")
            c = add_atom("C")
            rw.AddBond(atom_idx, o, Chem.BondType.SINGLE)
            rw.AddBond(o, c, Chem.BondType.SINGLE)
        else:
            return None
        return sanitize_product(rw.GetMol())
    except Exception:
        return None


def aromatic_additions(smiles: str) -> list[tuple[str, str]]:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return []
    out: list[tuple[str, str]] = []
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic() and atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() > 0:
            for label in ("F", "Cl", "Me", "OH", "NH2", "CN", "CHO", "OMe"):
                product = add_substituent(mol, atom.GetIdx(), label)
                if product:
                    out.append((product, f"add_aromatic_C{atom.GetIdx()}_{label}"))
    return out


def collect_branch(mol: Chem.Mol, start: int, blocked: int, max_atoms: int = 8) -> set[int]:
    seen = {blocked}
    stack = [start]
    branch: set[int] = set()
    while stack and len(branch) <= max_atoms:
        idx = stack.pop()
        if idx in seen:
            continue
        seen.add(idx)
        atom = mol.GetAtomWithIdx(idx)
        if atom.IsInRing():
            continue
        branch.add(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in seen:
                stack.append(nbr.GetIdx())
    return branch if len(branch) <= max_atoms else set()


def remove_atoms(smiles: str, atoms: set[int]) -> str | None:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    rw = Chem.RWMol(mol)
    for idx in sorted(atoms, reverse=True):
        rw.RemoveAtom(idx)
    return sanitize_product(rw.GetMol())


def deletions_and_shrinks(smiles: str, shrink: bool) -> list[tuple[str, str]]:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return []
    out: list[tuple[str, str]] = []
    for atom in mol.GetAtoms():
        if not atom.GetIsAromatic():
            continue
        for nbr in atom.GetNeighbors():
            if nbr.GetIsAromatic() or nbr.IsInRing():
                continue
            branch = collect_branch(mol, nbr.GetIdx(), atom.GetIdx())
            if not branch:
                continue
            deleted = remove_atoms(smiles, branch)
            if deleted:
                out.append((deleted, f"delete_C{atom.GetIdx()}_{len(branch)}ha"))
            if shrink:
                removed = remove_atoms(smiles, branch)
                if removed:
                    base = mol_from_smiles(removed)
                    if base is not None:
                        # atom index is stable only if removed atoms were after it; otherwise skip risky replacement.
                        if atom.GetIdx() < base.GetNumAtoms():
                            product = add_substituent(base, atom.GetIdx(), "Me")
                            if product:
                                out.append((product, f"shrink_C{atom.GetIdx()}_to_Me"))
    return out


def n_swaps(smiles: str) -> list[tuple[str, str]]:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return []
    out: list[tuple[str, str]] = []
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic() and atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() > 0:
            rw = Chem.RWMol(mol)
            rw.GetAtomWithIdx(atom.GetIdx()).SetAtomicNum(7)
            product = sanitize_product(rw.GetMol())
            if product:
                out.append((product, f"drastic_aromatic_C{atom.GetIdx()}_Nswap"))
    return out


def enumerate_edits(smiles: str, modes: tuple[str, ...]) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    if "add" in modes:
        out.extend(aromatic_additions(smiles))
    if "delete" in modes:
        out.extend(deletions_and_shrinks(smiles, shrink=False))
    if "shrink" in modes:
        out.extend(deletions_and_shrinks(smiles, shrink=True))
    if "drastic" in modes:
        out.extend(n_swaps(smiles))
        out.extend(deletions_and_shrinks(smiles, shrink=True))
        out.extend(aromatic_additions(smiles))
    dedup: dict[str, str] = {}
    for product, label in out:
        can = canonical(product)
        if can and can != canonical(smiles):
            dedup.setdefault(can, label)
    return [(product, label) for product, label in dedup.items()]


def parse_edit_modes(text: str) -> tuple[str, ...]:
    if text in ("", "default", "all"):
        return ("add", "delete", "shrink", "drastic")
    modes: list[str] = []
    aliases = {"del": "delete", "remove": "delete", "small": "shrink", "shirnk": "shrink"}
    for raw in re.split(r"[,; ]+", text):
        if not raw:
            continue
        mode = aliases.get(raw.lower(), raw.lower())
        if mode not in {"add", "delete", "shrink", "drastic"}:
            raise ValueError(f"unknown edit mode: {raw}")
        if mode not in modes:
            modes.append(mode)
    return tuple(modes)


def selected_edit_mode(args: argparse.Namespace) -> str:
    modes = []
    for flag, mode in (("add_only", "add"), ("delete_only", "delete"), ("shrink_only", "shrink"), ("drastic_only", "drastic")):
        if getattr(args, flag, False):
            modes.append(mode)
    return ",".join(modes) if modes else args.edit_mode


def row_identity(row: dict[str, str]) -> dict[str, str]:
    return {
        "canonical_smiles": row.get("canonical_smiles", ""),
    }


def find_existing(rows: list[dict[str, str]], ids: dict[str, str]) -> dict[str, str] | None:
    for row in rows:
        if same_molecule(row_identity(row), ids):
            return row
    return None


def dock_one(
    *,
    seq: str,
    smiles: str,
    input_smiles: str,
    nickname: str,
    ancestor_smiles: str,
    parent_smiles: str,
    edit_label: str,
    ledger_dir: Path,
    receptor: Path,
    receptor_atoms: list[dict[str, str | float]],
    box: VinaBox,
    vina: str,
    meeko: str,
    obabel: str,
    seed: int,
    cluster_cutoff: float,
) -> dict[str, str | int | float]:
    ids = smiles_identity(smiles, input_smiles)
    can = ids.get("canonical_smiles")
    base = {
        "seq_id": seq,
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "nickname": nickname,
        "ancestor_smiles": ancestor_smiles,
        "parent_smiles": parent_smiles,
        "edit_label": edit_label,
        **ids,
    }
    if not can:
        return {**base, "reason": "rdkit_failed"}
    for subdir in ("sdf", "pdbqt", "poses", "logs", "prep_logs"):
        (ledger_dir / subdir).mkdir(parents=True, exist_ok=True)
    sdf_path = ledger_dir / "sdf" / f"{seq}.sdf"
    pdbqt_path = ledger_dir / "pdbqt" / f"{seq}.pdbqt"
    pose_path = ledger_dir / "poses" / f"{seq}.pdbqt"
    log_path = ledger_dir / "logs" / f"{seq}.log"
    prep_log_path = ledger_dir / "prep_logs" / f"{seq}.log"
    reason = ""
    affinity: float | str = ""
    try:
        embed_sdf(can, seq, sdf_path, seed)
        prep_ok, prep_output = prepare_pdbqt(sdf_path, pdbqt_path, meeko, obabel)
        prep_log_path.write_text(prep_output, encoding="utf-8", errors="replace")
        if not prep_ok:
            reason = "pdbqt_preparation_failed"
        else:
            dock_ok, dock_affinity, dock_output = dock(vina, receptor, pdbqt_path, pose_path, log_path, box, seed)
            if dock_ok and dock_affinity is not None:
                affinity = dock_affinity
            else:
                if "No such file or directory" in dock_output or "not found" in dock_output:
                    reason = "vina_executable_not_found"
                else:
                    reason = "vina_failed_or_no_affinity"
    except Exception as exc:
        reason = str(exc)
    metrics = pose_metrics(pose_path, receptor_atoms, cluster_cutoff) if not reason else {}
    props = properties(can) if not reason else {}
    return {
        **base,
        "affinity_kcal_mol": affinity,
        **metrics,
        **props,
        "sdf_path": str(sdf_path),
        "pdbqt_path": str(pdbqt_path),
        "pose_path": str(pose_path),
        "log_path": str(log_path),
        "prep_log_path": str(prep_log_path),
        "reason": reason,
    }


def collect(args: argparse.Namespace) -> list[dict[str, str | int | float]]:
    ledger_dir = normalize_ledger_dir(getattr(args, "ledger_dir", None))
    args.ledger_dir = ledger_dir
    ledger_dir.mkdir(parents=True, exist_ok=True)
    save_config(ledger_dir, args)
    receptor_arg = arg_or_config(args, "receptor")
    config_arg = arg_or_config(args, "config")
    if not receptor_arg or not config_arg:
        raise ValueError("provide --receptor and --config once; later runs can read history_config.json")
    receptor = Path(str(receptor_arg)).expanduser().resolve()
    config = Path(str(config_arg)).expanduser().resolve()
    vina = resolve_executable(str(arg_or_config(args, "vina", "vina")), ledger_dir, "vina")
    meeko = resolve_executable(str(arg_or_config(args, "meeko", "mk_prepare_ligand.py")), ledger_dir, "mk_prepare_ligand.py")
    obabel = resolve_executable(str(arg_or_config(args, "obabel", "obabel")), ledger_dir, "obabel")
    box = parse_vina_config(config, args)
    print(f"vina_cpu={box.cpu} exhaustiveness={box.exhaustiveness} num_modes={box.num_modes}", flush=True)
    receptor_atoms = parse_receptor_atoms(receptor)
    ledger_file = ledger_path(ledger_dir)
    existing_rows = read_csv(ledger_file)
    seq_number = next_seq(existing_rows)
    rng = random.Random(args.seed)
    modes = parse_edit_modes(selected_edit_mode(args))

    input_smiles = clean_smiles(args.smiles)
    input_nickname = clean_smiles(getattr(args, "nickname", ""))
    input_ancestor_smiles = clean_smiles(getattr(args, "ancestor_smiles", ""))
    input_parent_smiles = clean_smiles(getattr(args, "parent_smiles", ""))
    input_edit_label = clean_smiles(getattr(args, "edit_label", "")) or "input_smiles"
    refinement_meta = refinement_meta_from_args(args)
    ids = smiles_identity(input_smiles, input_smiles)
    existing = find_existing(existing_rows, ids)
    new_rows: list[dict[str, str | int | float]] = []
    if existing and not args.redock_existing:
        print(f"input_exists=true seq_id={existing.get('seq_id')}", flush=True)
        if input_nickname:
            existing["nickname"] = input_nickname
        if input_ancestor_smiles and not existing.get("ancestor_smiles"):
            existing["ancestor_smiles"] = input_ancestor_smiles
        if input_parent_smiles and not existing.get("parent_smiles"):
            existing["parent_smiles"] = input_parent_smiles
        if input_edit_label and existing.get("edit_label") in ("", "input_smiles"):
            existing["edit_label"] = input_edit_label
        merge_refinement_meta(existing, refinement_meta)
        active = [existing]
    else:
        seq = seq_id(seq_number)
        seq_number += 1
        row = dock_one(
            seq=seq,
            smiles=ids.get("canonical_smiles") or input_smiles,
            input_smiles=input_smiles,
            nickname=input_nickname,
            ancestor_smiles=input_ancestor_smiles or ids.get("canonical_smiles") or input_smiles,
            parent_smiles=input_parent_smiles,
            edit_label=input_edit_label,
            ledger_dir=ledger_dir,
            receptor=receptor,
            receptor_atoms=receptor_atoms,
            box=box,
            vina=vina,
            meeko=meeko,
            obabel=obabel,
            seed=args.seed,
            cluster_cutoff=args.internal_cluster_rmsd_cutoff,
        )
        merge_refinement_meta(row, refinement_meta)
        print(f"[dock] {seq} ancestor affinity={row.get('affinity_kcal_mol')} reason={row.get('reason', '')}", flush=True)
        new_rows.append(row)
        active = [{key: str(value) for key, value in row.items()}]

    seen = [row_identity(row) for row in existing_rows]
    seen.extend(row_identity({key: str(value) for key, value in row.items()}) for row in new_rows)
    for generation in range(1, args.max_rounds + 1):
        next_active: list[dict[str, str]] = []
        docked = 0
        for parent in active:
            parent_smiles = parent.get("canonical_smiles") or parent.get("smiles") or ""
            candidates = enumerate_edits(parent_smiles, modes)
            if not args.deterministic_batch:
                rng.shuffle(candidates)
            for candidate_smiles, edit_label in candidates:
                cids = smiles_identity(candidate_smiles, candidate_smiles)
                if not cids.get("canonical_smiles") or any(same_molecule(cids, old) for old in seen):
                    continue
                rejection = analog_property_rejection(cids["canonical_smiles"], args)
                if rejection:
                    print(f"[skip] gen={generation} parent={parent.get('seq_id')} {edit_label} reason={rejection}", flush=True)
                    seen.append(cids)
                    continue
                seq = seq_id(seq_number)
                seq_number += 1
                row = dock_one(
                    seq=seq,
                    smiles=cids["canonical_smiles"],
                    input_smiles=candidate_smiles,
                    nickname="",
                    ancestor_smiles=parent.get("ancestor_smiles") or parent_smiles,
                    parent_smiles=parent_smiles,
                    edit_label=edit_label,
                    ledger_dir=ledger_dir,
                    receptor=receptor,
                    receptor_atoms=receptor_atoms,
                    box=box,
                    vina=vina,
                    meeko=meeko,
                    obabel=obabel,
                    seed=args.seed + seq_number,
                    cluster_cutoff=args.internal_cluster_rmsd_cutoff,
                )
                print(
                    f"[dock] {seq} gen={generation} parent={parent.get('seq_id')} {edit_label} "
                    f"affinity={row.get('affinity_kcal_mol')} inner={row.get('inner_rmsd', '')} whole={row.get('whole_rmsd', '')}",
                    flush=True,
                )
                new_rows.append(row)
                seen.append(row_identity({key: str(value) for key, value in row.items()}))
                if row.get("reason", "") == "":
                    next_active.append({key: str(value) for key, value in row.items()})
                docked += 1
                if docked >= args.batch_size:
                    break
            if docked >= args.batch_size:
                break
        active = next_active
        if not active:
            break
    write_csv(ledger_file, existing_rows + new_rows)
    print(f"wrote {ledger_file} new_rows={len(new_rows)} total_rows={len(existing_rows) + len(new_rows)}", flush=True)
    return new_rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smiles", required=True)
    parser.add_argument("--nickname", default="", help="Optional user-facing label to preserve in the ledger")
    parser.add_argument("--history-dir", "--ledger-dir", dest="ledger_dir", metavar="HISTORY_DIR", type=Path, default=DEFAULT_LEDGER_DIR)
    parser.add_argument("--receptor", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--vina")
    parser.add_argument("--meeko")
    parser.add_argument("--obabel")
    parser.add_argument("--max-rounds", type=int, default=0)
    parser.add_argument("--batch-size", type=int, default=24)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--exhaustiveness", type=int)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument("--energy-range", type=int)
    parser.add_argument("--internal-cluster-rmsd-cutoff", type=float, default=3.0)
    parser.add_argument("--edit-mode", default="default")
    parser.add_argument("--add", dest="add_only", action="store_true")
    parser.add_argument("--delete", dest="delete_only", action="store_true")
    parser.add_argument("--shrink", dest="shrink_only", action="store_true")
    parser.add_argument("--drastic", dest="drastic_only", action="store_true")
    parser.add_argument("--deterministic-batch", action="store_true")
    parser.add_argument("--redock-existing", action="store_true")
    parser.add_argument("--min-analog-qed", type=float, default=0.20)
    parser.add_argument("--max-analog-mw", type=float, default=650.0)
    parser.add_argument("--max-analog-rot-bonds", type=int, default=14)
    parser.add_argument("--max-analog-tpsa", type=float, default=180.0)
    return parser


def main() -> None:
    collect(build_parser().parse_args())


if __name__ == "__main__":
    main()
