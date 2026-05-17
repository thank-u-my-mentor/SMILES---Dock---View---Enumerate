#!/usr/bin/env python
"""Build Pure-SMILES and structural-interaction spaces, then learn binding score."""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import re
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger, rdBase
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, QED, rdFMCS, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.manifold import MDS
from sklearn.metrics import mean_absolute_error, pairwise_distances, r2_score
from sklearn.model_selection import cross_val_score, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")
rdBase.DisableLog("rdApp.*")


HISTORY_FIELDS = [
    "seq_id",
    "nickname",
    "input_smiles",
    "canonical_smiles",
    "official_binding_score",
    "affinity_kcal_mol",
    "inner_rmsd",
    "whole_rmsd",
    "hbond_count",
    "hydrophobic_count",
    "pi_contact_count",
    "ch_pi_count",
    "mw",
    "logp",
    "hbd",
    "hba",
    "tpsa",
    "rot_bonds",
    "druglike_refinement_score",
    "contact_residue_count_4a",
    "surface_contact_fraction_4a",
    "receptor_residue_coverage_4a",
    "ligand_atom_contact_fraction_4a",
    "min_ligand_receptor_distance",
    "pose_path",
]

POCKET_SUMMARY_FEATURES = [
    "pocket_atom_count_5a",
    "pocket_residue_count_5a",
    "pocket_rg_5a",
    "pocket_span_5a",
    "pocket_centroid_distance_5a",
    "pocket_hydrophobic_fraction_5a",
    "pocket_polar_fraction_5a",
    "pocket_charged_fraction_5a",
    "pocket_positive_fraction_5a",
    "pocket_negative_fraction_5a",
    "pocket_aromatic_fraction_5a",
]

POCKET_GEOMETRY_FIELDS = [f"pocket_geom_{i:02d}" for i in range(64)]
POCKET_PHYSCHEM_FIELDS = [f"pocket_physchem_{i:02d}" for i in range(24)]
POCKET_EMBEDDING_FIELDS = [*POCKET_GEOMETRY_FIELDS, *POCKET_PHYSCHEM_FIELDS]

STRUCTURAL_FEATURES = [
    "affinity_kcal_mol",
    "inner_rmsd",
    "whole_rmsd",
    "hbond_count",
    "hydrophobic_count",
    "pi_contact_count",
    "ch_pi_count",
    "logp",
    "hbd",
    "hba",
    "tpsa",
    "rot_bonds",
    "qed",
    "aromatic_rings",
    "chembl_scaffold_similarity",
    "contact_residue_count_4a",
    "surface_contact_fraction_4a",
    "receptor_residue_coverage_4a",
    "ligand_atom_contact_fraction_4a",
    "min_ligand_receptor_distance",
    *POCKET_SUMMARY_FEATURES,
]

AA_PROPERTIES = {
    # hydrophobic, polar, positive, negative, aromatic, size proxy
    "ALA": (0.62, 0.00, 0.00, 0.00, 0.00, 0.31),
    "ARG": (0.00, 0.70, 1.00, 0.00, 0.00, 1.00),
    "ASN": (0.00, 0.85, 0.00, 0.00, 0.00, 0.55),
    "ASP": (0.00, 0.90, 0.00, 1.00, 0.00, 0.54),
    "CYS": (0.68, 0.25, 0.00, 0.00, 0.00, 0.43),
    "GLN": (0.00, 0.82, 0.00, 0.00, 0.00, 0.68),
    "GLU": (0.00, 0.88, 0.00, 1.00, 0.00, 0.67),
    "GLY": (0.48, 0.00, 0.00, 0.00, 0.00, 0.20),
    "HIS": (0.20, 0.60, 0.45, 0.00, 1.00, 0.70),
    "ILE": (1.00, 0.00, 0.00, 0.00, 0.00, 0.72),
    "LEU": (0.98, 0.00, 0.00, 0.00, 0.00, 0.72),
    "LYS": (0.05, 0.65, 1.00, 0.00, 0.00, 0.82),
    "MET": (0.78, 0.10, 0.00, 0.00, 0.00, 0.74),
    "PHE": (0.88, 0.00, 0.00, 0.00, 1.00, 0.82),
    "PRO": (0.64, 0.00, 0.00, 0.00, 0.00, 0.57),
    "SER": (0.18, 0.75, 0.00, 0.00, 0.00, 0.38),
    "THR": (0.30, 0.70, 0.00, 0.00, 0.00, 0.49),
    "TRP": (0.85, 0.15, 0.00, 0.00, 1.00, 1.00),
    "TYR": (0.55, 0.45, 0.00, 0.00, 1.00, 0.88),
    "VAL": (0.86, 0.00, 0.00, 0.00, 0.00, 0.60),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    for encoding in ("utf-8-sig", "utf-8", "gb18030", "cp936", "latin1"):
        try:
            with path.open("r", encoding=encoding, newline="") as handle:
                return list(csv.DictReader(handle))
        except UnicodeDecodeError:
            continue
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def safe_float(value: object) -> float:
    try:
        if value in ("", None):
            return math.nan
        return float(value)
    except Exception:
        return math.nan


def finite(value: object) -> bool:
    parsed = safe_float(value)
    return not math.isnan(parsed) and math.isfinite(parsed)


def non_nan_values(rows: list[dict[str, object]], field: str) -> list[float]:
    values = [safe_float(row.get(field)) for row in rows]
    return [value for value in values if not math.isnan(value)]


def mol_from_smiles(smiles: str | None) -> Chem.Mol | None:
    if not smiles:
        return None
    try:
        return Chem.MolFromSmiles(str(smiles).strip())
    except Exception:
        return None


def canonical(smiles: str | None) -> str:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def resolve_maybe_path(value: object, base_dir: Path | None = None) -> Path | None:
    text = str(value or "").strip()
    if not text:
        return None
    path = Path(text).expanduser()
    candidates = [path]
    if base_dir is not None and not path.is_absolute():
        candidates.append(base_dir / path)
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return path.resolve() if path.is_absolute() else None


def parse_pdbqt_atom_line(line: str) -> dict[str, object] | None:
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except Exception:
        return None
    atom_name = line[12:16].strip()
    residue_name = line[17:20].strip().upper()
    chain_id = line[21:22].strip()
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
        "chain_id": chain_id,
    }


def read_pdbqt_atoms(path: Path | None) -> list[dict[str, object]]:
    if path is None or not path.exists():
        return []
    atoms: list[dict[str, object]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atom = parse_pdbqt_atom_line(line)
            if atom and str(atom.get("element", "")) != "H":
                atoms.append(atom)
    return atoms


def parse_ligand_best_mode(path: Path | None) -> list[dict[str, object]]:
    if path is None or not path.exists():
        return []
    modes: list[tuple[int, float | None, list[dict[str, object]]]] = []
    current_atoms: list[dict[str, object]] = []
    current_affinity: float | None = None
    current_index = 1
    saw_model = False
    affinity_re = re.compile(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)")

    def flush() -> None:
        nonlocal current_atoms, current_affinity
        if current_atoms:
            modes.append((current_index, current_affinity, current_atoms))
            current_atoms = []
            current_affinity = None

    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            saw_model = True
            flush()
            try:
                current_index = int(line.split()[1])
            except Exception:
                current_index = len(modes) + 1
        elif line.startswith("REMARK VINA RESULT:"):
            match = affinity_re.search(line)
            if match:
                current_affinity = float(match.group(1))
        elif line.startswith(("ATOM", "HETATM")):
            atom = parse_pdbqt_atom_line(line)
            if atom and str(atom.get("element", "")) != "H":
                current_atoms.append(atom)
        elif line.startswith("ENDMDL"):
            flush()
    if current_atoms:
        if not saw_model:
            current_index = 1
        flush()
    if not modes:
        return []
    best = min(modes, key=lambda item: item[1] if item[1] is not None else float("inf"))
    return best[2]


def atom_coord(atom: dict[str, object]) -> np.ndarray:
    return np.array([float(atom["x"]), float(atom["y"]), float(atom["z"])], dtype=float)


def residue_key(atom: dict[str, object]) -> tuple[str, str, str]:
    return (
        str(atom.get("chain_id", "")),
        str(atom.get("residue_name", "")),
        str(atom.get("residue_id", "")),
    )


def pocket_physchem_embedding(residue_names: list[str]) -> tuple[np.ndarray, dict[str, float]]:
    if not residue_names:
        return np.zeros((24,), dtype=float), {
            "hydrophobic": math.nan,
            "polar": math.nan,
            "charged": math.nan,
            "positive": math.nan,
            "negative": math.nan,
            "aromatic": math.nan,
        }
    props = np.array([AA_PROPERTIES.get(name.upper(), (0.0, 0.0, 0.0, 0.0, 0.0, 0.5)) for name in residue_names], dtype=float)
    mean = props.mean(axis=0)
    std = props.std(axis=0)
    minv = props.min(axis=0)
    maxv = props.max(axis=0)
    embedding = np.concatenate([mean, std, minv, maxv])
    return embedding[:24], {
        "hydrophobic": float(mean[0]),
        "polar": float(mean[1]),
        "charged": float(mean[2] + mean[3]),
        "positive": float(mean[2]),
        "negative": float(mean[3]),
        "aromatic": float(mean[4]),
    }


def pocket_geometry_embedding(receptor_atoms: list[dict[str, object]], ligand_atoms: list[dict[str, object]]) -> tuple[np.ndarray, dict[str, float]]:
    if not receptor_atoms:
        return np.zeros((64,), dtype=float), {
            "rg": math.nan,
            "span": math.nan,
            "centroid_distance": math.nan,
        }
    coords = np.array([atom_coord(atom) for atom in receptor_atoms], dtype=float)
    ligand_coords = np.array([atom_coord(atom) for atom in ligand_atoms], dtype=float) if ligand_atoms else np.empty((0, 3))
    centroid = coords.mean(axis=0)
    centered = coords - centroid
    var = centered.var(axis=0)
    cov = np.cov(centered.T) if len(coords) >= 2 else np.zeros((3, 3), dtype=float)
    cov_upper = cov[np.triu_indices(3)]
    rg = float(np.sqrt(np.mean(np.sum(centered**2, axis=1)))) if len(coords) else math.nan
    span_vec = coords.max(axis=0) - coords.min(axis=0)
    span = float(np.linalg.norm(span_vec))
    skew = np.mean(centered**3, axis=0) / np.maximum(np.std(centered, axis=0) ** 3, 1e-6)
    ligand_centroid = ligand_coords.mean(axis=0) if len(ligand_coords) else centroid
    centroid_distance = float(np.linalg.norm(centroid - ligand_centroid))
    radial = np.linalg.norm(centered, axis=1)
    quantiles = np.quantile(radial, [0.1, 0.25, 0.5, 0.75, 0.9]) if len(radial) else np.zeros((5,), dtype=float)
    base = np.concatenate([
        centroid,
        var,
        cov_upper,
        skew,
        span_vec,
        np.array([rg, span, len(coords), centroid_distance], dtype=float),
        quantiles,
    ])
    embedding = np.zeros((64,), dtype=float)
    n = min(len(base), 64)
    embedding[:n] = base[:n]
    return embedding, {"rg": rg, "span": span, "centroid_distance": centroid_distance}


def pocket_embedding_features(
    ligand_atoms: list[dict[str, object]],
    receptor_atoms: list[dict[str, object]],
    cutoff: float,
) -> dict[str, float]:
    out: dict[str, float] = {field: math.nan for field in POCKET_SUMMARY_FEATURES}
    out.update({field: math.nan for field in POCKET_EMBEDDING_FIELDS})
    if not ligand_atoms or not receptor_atoms:
        return out
    ligand_coords = np.array([atom_coord(atom) for atom in ligand_atoms], dtype=float)
    nearby_atoms: list[dict[str, object]] = []
    residues: dict[tuple[str, str, str], str] = {}
    cutoff2 = cutoff * cutoff
    for atom in receptor_atoms:
        coord = atom_coord(atom)
        d2 = np.min(np.sum((ligand_coords - coord) ** 2, axis=1))
        if float(d2) <= cutoff2:
            nearby_atoms.append(atom)
            key = residue_key(atom)
            residues[key] = str(atom.get("residue_name", "")).upper()
    residue_names = list(residues.values())
    geom, geom_summary = pocket_geometry_embedding(nearby_atoms, ligand_atoms)
    phys, phys_summary = pocket_physchem_embedding(residue_names)
    out["pocket_atom_count_5a"] = float(len(nearby_atoms))
    out["pocket_residue_count_5a"] = float(len(residue_names))
    out["pocket_rg_5a"] = round(float(geom_summary["rg"]), 6) if not math.isnan(geom_summary["rg"]) else math.nan
    out["pocket_span_5a"] = round(float(geom_summary["span"]), 6) if not math.isnan(geom_summary["span"]) else math.nan
    out["pocket_centroid_distance_5a"] = round(float(geom_summary["centroid_distance"]), 6) if not math.isnan(geom_summary["centroid_distance"]) else math.nan
    out["pocket_hydrophobic_fraction_5a"] = round(phys_summary["hydrophobic"], 6) if not math.isnan(phys_summary["hydrophobic"]) else math.nan
    out["pocket_polar_fraction_5a"] = round(phys_summary["polar"], 6) if not math.isnan(phys_summary["polar"]) else math.nan
    out["pocket_charged_fraction_5a"] = round(phys_summary["charged"], 6) if not math.isnan(phys_summary["charged"]) else math.nan
    out["pocket_positive_fraction_5a"] = round(phys_summary["positive"], 6) if not math.isnan(phys_summary["positive"]) else math.nan
    out["pocket_negative_fraction_5a"] = round(phys_summary["negative"], 6) if not math.isnan(phys_summary["negative"]) else math.nan
    out["pocket_aromatic_fraction_5a"] = round(phys_summary["aromatic"], 6) if not math.isnan(phys_summary["aromatic"]) else math.nan
    for field, value in zip(POCKET_GEOMETRY_FIELDS, geom):
        out[field] = round(float(value), 6)
    for field, value in zip(POCKET_PHYSCHEM_FIELDS, phys):
        out[field] = round(float(value), 6)
    return out


def morgan_vector(smiles: str, n_bits: int = 2048) -> np.ndarray:
    mol = mol_from_smiles(smiles)
    arr = np.zeros((n_bits,), dtype=np.float32)
    if mol is None:
        return arr
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def morgan_fp(smiles: str):
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def tanimoto_similarity(smiles_a: str, smiles_b: str) -> float:
    fp_a = morgan_fp(smiles_a)
    fp_b = morgan_fp(smiles_b)
    if fp_a is None or fp_b is None:
        return math.nan
    return float(DataStructs.TanimotoSimilarity(fp_a, fp_b))


def coordinate_distance(row_a: dict[str, object], row_b: dict[str, object], xfield: str, yfield: str) -> float:
    ax, ay = safe_float(row_a.get(xfield)), safe_float(row_a.get(yfield))
    bx, by = safe_float(row_b.get(xfield)), safe_float(row_b.get(yfield))
    if any(math.isnan(v) for v in (ax, ay, bx, by)):
        return math.nan
    return math.sqrt((ax - bx) ** 2 + (ay - by) ** 2)


def mcs_summary(smiles_a: str, smiles_b: str) -> dict[str, object]:
    mol_a = mol_from_smiles(smiles_a)
    mol_b = mol_from_smiles(smiles_b)
    if mol_a is None or mol_b is None:
        return {"mcs_atoms": "", "mcs_bonds": "", "mcs_fraction_min": "", "mcs_fraction_max": ""}
    try:
        mcs = rdFMCS.FindMCS([mol_a, mol_b], timeout=5, ringMatchesRingOnly=True, completeRingsOnly=True)
        min_heavy = max(1, min(mol_a.GetNumHeavyAtoms(), mol_b.GetNumHeavyAtoms()))
        max_heavy = max(1, max(mol_a.GetNumHeavyAtoms(), mol_b.GetNumHeavyAtoms()))
        return {
            "mcs_atoms": int(mcs.numAtoms),
            "mcs_bonds": int(mcs.numBonds),
            "mcs_fraction_min": round(float(mcs.numAtoms) / min_heavy, 4),
            "mcs_fraction_max": round(float(mcs.numAtoms) / max_heavy, 4),
        }
    except Exception:
        return {"mcs_atoms": "", "mcs_bonds": "", "mcs_fraction_min": "", "mcs_fraction_max": ""}


def murcko_scaffold(smiles: str) -> str:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return ""
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold is None or scaffold.GetNumHeavyAtoms() == 0:
            return ""
        return Chem.MolToSmiles(scaffold, canonical=True)
    except Exception:
        return ""


def scaffold_fp(smiles: str):
    scaffold = murcko_scaffold(smiles)
    mol = mol_from_smiles(scaffold)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def molecule_props(smiles: str) -> dict[str, float | int | str]:
    mol = mol_from_smiles(smiles)
    scaffold = murcko_scaffold(smiles)
    if mol is None:
        return {
            "qed": math.nan,
            "aromatic_rings": math.nan,
            "murcko_scaffold": scaffold,
            "murcko_heavy_atoms": math.nan,
        }
    scaffold_mol = mol_from_smiles(scaffold)
    try:
        qed = float(QED.qed(mol))
    except Exception:
        qed = math.nan
    return {
        "qed": round(qed, 4) if not math.isnan(qed) else math.nan,
        "aromatic_rings": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "murcko_scaffold": scaffold,
        "murcko_heavy_atoms": scaffold_mol.GetNumHeavyAtoms() if scaffold_mol is not None else 0,
        "mw": round(float(Descriptors.MolWt(mol)), 3),
        "logp": round(float(Crippen.MolLogP(mol)), 3),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(mol)), 3),
        "rot_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
    }


def load_reference_scaffold_fps(path: Path | None, max_rows: int) -> list[object]:
    if path is None or not path.expanduser().exists():
        return []
    refs: list[object] = []
    seen: set[str] = set()
    for row in read_csv(path.expanduser())[:max_rows]:
        smiles = row.get("canonical_smiles") or row.get("smiles") or row.get("SMILES")
        scaffold = murcko_scaffold(canonical(smiles))
        if not scaffold or scaffold in seen:
            continue
        fp = scaffold_fp(scaffold)
        if fp is None:
            continue
        seen.add(scaffold)
        refs.append(fp)
    return refs


def max_scaffold_similarity(smiles: str, reference_fps: list[object]) -> float:
    fp = scaffold_fp(smiles)
    if fp is None or not reference_fps:
        return 0.0
    return round(max(DataStructs.BulkTanimotoSimilarity(fp, reference_fps)), 4)


def merge_analysis_scores(history_rows: list[dict[str, str]], analysis_rows: list[dict[str, str]]) -> None:
    by_can: dict[str, str] = {}
    for row in analysis_rows:
        smiles = row.get("canonical_smiles") or row.get("smiles") or row.get("mol_smiles")
        score = row.get("official_binding_score") or row.get("binding_score") or row.get("binding_scores")
        can = canonical(smiles)
        if can and score not in ("", None):
            by_can[can] = str(score)
    for row in history_rows:
        if row.get("official_binding_score") not in ("", None):
            continue
        can = canonical(row.get("canonical_smiles") or row.get("smiles") or row.get("input_smiles"))
        if can in by_can:
            row["official_binding_score"] = by_can[can]


def build_rows(args: argparse.Namespace) -> list[dict[str, object]]:
    history_csv = args.history_csv.expanduser()
    history_rows = read_csv(history_csv)
    if args.analysis_csv:
        merge_analysis_scores(history_rows, read_csv(args.analysis_csv.expanduser()))
    reference_fps = load_reference_scaffold_fps(args.reference_smiles_csv, args.max_reference_rows)
    receptor = args.receptor.expanduser() if args.receptor else None
    receptor_atoms = read_pdbqt_atoms(receptor) if receptor else []
    history_dir = history_csv.parent
    out: list[dict[str, object]] = []
    seen: set[str] = set()
    for row in history_rows:
        smiles = canonical(row.get("canonical_smiles") or row.get("smiles") or row.get("input_smiles"))
        if not smiles or smiles in seen:
            continue
        seen.add(smiles)
        props = molecule_props(smiles)
        record: dict[str, object] = {field: row.get(field, "") for field in HISTORY_FIELDS}
        record.update(props)
        record["canonical_smiles"] = smiles
        record["official_binding_score"] = safe_float(row.get("official_binding_score"))
        record["chembl_scaffold_similarity"] = max_scaffold_similarity(smiles, reference_fps)
        if receptor_atoms:
            pose_path = resolve_maybe_path(row.get("pose_path"), history_dir)
            ligand_atoms = parse_ligand_best_mode(pose_path)
            record.update(pocket_embedding_features(ligand_atoms, receptor_atoms, args.pocket_embedding_cutoff))
        for field in STRUCTURAL_FEATURES:
            if field in props or field == "chembl_scaffold_similarity":
                continue
            record[field] = safe_float(record.get(field))
        out.append(record)
    return out


def pca_coords(matrix: np.ndarray) -> np.ndarray:
    if len(matrix) == 0:
        return np.empty((0, 2))
    if len(matrix) == 1:
        return np.zeros((1, 2))
    return PCA(n_components=2, random_state=13).fit_transform(matrix)


def cosine_embedding_coords(matrix: np.ndarray, seed: int = 13) -> tuple[np.ndarray, str]:
    if len(matrix) == 0:
        return np.empty((0, 2)), "empty"
    if len(matrix) == 1:
        return np.zeros((1, 2)), "single_point"
    distances = pairwise_distances(matrix, metric="cosine")
    distances = np.nan_to_num(distances, nan=1.0, posinf=1.0, neginf=1.0)
    np.fill_diagonal(distances, 0.0)
    try:
        import umap  # type: ignore

        reducer = umap.UMAP(n_components=2, metric="precomputed", random_state=seed, min_dist=0.08, n_neighbors=min(25, max(2, len(matrix) - 1)))
        return reducer.fit_transform(distances), "cosine_umap"
    except Exception:
        try:
            reducer = MDS(n_components=2, dissimilarity="precomputed", random_state=seed, normalized_stress="auto")
        except TypeError:
            reducer = MDS(n_components=2, dissimilarity="precomputed", random_state=seed)
        return reducer.fit_transform(distances), "cosine_mds"


def add_spaces(rows: list[dict[str, object]], args: argparse.Namespace) -> dict[str, object]:
    smiles_matrix = np.vstack([morgan_vector(str(row["canonical_smiles"])) for row in rows])
    smiles_xy = pca_coords(smiles_matrix)
    cosine_xy, cosine_method = cosine_embedding_coords(smiles_matrix, args.seed)
    structural_fields = [
        field for field in STRUCTURAL_FEATURES
        if any(not math.isnan(safe_float(row.get(field))) for row in rows)
    ]
    structural_matrix = np.array([[safe_float(row.get(field)) for field in structural_fields] for row in rows], dtype=float)
    structural_xy = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), PCA(n_components=2, random_state=13)).fit_transform(structural_matrix)
    for row, xy1, xy_cos, xy2 in zip(rows, smiles_xy, cosine_xy, structural_xy):
        row["pure_smiles_pc1"] = round(float(xy1[0]), 6)
        row["pure_smiles_pc2"] = round(float(xy1[1]), 6)
        row["pure_smiles_cosine1"] = round(float(xy_cos[0]), 6)
        row["pure_smiles_cosine2"] = round(float(xy_cos[1]), 6)
        row["structural_interaction_pc1"] = round(float(xy2[0]), 6)
        row["structural_interaction_pc2"] = round(float(xy2[1]), 6)
    return {
        "pure_smiles_space_method": cosine_method,
        "pure_smiles_space_note": "Morgan fingerprint visualized from a cosine distance matrix; PCA columns are retained only for backward compatibility.",
        "structural_space_method": "standardized_raw_features_pca",
        "structural_space_features": structural_fields,
    }


def train_predict(rows: list[dict[str, object]], args: argparse.Namespace) -> dict[str, object]:
    scored = [row for row in rows if not math.isnan(safe_float(row.get("official_binding_score")))]
    rng = random.Random(args.seed)
    rng.shuffle(scored)
    train_n = min(args.train_size, max(1, len(scored) - args.test_size))
    test_n = min(args.test_size, max(0, len(scored) - train_n))
    train_rows = scored[:train_n]
    test_rows = scored[train_n : train_n + test_n]
    # use raw structural features only; PCA coordinates are for visualization, not model input
    feature_fields = [
        field for field in STRUCTURAL_FEATURES
        if any(not math.isnan(safe_float(row.get(field))) for row in train_rows)
    ]
    x_train = np.array([[safe_float(row.get(field)) for field in feature_fields] for row in train_rows], dtype=float)
    y_train = np.array([safe_float(row.get("official_binding_score")) for row in train_rows], dtype=float)
    model = make_pipeline(
        SimpleImputer(strategy="median"),
        RandomForestRegressor(n_estimators=args.trees, random_state=args.seed, min_samples_leaf=5),
    )
    model.fit(x_train, y_train)
    # 5-fold CV on training set
    cv = KFold(n_splits=5, shuffle=True, random_state=args.seed)
    cv_r2 = cross_val_score(model, x_train, y_train, cv=cv, scoring="r2")
    cv_mae = -cross_val_score(model, x_train, y_train, cv=cv, scoring="neg_mean_absolute_error")
    x_all = np.array([[safe_float(row.get(field)) for field in feature_fields] for row in rows], dtype=float)
    preds = model.predict(x_all)
    for row, pred in zip(rows, preds):
        row["predicted_binding_score"] = round(float(pred), 6)
        row["score_set"] = "unscored"
    for row in train_rows:
        row["score_set"] = "train"
    for row in test_rows:
        row["score_set"] = "test"
    # adaptive score bins from scored data distribution
    official_scores = [safe_float(r.get("official_binding_score")) for r in scored]
    valid_scores = sorted([v for v in official_scores if not math.isnan(v)])
    if len(valid_scores) >= 5:
        n = len(valid_scores)
        score_bins = [round(valid_scores[int(n * p)], 4) for p in [0.2, 0.4, 0.6, 0.8]]
    else:
        score_bins = [0.30, 0.35, 0.38, 0.40]
    metrics: dict[str, object] = {
        "scored_rows": len(scored),
        "train_rows": len(train_rows),
        "test_rows": len(test_rows),
        "feature_fields": feature_fields,
        "score_bins": score_bins,
        "cv_r2_mean": round(float(np.mean(cv_r2)), 6),
        "cv_r2_std": round(float(np.std(cv_r2)), 6),
        "cv_mae_mean": round(float(np.mean(cv_mae)), 6),
        "cv_mae_std": round(float(np.std(cv_mae)), 6),
    }
    if test_rows:
        x_test = np.array([[safe_float(row.get(field)) for field in feature_fields] for row in test_rows], dtype=float)
        y_test = np.array([safe_float(row.get("official_binding_score")) for row in test_rows], dtype=float)
        pred = model.predict(x_test)
        metrics["test_mae"] = round(float(mean_absolute_error(y_test, pred)), 6)
        metrics["test_r2"] = round(float(r2_score(y_test, pred)), 6) if len(test_rows) >= 2 else ""
        for row, value in zip(test_rows, pred):
            row["prediction_residual"] = round(safe_float(row.get("official_binding_score")) - float(value), 6)
    rf = model.named_steps.get("randomforestregressor")
    if rf is not None and hasattr(rf, "feature_importances_"):
        importances = [
            {"feature": field, "importance": round(float(value), 6)}
            for field, value in zip(feature_fields, rf.feature_importances_)
        ]
        metrics["feature_importance"] = sorted(importances, key=lambda item: item["importance"], reverse=True)
        try:
            import shap  # type: ignore

            imputer = model.named_steps.get("simpleimputer")
            if imputer is None:
                raise RuntimeError("pipeline imputer not found")
            sample_n = min(max(1, args.shap_sample_size), len(x_train))
            x_shap = imputer.transform(x_train[:sample_n])
            explainer = shap.TreeExplainer(rf)
            shap_values = explainer.shap_values(x_shap)
            shap_arr = np.asarray(shap_values)
            if shap_arr.ndim == 3:
                shap_arr = shap_arr[0]
            mean_abs = np.mean(np.abs(shap_arr), axis=0)
            shap_rows = [
                {"feature": field, "mean_abs_shap": round(float(value), 6)}
                for field, value in zip(feature_fields, mean_abs)
            ]
            metrics["shap_status"] = "ok"
            metrics["shap_sample_rows"] = sample_n
            metrics["shap_mean_abs"] = sorted(shap_rows, key=lambda item: item["mean_abs_shap"], reverse=True)
        except Exception as exc:
            metrics["shap_status"] = f"unavailable: {exc.__class__.__name__}"
            metrics["shap_install_hint"] = "/mnt/l/WSL/softwares/conda_envs/md/bin/python -m pip install shap"
            if args.require_shap:
                raise RuntimeError(
                    "--require-shap was set, but SHAP could not run. Install it with: "
                    "/mnt/l/WSL/softwares/conda_envs/md/bin/python -m pip install shap"
                ) from exc
    test_r2 = safe_float(metrics.get("test_r2"))
    cv_r2 = safe_float(metrics.get("cv_r2_mean"))
    if math.isnan(test_r2) and math.isnan(cv_r2):
        interpretation = "not_enough_test_rows"
    elif (math.isnan(test_r2) or test_r2 < 0.3) and (math.isnan(cv_r2) or cv_r2 < 0.3):
        interpretation = "poor_predictor_use_for_diagnosis_only"
    elif (math.isnan(test_r2) or test_r2 < 0.6) and (math.isnan(cv_r2) or cv_r2 < 0.6):
        interpretation = "weak_predictor_hypothesis_only"
    else:
        interpretation = "usable_screening_signal_needs_external_validation"
    metrics["model_interpretation"] = interpretation
    metrics["prediction_warning"] = (
        "Do not treat predicted_binding_score as authoritative when test_r2 or cv_r2 is low. "
        "Use the dashboard to inspect chemical space, pocket-contact space, QED/descriptors, "
        "ChEMBL scaffold similarity, and official-score neighborhoods."
    )
    return metrics


def descriptor_delta(row_a: dict[str, object], row_b: dict[str, object], field: str) -> float:
    a = safe_float(row_a.get(field))
    b = safe_float(row_b.get(field))
    if math.isnan(a) or math.isnan(b):
        return math.nan
    return b - a


def penalty_hypotheses(row_a: dict[str, object], row_b: dict[str, object]) -> tuple[str, str]:
    score_a = safe_float(row_a.get("official_binding_score"))
    score_b = safe_float(row_b.get("official_binding_score"))
    if math.isnan(score_a) or math.isnan(score_b):
        return "", ""
    high, low = (row_a, row_b) if score_a >= score_b else (row_b, row_a)
    reasons: list[str] = []
    evidence: list[str] = []

    affinity_high = safe_float(high.get("affinity_kcal_mol"))
    affinity_low = safe_float(low.get("affinity_kcal_mol"))
    if not math.isnan(affinity_high) and not math.isnan(affinity_low):
        if affinity_low < affinity_high - 0.5:
            reasons.append("vina_false_positive_or_missing_physics")
            evidence.append(f"low_official_has_better_vina({affinity_low:.2f}<{affinity_high:.2f})")
        elif affinity_low > affinity_high + 1.0:
            reasons.append("vina_affinity_consistent_with_cliff")
            evidence.append(f"low_official_has_weaker_vina({affinity_low:.2f}>{affinity_high:.2f})")

    hbond_delta = descriptor_delta(high, low, "hbond_count")
    if not math.isnan(hbond_delta) and hbond_delta <= -3:
        reasons.append("lost_hbond_network")
        evidence.append(f"hbond_delta_low_minus_high={hbond_delta:.1f}")
    if not math.isnan(hbond_delta) and hbond_delta >= 3:
        reasons.append("more_hbonds_but_score_lower_possible_bad_geometry_or_desolvation")
        evidence.append(f"hbond_delta_low_minus_high=+{hbond_delta:.1f}")

    hydrophobic_delta = descriptor_delta(high, low, "hydrophobic_count")
    if not math.isnan(hydrophobic_delta) and hydrophobic_delta <= -5:
        reasons.append("lost_hydrophobic_contacts")
        evidence.append(f"hydrophobic_delta_low_minus_high={hydrophobic_delta:.1f}")
    if not math.isnan(hydrophobic_delta) and hydrophobic_delta >= 5:
        reasons.append("hydrophobic_gain_not_translating_possible_wrong_region")
        evidence.append(f"hydrophobic_delta_low_minus_high=+{hydrophobic_delta:.1f}")

    vdw_delta = descriptor_delta(high, low, "vdw_contact_count")
    if not math.isnan(vdw_delta) and vdw_delta >= 35:
        reasons.append("more_vdw_but_score_lower_possible_overpacking_or_pose_artifact")
        evidence.append(f"vdw_delta_low_minus_high=+{vdw_delta:.1f}")
    if not math.isnan(vdw_delta) and vdw_delta <= -25:
        reasons.append("lost_vdw_surface_contact")
        evidence.append(f"vdw_delta_low_minus_high={vdw_delta:.1f}")

    inner_delta = descriptor_delta(high, low, "inner_rmsd")
    whole_delta = descriptor_delta(high, low, "whole_rmsd")
    if (not math.isnan(inner_delta) and inner_delta >= 1.5) or (not math.isnan(whole_delta) and whole_delta >= 1.5):
        reasons.append("less_consistent_pose_family")
        evidence.append(f"inner_rmsd_delta={inner_delta:.2f};whole_rmsd_delta={whole_delta:.2f}")

    tpsa_delta = descriptor_delta(high, low, "tpsa")
    logp_delta = descriptor_delta(high, low, "logp")
    rot_delta = descriptor_delta(high, low, "rot_bonds")
    if not math.isnan(tpsa_delta) and abs(tpsa_delta) >= 25:
        reasons.append("polar_surface_or_desolvation_shift")
        evidence.append(f"tpsa_delta_low_minus_high={tpsa_delta:.1f}")
    if not math.isnan(logp_delta) and abs(logp_delta) >= 0.8:
        reasons.append("lipophilicity_shift")
        evidence.append(f"logp_delta_low_minus_high={logp_delta:.2f}")
    if not math.isnan(rot_delta) and rot_delta >= 3:
        reasons.append("higher_conformational_entropy_penalty")
        evidence.append(f"rot_bonds_delta_low_minus_high=+{rot_delta:.1f}")

    if not reasons:
        reasons.append("unexplained_by_current_features_check_pose_or_assay")
        evidence.append("current_descriptor_contact_features_do_not_explain_cliff")
    return ";".join(dict.fromkeys(reasons)), ";".join(evidence)


def activity_cliff_report(rows: list[dict[str, object]], args: argparse.Namespace) -> list[dict[str, object]]:
    scored = [row for row in rows if finite(row.get("official_binding_score")) and row.get("canonical_smiles")]
    report: list[dict[str, object]] = []
    for i, row_a in enumerate(scored):
        for row_b in scored[i + 1:]:
            score_a = safe_float(row_a.get("official_binding_score"))
            score_b = safe_float(row_b.get("official_binding_score"))
            score_delta = abs(score_a - score_b)
            if score_delta < args.activity_cliff_min_score_delta:
                continue
            smiles_a = str(row_a.get("canonical_smiles", ""))
            smiles_b = str(row_b.get("canonical_smiles", ""))
            tanimoto = tanimoto_similarity(smiles_a, smiles_b)
            pure_dist = coordinate_distance(row_a, row_b, "pure_smiles_cosine1", "pure_smiles_cosine2")
            structural_dist = coordinate_distance(row_a, row_b, "structural_interaction_pc1", "structural_interaction_pc2")
            mcs = mcs_summary(smiles_a, smiles_b)
            mcs_fraction = safe_float(mcs.get("mcs_fraction_min"))
            chemically_close = (
                (not math.isnan(tanimoto) and tanimoto >= args.activity_cliff_min_tanimoto)
                or (not math.isnan(mcs_fraction) and mcs_fraction >= args.activity_cliff_min_mcs_fraction)
                or (not math.isnan(pure_dist) and pure_dist <= args.activity_cliff_max_pure_distance)
            )
            space_close = not math.isnan(structural_dist) and structural_dist <= args.activity_cliff_max_structural_distance
            if not chemically_close and not space_close:
                continue
            affinity_a = safe_float(row_a.get("affinity_kcal_mol"))
            affinity_b = safe_float(row_b.get("affinity_kcal_mol"))
            affinity_delta = abs(affinity_a - affinity_b) if not math.isnan(affinity_a) and not math.isnan(affinity_b) else math.nan
            high = row_a if score_a >= score_b else row_b
            low = row_b if high is row_a else row_a
            high_affinity = safe_float(high.get("affinity_kcal_mol"))
            low_affinity = safe_float(low.get("affinity_kcal_mol"))
            vina_direction = "unknown"
            if not math.isnan(high_affinity) and not math.isnan(low_affinity):
                if low_affinity < high_affinity - 0.5:
                    vina_direction = "opposite_low_official_better_vina"
                elif low_affinity > high_affinity + 0.5:
                    vina_direction = "same_low_official_weaker_vina"
                else:
                    vina_direction = "similar_affinity"
            if (
                args.activity_cliff_require_vina_signal
                and vina_direction not in {"opposite_low_official_better_vina", "similar_affinity"}
                and (math.isnan(affinity_delta) or affinity_delta > args.activity_cliff_max_affinity_delta)
            ):
                continue
            hypotheses, evidence = penalty_hypotheses(row_a, row_b)
            report.append(
                {
                    "seq_id_a": row_a.get("seq_id", ""),
                    "seq_id_b": row_b.get("seq_id", ""),
                    "higher_score_seq_id": high.get("seq_id", ""),
                    "lower_score_seq_id": low.get("seq_id", ""),
                    "official_score_a": round(score_a, 6),
                    "official_score_b": round(score_b, 6),
                    "official_score_delta": round(score_delta, 6),
                    "affinity_a": row_a.get("affinity_kcal_mol", ""),
                    "affinity_b": row_b.get("affinity_kcal_mol", ""),
                    "affinity_delta": round(float(affinity_delta), 4) if not math.isnan(affinity_delta) else "",
                    "vina_direction": vina_direction,
                    "morgan_tanimoto": round(float(tanimoto), 4) if not math.isnan(tanimoto) else "",
                    **mcs,
                    "pure_smiles_distance": round(float(pure_dist), 4) if not math.isnan(pure_dist) else "",
                    "structural_interaction_distance": round(float(structural_dist), 4) if not math.isnan(structural_dist) else "",
                    "penalty_hypotheses": hypotheses,
                    "evidence": evidence,
                    "smiles_a": smiles_a,
                    "smiles_b": smiles_b,
                }
            )
    report.sort(
        key=lambda row: (
            -safe_float(row.get("official_score_delta")),
            safe_float(row.get("pure_smiles_distance")) if finite(row.get("pure_smiles_distance")) else 999.0,
            safe_float(row.get("structural_interaction_distance")) if finite(row.get("structural_interaction_distance")) else 999.0,
        )
    )
    return report[: args.activity_cliff_top_n]


def group_stats(rows: list[dict[str, object]], feature_fields: list[str], prefix: str) -> dict[str, object]:
    out: dict[str, object] = {"group": prefix, "count": len(rows)}
    for field in feature_fields:
        values = non_nan_values(rows, field)
        if not values:
            continue
        out[f"{field}_mean"] = round(float(np.mean(values)), 6)
        out[f"{field}_std"] = round(float(np.std(values)), 6)
        out[f"{field}_median"] = round(float(np.median(values)), 6)
    return out


def compare_high_low_groups(rows: list[dict[str, object]], args: argparse.Namespace) -> tuple[list[dict[str, object]], dict[str, object], list[dict[str, object]]]:
    scored = [row for row in rows if finite(row.get("official_binding_score")) and row.get("canonical_smiles")]
    high = [row for row in scored if safe_float(row.get("official_binding_score")) > args.high_score_threshold]
    low = [row for row in scored if 0.0 <= safe_float(row.get("official_binding_score")) <= args.low_score_max]
    compare_fields = [
        field for field in [
            "affinity_kcal_mol",
            "inner_rmsd",
            "whole_rmsd",
            "hbond_count",
            "hydrophobic_count",
            "pi_contact_count",
            "ch_pi_count",
            "logp",
            "hbd",
            "hba",
            "tpsa",
            "rot_bonds",
            "qed",
            "aromatic_rings",
            "chembl_scaffold_similarity",
            "contact_residue_count_4a",
            "surface_contact_fraction_4a",
            "receptor_residue_coverage_4a",
            "ligand_atom_contact_fraction_4a",
            "min_ligand_receptor_distance",
            "pocket_atom_count_5a",
            "pocket_residue_count_5a",
            "pocket_rg_5a",
            "pocket_span_5a",
            "pocket_centroid_distance_5a",
            "pocket_hydrophobic_fraction_5a",
            "pocket_polar_fraction_5a",
            "pocket_charged_fraction_5a",
            "pocket_aromatic_fraction_5a",
            "pure_smiles_cosine1",
            "pure_smiles_cosine2",
            "structural_interaction_pc1",
            "structural_interaction_pc2",
        ]
        if any(not math.isnan(safe_float(row.get(field))) for row in scored)
    ]
    high_stats = group_stats(high, compare_fields, "high")
    low_stats = group_stats(low, compare_fields, "low")

    rows_out: list[dict[str, object]] = []
    for field in compare_fields:
        high_vals = non_nan_values(high, field)
        low_vals = non_nan_values(low, field)
        if not high_vals or not low_vals:
            continue
        high_mean = float(np.mean(high_vals))
        low_mean = float(np.mean(low_vals))
        pooled = float(np.sqrt((np.var(high_vals) + np.var(low_vals)) / 2.0)) if len(high_vals) > 1 and len(low_vals) > 1 else math.nan
        effect = (high_mean - low_mean) / pooled if pooled and not math.isnan(pooled) and pooled > 1e-12 else math.nan
        rows_out.append(
            {
                "feature": field,
                "high_mean": round(high_mean, 6),
                "low_mean": round(low_mean, 6),
                "delta_high_minus_low": round(high_mean - low_mean, 6),
                "abs_delta": round(abs(high_mean - low_mean), 6),
                "high_std": round(float(np.std(high_vals)), 6),
                "low_std": round(float(np.std(low_vals)), 6),
                "effect_size_approx": round(float(effect), 6) if not math.isnan(effect) else "",
                "high_non_nan": len(high_vals),
                "low_non_nan": len(low_vals),
            }
        )
    rows_out.sort(key=lambda item: item["abs_delta"], reverse=True)

    group_rows: list[dict[str, object]] = []
    for row in high:
        group_rows.append(
            {
                "group": "high",
                "seq_id": row.get("seq_id", ""),
                "nickname": row.get("nickname", ""),
                "official_binding_score": row.get("official_binding_score", ""),
                "predicted_binding_score": row.get("predicted_binding_score", ""),
                "canonical_smiles": row.get("canonical_smiles", ""),
                "pure_smiles_cosine1": row.get("pure_smiles_cosine1", ""),
                "pure_smiles_cosine2": row.get("pure_smiles_cosine2", ""),
                "structural_interaction_pc1": row.get("structural_interaction_pc1", ""),
                "structural_interaction_pc2": row.get("structural_interaction_pc2", ""),
                "affinity_kcal_mol": row.get("affinity_kcal_mol", ""),
            }
        )
    for row in low:
        group_rows.append(
            {
                "group": "low",
                "seq_id": row.get("seq_id", ""),
                "nickname": row.get("nickname", ""),
                "official_binding_score": row.get("official_binding_score", ""),
                "predicted_binding_score": row.get("predicted_binding_score", ""),
                "canonical_smiles": row.get("canonical_smiles", ""),
                "pure_smiles_cosine1": row.get("pure_smiles_cosine1", ""),
                "pure_smiles_cosine2": row.get("pure_smiles_cosine2", ""),
                "structural_interaction_pc1": row.get("structural_interaction_pc1", ""),
                "structural_interaction_pc2": row.get("structural_interaction_pc2", ""),
                "affinity_kcal_mol": row.get("affinity_kcal_mol", ""),
            }
        )

    center_distance = math.nan
    if high and low:
        high_center = np.array([
            float(np.mean([safe_float(row.get("pure_smiles_cosine1")) for row in high if not math.isnan(safe_float(row.get("pure_smiles_cosine1")))])) if any(not math.isnan(safe_float(row.get("pure_smiles_cosine1"))) for row in high) else 0.0,
            float(np.mean([safe_float(row.get("pure_smiles_cosine2")) for row in high if not math.isnan(safe_float(row.get("pure_smiles_cosine2")))])) if any(not math.isnan(safe_float(row.get("pure_smiles_cosine2"))) for row in high) else 0.0,
            float(np.mean([safe_float(row.get("structural_interaction_pc1")) for row in high if not math.isnan(safe_float(row.get("structural_interaction_pc1")))])) if any(not math.isnan(safe_float(row.get("structural_interaction_pc1"))) for row in high) else 0.0,
            float(np.mean([safe_float(row.get("structural_interaction_pc2")) for row in high if not math.isnan(safe_float(row.get("structural_interaction_pc2")))])) if any(not math.isnan(safe_float(row.get("structural_interaction_pc2"))) for row in high) else 0.0,
        ], dtype=float)
        low_center = np.array([
            float(np.mean([safe_float(row.get("pure_smiles_cosine1")) for row in low if not math.isnan(safe_float(row.get("pure_smiles_cosine1")))])) if any(not math.isnan(safe_float(row.get("pure_smiles_cosine1"))) for row in low) else 0.0,
            float(np.mean([safe_float(row.get("pure_smiles_cosine2")) for row in low if not math.isnan(safe_float(row.get("pure_smiles_cosine2")))])) if any(not math.isnan(safe_float(row.get("pure_smiles_cosine2"))) for row in low) else 0.0,
            float(np.mean([safe_float(row.get("structural_interaction_pc1")) for row in low if not math.isnan(safe_float(row.get("structural_interaction_pc1")))])) if any(not math.isnan(safe_float(row.get("structural_interaction_pc1"))) for row in low) else 0.0,
            float(np.mean([safe_float(row.get("structural_interaction_pc2")) for row in low if not math.isnan(safe_float(row.get("structural_interaction_pc2")))])) if any(not math.isnan(safe_float(row.get("structural_interaction_pc2"))) for row in low) else 0.0,
        ], dtype=float)
        center_distance = float(np.linalg.norm(high_center - low_center))

    summary = {
        "high_threshold": args.high_score_threshold,
        "low_max": args.low_score_max,
        "high_count": len(high),
        "low_count": len(low),
        "feature_fields_compared": compare_fields,
        "embedding_center_distance_4d": round(center_distance, 6) if not math.isnan(center_distance) else "",
        "high_stats": high_stats,
        "low_stats": low_stats,
        "note": "High/low group comparison is descriptive, not causal. Use it to inspect which pocket/pose/descriptor axes separate strong and weak binders.",
    }
    return rows_out, summary, group_rows


def score_color(value: object, bins: list[float] | None = None) -> str:
    score = safe_float(value)
    if math.isnan(score):
        return "#9aa0a6"
    if bins is None or len(bins) < 4:
        bins = [0.30, 0.35, 0.38, 0.40]
    colors = ["#f28482", "#f6bd60", "#84a59d", "#5dade2", "#7bc96f"]
    for threshold, color in zip(bins, colors):
        if score < threshold:
            return color
    return colors[-1]


def render_space_html(rows: list[dict[str, object]], xfield: str, yfield: str, title: str, path: Path, score_bins: list[float] | None = None) -> None:
    xs = [safe_float(row.get(xfield)) for row in rows]
    ys = [safe_float(row.get(yfield)) for row in rows]
    finite = [(x, y) for x, y in zip(xs, ys) if not math.isnan(x) and not math.isnan(y)]
    if not finite:
        path.write_text("<html><body>No points.</body></html>", encoding="utf-8")
        return
    min_x, max_x = min(x for x, _ in finite), max(x for x, _ in finite)
    min_y, max_y = min(y for _, y in finite), max(y for _, y in finite)
    width, height, pad = 980, 700, 48
    def sx(x: float) -> float:
        return pad + (x - min_x) / (max_x - min_x or 1.0) * (width - 2 * pad)
    def sy(y: float) -> float:
        return height - pad - (y - min_y) / (max_y - min_y or 1.0) * (height - 2 * pad)
    circles: list[str] = []
    for row in rows:
        x, y = safe_float(row.get(xfield)), safe_float(row.get(yfield))
        if math.isnan(x) or math.isnan(y):
            continue
        scored = not math.isnan(safe_float(row.get("official_binding_score")))
        r = 6 if scored else 3
        stroke = "#111827" if scored else "none"
        title_text = (
            f"{row.get('seq_id','')} score={row.get('official_binding_score','')} "
            f"pred={row.get('predicted_binding_score','')} scaffold={row.get('murcko_scaffold','')} "
            f"smiles={row.get('canonical_smiles','')}"
        )
        circles.append(
            f"<circle cx='{sx(x):.2f}' cy='{sy(y):.2f}' r='{r}' fill='{score_color(row.get('official_binding_score'), score_bins)}' "
            f"stroke='{stroke}' stroke-width='1.2' opacity='0.82'><title>{title_text}</title></circle>"
        )
    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{title}</title>
<style>body{{font-family:Arial,sans-serif;margin:24px;color:#111827}} svg{{border:1px solid #d1d5db;background:#fff}} .note{{color:#4b5563}}</style>
</head><body>
<h1>{title}</h1>
<p class="note">Color = official binding score bin; gray = unscored. Large outlined points have teacher binding scores.</p>
<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>
{''.join(circles)}
</svg>
</body></html>"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--history-csv", type=Path, default=Path("~/vina_task2/dock_history/dock_history.csv"))
    parser.add_argument("--analysis-csv", type=Path, help="Optional binding_score_analysis.csv to merge scores")
    parser.add_argument("--reference-smiles-csv", type=Path, default=Path("~/vina_task2/references/chembl_kinase_smiles.csv"))
    parser.add_argument("--receptor", type=Path, default=Path("~/vina_task2/vina_bin/target.pdbqt"), help="Optional receptor PDBQT/PDB used to build ligand-local 5A pocket embeddings")
    parser.add_argument("--pocket-embedding-cutoff", type=float, default=5.0, help="Receptor atoms within this Angstrom cutoff of the best ligand pose define the local pocket embedding")
    parser.add_argument("--outdir", type=Path, default=Path("~/vina_task2/score_space_model"))
    parser.add_argument("--train-size", type=int, default=60)
    parser.add_argument("--test-size", type=int, default=20)
    parser.add_argument("--trees", type=int, default=400)
    parser.add_argument("--max-reference-rows", type=int, default=5000)
    parser.add_argument("--shap-sample-size", type=int, default=80, help="Rows sampled from the training set for optional SHAP TreeExplainer output")
    parser.add_argument("--require-shap", action="store_true", help="Fail if SHAP is not installed or cannot explain the trained model")
    parser.add_argument("--activity-cliff-top-n", type=int, default=200, help="Maximum activity-cliff molecule pairs to write")
    parser.add_argument("--activity-cliff-min-score-delta", type=float, default=0.15, help="Minimum official binding score delta for activity-cliff pairs")
    parser.add_argument("--activity-cliff-min-tanimoto", type=float, default=0.55, help="Morgan Tanimoto threshold for chemical-neighbor cliffs")
    parser.add_argument("--activity-cliff-min-mcs-fraction", type=float, default=0.80, help="MCS fraction threshold for scaffold-neighbor cliffs")
    parser.add_argument("--activity-cliff-max-pure-distance", type=float, default=0.45, help="Maximum Pure-SMILES cosine-space distance for visual-neighbor cliffs")
    parser.add_argument("--activity-cliff-max-structural-distance", type=float, default=0.50, help="Maximum structural-interaction PCA distance for interaction-neighbor cliffs")
    parser.add_argument("--activity-cliff-max-affinity-delta", type=float, default=1.0, help="Vina affinity delta treated as similar when filtering cliffs")
    parser.add_argument("--activity-cliff-require-vina-signal", action="store_true", help="Only report cliffs where Vina affinity is similar or points in the opposite direction")
    parser.add_argument("--high-score-threshold", type=float, default=0.4, help="High official binding score threshold for group comparison")
    parser.add_argument("--low-score-max", type=float, default=0.25, help="Low official binding score upper bound for group comparison")
    parser.add_argument("--seed", type=int, default=13)
    args = parser.parse_args()

    rows = build_rows(args)
    space_metrics = add_spaces(rows, args)
    metrics = train_predict(rows, args)
    metrics.update(space_metrics)
    metrics["pocket_embedding"] = {
        "receptor": str(args.receptor.expanduser()) if args.receptor else "",
        "cutoff_angstrom": args.pocket_embedding_cutoff,
        "summary_features": POCKET_SUMMARY_FEATURES,
        "geometry_embedding_fields": POCKET_GEOMETRY_FIELDS,
        "physchem_embedding_fields": POCKET_PHYSCHEM_FIELDS,
    }
    outdir = args.outdir.expanduser().resolve()
    fields = [
        "seq_id", "nickname", "score_set", "official_binding_score", "predicted_binding_score",
        "prediction_residual",
        "canonical_smiles", "murcko_scaffold", "chembl_scaffold_similarity",
        "pure_smiles_cosine1", "pure_smiles_cosine2", "pure_smiles_pc1", "pure_smiles_pc2",
        "structural_interaction_pc1", "structural_interaction_pc2",
        *STRUCTURAL_FEATURES, *POCKET_EMBEDDING_FIELDS,
    ]
    write_csv(outdir / "binding_score_feature_matrix.csv", rows, fields)
    shap_rows = metrics.get("shap_mean_abs")
    if isinstance(shap_rows, list) and shap_rows:
        write_csv(outdir / "shap_feature_importance.csv", shap_rows, ["feature", "mean_abs_shap"])
    cliff_rows = activity_cliff_report(rows, args)
    cliff_fields = [
        "seq_id_a", "seq_id_b", "higher_score_seq_id", "lower_score_seq_id",
        "official_score_a", "official_score_b", "official_score_delta",
        "affinity_a", "affinity_b", "affinity_delta", "vina_direction",
        "morgan_tanimoto", "mcs_atoms", "mcs_bonds", "mcs_fraction_min", "mcs_fraction_max",
        "pure_smiles_distance", "structural_interaction_distance",
        "penalty_hypotheses", "evidence", "smiles_a", "smiles_b",
    ]
    write_csv(outdir / "activity_cliff_report.csv", cliff_rows, cliff_fields)
    metrics["activity_cliff_rows"] = len(cliff_rows)
    metrics["activity_cliff_note"] = (
        "Rule-based hypotheses only. Use these pairs to inspect pose families, local contacts, "
        "conformational/desolvation penalties, or assay-specific effects."
    )
    group_compare_rows, group_summary, group_members = compare_high_low_groups(rows, args)
    group_fields = [
        "feature", "high_mean", "low_mean", "delta_high_minus_low", "abs_delta",
        "high_std", "low_std", "effect_size_approx", "high_non_nan", "low_non_nan",
    ]
    write_csv(outdir / "high_low_group_feature_comparison.csv", group_compare_rows, group_fields)
    write_csv(
        outdir / "high_low_group_members.csv",
        group_members,
        [
            "group", "seq_id", "nickname", "official_binding_score", "predicted_binding_score",
            "affinity_kcal_mol", "pure_smiles_cosine1", "pure_smiles_cosine2",
            "structural_interaction_pc1", "structural_interaction_pc2", "canonical_smiles",
        ],
    )
    with (outdir / "high_low_group_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(group_summary, handle, indent=2)
    metrics["high_low_group_summary"] = {
        "high_count": group_summary.get("high_count", 0),
        "low_count": group_summary.get("low_count", 0),
        "embedding_center_distance_4d": group_summary.get("embedding_center_distance_4d", ""),
    }
    with (outdir / "binding_score_model_metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)
    score_bins = metrics.get("score_bins")
    render_space_html(rows, "pure_smiles_cosine1", "pure_smiles_cosine2", "Pure-SMILES-space (cosine UMAP/MDS)", outdir / "pure_smiles_space.html", score_bins)
    render_space_html(rows, "structural_interaction_pc1", "structural_interaction_pc2", "structural-interaction-space", outdir / "structural_interaction_space.html", score_bins)
    print(
        f"rows={len(rows)} scored={metrics['scored_rows']} train={metrics['train_rows']} test={metrics['test_rows']} "
        f"test_mae={metrics.get('test_mae', '')} test_r2={metrics.get('test_r2', '')} "
        f"activity_cliffs={len(cliff_rows)} high={group_summary.get('high_count', 0)} "
        f"low={group_summary.get('low_count', 0)} wrote={outdir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
