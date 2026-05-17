#!/usr/bin/env python
"""Build Pure-SMILES (cosine-UMAP) and structural-interaction spaces, then learn binding score.

Improvements over original:
- Pure-SMILES space uses cosine-distance UMAP instead of PCA
- Pocket geometry + physicochemical embeddings extracted from docked pose + receptor
- Receptor pocket features added to RF model and high/low group comparison
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger, rdBase
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, QED, rdFMCS, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import cross_val_score, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

try:
    import umap
    UMAP_AVAILABLE = True
except Exception:
    UMAP_AVAILABLE = False

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
]

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
    # Pocket geometry + physicochemical embeddings from S1 workflow
    "pocket_geom_pc1",
    "pocket_geom_pc2",
    "pocket_phys_pc1",
    "pocket_phys_pc2",
]

# ====================== S1 Pocket Embedding Constants ======================
AA_PHYSCHEM = {
    'A': [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    'R': [0.0, 1.0, 1.0, 0.8, 0.8, 1.0],
    'N': [0.0, 0.0, 1.0, 0.6, 0.6, 1.0],
    'D': [0.0, -1.0, 1.0, 0.5, 0.4, 1.0],
    'C': [0.5, 0.0, 0.5, 0.7, 0.3, 0.5],
    'Q': [0.0, 0.0, 1.0, 0.7, 0.7, 1.0],
    'E': [0.0, -1.0, 1.0, 0.6, 0.5, 1.0],
    'G': [0.0, 0.0, 0.0, 0.3, 1.0, 0.0],
    'H': [0.0, 0.5, 1.0, 0.7, 0.6, 1.0],
    'I': [1.0, 0.0, 0.0, 0.9, 0.2, 0.0],
    'L': [1.0, 0.0, 0.0, 0.9, 0.2, 0.0],
    'K': [0.0, 1.0, 1.0, 0.8, 0.7, 1.0],
    'M': [0.8, 0.0, 0.2, 0.8, 0.4, 0.3],
    'F': [0.9, 0.0, 0.1, 0.9, 0.3, 0.2],
    'P': [0.2, 0.0, 0.3, 0.6, 0.9, 0.1],
    'S': [0.1, 0.0, 0.8, 0.4, 0.5, 0.8],
    'T': [0.2, 0.0, 0.7, 0.5, 0.4, 0.7],
    'W': [0.7, 0.0, 0.3, 1.0, 0.5, 0.4],
    'Y': [0.6, 0.0, 0.6, 0.9, 0.4, 0.8],
    'V': [0.9, 0.0, 0.0, 0.8, 0.2, 0.0]
}


def parse_pdbqt_atoms(path: Path):
    """Parse ATOM/HETATM lines from a PDBQT file."""
    atoms = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                resn = line[17:20].strip()
                chain = line[21].strip() or "A"
                resi = line[22:26].strip()
                atom_name = line[12:16].strip()
                element = (line[76:78].strip() if len(line) > 78 else "")
                if not element:
                    element = atom_name[0] if atom_name else "C"
                atoms.append({
                    "coord": np.array([x, y, z], dtype=np.float64),
                    "resn": resn,
                    "chain": chain,
                    "resi": resi,
                    "atom_name": atom_name,
                    "element": element,
                })
            except Exception:
                continue
    return atoms


def extract_pose_pocket_residues(ligand_atoms: list, receptor_atoms: list, cutoff: float = 5.0):
    """Find receptor residues within cutoff of any ligand atom."""
    if not ligand_atoms or not receptor_atoms:
        return []
    lig_coords = np.vstack([a["coord"] for a in ligand_atoms])
    res_map = {}
    for ra in receptor_atoms:
        key = (ra["chain"], ra["resi"], ra["resn"])
        res_map.setdefault(key, []).append(ra["coord"])
    pocket_residues = []
    for key, coords in res_map.items():
        center = np.mean(coords, axis=0)
        dists = np.linalg.norm(lig_coords - center, axis=1)
        if np.min(dists) <= cutoff:
            pocket_residues.append({"chain": key[0], "resi": key[1], "resn": key[2]})
    return pocket_residues


def pocket_geometry_embedding(residues: list, waters: list | None = None, target_dim: int = 64, water_weight: float = 3.0):
    """Weighted geometric embedding of pocket residues + waters (from S1)."""
    atom_coords = []
    atom_weights = []
    for r in residues:
        atom_coords.append(r.get("coord", np.zeros(3)))
        atom_weights.append(1.0)
    if waters:
        for w in waters:
            if isinstance(w, dict) and "coord" in w:
                atom_coords.append(w["coord"])
                atom_weights.append(water_weight)
            elif isinstance(w, (list, tuple)) and len(w) == 3:
                atom_coords.append(np.array(w))
                atom_weights.append(water_weight)
    if not atom_coords:
        return np.zeros(target_dim, dtype=np.float32)
    points = np.array(atom_coords, dtype=np.float64)
    weights = np.array(atom_weights, dtype=np.float64)
    weights = weights / (weights.sum() + 1e-12)
    centroid = np.average(points, axis=0, weights=weights)
    centered = points - centroid
    variance = np.average(centered ** 2, axis=0, weights=weights)
    variance = np.nan_to_num(variance, nan=0.0, posinf=1e5, neginf=-1e5)
    if len(points) >= 4:
        cov = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                cov[i, j] = np.average(centered[:, i] * centered[:, j], weights=weights)
        cov_upper = cov[np.triu_indices(3)]
    else:
        cov_upper = np.zeros(6)
    skew = np.average(centered ** 3, axis=0, weights=weights)
    rg = np.sqrt(np.average(np.sum(centered ** 2, axis=1), weights=weights))
    span = float(np.ptp(points, axis=0).max()) if len(points) > 1 else 0.0
    num_atoms = sum(1 for w in atom_weights if abs(w - 1.0) < 1e-6)
    num_waters = len(atom_weights) - num_atoms
    water_spread = 0.0
    if num_waters >= 2:
        water_points = points[np.array(atom_weights) == water_weight]
        water_centroid = water_points.mean(axis=0)
        water_spread = float(np.std(np.linalg.norm(water_points - water_centroid, axis=1)))
    geom_vec = np.concatenate([centroid, variance, cov_upper, skew, [rg, span], [num_atoms, num_waters, water_spread]])
    if len(geom_vec) < target_dim:
        padded = np.zeros(target_dim, dtype=np.float32)
        padded[:len(geom_vec)] = geom_vec
        geom_vec = padded
    else:
        geom_vec = geom_vec[:target_dim].astype(np.float32)
    geom_vec = np.nan_to_num(geom_vec, nan=0.0, posinf=1e6, neginf=-1e6)
    std = geom_vec.std()
    if std < 1e-8:
        return np.zeros(target_dim, dtype=np.float32)
    return (geom_vec - geom_vec.mean()) / (std + 1e-8)


def physicochemical_embedding(residues: list, target_dim: int = 24):
    """Physicochemical embedding of pocket residues (from S1)."""
    if not residues:
        return np.zeros(target_dim, dtype=np.float32)
    props = []
    for r in residues:
        aa = str(r.get("resn", "A")).upper()
        vec = AA_PHYSCHEM.get(aa, AA_PHYSCHEM['A'])
        props.append(vec)
    props = np.array(props, dtype=np.float32)
    mean = props.mean(axis=0)
    std = props.std(axis=0)
    minv = props.min(axis=0)
    maxv = props.max(axis=0)
    hydrophobic = sum(1 for r in residues if str(r.get("resn", "A")).upper() in "ILVFWM") / max(len(residues), 1)
    charged = sum(1 for r in residues if str(r.get("resn", "A")).upper() in "DEKRH") / max(len(residues), 1)
    phys_vec = np.concatenate([mean, std, minv, maxv, [hydrophobic, charged]])
    phys_vec = np.nan_to_num(phys_vec, nan=0.0)
    if phys_vec.std() > 1e-8:
        phys_vec = (phys_vec - phys_vec.mean()) / (phys_vec.std() + 1e-8)
    return phys_vec[:target_dim]


def cosine_umap(matrix: np.ndarray, n_components: int = 2, random_state: int = 13) -> np.ndarray:
    """UMAP on cosine-distance matrix instead of raw PCA."""
    if len(matrix) == 0:
        return np.empty((0, n_components))
    if len(matrix) == 1:
        return np.zeros((1, n_components))
    if not UMAP_AVAILABLE:
        return PCA(n_components=n_components, random_state=random_state).fit_transform(matrix)
    reducer = umap.UMAP(
        n_components=n_components,
        metric="cosine",
        random_state=random_state,
        n_neighbors=min(15, len(matrix) - 1),
        min_dist=0.1,
    )
    return reducer.fit_transform(matrix)


def build_pose_pocket_features(pose_path: str, receptor_path: str):
    """Extract pocket geometry + physicochemical embeddings from a docked pose."""
    pose_file = Path(pose_path).expanduser()
    receptor_file = Path(receptor_path).expanduser()
    if not pose_file.exists() or not receptor_file.exists():
        return {
            "pocket_geom_pc1": math.nan, "pocket_geom_pc2": math.nan,
            "pocket_phys_pc1": math.nan, "pocket_phys_pc2": math.nan,
        }
    ligand_atoms = parse_pdbqt_atoms(pose_file)
    receptor_atoms = parse_pdbqt_atoms(receptor_file)
    pocket_residues = extract_pose_pocket_residues(ligand_atoms, receptor_atoms, cutoff=5.0)
    if not pocket_residues:
        return {
            "pocket_geom_pc1": math.nan, "pocket_geom_pc2": math.nan,
            "pocket_phys_pc1": math.nan, "pocket_phys_pc2": math.nan,
        }
    pocket_atoms = []
    for r in pocket_residues:
        for ra in receptor_atoms:
            if ra["chain"] == r["chain"] and ra["resi"] == r["resi"] and ra["resn"] == r["resn"]:
                pocket_atoms.append({"coord": ra["coord"], "resn": r["resn"]})
                break
    geom_emb = pocket_geometry_embedding(pocket_atoms, waters=None, target_dim=64)
    phys_emb = physicochemical_embedding(pocket_atoms, target_dim=24)
    geom_pc = PCA(n_components=2, random_state=13).fit_transform(geom_emb.reshape(1, -1))[0]
    phys_pc = PCA(n_components=2, random_state=13).fit_transform(phys_emb.reshape(1, -1))[0]
    return {
        "pocket_geom_pc1": round(float(geom_pc[0]), 6),
        "pocket_geom_pc2": round(float(geom_pc[1]), 6),
        "pocket_phys_pc1": round(float(phys_pc[0]), 6),
        "pocket_phys_pc2": round(float(phys_pc[1]), 6),
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
    history_rows = read_csv(args.history_csv.expanduser())
    if args.analysis_csv:
        merge_analysis_scores(history_rows, read_csv(args.analysis_csv.expanduser()))
    reference_fps = load_reference_scaffold_fps(args.reference_smiles_csv, args.max_reference_rows)

    receptor_path = args.receptor
    if not receptor_path:
        config_path = args.history_csv.expanduser().parent / "history_config.json"
        if config_path.exists():
            try:
                config = json.loads(config_path.read_text(encoding="utf-8"))
                receptor_path = config.get("receptor", "")
            except Exception:
                pass

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
        for field in STRUCTURAL_FEATURES:
            if field in props or field == "chembl_scaffold_similarity":
                continue
            record[field] = safe_float(record.get(field))
        pose_path = row.get("pose_path", "")
        if receptor_path and pose_path:
            pocket_feats = build_pose_pocket_features(pose_path, receptor_path)
            record.update(pocket_feats)
        out.append(record)
    return out


def pca_coords(matrix: np.ndarray) -> np.ndarray:
    if len(matrix) == 0:
        return np.empty((0, 2))
    if len(matrix) == 1:
        return np.zeros((1, 2))
    return PCA(n_components=2, random_state=13).fit_transform(matrix)


def add_spaces(rows: list[dict[str, object]]) -> None:
    smiles_matrix = np.vstack([morgan_vector(str(row["canonical_smiles"])) for row in rows])
    smiles_xy = cosine_umap(smiles_matrix)
    structural_fields = [
        field for field in STRUCTURAL_FEATURES
        if any(not math.isnan(safe_float(row.get(field))) for row in rows)
    ]
    structural_matrix = np.array([[safe_float(row.get(field)) for field in structural_fields] for row in rows], dtype=float)
    structural_xy = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), PCA(n_components=2, random_state=13)).fit_transform(structural_matrix)
    for row, xy1, xy2 in zip(rows, smiles_xy, structural_xy):
        row["pure_smiles_pc1"] = round(float(xy1[0]), 6)
        row["pure_smiles_pc2"] = round(float(xy1[1]), 6)
        row["structural_interaction_pc1"] = round(float(xy2[0]), 6)
        row["structural_interaction_pc2"] = round(float(xy2[1]), 6)

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
            pure_dist = coordinate_distance(row_a, row_b, "pure_smiles_pc1", "pure_smiles_pc2")
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
            # S1 pocket embeddings
            "pocket_geom_pc1",
            "pocket_geom_pc2",
            "pocket_phys_pc1",
            "pocket_phys_pc2",
            "pure_smiles_pc1",
            "pure_smiles_pc2",
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
                "pure_smiles_pc1": row.get("pure_smiles_pc1", ""),
                "pure_smiles_pc2": row.get("pure_smiles_pc2", ""),
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
                "pure_smiles_pc1": row.get("pure_smiles_pc1", ""),
                "pure_smiles_pc2": row.get("pure_smiles_pc2", ""),
                "structural_interaction_pc1": row.get("structural_interaction_pc1", ""),
                "structural_interaction_pc2": row.get("structural_interaction_pc2", ""),
                "affinity_kcal_mol": row.get("affinity_kcal_mol", ""),
            }
        )

    center_distance = math.nan
    if high and low:
        high_center = np.array([
            float(np.mean([safe_float(row.get("pure_smiles_pc1")) for row in high if not math.isnan(safe_float(row.get("pure_smiles_pc1")))])) if any(not math.isnan(safe_float(row.get("pure_smiles_pc1"))) for row in high) else 0.0,
            float(np.mean([safe_float(row.get("pure_smiles_pc2")) for row in high if not math.isnan(safe_float(row.get("pure_smiles_pc2")))])) if any(not math.isnan(safe_float(row.get("pure_smiles_pc2"))) for row in high) else 0.0,
            float(np.mean([safe_float(row.get("structural_interaction_pc1")) for row in high if not math.isnan(safe_float(row.get("structural_interaction_pc1")))])) if any(not math.isnan(safe_float(row.get("structural_interaction_pc1"))) for row in high) else 0.0,
            float(np.mean([safe_float(row.get("structural_interaction_pc2")) for row in high if not math.isnan(safe_float(row.get("structural_interaction_pc2")))])) if any(not math.isnan(safe_float(row.get("structural_interaction_pc2"))) for row in high) else 0.0,
        ], dtype=float)
        low_center = np.array([
            float(np.mean([safe_float(row.get("pure_smiles_pc1")) for row in low if not math.isnan(safe_float(row.get("pure_smiles_pc1")))])) if any(not math.isnan(safe_float(row.get("pure_smiles_pc1"))) for row in low) else 0.0,
            float(np.mean([safe_float(row.get("pure_smiles_pc2")) for row in low if not math.isnan(safe_float(row.get("pure_smiles_pc2")))])) if any(not math.isnan(safe_float(row.get("pure_smiles_pc2"))) for row in low) else 0.0,
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
    parser.add_argument("--receptor", type=Path, default=None, help="Receptor PDBQT path (auto-detected from history_config.json if omitted)")
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
    parser.add_argument("--activity-cliff-max-pure-distance", type=float, default=0.45, help="Maximum Pure-SMILES PCA distance for visual-neighbor cliffs")
    parser.add_argument("--activity-cliff-max-structural-distance", type=float, default=0.50, help="Maximum structural-interaction PCA distance for interaction-neighbor cliffs")
    parser.add_argument("--activity-cliff-max-affinity-delta", type=float, default=1.0, help="Vina affinity delta treated as similar when filtering cliffs")
    parser.add_argument("--activity-cliff-require-vina-signal", action="store_true", help="Only report cliffs where Vina affinity is similar or points in the opposite direction")
    parser.add_argument("--high-score-threshold", type=float, default=0.4, help="High official binding score threshold for group comparison")
    parser.add_argument("--low-score-max", type=float, default=0.25, help="Low official binding score upper bound for group comparison")
    parser.add_argument("--seed", type=int, default=13)
    args = parser.parse_args()

    rows = build_rows(args)
    add_spaces(rows)
    metrics = train_predict(rows, args)
    outdir = args.outdir.expanduser().resolve()
    fields = [
        "seq_id", "nickname", "score_set", "official_binding_score", "predicted_binding_score",
        "prediction_residual",
        "canonical_smiles", "murcko_scaffold", "chembl_scaffold_similarity",
        "pure_smiles_pc1", "pure_smiles_pc2", "structural_interaction_pc1", "structural_interaction_pc2",
        *STRUCTURAL_FEATURES,
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
            "affinity_kcal_mol", "pure_smiles_pc1", "pure_smiles_pc2",
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
    render_space_html(rows, "pure_smiles_pc1", "pure_smiles_pc2", "Pure-SMILES-space (cosine-UMAP)", outdir / "pure_smiles_space.html", score_bins)
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
