#!/usr/bin/env python
"""Generate more drug-like SMILES candidates from expert/docking anchors."""

from __future__ import annotations

import argparse
import os
import csv
import importlib.util
import math
import random
import sys
from argparse import Namespace
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, QED, rdChemReactions, rdMolDescriptors

SMILES_SKILL = Path(__file__).resolve().parents[2] / "smiles-to-vina-docking" / "scripts"
sys.path.insert(0, str(SMILES_SKILL))
import dock_utils as du  # noqa: E402

DEFAULT_REFINEMENT_DIR = Path("~/vina_task2/druglike_refinement")
_WORKER_REFERENCE_FPS: list[object] = []
_WORKER_REFERENCE_PATTERN_FPS: list[object] = []
_WORKER_EXPERT_FPS: list[object] = []


def default_refinement_workers() -> int:
    logical = os.cpu_count() or 1
    if logical <= 2:
        return logical
    return max(1, min(8, logical // 2))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def clean_float(value: object, default: float = 0.0) -> float:
    try:
        if value in ("", None):
            return default
        return float(value)
    except Exception:
        return default


def mol(smiles: str | None) -> Chem.Mol | None:
    return du.mol_from_smiles(smiles)


def can(smiles: str | None) -> str:
    return du.canonical(smiles) or ""


@lru_cache(maxsize=200000)
def fingerprint(smiles: str) -> object | None:
    m = mol(smiles)
    if m is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)


@lru_cache(maxsize=200000)
def pattern_fingerprint(smiles: str) -> object | None:
    m = mol(smiles)
    if m is None:
        return None
    return Chem.PatternFingerprint(m)


def tanimoto(a: object | None, b: object | None) -> float:
    if a is None or b is None:
        return 0.0
    return float(DataStructs.TanimotoSimilarity(a, b))


@lru_cache(maxsize=200000)
def qed_score(smiles: str) -> float:
    m = mol(smiles)
    if m is None:
        return 0.0
    try:
        return float(QED.qed(m))
    except Exception:
        return 0.0


@lru_cache(maxsize=200000)
def mol_properties(smiles: str) -> dict[str, float | int]:
    m = mol(smiles)
    if m is None:
        return {}
    return {
        "mw": round(float(Descriptors.MolWt(m)), 3),
        "logp": round(float(Crippen.MolLogP(m)), 3),
        "hbd": int(Lipinski.NumHDonors(m)),
        "hba": int(Lipinski.NumHAcceptors(m)),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(m)), 3),
        "rot_bonds": int(Lipinski.NumRotatableBonds(m)),
        "heavy_atoms": int(m.GetNumHeavyAtoms()),
        "aromatic_rings": int(rdMolDescriptors.CalcNumAromaticRings(m)),
        "formal_charge": int(sum(atom.GetFormalCharge() for atom in m.GetAtoms())),
    }


def window_score(props: dict[str, float | int]) -> float:
    if not props:
        return 0.0
    score = 1.0
    mw = float(props.get("mw", 999))
    logp = float(props.get("logp", 99))
    tpsa = float(props.get("tpsa", 999))
    rot = float(props.get("rot_bonds", 99))
    hbd = float(props.get("hbd", 99))
    hba = float(props.get("hba", 99))
    if mw > 900:
        score -= min(0.06, (mw - 900) / 1500)
    elif mw > 500:
        score -= min(0.03, (mw - 500) / 1000)
    if mw < 250:
        score -= min(0.15, (250 - mw) / 250)
    if logp > 5:
        score -= min(0.25, (logp - 5) / 5)
    if logp < -1:
        score -= min(0.15, (-1 - logp) / 4)
    if tpsa > 140:
        score -= min(0.25, (tpsa - 140) / 200)
    if rot > 10:
        score -= min(0.20, (rot - 10) / 20)
    if hbd > 5:
        score -= min(0.15, (hbd - 5) / 10)
    if hba > 10:
        score -= min(0.15, (hba - 10) / 15)
    return round(max(0.0, min(1.0, score)), 4)


def clamp01(value: float) -> float:
    return max(0.0, min(1.0, value))


def load_pose_analyzer_module():
    pose_script = Path(__file__).resolve().parents[2] / "pose-analyzer" / "scripts" / "infer_binding_mode_families.py"
    if not pose_script.exists():
        return None
    spec = importlib.util.spec_from_file_location("pose_analyzer_infer", pose_script)
    if spec is None or spec.loader is None:
        return None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@lru_cache(maxsize=200000)
def synthetic_score_proxy(smiles: str) -> float:
    """A lightweight 0-1 proxy; higher means easier-looking, not true retrosynthesis."""
    m = mol(smiles)
    if m is None:
        return 0.0
    props = mol_properties(smiles)
    heavy = float(props.get("heavy_atoms", 99))
    rot = float(props.get("rot_bonds", 99))
    rings = rdMolDescriptors.CalcNumRings(m)
    aromatic = float(props.get("aromatic_rings", 0))
    stereo = len(Chem.FindMolChiralCenters(m, includeUnassigned=True))
    bertz = float(Descriptors.BertzCT(m))
    macrocycles = sum(1 for ring in m.GetRingInfo().AtomRings() if len(ring) >= 8)
    risk = 0.0
    risk += max(0.0, heavy - 55) / 180
    risk += max(0.0, rot - 8) / 20
    risk += max(0.0, rings - 5) / 10
    risk += max(0.0, aromatic - 4) / 8
    risk += stereo / 12
    risk += max(0.0, bertz - 800) / 1800
    risk += macrocycles * 0.2
    return round(max(0.0, min(1.0, 1.0 - risk)), 4)


def structural_alert_penalty(smiles: str, edit_label: str = "") -> float:
    m = mol(smiles)
    if m is None:
        return 1.0
    alerts = [
        ("quinone", "[#6]1(=O)[#6]=[#6][#6](=O)[#6]=[#6]1", 0.25),
        ("catechol", "c1cc(O)c(O)cc1", 0.15),
        ("michael_acceptor", "C=CC=O", 0.15),
        ("thiol", "[SH]", 0.12),
    ]
    penalty = 0.0
    for _, smarts, value in alerts:
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None and m.HasSubstructMatch(patt):
            penalty += value
    if "quinone_alert" in edit_label:
        penalty += 0.20
    return round(min(0.6, penalty), 4)


def row_contact_score(row: dict[str, str]) -> float:
    hbond = clean_float(row.get("hbond_count"))
    hydrophobic = clean_float(row.get("hydrophobic_count") or row.get("hydrophobic_count ç–æ°´"))
    vdw = clean_float(row.get("vdw_contact_count"))
    pi = clean_float(row.get("pi_contact_count"))
    ch_pi = clean_float(row.get("ch_pi_count"))
    return 3.0 * hbond + 0.7 * hydrophobic + 0.05 * vdw + 2.0 * pi + 2.0 * ch_pi


def structural_footprint_score(row: dict[str, str]) -> float:
    raw_contact = row_contact_score(row)
    contact_score = clamp01(raw_contact / 70.0)
    residue_score = clamp01(clean_float(row.get("contact_residue_count_4a")) / 18.0)
    receptor_cov = clamp01(clean_float(row.get("receptor_residue_coverage_4a")) / 0.12)
    surface_fraction = clamp01(clean_float(row.get("surface_contact_fraction_4a")))
    ligand_fraction = clamp01(clean_float(row.get("ligand_atom_contact_fraction_4a")))
    return round(
        0.30 * contact_score
        + 0.25 * residue_score
        + 0.20 * receptor_cov
        + 0.15 * surface_fraction
        + 0.10 * ligand_fraction,
        4,
    )


def parent_context_from_row(row: dict[str, str]) -> dict[str, object]:
    smiles = can(row.get("canonical_smiles") or row.get("smiles") or row.get("input_smiles"))
    props = mol_properties(smiles)
    q = qed_score(smiles) if smiles else 0.0
    synth = synthetic_score_proxy(smiles) if smiles else 0.0
    affinity = clean_float(row.get("affinity_kcal_mol"), default=0.0)
    contact = row_contact_score(row)
    heavy = max(1.0, float(props.get("heavy_atoms", 1) or 1))
    efficiency = contact / heavy
    footprint = structural_footprint_score(row)
    prop_window = window_score(props)
    anchor_score = round(
        0.24 * q
        + 0.18 * synth
        + 0.14 * prop_window
        + 0.26 * footprint
        + 0.08 * efficiency
        + 0.10 * max(0, -affinity / 12),
        4,
    )
    return {
        "seq_id": row.get("seq_id", ""),
        "nickname": row.get("nickname", ""),
        "ancestor_smiles": row.get("ancestor_smiles", "") or smiles,
        "canonical_smiles": smiles,
        "anchor_score": anchor_score,
        "structural_footprint_score": footprint,
    }


def best_history_row_by_smiles(history_rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    best: dict[str, dict[str, str]] = {}
    for row in history_rows:
        canonical = can(row.get("canonical_smiles") or row.get("smiles") or row.get("input_smiles"))
        if not canonical:
            continue
        current = best.get(canonical)
        if current is None:
            best[canonical] = row
            continue
        row_key = (
            structural_footprint_score(row),
            row_contact_score(row),
            clean_float(row.get("affinity_kcal_mol"), default=-999.0),
        )
        current_key = (
            structural_footprint_score(current),
            row_contact_score(current),
            clean_float(current.get("affinity_kcal_mol"), default=-999.0),
        )
        if row_key > current_key:
            best[canonical] = row
    return best


def enrich_history_with_pose_context(history_rows: list[dict[str, str]], receptor_path: Path | None) -> None:
    if receptor_path is None:
        return
    receptor_path = receptor_path.expanduser()
    if not receptor_path.exists():
        return
    pose_module = load_pose_analyzer_module()
    if pose_module is None:
        return
    receptor_atoms = pose_module.read_receptor_atoms(receptor_path)
    receptor_residues = pose_module.receptor_residue_index(receptor_atoms)
    receptor_residue_count = len(receptor_residues)
    try:
        _, surface_labels, _ = pose_module.pymol_surface_residue_table(receptor_path, receptor_atoms, 0.30)
    except Exception:
        _, surface_labels, _ = pose_module.surface_residue_table(receptor_atoms)
    cache: dict[str, dict[str, object]] = {}
    for row in history_rows:
        pose_path_text = row.get("pose_path") or row.get("pdbqt_path")
        if not pose_path_text:
            continue
        pose_path = Path(pose_path_text).expanduser()
        if not pose_path.exists():
            continue
        cache_key = str(pose_path)
        if cache_key not in cache:
            modes = pose_module.parse_pdbqt_modes(pose_path)
            if not modes:
                cache[cache_key] = {}
            else:
                best_mode = min(
                    modes,
                    key=lambda mode: clean_float(mode.get("affinity"), default=999.0),
                )
                cache[cache_key] = pose_module.contact_features(
                    best_mode["atoms"],
                    receptor_atoms,
                    cutoff=4.0,
                    surface_labels=surface_labels,
                    binding_surface_labels=surface_labels,
                    receptor_residue_count=receptor_residue_count,
                )
        features = cache[cache_key]
        if not features:
            continue
        for key in (
            "contact_residue_count_4a",
            "contact_residue_count_5a",
            "contact_residue_count_6a",
            "surface_contact_residue_count_4a",
            "surface_contact_residue_count_5a",
            "surface_contact_fraction_4a",
            "binding_surface_contact_count_4a",
            "binding_surface_contact_count_5a",
            "binding_surface_coverage_4a",
            "binding_surface_coverage_5a",
            "binding_surface_total_residues",
            "receptor_residue_coverage_4a",
            "receptor_residue_coverage_5a",
            "ligand_atom_contact_fraction_4a",
            "min_ligand_receptor_distance",
            "closest_residues",
        ):
            value = features.get(key, "")
            row[key] = str(value) if value != "" else ""


def anchor_rows(history_rows: list[dict[str, str]], limit: int) -> list[dict[str, object]]:
    anchors: list[dict[str, object]] = []
    best_rows = best_history_row_by_smiles(history_rows)
    for row in best_rows.values():
        smiles = row.get("canonical_smiles") or row.get("smiles") or row.get("input_smiles")
        canonical = can(smiles)
        if not canonical:
            continue
        props = mol_properties(canonical)
        heavy = max(1.0, float(props.get("heavy_atoms", 1)))
        contact = row_contact_score(row)
        affinity = clean_float(row.get("affinity_kcal_mol"), default=0.0)
        q = qed_score(canonical)
        synth = synthetic_score_proxy(canonical)
        efficiency = contact / heavy
        footprint = structural_footprint_score(row)
        prop_window = window_score(props)
        anchors.append(
            {
                "seq_id": row.get("seq_id", ""),
                "nickname": row.get("nickname", ""),
                "ancestor_smiles": row.get("ancestor_smiles", "") or canonical,
                "canonical_smiles": canonical,
                "affinity_kcal_mol": affinity,
                "contact_score": round(contact, 4),
                "contact_efficiency": round(efficiency, 4),
                "structural_footprint_score": footprint,
                "contact_residue_count_4a": row.get("contact_residue_count_4a", ""),
                "surface_contact_fraction_4a": row.get("surface_contact_fraction_4a", ""),
                "receptor_residue_coverage_4a": row.get("receptor_residue_coverage_4a", ""),
                "ligand_atom_contact_fraction_4a": row.get("ligand_atom_contact_fraction_4a", ""),
                "qed": round(q, 4),
                "synthetic_score_proxy": synth,
                "anchor_score": round(
                    0.24 * q
                    + 0.18 * synth
                    + 0.14 * prop_window
                    + 0.26 * footprint
                    + 0.08 * efficiency
                    + 0.10 * max(0, -affinity / 12),
                    4,
                ),
                **props,
            }
        )
    anchors.sort(key=lambda row: float(row["anchor_score"]), reverse=True)
    return anchors[:limit]


def load_reference_fps(path: Path | None, max_rows: int) -> list[object]:
    if path is None or not path.expanduser().exists():
        return []
    rows = du.read_csv(path.expanduser())
    fps: list[object] = []
    for row in rows[:max_rows]:
        smiles = row.get("smiles") or row.get("canonical_smiles") or row.get("canonical_smiles ") or row.get("SMILES")
        fp = fingerprint(can(smiles))
        if fp is not None:
            fps.append(fp)
    return fps


def load_reference_pattern_fps(path: Path | None, max_rows: int) -> list[object]:
    if path is None or not path.expanduser().exists():
        return []
    rows = du.read_csv(path.expanduser())
    fps: list[object] = []
    for row in rows[:max_rows]:
        smiles = row.get("smiles") or row.get("canonical_smiles") or row.get("canonical_smiles ") or row.get("SMILES")
        fp = pattern_fingerprint(can(smiles))
        if fp is not None:
            fps.append(fp)
    return fps


def max_reference_similarity(smiles: str, reference_fps: list[object]) -> float:
    fp = fingerprint(smiles)
    if fp is None or not reference_fps:
        return 0.0
    return round(max(DataStructs.BulkTanimotoSimilarity(fp, reference_fps)), 4)


def max_reference_pattern_similarity(smiles: str, reference_pattern_fps: list[object]) -> float:
    fp = pattern_fingerprint(smiles)
    if fp is None or not reference_pattern_fps:
        return 0.0
    return round(max(DataStructs.BulkTanimotoSimilarity(fp, reference_pattern_fps)), 4)


REACTION_SMARTS: list[tuple[str, str]] = [
    ("amide_to_ketone", "[C:1](=[O:2])[N:3]>>[C:1](=[O:2])[C:3]"),
    ("amide_to_amine_linker", "[C:1](=[O:2])[N:3]>>[C:1][N:3]"),
    ("amide_to_ether_linker", "[C:1](=[O:2])[N:3]>>[C:1][O:3]"),
    ("amide_to_thioether_linker", "[C:1](=[O:2])[N:3]>>[C:1][S:3]"),
    ("methoxy_to_fluoro", "[c:1][O:2][CH3:3]>>[c:1][F:2]"),
    ("methoxy_to_hydroxy", "[c:1][O:2][CH3:3]>>[c:1][O:2]"),
    ("ethyl_to_methyl", "[c:1][CH2:2][CH3:3]>>[c:1][CH3:2]"),
    ("propyl_alcohol_to_methoxy", "[c:1][CH2:2][CH2:3][CH2:4][O:5]>>[c:1][O:2][CH3:3]"),
    ("propyl_alcohol_to_hydroxy", "[c:1][CH2:2][CH2:3][CH2:4][O:5]>>[c:1][O:2]"),
    ("phenyl_CH_to_pyridyl_N", "[cH:1]1[c:2][c:3][c:4][c:5][c:6]1>>[n:1]1[c:2][c:3][c:4][c:5][c:6]1"),
    ("phenyl_diCH_to_pyrimidyl_diN", "[cH:1]1[c:2][cH:3][c:4][c:5][c:6]1>>[n:1]1[c:2][n:3][c:4][c:5][c:6]1"),
    ("phenol_to_benzoquinone_alert", "[c:1]1[cH:2][cH:3][c:4]([O:7])[cH:5][cH:6]1>>[O:7]=[C:4]1[C:3]=[C:2][C:1]=[C:6][C:5]=1"),
]


def reaction_edits(smiles: str, max_products: int = 80) -> list[tuple[str, str]]:
    reactant = mol(smiles)
    if reactant is None:
        return []
    out: list[tuple[str, str]] = []
    seen: set[str] = set()
    for label, smarts in REACTION_SMARTS:
        try:
            reaction = rdChemReactions.ReactionFromSmarts(smarts)
            product_sets = reaction.RunReactants((reactant,))
        except Exception:
            continue
        for products in product_sets:
            if not products:
                continue
            product = products[0]
            try:
                Chem.SanitizeMol(product)
            except Exception:
                continue
            canonical = Chem.MolToSmiles(product, canonical=True)
            if canonical and canonical != can(smiles):
                if canonical in seen:
                    continue
                seen.add(canonical)
                out.append((canonical, f"bioisostere_{label}"))
                if len(out) >= max_products:
                    return out
    return out


FIVE_MEMBER_HETEROARYL_TEMPLATES: list[tuple[str, tuple[str, ...], int]] = [
    ("phenyl_to_furyl_manual", ("O", "C", "C", "C", "C"), 1),
    ("phenyl_to_thienyl_manual", ("S", "C", "C", "C", "C"), 1),
    ("phenyl_to_oxazolyl_manual", ("O", "C", "N", "C", "C"), 1),
    ("phenyl_to_thiazolyl_manual", ("S", "C", "N", "C", "C"), 1),
]


def five_member_heteroaryl_replacements(smiles: str, max_products: int = 40) -> list[tuple[str, str]]:
    """Conservative 6->5 aromatic replacement for mono-substituted phenyl rings.

    General 6-to-5 ring contraction is topology-changing, so reaction SMARTS are
    unsafe here. This manual version only handles phenyl rings with one external
    attachment and drops an unsubstituted ring carbon.
    """
    base = mol(smiles)
    if base is None:
        return []
    out: list[tuple[str, str]] = []
    ring_info = base.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        ring_set = set(ring)
        if not all(base.GetAtomWithIdx(idx).GetIsAromatic() and base.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue
        external: list[tuple[int, int, Chem.BondType]] = []
        for idx in ring:
            atom = base.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in ring_set:
                    continue
                bond = base.GetBondBetweenAtoms(idx, nbr_idx)
                if bond is not None:
                    external.append((idx, nbr_idx, bond.GetBondType()))
        if len(external) != 1:
            continue
        attached_ring_idx, external_idx, external_bond_type = external[0]
        ordered = list(ring)
        if attached_ring_idx not in ordered:
            continue
        attach_pos = ordered.index(attached_ring_idx)
        for drop_offset in (2, 3, 4):
            drop_pos = (attach_pos + drop_offset) % 6
            drop_idx = ordered[drop_pos]
            drop_atom = base.GetAtomWithIdx(drop_idx)
            if drop_atom.GetTotalNumHs() < 1:
                continue
            kept = [idx for i, idx in enumerate(ordered) if i != drop_pos]
            mapped_attach_pos = kept.index(attached_ring_idx)
            for label, atom_symbols, template_attach_pos in FIVE_MEMBER_HETEROARYL_TEMPLATES:
                rw = Chem.RWMol(base)
                new_indices: list[int] = []
                for symbol in atom_symbols:
                    atom = Chem.Atom(symbol)
                    atom.SetIsAromatic(True)
                    new_indices.append(rw.AddAtom(atom))
                for i in range(5):
                    rw.AddBond(new_indices[i], new_indices[(i + 1) % 5], Chem.BondType.AROMATIC)
                # Rotate template so the original substituent lands on a carbon attachment site.
                attach_new_idx = new_indices[template_attach_pos]
                if atom_symbols[template_attach_pos] != "C":
                    attach_new_idx = new_indices[1]
                try:
                    rw.AddBond(external_idx, attach_new_idx, external_bond_type)
                    for old_idx in sorted(ring_set, reverse=True):
                        rw.RemoveAtom(old_idx)
                    product = rw.GetMol()
                    Chem.SanitizeMol(product)
                    canonical = Chem.MolToSmiles(product, canonical=True)
                except Exception:
                    continue
                if canonical and canonical != can(smiles):
                    out.append((canonical, f"bioisostere_{label}"))
                    if len(out) >= max_products:
                        return out
    dedup: dict[str, str] = {}
    for product, label in out:
        dedup.setdefault(product, label)
    return [(product, label) for product, label in dedup.items()]


def ring_system_size(atom: Chem.Atom) -> int:
    info = atom.GetOwningMol().GetRingInfo()
    rings = [set(ring) for ring in info.AtomRings() if atom.GetIdx() in ring]
    if not rings:
        return 0
    merged = set().union(*rings)
    changed = True
    while changed:
        changed = False
        for ring in rings:
            if merged & ring and not ring <= merged:
                merged |= ring
                changed = True
    return len(merged)


def intramolecular_cyclizations(smiles: str, max_products: int = 50) -> list[tuple[str, str]]:
    base = mol(smiles)
    if base is None:
        return []
    candidates: list[tuple[int, int, int]] = []
    for a in base.GetAtoms():
        if a.GetAtomicNum() not in {6, 7, 8, 16}:
            continue
        if a.GetTotalNumHs() < 1:
            continue
        for b in base.GetAtoms():
            if b.GetIdx() <= a.GetIdx():
                continue
            if b.GetAtomicNum() not in {6, 7, 8, 16} or b.GetTotalNumHs() < 1:
                continue
            if a.IsInRing() or b.IsInRing():
                continue
            if base.GetBondBetweenAtoms(a.GetIdx(), b.GetIdx()) is not None:
                continue
            path = Chem.rdmolops.GetShortestPath(base, a.GetIdx(), b.GetIdx())
            if not 4 <= len(path) <= 6:
                continue
            if any(base.GetAtomWithIdx(idx).IsInRing() for idx in path[1:-1]):
                continue
            candidates.append((a.GetIdx(), b.GetIdx(), len(path)))
    candidates.sort(key=lambda item: item[2])
    out: list[tuple[str, str]] = []
    for a_idx, b_idx, path_len in candidates[: max_products * 2]:
        rw = Chem.RWMol(base)
        try:
            rw.AddBond(a_idx, b_idx, Chem.BondType.SINGLE)
            product = rw.GetMol()
            Chem.SanitizeMol(product)
            canonical = Chem.MolToSmiles(product, canonical=True)
        except Exception:
            continue
        if canonical and canonical != can(smiles):
            out.append((canonical, f"cyclize_path{path_len}_A{a_idx}_B{b_idx}"))
            if len(out) >= max_products:
                break
    return out


def qve_delta(parent_smiles: str, candidate_smiles: str) -> float:
    parent = mol_properties(parent_smiles)
    child = mol_properties(candidate_smiles)
    if not parent or not child:
        return -999.0
    return round(
        0.55 * (qed_score(candidate_smiles) - qed_score(parent_smiles))
        + 0.20 * (window_score(child) - window_score(parent))
        + 0.15 * (synthetic_score_proxy(candidate_smiles) - synthetic_score_proxy(parent_smiles))
        + 0.10 * max(0.0, float(parent.get("rot_bonds", 0)) - float(child.get("rot_bonds", 0))) / 10.0,
        4,
    )


def refinement_edits(smiles: str, allow_add: bool, scaffold_edits: bool, require_qve_gain: bool) -> list[tuple[str, str]]:
    modes = ("delete", "shrink", "drastic", "add") if allow_add else ("delete", "shrink", "drastic")
    edits = list(du.enumerate_edits(smiles, modes))
    if scaffold_edits:
        edits.extend(reaction_edits(smiles))
        edits.extend(five_member_heteroaryl_replacements(smiles))
        edits.extend(intramolecular_cyclizations(smiles))
    dedup: dict[str, str] = {}
    for product, label in edits:
        if not allow_add and label.startswith("add_"):
            continue
        canonical = can(product)
        if not canonical or canonical == can(smiles):
            continue
        if "." in canonical:
            continue
        if require_qve_gain and qve_delta(smiles, canonical) <= 0:
            continue
        dedup.setdefault(canonical, label)
    return [(product, label) for product, label in dedup.items()]


def generate_candidates(
    parents: list[dict[str, object]],
    max_rounds: int,
    batch_size: int,
    seed: int,
    allow_add: bool,
    scaffold_edits: bool,
    require_qve_gain: bool,
    score_fn=None,
    target_score: float | None = None,
    max_total_candidates: int | None = None,
    progress_interval: int = 25,
) -> list[dict[str, object]]:
    rng = random.Random(seed)
    seen = {str(parent["canonical_smiles"]) for parent in parents}
    frontier = list(parents)
    generated: list[dict[str, object]] = []
    best_score = float("-inf")
    stop_reason = "max_rounds"
    for round_index in range(1, max_rounds + 1):
        next_frontier: list[dict[str, object]] = []
        round_count = 0
        stop_now = False
        rng.shuffle(frontier)
        print(
            f"[refine] round={round_index}/{max_rounds} parents={len(frontier)} "
            f"batch_size={batch_size} generated_total={len(generated)} best_score={best_score if best_score > -999 else ''}",
            flush=True,
        )
        for parent in frontier:
            parent_smiles = str(parent["canonical_smiles"])
            edits = refinement_edits(parent_smiles, allow_add, scaffold_edits, require_qve_gain)
            rng.shuffle(edits)
            for product, label in edits:
                canonical = can(product)
                if not canonical or canonical in seen:
                    continue
                seen.add(canonical)
                row = {
                    "source_parent": parent.get("seq_id", "") or parent.get("source_parent", "expert"),
                    "ancestor_smiles": parent.get("ancestor_smiles", "") or parent_smiles,
                    "parent_smiles": parent_smiles,
                    "parent_anchor_score": parent.get("anchor_score", ""),
                    "parent_structural_footprint_score": parent.get("structural_footprint_score", ""),
                    "candidate_smiles": canonical,
                    "edit_label": label,
                    "qve_delta": qve_delta(parent_smiles, canonical),
                    "generation": round_index,
                }
                generated.append(row)
                row_score = float("-inf")
                if score_fn is not None:
                    row_score = clean_float(score_fn(row).get("druglike_refinement_score"), default=float("-inf"))
                    best_score = max(best_score, row_score)
                next_frontier.append(
                    {
                        "canonical_smiles": canonical,
                        "seq_id": row["source_parent"],
                        "ancestor_smiles": row["ancestor_smiles"],
                        "anchor_score": row["parent_anchor_score"],
                        "structural_footprint_score": row["parent_structural_footprint_score"],
                        "_score": row_score,
                    }
                )
                round_count += 1
                if progress_interval > 0 and round_count % progress_interval == 0:
                    print(
                        f"[refine] round={round_index} {round_count}/{batch_size} "
                        f"generated_total={len(generated)} best_score={best_score:.4f}",
                        flush=True,
                    )
                if max_total_candidates is not None and len(generated) >= max_total_candidates:
                    stop_reason = "max_total_candidates"
                    print(f"[refine] stop reason={stop_reason} generated_total={len(generated)} best_score={best_score:.4f}", flush=True)
                    return generated
                if target_score is not None and best_score >= target_score:
                    stop_reason = "target_score"
                    print(f"[refine] stop reason={stop_reason} best_score={best_score:.4f}", flush=True)
                    return generated
                if round_count >= batch_size:
                    stop_now = True
                    break
            if stop_now:
                break
        if score_fn is not None:
            next_frontier.sort(key=lambda parent: clean_float(parent.get("_score"), default=float("-inf")), reverse=True)
        print(
            f"[refine] round={round_index} done generated_round={round_count} "
            f"next_parents={len(next_frontier)} best_score={best_score:.4f}",
            flush=True,
        )
        frontier = next_frontier[:batch_size]
        if not frontier:
            stop_reason = "no_next_parents"
            print(f"[refine] stop reason={stop_reason} generated_total={len(generated)} best_score={best_score:.4f}", flush=True)
            break
    else:
        print(f"[refine] stop reason={stop_reason} generated_total={len(generated)} best_score={best_score:.4f}", flush=True)
    return generated


def score_candidate(
    row: dict[str, object],
    reference_fps: list[object],
    reference_pattern_fps: list[object],
    expert_fps: list[object],
) -> dict[str, object]:
    smiles = str(row["candidate_smiles"])
    props = mol_properties(smiles)
    q = qed_score(smiles)
    synth = synthetic_score_proxy(smiles)
    prop_score = window_score(props)
    ref_sim = max_reference_similarity(smiles, reference_fps)
    ref_partial_sim = max_reference_pattern_similarity(smiles, reference_pattern_fps)
    fp = fingerprint(smiles)
    expert_sim = round(max([tanimoto(fp, expert_fp) for expert_fp in expert_fps] or [0.0]), 4)
    parent_sim = tanimoto(fp, fingerprint(str(row.get("parent_smiles", ""))))
    alert_penalty = structural_alert_penalty(smiles, str(row.get("edit_label", "")))
    parent_anchor_score = clean_float(row.get("parent_anchor_score"))
    parent_footprint = clean_float(row.get("parent_structural_footprint_score"))
    inherited_structure = parent_sim * max(parent_anchor_score, parent_footprint)
    qve_component = clamp01(max(0.0, float(row.get("qve_delta", 0.0))) * 8.0)
    components = {
        "qed_component_score": 10.0 * q,
        "synthetic_component_score": 10.0 * synth,
        "property_component_score": 10.0 * prop_score,
        "parent_similarity_score": 10.0 * parent_sim,
        "expert_similarity_score": 10.0 * expert_sim,
        "reference_similarity_score": 10.0 * ref_sim,
        "reference_partial_similarity_score": 10.0 * ref_partial_sim,
        "qve_component_score": 10.0 * qve_component,
        "inherited_structure_score": 10.0 * inherited_structure,
        "structural_alert_penalty_score": 10.0 * alert_penalty,
    }
    score = (
        0.16 * components["qed_component_score"]
        + 0.14 * components["synthetic_component_score"]
        + 0.10 * components["property_component_score"]
        + 0.12 * components["parent_similarity_score"]
        + 0.10 * components["expert_similarity_score"]
        + 0.12 * components["reference_similarity_score"]
        + 0.10 * components["reference_partial_similarity_score"]
        + 0.08 * components["qve_component_score"]
        + 0.12 * components["inherited_structure_score"]
        - 0.18 * components["structural_alert_penalty_score"]
    )
    return {
        **row,
        "druglike_refinement_score": round(score, 4),
        **{key: round(value, 4) for key, value in components.items()},
        "qed": round(q, 4),
        "synthetic_score_proxy": synth,
        "property_window_score": prop_score,
        "parent_similarity": round(parent_sim, 4),
        "expert_similarity": expert_sim,
        "kinase_reference_similarity": ref_sim,
        "kinase_reference_partial_similarity": ref_partial_sim,
        "parent_anchor_score": round(parent_anchor_score, 4),
        "parent_structural_footprint_score": round(parent_footprint, 4),
        "inherited_structure_raw": round(inherited_structure, 4),
        "qve_delta": row.get("qve_delta", ""),
        "structural_alert_penalty": alert_penalty,
        **props,
    }


def init_score_worker(reference_fps: list[object], reference_pattern_fps: list[object], expert_fps: list[object]) -> None:
    global _WORKER_REFERENCE_FPS, _WORKER_REFERENCE_PATTERN_FPS, _WORKER_EXPERT_FPS
    _WORKER_REFERENCE_FPS = reference_fps
    _WORKER_REFERENCE_PATTERN_FPS = reference_pattern_fps
    _WORKER_EXPERT_FPS = expert_fps


def score_candidate_worker(row: dict[str, object]) -> dict[str, object]:
    return score_candidate(row, _WORKER_REFERENCE_FPS, _WORKER_REFERENCE_PATTERN_FPS, _WORKER_EXPERT_FPS)


def score_rows(
    rows: list[dict[str, object]],
    reference_fps: list[object],
    reference_pattern_fps: list[object],
    expert_fps: list[object],
    workers: int,
) -> list[dict[str, object]]:
    if workers <= 1 or len(rows) < 2:
        return [score_candidate(row, reference_fps, reference_pattern_fps, expert_fps) for row in rows]
    worker_count = max(1, min(workers, len(rows)))
    chunksize = max(1, len(rows) // (worker_count * 4))
    try:
        with ProcessPoolExecutor(
            max_workers=worker_count,
            initializer=init_score_worker,
            initargs=(reference_fps, reference_pattern_fps, expert_fps),
        ) as pool:
            return list(pool.map(score_candidate_worker, rows, chunksize=chunksize))
    except Exception as exc:
        print(f"[refine] parallel scoring failed; falling back to serial scoring: {exc}", flush=True)
        return [score_candidate(row, reference_fps, reference_pattern_fps, expert_fps) for row in rows]


def refinement_meta_for_row(row: dict[str, object]) -> dict[str, str]:
    return {
        "refinement_source": "druglike-pocket-refiner",
        "druglike_refinement_score": str(row.get("druglike_refinement_score", "")),
        "refinement_qed": str(row.get("qed", "")),
        "synthetic_score_proxy": str(row.get("synthetic_score_proxy", "")),
        "property_window_score": str(row.get("property_window_score", "")),
        "qed_component_score": str(row.get("qed_component_score", "")),
        "synthetic_component_score": str(row.get("synthetic_component_score", "")),
        "property_component_score": str(row.get("property_component_score", "")),
        "parent_similarity_score": str(row.get("parent_similarity_score", "")),
        "expert_similarity_score": str(row.get("expert_similarity_score", "")),
        "reference_similarity_score": str(row.get("reference_similarity_score", "")),
        "reference_partial_similarity_score": str(row.get("reference_partial_similarity_score", "")),
        "qve_component_score": str(row.get("qve_component_score", "")),
        "inherited_structure_score": str(row.get("inherited_structure_score", "")),
        "structural_alert_penalty_score": str(row.get("structural_alert_penalty_score", "")),
        "refinement_parent_similarity": str(row.get("parent_similarity", "")),
        "refinement_expert_similarity": str(row.get("expert_similarity", "")),
        "kinase_reference_similarity": str(row.get("kinase_reference_similarity", "")),
        "kinase_reference_partial_similarity": str(row.get("kinase_reference_partial_similarity", "")),
        "qve_delta": str(row.get("qve_delta", "")),
        "structural_alert_penalty": str(row.get("structural_alert_penalty", "")),
        "refinement_generation": str(row.get("generation", "")),
    }


def dock_refined_worker(payload: dict[str, object]) -> dict[str, str | int | float]:
    row = du.dock_one(
        seq=str(payload["seq"]),
        smiles=str(payload["smiles"]),
        input_smiles=str(payload["input_smiles"]),
        nickname=str(payload["nickname"]),
        ancestor_smiles=str(payload["ancestor_smiles"]),
        parent_smiles=str(payload["parent_smiles"]),
        edit_label=str(payload["edit_label"]),
        ledger_dir=Path(str(payload["ledger_dir"])),
        receptor=Path(str(payload["receptor"])),
        receptor_atoms=payload["receptor_atoms"],  # type: ignore[arg-type]
        box=payload["box"],  # type: ignore[arg-type]
        vina=str(payload["vina"]),
        meeko=str(payload["meeko"]),
        obabel=str(payload["obabel"]),
        seed=int(payload["seed"]),
        cluster_cutoff=float(payload["cluster_cutoff"]),
    )
    du.merge_refinement_meta(row, payload["refinement_meta"])  # type: ignore[arg-type]
    return row


def dock_refined_candidates(rows: list[dict[str, object]], args: argparse.Namespace) -> None:
    if args.dock_top_candidates <= 0:
        return
    if args.receptor is None and not (args.history_dir / du.CONFIG_FILENAME).exists():
        raise ValueError("--dock-top-candidates needs --receptor/--config once, or an existing history_config.json")
    dock_args = Namespace(**vars(args), ledger_dir=args.history_dir)
    du.save_config(args.history_dir, dock_args)
    receptor_arg = du.arg_or_config(dock_args, "receptor")
    config_arg = du.arg_or_config(dock_args, "config")
    if not receptor_arg or not config_arg:
        raise ValueError("--dock-top-candidates needs --receptor and --config once")
    receptor = Path(str(receptor_arg)).expanduser().resolve()
    config = Path(str(config_arg)).expanduser().resolve()
    vina = du.resolve_executable(str(du.arg_or_config(dock_args, "vina", "vina")), args.history_dir, "vina")
    meeko = du.resolve_executable(str(du.arg_or_config(dock_args, "meeko", "mk_prepare_ligand.py")), args.history_dir, "mk_prepare_ligand.py")
    obabel = du.resolve_executable(str(du.arg_or_config(dock_args, "obabel", "obabel")), args.history_dir, "obabel")
    box = du.parse_vina_config(config, dock_args)
    receptor_atoms = du.parse_receptor_atoms(receptor)
    ledger_file = du.ledger_path(args.history_dir)
    existing_rows = du.read_csv(ledger_file)
    seq_number = du.next_seq(existing_rows)
    seen = [du.row_identity(row) for row in existing_rows]
    jobs: list[dict[str, object]] = []
    new_rows: list[dict[str, str | int | float]] = []
    updated_existing = 0
    for index, row in enumerate(rows[: args.dock_top_candidates], start=1):
        smiles = str(row.get("candidate_smiles", ""))
        if not smiles:
            continue
        ids = du.smiles_identity(smiles, smiles)
        if not ids.get("canonical_smiles"):
            continue
        meta = refinement_meta_for_row(row)
        existing = du.find_existing(existing_rows, ids)
        if existing and not args.redock_refined:
            du.merge_refinement_meta(existing, meta)
            if not existing.get("nickname"):
                existing["nickname"] = f"REFINE_{index:04d}"
            updated_existing += 1
            print(f"[dock-refine] existing seq_id={existing.get('seq_id')} updated_refinement_score={row.get('druglike_refinement_score', '')}", flush=True)
            continue
        if any(du.same_molecule(ids, old) for old in seen):
            continue
        seq = du.seq_id(seq_number)
        seq_number += 1
        seen.append(ids)
        jobs.append(
            {
                "seq": seq,
                "smiles": ids["canonical_smiles"],
                "input_smiles": smiles,
                "nickname": f"REFINE_{index:04d}",
                "ledger_dir": args.history_dir,
                "receptor": receptor,
                "receptor_atoms": receptor_atoms,
                "box": box,
                "vina": vina,
                "meeko": meeko,
                "obabel": obabel,
                "ancestor_smiles": row.get("ancestor_smiles", "") or row.get("parent_smiles", ""),
                "parent_smiles": row.get("parent_smiles", ""),
                "edit_label": f"druglike_refine_g{row.get('generation', '')}_{row.get('edit_label', '')}",
                "refinement_meta": meta,
                "seed": args.seed + index,
                "cluster_cutoff": args.internal_cluster_rmsd_cutoff,
            }
        )
    print(
        f"[dock-refine] new_jobs={len(jobs)} existing_updates={updated_existing} "
        f"dock_workers={args.dock_workers} vina_cpu={box.cpu}",
        flush=True,
    )
    if args.dock_workers > 1 and len(jobs) > 1:
        worker_count = max(1, min(args.dock_workers, len(jobs)))
        try:
            with ProcessPoolExecutor(max_workers=worker_count) as pool:
                for docked in pool.map(dock_refined_worker, jobs):
                    print(
                        f"[dock-refine] {docked.get('seq_id')} affinity={docked.get('affinity_kcal_mol')} "
                        f"reason={docked.get('reason', '')}",
                        flush=True,
                    )
                    new_rows.append(docked)
        except Exception as exc:
            print(f"[dock-refine] parallel docking failed; falling back to serial docking: {exc}", flush=True)
            for job in jobs:
                docked = dock_refined_worker(job)
                print(
                    f"[dock-refine] {docked.get('seq_id')} affinity={docked.get('affinity_kcal_mol')} "
                    f"reason={docked.get('reason', '')}",
                    flush=True,
                )
                new_rows.append(docked)
    else:
        for job in jobs:
            docked = dock_refined_worker(job)
            print(
                f"[dock-refine] {docked.get('seq_id')} affinity={docked.get('affinity_kcal_mol')} "
                f"reason={docked.get('reason', '')}",
                flush=True,
            )
            new_rows.append(docked)
    du.write_csv(ledger_file, existing_rows + new_rows)
    print(f"[dock-refine] wrote {ledger_file} new_rows={len(new_rows)} total_rows={len(existing_rows) + len(new_rows)}", flush=True)


def add_docked_scores(
    rows: list[dict[str, object]],
    history_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    best = best_history_row_by_smiles(history_rows)
    out: list[dict[str, object]] = []
    for row in rows:
        canonical = can(str(row.get("candidate_smiles", "")))
        docked = best.get(canonical)
        merged = dict(row)
        docked_reason = docked.get("reason", "") if docked else ""
        if docked and not docked_reason and docked.get("affinity_kcal_mol") not in ("", None):
            footprint = structural_footprint_score(docked)
            affinity = clean_float(docked.get("affinity_kcal_mol"), default=0.0)
            affinity_component = clamp01(max(0.0, -affinity) / 14.0)
            contact = row_contact_score(docked)
            docked_score = (
                0.58 * clean_float(row.get("druglike_refinement_score"))
                + 0.24 * (10.0 * footprint)
                + 0.12 * (10.0 * affinity_component)
                + 0.06 * (10.0 * clamp01(contact / 70.0))
            )
            merged.update(
                {
                    "docked_seq_id": docked.get("seq_id", ""),
                    "docked_reason": "",
                    "docked_affinity_kcal_mol": docked.get("affinity_kcal_mol", ""),
                    "docked_pose_path": docked.get("pose_path", ""),
                    "docked_contact_score": round(contact, 4),
                    "docked_structural_footprint_score": footprint,
                    "docked_contact_residue_count_4a": docked.get("contact_residue_count_4a", ""),
                    "docked_surface_contact_fraction_4a": docked.get("surface_contact_fraction_4a", ""),
                    "docked_receptor_residue_coverage_4a": docked.get("receptor_residue_coverage_4a", ""),
                    "docked_ligand_atom_contact_fraction_4a": docked.get("ligand_atom_contact_fraction_4a", ""),
                    "docked_refinement_score": round(docked_score, 4),
                }
            )
        else:
            merged.update(
                {
                    "docked_seq_id": docked.get("seq_id", "") if docked else "",
                    "docked_reason": docked_reason,
                    "docked_affinity_kcal_mol": "",
                    "docked_pose_path": "",
                    "docked_contact_score": "",
                    "docked_structural_footprint_score": "",
                    "docked_contact_residue_count_4a": "",
                    "docked_surface_contact_fraction_4a": "",
                    "docked_receptor_residue_coverage_4a": "",
                    "docked_ligand_atom_contact_fraction_4a": "",
                    "docked_refinement_score": "",
                }
            )
        out.append(merged)
    out.sort(key=lambda item: clean_float(item.get("docked_refinement_score"), default=-999.0), reverse=True)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--history-dir", type=Path, default=du.DEFAULT_LEDGER_DIR)
    parser.add_argument("--expert-smiles", action="append", default=[], help="High-confidence expert pocket-filling SMILES; may be repeated")
    parser.add_argument("--reference-smiles-csv", type=Path, help="Optional kinase-like reference SMILES CSV")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_REFINEMENT_DIR)
    parser.add_argument("--top-anchors", type=int, default=12)
    parser.add_argument("--include-history-anchors", action="store_true", help="Also use top historical anchors when expert SMILES are provided")
    parser.add_argument("--max-rounds", "--max-iterations", "--rounds", dest="max_rounds", type=int, default=10, help="Maximum refinement generations before stopping")
    parser.add_argument("--target-refinement-score", type=float, default=5.0, help="Stop early once any candidate reaches this 0-10 score")
    parser.add_argument("--batch-size", "--max-candidates", dest="batch_size", type=int, default=300, help="Maximum new candidates generated per refinement round")
    parser.add_argument("--progress-interval", type=int, default=25, help="Print refinement progress every N generated candidates within a round; 0 disables intra-round progress")
    parser.add_argument("--top-to-dock", type=int, default=50)
    parser.add_argument("--workers", type=int, default=default_refinement_workers(), help="Parallel workers for candidate scoring")
    parser.add_argument("--min-qed", type=float, default=0.0)
    parser.add_argument("--max-mw", type=float, default=1200.0)
    parser.add_argument("--max-rot-bonds", type=int, default=40)
    parser.add_argument("--receptor", type=Path, help="Optional receptor PDB/PDBQT for pose-context enrichment of historical anchors")
    parser.add_argument("--config", type=Path, help="Vina config txt, required when --dock-top-candidates cannot use saved history config")
    parser.add_argument("--vina")
    parser.add_argument("--meeko")
    parser.add_argument("--obabel")
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--exhaustiveness", type=int)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument("--energy-range", type=int)
    parser.add_argument("--internal-cluster-rmsd-cutoff", type=float, default=3.0)
    parser.add_argument("--dock-top-candidates", type=int, default=10, help="Dock the current top N refined candidates, then write docked re-ranked output")
    parser.add_argument("--dock-workers", type=int, default=1, help="Parallel workers for docking refined candidates")
    parser.add_argument("--redock-refined", action="store_true", help="Dock refined candidates even if the same canonical SMILES already exists in history")
    parser.add_argument("--allow-add", action="store_true", help="Allow growth edits during refinement; default is trim/shrink/drastic only")
    parser.add_argument("--no-scaffold-edits", dest="scaffold_edits", action="store_false", help="Disable bioisostere, ring/side-chain replacement, and cyclization edits")
    parser.add_argument("--allow-qve-loss", dest="require_qve_gain", action="store_false", help="Allow scaffold edits even when QED/window/synthetic proxy does not improve")
    parser.add_argument("--min-reference-similarity", type=float, default=0.0, help="Optional hard filter against kinase reference SMILES fingerprints")
    parser.add_argument("--seed", type=int, default=42)
    parser.set_defaults(scaffold_edits=True, require_qve_gain=True)
    args = parser.parse_args()

    args.history_dir = du.normalize_ledger_dir(args.history_dir)
    history_file = args.history_dir / du.LEDGER_FILENAME
    history_rows = du.read_csv(history_file)
    enrich_history_with_pose_context(history_rows, args.receptor)
    anchors = anchor_rows(history_rows, args.top_anchors)
    history_best = best_history_row_by_smiles(history_rows)
    expert_parents: list[dict[str, object]] = []
    for i, smiles in enumerate(args.expert_smiles, start=1):
        canonical = can(smiles)
        if not canonical:
            raise ValueError(f"--expert-smiles #{i} is not a valid RDKit SMILES; refinement would otherwise ignore it")
        matched = history_best.get(canonical, {})
        if matched:
            expert_parent = parent_context_from_row(matched)
            expert_parent["seq_id"] = expert_parent.get("seq_id") or f"EXPERT_{i:02d}"
            expert_parent["nickname"] = expert_parent.get("nickname") or "expert_anchor"
            expert_parents.append(expert_parent)
        else:
            expert_parents.append(
                {
                        "seq_id": f"EXPERT_{i:02d}",
                        "nickname": "expert_anchor",
                        "ancestor_smiles": canonical,
                        "canonical_smiles": canonical,
                        "anchor_score": 0.0,
                        "structural_footprint_score": 0.0,
                }
            )
    parent_map: dict[str, dict[str, object]] = {}
    if expert_parents and not args.include_history_anchors:
        candidate_parents = expert_parents
    else:
        candidate_parents = expert_parents + anchors
    for parent in candidate_parents:
        canonical = str(parent.get("canonical_smiles", ""))
        if canonical and canonical not in parent_map:
            parent_map[canonical] = parent
    parents = list(parent_map.values())
    if not parents:
        raise ValueError("no valid expert or history anchors found")

    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    anchor_fields = [
        "seq_id", "nickname", "canonical_smiles", "anchor_score", "contact_score", "contact_efficiency",
        "structural_footprint_score", "contact_residue_count_4a", "surface_contact_fraction_4a",
        "receptor_residue_coverage_4a", "ligand_atom_contact_fraction_4a",
        "affinity_kcal_mol", "qed", "synthetic_score_proxy", "mw", "logp", "hbd", "hba",
        "tpsa", "rot_bonds", "heavy_atoms", "aromatic_rings", "formal_charge",
    ]
    write_csv(outdir / "anchor_candidates.csv", anchors, anchor_fields)

    reference_fps = load_reference_fps(args.reference_smiles_csv, max_rows=5000)
    reference_pattern_fps = load_reference_pattern_fps(args.reference_smiles_csv, max_rows=5000)
    expert_fps = [fingerprint(str(parent["canonical_smiles"])) for parent in expert_parents]
    expert_fps = [fp for fp in expert_fps if fp is not None]
    def live_score(row: dict[str, object]) -> dict[str, object]:
        return score_candidate(row, reference_fps, reference_pattern_fps, expert_fps)

    generated = generate_candidates(
        parents,
        args.max_rounds,
        args.batch_size,
        args.seed,
        args.allow_add,
        args.scaffold_edits,
        args.require_qve_gain,
        score_fn=live_score,
        target_score=args.target_refinement_score,
        progress_interval=args.progress_interval,
    )
    seed_rows = []
    for parent in expert_parents:
        seed_rows.append(
            {
                "source_parent": parent.get("seq_id", "expert") or "expert",
                "ancestor_smiles": parent.get("canonical_smiles", ""),
                "parent_smiles": parent.get("canonical_smiles", ""),
                "parent_anchor_score": parent.get("anchor_score", ""),
                "parent_structural_footprint_score": parent.get("structural_footprint_score", ""),
                "candidate_smiles": parent.get("canonical_smiles", ""),
                "edit_label": "expert_input",
                "qve_delta": 0.0,
                "generation": 0,
            }
        )
    scored = score_rows(generated + seed_rows, reference_fps, reference_pattern_fps, expert_fps, args.workers)
    scored = [
        row for row in scored
        if float(row.get("qed", 0)) >= args.min_qed
        and float(row.get("mw", math.inf)) <= args.max_mw
        and int(row.get("rot_bonds", 99)) <= args.max_rot_bonds
        and (not reference_fps or float(row.get("kinase_reference_similarity", 0)) >= args.min_reference_similarity)
    ]
    scored.sort(key=lambda row: float(row["druglike_refinement_score"]), reverse=True)
    fields = [
        "druglike_refinement_score", "qed", "synthetic_score_proxy", "property_window_score",
        "qed_component_score", "synthetic_component_score", "property_component_score",
        "parent_similarity_score", "expert_similarity_score", "reference_similarity_score",
        "reference_partial_similarity_score", "qve_component_score", "inherited_structure_score",
        "structural_alert_penalty_score",
        "parent_similarity", "expert_similarity", "kinase_reference_similarity", "kinase_reference_partial_similarity",
        "parent_anchor_score", "parent_structural_footprint_score", "inherited_structure_raw",
        "qve_delta", "structural_alert_penalty", "generation", "source_parent", "edit_label", "candidate_smiles", "parent_smiles", "ancestor_smiles",
        "mw", "logp", "hbd", "hba", "tpsa", "rot_bonds", "heavy_atoms", "aromatic_rings", "formal_charge",
    ]
    write_csv(outdir / "druglike_refinement_candidates.csv", scored, fields)
    write_csv(outdir / "druglike_refinement_ranked.csv", scored[: args.top_to_dock], fields)
    with (outdir / "druglike_refinement_to_dock.smi").open("w", encoding="utf-8") as handle:
        for i, row in enumerate(scored[: args.top_to_dock], start=1):
            handle.write(f"{row['candidate_smiles']} DRUGLIKE_{i:04d}\n")
    if args.dock_top_candidates > 0:
        dock_refined_candidates(scored, args)
        refreshed_history_rows = du.read_csv(history_file)
        enrich_history_with_pose_context(refreshed_history_rows, args.receptor)
        docked_scored = add_docked_scores(scored, refreshed_history_rows)
        docked_fields = [
            "docked_refinement_score", "docked_seq_id", "docked_reason", "docked_affinity_kcal_mol",
            "docked_structural_footprint_score", "docked_contact_score",
            "docked_contact_residue_count_4a", "docked_surface_contact_fraction_4a",
            "docked_receptor_residue_coverage_4a", "docked_ligand_atom_contact_fraction_4a",
            "docked_pose_path",
            *fields,
        ]
        write_csv(outdir / "druglike_refinement_ranked.csv", docked_scored[: args.top_to_dock], docked_fields)
        write_csv(outdir / "druglike_refinement_docked_ranked.csv", docked_scored[: args.top_to_dock], docked_fields)
    print(
        f"parents={len(parents)} max_rounds={args.max_rounds} batch_size={args.batch_size} "
        f"target_score={args.target_refinement_score} generated={len(generated)} kept={len(scored)} "
        f"completed_rounds={max([int(row.get('generation', 0)) for row in generated] or [0])} "
        f"wrote={outdir / 'druglike_refinement_ranked.csv'}",
        flush=True,
    )


if __name__ == "__main__":
    main()
