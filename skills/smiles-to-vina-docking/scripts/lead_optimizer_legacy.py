#!/usr/bin/env python
"""Bounded SMILES analog enumeration and AutoDock Vina screening.

This script is intentionally conservative: generate small RDKit-valid analogs,
prepare 3D/PDBQT ligands, dock them with a fixed Vina box, and record every
accepted or failed candidate for iterative lead optimization.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError:  # pragma: no cover - depends on RDKit build
    rdMolStandardize = None  # type: ignore[assignment]


RDLogger.DisableLog("rdApp.warning")
RDLogger.DisableLog("rdApp.error")
AFFINITY_LOG_RE = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+", re.MULTILINE)
AFFINITY_PDBQT_RE = re.compile(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)")
ATTEMPTED_SMILES_FILENAME = "attempted_smiles.smi"
SCREENING_FIELDNAMES = [
    "event_id",
    "nickname",
    "parent",
    "edit_label",
    "input_smiles",
    "smiles",
    "canonical_smiles",
    "standard_smiles",
    "tautomer_smiles",
    "inchikey",
    "affinity_kcal_mol",
    "live_min_pose_rmsd_to_lead",
    "live_near_lead_mode_count",
    "live_mode_count",
    "live_inner_rmsd",
    "live_inner_cluster_fraction",
    "live_whole_rmsd",
    "pdbqt_path",
    "pose_path",
    "log_path",
    "mw",
    "logp",
    "hbd",
    "hba",
    "tpsa",
    "rot_bonds",
    "heavy_atoms",
    "formal_charge",
]
REJECTED_FIELDNAMES = [
    "round",
    "parent_id",
    "edit_label",
    "input_smiles",
    "canonical_smiles",
    "reason",
    "mw",
    "logp",
    "hbd",
    "hba",
    "tpsa",
    "rot_bonds",
    "heavy_atoms",
    "formal_charge",
]
POSE_RMSD_FIELDNAMES = [
    "candidate_id",
    "status",
    "affinity_kcal_mol",
    "delta_vs_baseline",
    "pose_path",
    "mode_count",
    "matched_heavy_atoms",
    "lead_heavy_atoms",
    "candidate_mode_heavy_atoms_min",
    "min_pose_rmsd_to_lead",
    "best_rmsd_mode",
    "near_lead_mode_count",
    "near_lead_mode_fraction",
    "best_internal_cluster_size",
    "best_internal_cluster_fraction",
    "best_internal_cluster_representative_mode",
    "best_internal_cluster_mean_rmsd",
    "rmsd_q25_cutoff",
    "rmsd_q25_keep",
    "smiles",
]
POSE_MODE_RMSD_FIELDNAMES = [
    "candidate_id",
    "mode_index",
    "mode_affinity_kcal_mol",
    "mode_rmsd_to_lead",
    "within_rmsd_cutoff",
    "matched_heavy_atoms",
    "mode_heavy_atoms",
    "lead_heavy_atoms",
    "pose_path",
    "smiles",
]
FINAL_RANK_FIELDNAMES = [
    "rank",
    "mol_id",
    "nickname",
    "event_id",
    "parent",
    "edit_label",
    "input_smiles",
    "smiles",
    "canonical_smiles",
    "standard_smiles",
    "tautomer_smiles",
    "inchikey",
    "final_score",
    "affinity_score",
    "inner_score",
    "lead_family_score",
    "whole_score",
    "interaction_score",
    "affinity_kcal_mol",
    "classification",
    "lead_like_score",
    "alternative_score",
    "lead_difference_rmsd",
    "inner_rmsd",
    "whole_rmsd",
    "inner_cluster_fraction",
    "hbond_count",
    "hydrophobic_count",
    "vdw_contact_count",
    "pi_contact_count",
    "ch_pi_count",
    "near_lead_mode_count",
    "mode_count",
    "best_affinity_mode",
    "pose_path",
    "log_path",
]
COMBINED_SUMMARY_FIELDNAMES = FINAL_RANK_FIELDNAMES
COMPACT_MODE_FIELDNAMES = [
    "mol_id",
    "candidate_id",
    "mode_index",
    "mode_affinity_kcal_mol",
    "mode_rmsd_to_lead",
    "within_rmsd_cutoff",
    "pose_path",
]
NEXT_PARENT_FIELDNAMES = [
    "smiles",
    "nickname",
    "event_id",
    "parent_rank",
    "candidate_id",
    "selection_reason",
    "lead_difference_rmsd",
    "inner_rmsd",
    "affinity_kcal_mol",
]
LEAD_FIELDNAMES = [
    "lead_index",
    "first_seen_at",
    "event_id",
    "candidate_id",
    "input_smiles",
    "smiles",
    "canonical_smiles",
    "standard_smiles",
    "tautomer_smiles",
    "inchikey",
    "affinity_kcal_mol",
    "live_min_pose_rmsd_to_lead",
    "live_near_lead_mode_count",
    "live_mode_count",
    "live_inner_rmsd",
    "live_inner_cluster_fraction",
    "live_whole_rmsd",
    "pose_path",
    "log_path",
    "mw",
    "logp",
    "hbd",
    "hba",
    "tpsa",
    "rot_bonds",
    "heavy_atoms",
    "formal_charge",
]
ALTERNATIVE_MODE_FIELDNAMES = [
    "rank",
    "candidate_id",
    "alternative_score",
    "affinity_kcal_mol",
    "min_pose_rmsd_to_lead",
    "near_lead_mode_count",
    "mode_count",
    "best_internal_cluster_size",
    "best_internal_cluster_fraction",
    "best_internal_cluster_representative_mode",
    "best_internal_cluster_mean_rmsd",
    "pose_path",
    "smiles",
]
HISTORY_SUMMARY_FIELDNAMES = [
    "mol_id",
    "nickname",
    "event_id",
    "parent",
    "edit_label",
    "input_smiles",
    "smiles",
    "canonical_smiles",
    "standard_smiles",
    "tautomer_smiles",
    "inchikey",
    "final_score",
    "affinity_score",
    "inner_score",
    "lead_family_score",
    "whole_score",
    "interaction_score",
    "best_affinity_kcal_mol",
    "best_rmsd_to_lead",
    "inner_rmsd",
    "whole_rmsd",
    "inner_cluster_fraction",
    "hbond_count",
    "hydrophobic_count",
    "vdw_contact_count",
    "pi_contact_count",
    "ch_pi_count",
    "classification",
    "pose_path",
    "log_path",
    "source_file",
]


@dataclass(frozen=True)
class ConfigBox:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    exhaustiveness: int = 8
    cpu: int = 1
    energy_range: int = 3


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    parent_id: str
    smiles: str
    edit_label: str
    round_index: int


@dataclass(frozen=True)
class PdbqtMode:
    mode_index: int
    affinity: float | None
    coords: list[tuple[float, float, float]]
    atoms: list[dict[str, str | float]]


PDBQT_ELEMENT_FALLBACKS = {
    "A": "C",
    "C": "C",
    "HD": "H",
    "HS": "H",
    "N": "N",
    "NA": "N",
    "NS": "N",
    "OA": "O",
    "OS": "O",
    "O": "O",
    "S": "S",
    "SA": "S",
    "P": "P",
    "F": "F",
    "CL": "Cl",
    "BR": "Br",
    "I": "I",
}
HBOND_ELEMENTS = {"N", "O", "S"}
HYDROPHOBIC_ELEMENTS = {"C", "Cl", "Br", "I", "F"}
AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS"}
SMILES_IDENTITY_CACHE: dict[tuple[str, str], dict[str, str]] = {}


def run(cmd: list[str], cwd: Path | None = None) -> tuple[int, str]:
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
    except FileNotFoundError as exc:
        return 127, str(exc)
    return proc.returncode, proc.stdout


def command_exists(command: str) -> bool:
    if str(Path(command)).count("/") or str(Path(command)).count("\\"):
        return Path(command).exists()
    return shutil.which(command) is not None


def canonical(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def clean_smiles_text(smiles: str | None) -> str:
    if smiles is None:
        return ""
    return re.sub(r"\s+", "", str(smiles).strip())


def standardized_mol(mol: Chem.Mol) -> Chem.Mol:
    out = Chem.Mol(mol)
    if rdMolStandardize is None:
        return out
    try:
        out = rdMolStandardize.Cleanup(out)
    except Exception:
        out = Chem.Mol(mol)
    try:
        out = rdMolStandardize.FragmentParent(out)
    except Exception:
        pass
    try:
        out = rdMolStandardize.Uncharger().uncharge(out)
    except Exception:
        pass
    try:
        Chem.SanitizeMol(out)
    except Exception:
        out = Chem.Mol(mol)
    return out


def tautomer_canonical_smiles(mol: Chem.Mol) -> str:
    if rdMolStandardize is None:
        return Chem.MolToSmiles(mol, canonical=True)
    try:
        tautomer = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
        return Chem.MolToSmiles(tautomer, canonical=True)
    except Exception:
        return Chem.MolToSmiles(mol, canonical=True)


def smiles_identity_fields(smiles: str | None, input_smiles: str | None = None) -> dict[str, str]:
    text = clean_smiles_text(smiles)
    original = clean_smiles_text(input_smiles) or text
    cache_key = (text, original)
    cached = SMILES_IDENTITY_CACHE.get(cache_key)
    if cached is not None:
        return dict(cached)
    if not text:
        fields = {
            "input_smiles": original,
            "smiles": "",
            "canonical_smiles": "",
            "standard_smiles": "",
            "tautomer_smiles": "",
            "inchikey": "",
        }
        SMILES_IDENTITY_CACHE[cache_key] = fields
        return dict(fields)
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        fields = {
            "input_smiles": original,
            "smiles": text,
            "canonical_smiles": "",
            "standard_smiles": "",
            "tautomer_smiles": "",
            "inchikey": "",
        }
        SMILES_IDENTITY_CACHE[cache_key] = fields
        return dict(fields)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        fields = {
            "input_smiles": original,
            "smiles": text,
            "canonical_smiles": "",
            "standard_smiles": "",
            "tautomer_smiles": "",
            "inchikey": "",
        }
        SMILES_IDENTITY_CACHE[cache_key] = fields
        return dict(fields)
    can = Chem.MolToSmiles(mol, canonical=True)
    standard = standardized_mol(mol)
    standard_can = Chem.MolToSmiles(standard, canonical=True)
    tautomer_can = tautomer_canonical_smiles(standard)
    try:
        inchikey = Chem.MolToInchiKey(standard)
    except Exception:
        inchikey = ""
    fields = {
        "input_smiles": original,
        "smiles": can,
        "canonical_smiles": can,
        "standard_smiles": standard_can,
        "tautomer_smiles": tautomer_can,
        "inchikey": inchikey,
    }
    SMILES_IDENTITY_CACHE[cache_key] = fields
    return dict(fields)


def ensure_light_smiles_identity(row: dict[str, str | float | int]) -> dict[str, str | float | int]:
    enriched = dict(row)
    smiles = row_smiles_text(enriched)
    if not smiles:
        return enriched
    can = str(enriched.get("canonical_smiles") or enriched.get("smiles") or "")
    if not can:
        can = canonical(smiles) or smiles
    enriched.setdefault("input_smiles", str(row.get("input_smiles", "") or smiles))
    enriched["smiles"] = can
    enriched["canonical_smiles"] = str(enriched.get("canonical_smiles") or can)
    enriched.setdefault("standard_smiles", "")
    enriched.setdefault("tautomer_smiles", "")
    enriched.setdefault("inchikey", "")
    return enriched


def row_smiles_text(row: dict[str, str | float | int]) -> str:
    for field in ("smiles", "canonical_smiles", "standard_smiles", "tautomer_smiles", "input_smiles", "mol_smiles"):
        value = clean_smiles_text(str(row.get(field, "") or ""))
        if value:
            return value
    return ""


def nickname_is_machine_alias(row: dict[str, str | float | int], nickname: str) -> bool:
    if not nickname:
        return True
    machine_aliases = {
        str(row.get("mol_id", "") or ""),
        str(row.get("candidate_id", "") or ""),
        str(row.get("event_id", "") or ""),
        row_event_id(row),
    }
    parent_id = str(row.get("parent_id", "") or row.get("parent", "") or "")
    edit_label = str(row.get("edit_label", "") or "")
    if parent_id and edit_label:
        machine_aliases.add(f"{parent_id}+{edit_label}")
    return nickname in {alias for alias in machine_aliases if alias}


def load_history(history_dirs: list[Path]) -> tuple[set[str], set[str], list[str]]:
    seen_smiles: set[str] = set()
    seen_edit_labels: set[str] = set()
    notes: list[str] = []
    for history_dir in history_dirs:
        directory = history_dir.expanduser().resolve()
        if not directory.exists():
            notes.append(f"history_missing={directory}")
            continue
        file_count = 0
        row_count = 0
        attempted_path = directory / ATTEMPTED_SMILES_FILENAME
        if attempted_path.exists():
            file_count += 1
            try:
                for line in attempted_path.read_text(encoding="utf-8", errors="replace").splitlines():
                    parts = line.strip().split()
                    if not parts:
                        continue
                    can = canonical(parts[0])
                    if can:
                        row_count += 1
                        seen_smiles.add(can)
            except Exception as exc:
                notes.append(f"attempted_smiles_read_failed={attempted_path}: {exc}")
        for csv_path in directory.rglob("*.csv"):
            file_count += 1
            try:
                with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
                    reader = csv.DictReader(handle)
                    for row in reader:
                        row_count += 1
                        smiles = row.get("smiles") or row.get("canonical_smiles") or row.get("input_smiles") or ""
                        can = canonical(smiles) if smiles else None
                        if can:
                            seen_smiles.add(can)
                        edit_label = row.get("edit_label") or ""
                        if edit_label:
                            seen_edit_labels.add(edit_label)
            except Exception as exc:
                notes.append(f"history_read_failed={csv_path}: {exc}")
        notes.append(
            f"history_loaded={directory} files={file_count} rows={row_count} unique_smiles_total={len(seen_smiles)} unique_edit_labels_total={len(seen_edit_labels)}"
        )
    return seen_smiles, seen_edit_labels, notes


def rejection_reason_counts(rows: list[dict[str, str | int | float]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        reason = str(row.get("reason", "unknown") or "unknown")
        counts[reason] = counts.get(reason, 0) + 1
    return dict(sorted(counts.items(), key=lambda item: (-item[1], item[0])))


def safe_float(value: str | int | float | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def host_path(path_text: str | Path | None) -> Path | None:
    if path_text is None:
        return None
    raw = str(path_text).strip()
    if not raw:
        return None
    direct = Path(raw).expanduser()
    if direct.exists():
        return direct
    wsl_match = re.match(r"^/mnt/([A-Za-z])/(.*)$", raw)
    if wsl_match:
        drive = wsl_match.group(1).upper()
        rest = wsl_match.group(2).replace("/", "\\")
        windows_path = Path(f"{drive}:\\{rest}")
        if windows_path.exists():
            return windows_path
    windows_match = re.match(r"^([A-Za-z]):\\(.*)$", raw)
    if windows_match:
        drive = windows_match.group(1).lower()
        rest = windows_match.group(2).replace("\\", "/")
        wsl_path = Path(f"/mnt/{drive}/{rest}")
        if wsl_path.exists():
            return wsl_path
    return direct


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def mol_id_for_index(index: int) -> str:
    return f"M{index:05d}"


def row_event_id(row: dict[str, str | float | int]) -> str:
    existing = str(row.get("event_id", "") or "")
    if existing:
        return existing
    candidate_id = str(row.get("candidate_id", "") or "")
    run_id = str(row.get("run_id", "") or "")
    if run_id and candidate_id:
        return f"{run_id}:{candidate_id}"
    return candidate_id or run_id


def split_event_id(event_id: str) -> tuple[str, str]:
    if ":" in event_id:
        run_id, candidate_id = event_id.rsplit(":", 1)
        return run_id, candidate_id
    return "", event_id


def row_nickname(row: dict[str, str | float | int]) -> str:
    existing = str(row.get("nickname", "") or "")
    if nickname_is_machine_alias(row, existing):
        return ""
    return existing


def enrich_display_row(row: dict[str, str | float | int]) -> dict[str, str | float | int]:
    enriched = ensure_light_smiles_identity(row)
    enriched["event_id"] = row_event_id(enriched)
    enriched["nickname"] = row_nickname(enriched)
    enriched["parent"] = enriched.get("parent") or enriched.get("parent_id", "")
    return enriched


def collect_user_nicknames_from_rows(rows: list[dict[str, str | float | int]]) -> dict[str, str]:
    nicknames: dict[str, str] = {}
    for row in rows:
        nickname = row_nickname(row)
        if not nickname:
            continue
        can = canonical(row_smiles_text(row))
        if can and can not in nicknames:
            nicknames[can] = nickname
    return nicknames


def load_user_nickname_map(outdir: Path) -> dict[str, str]:
    nicknames: dict[str, str] = {}
    for csv_name in (
        "final_ranked_candidates.csv",
        "history_summary.csv",
        "all_screening_results.csv",
        "lead.csv",
        "next_round_parents.csv",
    ):
        csv_path = outdir / csv_name
        if not csv_path.exists():
            continue
        try:
            rows = read_csv_rows(csv_path)
        except Exception:
            continue
        for can, nickname in collect_user_nicknames_from_rows(rows).items():
            nicknames.setdefault(can, nickname)
    return nicknames


def apply_user_nicknames(
    rows: list[dict[str, str | float | int]],
    nickname_map: dict[str, str],
) -> list[dict[str, str | float | int]]:
    out: list[dict[str, str | float | int]] = []
    for row in rows:
        updated = dict(row)
        can = canonical(row_smiles_text(updated))
        existing = row_nickname(updated)
        if existing:
            if can:
                nickname_map.setdefault(can, existing)
            updated["nickname"] = existing
        elif can and can in nickname_map:
            updated["nickname"] = nickname_map[can]
        else:
            updated["nickname"] = ""
        out.append(updated)
    return out


def compact_molecule_rows(
    combined_rows: list[dict[str, str | float | int]],
    all_screening_rows: list[dict[str, str | float | int]],
) -> list[dict[str, str | float | int]]:
    screening_by_key: dict[tuple[str, str], dict[str, str | float | int]] = {}
    screening_by_smiles: dict[str, dict[str, str | float | int]] = {}
    for row in all_screening_rows:
        smiles = row_smiles_text(row)
        can = canonical(smiles) if smiles else None
        if not can:
            continue
        key = (str(row.get("candidate_id", "")), can)
        screening_by_key[key] = row
        screening_by_smiles.setdefault(can, row)

    best_by_smiles: dict[str, dict[str, str | float | int]] = {}
    for row in combined_rows:
        smiles = row_smiles_text(row)
        can = canonical(smiles) if smiles else None
        if not can:
            continue
        current = best_by_smiles.get(can)
        if current is None:
            best_by_smiles[can] = dict(row)
            continue
        current_aff = safe_float(current.get("affinity_kcal_mol"))
        new_aff = safe_float(row.get("affinity_kcal_mol"))
        current_score = safe_float(current.get("lead_like_score")) or safe_float(current.get("alternative_score")) or -999.0
        new_score = safe_float(row.get("lead_like_score")) or safe_float(row.get("alternative_score")) or -999.0
        current_rmsd = safe_float(current.get("lead_difference_rmsd"))
        new_rmsd = safe_float(row.get("lead_difference_rmsd"))
        if (
            new_aff if new_aff is not None else 999.0,
            new_rmsd if new_rmsd is not None else 999.0,
            -new_score,
        ) < (
            current_aff if current_aff is not None else 999.0,
            current_rmsd if current_rmsd is not None else 999.0,
            -current_score,
        ):
            best_by_smiles[can] = dict(row)

    for row in all_screening_rows:
        if row.get("status") not in {"ok", "reference", ""}:
            continue
        smiles = row_smiles_text(row)
        can = canonical(smiles) if smiles else None
        if not can:
            continue
        identity = smiles_identity_fields(can, str(row.get("input_smiles", "") or ""))
        fallback_row = {
            "candidate_id": row.get("candidate_id", ""),
            "classification": "history_unscored_pose",
            "final_score": "",
                "affinity_score": "",
                "inner_score": "",
                "lead_family_score": "",
                "whole_score": "",
                "interaction_score": "",
                "lead_like_score": "",
            "alternative_score": "",
            "affinity_kcal_mol": row.get("affinity_kcal_mol") or row.get("best_affinity_kcal_mol") or "",
            "lead_difference_rmsd": row.get("lead_difference_rmsd") or row.get("min_pose_rmsd_to_lead") or row.get("live_min_pose_rmsd_to_lead") or row.get("best_rmsd_to_lead") or "",
            "inner_rmsd": row.get("inner_rmsd") or row.get("best_internal_cluster_mean_rmsd") or row.get("live_inner_rmsd") or "",
            "whole_rmsd": row.get("whole_rmsd") or row.get("live_whole_rmsd") or "",
            "hbond_count": row.get("hbond_count", ""),
            "hydrophobic_count": row.get("hydrophobic_count", ""),
            "vdw_contact_count": row.get("vdw_contact_count", ""),
            "pi_contact_count": row.get("pi_contact_count", ""),
            "ch_pi_count": row.get("ch_pi_count", ""),
            "inner_cluster_fraction": row.get("inner_cluster_fraction") or row.get("best_internal_cluster_fraction") or row.get("live_inner_cluster_fraction") or "",
            "near_lead_mode_count": row.get("near_lead_mode_count") or row.get("live_near_lead_mode_count") or "",
            "mode_count": row.get("mode_count") or row.get("live_mode_count") or "",
            "best_affinity_mode": row.get("best_affinity_mode", ""),
            "pose_path": row.get("pose_path", ""),
            "smiles": can,
            **identity,
        }
        current = best_by_smiles.get(can)
        if current is None:
            best_by_smiles[can] = fallback_row
            continue
        if safe_float(current.get("final_score")) is not None:
            continue
        current_aff = safe_float(current.get("affinity_kcal_mol"))
        fallback_aff = safe_float(fallback_row.get("affinity_kcal_mol"))
        if fallback_aff is not None and (current_aff is None or fallback_aff < current_aff):
            for key in ("affinity_kcal_mol", "candidate_id", "pose_path", "lead_difference_rmsd", "inner_rmsd", "whole_rmsd", "inner_cluster_fraction", "near_lead_mode_count", "mode_count"):
                if fallback_row.get(key) != "":
                    current[key] = fallback_row[key]
            if current.get("classification") == "":
                current["classification"] = fallback_row["classification"]

    rows: list[dict[str, str | float | int]] = []
    for index, (can, row) in enumerate(
        sorted(
            best_by_smiles.items(),
            key=lambda item: (
                safe_float(item[1].get("final_score")) is None,
                -(safe_float(item[1].get("final_score")) if safe_float(item[1].get("final_score")) is not None else -1.0),
                safe_float(item[1].get("affinity_kcal_mol")) if safe_float(item[1].get("affinity_kcal_mol")) is not None else 999.0,
                safe_float(item[1].get("lead_difference_rmsd")) if safe_float(item[1].get("lead_difference_rmsd")) is not None else 999.0,
            ),
        ),
        start=1,
    ):
        screening = screening_by_key.get((str(row.get("candidate_id", "")), can), screening_by_smiles.get(can, {}))
        identity = smiles_identity_fields(can, str(row.get("input_smiles", "") or screening.get("input_smiles", "") or ""))
        rows.append(
            {
                "rank": index,
                "mol_id": mol_id_for_index(index),
                **identity,
                "smiles": can,
                "final_score": row.get("final_score", ""),
                "affinity_score": row.get("affinity_score", ""),
                "inner_score": row.get("inner_score", ""),
                "lead_family_score": row.get("lead_family_score", ""),
                "whole_score": row.get("whole_score", ""),
                "interaction_score": row.get("interaction_score", ""),
                "affinity_kcal_mol": row.get("affinity_kcal_mol", ""),
                "classification": row.get("classification", ""),
                "lead_like_score": row.get("lead_like_score", ""),
                "alternative_score": row.get("alternative_score", ""),
                "lead_difference_rmsd": row.get("lead_difference_rmsd", ""),
                "inner_rmsd": row.get("inner_rmsd", ""),
                "whole_rmsd": row.get("whole_rmsd", ""),
                "inner_cluster_fraction": row.get("inner_cluster_fraction", ""),
                "hbond_count": row.get("hbond_count", ""),
                "hydrophobic_count": row.get("hydrophobic_count", ""),
                "vdw_contact_count": row.get("vdw_contact_count", ""),
                "pi_contact_count": row.get("pi_contact_count", ""),
                "ch_pi_count": row.get("ch_pi_count", ""),
                "near_lead_mode_count": row.get("near_lead_mode_count", ""),
                "mode_count": row.get("mode_count", ""),
                "best_affinity_mode": row.get("best_affinity_mode", ""),
                "candidate_id": row.get("candidate_id", ""),
                "run_id": screening.get("run_id", ""),
                "parent_id": screening.get("parent_id", ""),
                "edit_label": screening.get("edit_label", ""),
                "pose_path": row.get("pose_path", ""),
                "log_path": screening.get("log_path", ""),
            }
        )
    return rows


def molecule_rows_to_history_rows(molecule_rows: list[dict[str, str | float | int]]) -> list[dict[str, str | float | int]]:
    return [
        {
            "mol_id": row.get("mol_id", ""),
            "smiles": row.get("smiles", ""),
            "input_smiles": row.get("input_smiles", ""),
            "canonical_smiles": row.get("canonical_smiles", ""),
            "standard_smiles": row.get("standard_smiles", ""),
            "tautomer_smiles": row.get("tautomer_smiles", ""),
            "inchikey": row.get("inchikey", ""),
            "final_score": row.get("final_score", ""),
            "affinity_score": row.get("affinity_score", ""),
            "inner_score": row.get("inner_score", ""),
            "lead_family_score": row.get("lead_family_score", ""),
            "whole_score": row.get("whole_score", ""),
            "interaction_score": row.get("interaction_score", ""),
            "best_affinity_kcal_mol": row.get("affinity_kcal_mol", ""),
            "best_rmsd_to_lead": row.get("lead_difference_rmsd", ""),
            "inner_rmsd": row.get("inner_rmsd", ""),
            "whole_rmsd": row.get("whole_rmsd", ""),
            "inner_cluster_fraction": row.get("inner_cluster_fraction", ""),
            "hbond_count": row.get("hbond_count", ""),
            "hydrophobic_count": row.get("hydrophobic_count", ""),
            "vdw_contact_count": row.get("vdw_contact_count", ""),
            "pi_contact_count": row.get("pi_contact_count", ""),
            "ch_pi_count": row.get("ch_pi_count", ""),
            "classification": row.get("classification", ""),
            "candidate_id": row.get("candidate_id", ""),
            "run_id": row.get("run_id", ""),
            "parent_id": row.get("parent_id", ""),
            "edit_label": row.get("edit_label", ""),
            "pose_path": row.get("pose_path", ""),
            "log_path": row.get("log_path", ""),
            "source_file": "final_ranked_candidates.csv",
        }
        for row in molecule_rows
    ]


def compact_mode_rows(
    mode_rows: list[dict[str, str | float | int]],
    molecule_rows: list[dict[str, str | float | int]],
) -> list[dict[str, str | float | int]]:
    selected_by_smiles: dict[str, tuple[str, str, str]] = {}
    for row in molecule_rows:
        smiles = row_smiles_text(row)
        if smiles:
            selected_by_smiles[smiles] = (
                str(row.get("mol_id", "")),
                str(row.get("candidate_id", "")),
                str(row.get("pose_path", "")),
            )
    compact: list[dict[str, str | float | int]] = []
    seen: set[tuple[str, str, str]] = set()
    for row in mode_rows:
        smiles = canonical(row_smiles_text(row)) or row_smiles_text(row)
        selected = selected_by_smiles.get(smiles)
        if not selected:
            continue
        mol_id, selected_candidate_id, selected_pose_path = selected
        if str(row.get("candidate_id", "")) != selected_candidate_id or str(row.get("pose_path", "")) != selected_pose_path:
            continue
        key = (mol_id, str(row.get("pose_path", "")), str(row.get("mode_index", "")))
        if key in seen:
            continue
        seen.add(key)
        compact.append(
            {
                "mol_id": mol_id,
                "candidate_id": row.get("candidate_id", ""),
                "mode_index": row.get("mode_index", ""),
                "mode_affinity_kcal_mol": row.get("mode_affinity_kcal_mol", ""),
                "mode_rmsd_to_lead": row.get("mode_rmsd_to_lead", ""),
                "within_rmsd_cutoff": row.get("within_rmsd_cutoff", ""),
                "pose_path": row.get("pose_path", ""),
            }
        )
    compact.sort(key=lambda row: (str(row.get("mol_id", "")), int(row.get("mode_index", 0) or 0)))
    return compact


def summarize_history(history_dirs: list[Path]) -> list[dict[str, str | float | int]]:
    summary: dict[str, dict[str, str | float | int]] = {}
    for history_dir in history_dirs:
        directory = history_dir.expanduser().resolve()
        if not directory.exists():
            continue
        csv_paths = (
            list(directory.rglob("history_summary.csv"))
            + list(directory.rglob("all_screening_results.csv"))
            + list(directory.rglob("candidate_pose_summary.csv"))
            + list(directory.rglob("pose_rmsd_analysis.csv"))
            + list(directory.rglob("final_ranked_candidates.csv"))
            + list(directory.rglob("screening_results.csv"))
        )
        seen_paths: set[Path] = set()
        ordered_csv_paths: list[Path] = []
        for csv_path in csv_paths:
            resolved = csv_path.resolve()
            if resolved not in seen_paths:
                seen_paths.add(resolved)
                ordered_csv_paths.append(csv_path)
        for csv_path in ordered_csv_paths:
            try:
                with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
                    for row in csv.DictReader(handle):
                        smiles = row_smiles_text(row)
                        can = canonical(smiles) if smiles else None
                        if not can:
                            continue
                        edit_label = row.get("edit_label") or row.get("candidate_id") or ""
                        key = can
                        current = summary.get(
                            key,
                            {
                                "mol_id": row.get("mol_id", ""),
                                "candidate_id": row.get("candidate_id", ""),
                                "run_id": row.get("run_id", ""),
                                "parent_id": row.get("parent_id", ""),
                                "edit_label": edit_label,
                                "smiles": can,
                                "best_affinity_kcal_mol": "",
                                "best_rmsd_to_lead": "",
                                "inner_rmsd": row.get("inner_rmsd", ""),
                                "inner_cluster_fraction": row.get("inner_cluster_fraction", ""),
                                "classification": row.get("classification", ""),
                                "pose_path": row.get("pose_path", ""),
                                "log_path": row.get("log_path", ""),
                                "source_file": str(csv_path),
                            },
                        )
                        for field in (
                            "mol_id",
                            "candidate_id",
                            "run_id",
                            "parent_id",
                            "edit_label",
                            "classification",
                            "pose_path",
                            "log_path",
                            "inner_rmsd",
                            "inner_cluster_fraction",
                        ):
                            value = row.get(field, "")
                            if value and not current.get(field):
                                current[field] = value
                        affinity = safe_float(row.get("affinity_kcal_mol") or row.get("best_affinity_kcal_mol"))
                        current_affinity = safe_float(current.get("best_affinity_kcal_mol"))
                        if affinity is not None and (current_affinity is None or affinity < current_affinity):
                            current["best_affinity_kcal_mol"] = affinity
                            current["source_file"] = str(csv_path)
                        rmsd = safe_float(row.get("lead_difference_rmsd") or row.get("min_pose_rmsd_to_lead") or row.get("live_min_pose_rmsd_to_lead") or row.get("best_rmsd_to_lead"))
                        current_rmsd = safe_float(current.get("best_rmsd_to_lead"))
                        if rmsd is not None and (current_rmsd is None or rmsd < current_rmsd):
                            current["best_rmsd_to_lead"] = rmsd
                        summary[key] = current
            except Exception:
                continue
    rows = list(summary.values())
    rows.sort(
        key=lambda row: (
            safe_float(row.get("best_rmsd_to_lead")) if safe_float(row.get("best_rmsd_to_lead")) is not None else 999.0,
            safe_float(row.get("best_affinity_kcal_mol")) if safe_float(row.get("best_affinity_kcal_mol")) is not None else 999.0,
        )
    )
    for index, row in enumerate(rows, start=1):
        if not row.get("mol_id"):
            row["mol_id"] = mol_id_for_index(index)
    return rows


def load_history_screening_rows(history_dirs: list[Path]) -> list[dict[str, str | float | int]]:
    loaded: list[dict[str, str | float | int]] = []
    seen_keys: set[tuple[str, str, str]] = set()
    for history_dir in history_dirs:
        directory = history_dir.expanduser().resolve()
        if not directory.exists():
            continue
        csv_paths = (
            list(directory.rglob("all_screening_results.csv"))
            + list(directory.rglob("history_summary.csv"))
            + list(directory.rglob("screening_results.csv"))
        )
        for csv_path in csv_paths:
            try:
                with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
                    for row in csv.DictReader(handle):
                        status = row.get("status") or ("ok" if row.get("smiles") else "")
                        if status != "ok":
                            continue
                        pose_path = row.get("pose_path") or ""
                        resolved_pose_path = host_path(pose_path)
                        event_id = row.get("event_id") or ""
                        event_run_id, event_candidate_id = split_event_id(event_id)
                        candidate_id = row.get("candidate_id") or event_candidate_id or row.get("nickname") or ""
                        smiles = row.get("smiles") or ""
                        can = canonical(smiles) if smiles else None
                        if not can:
                            continue
                        run_id = row.get("run_id") or event_run_id or ""
                        key = (candidate_id, run_id, can)
                        if key in seen_keys:
                            continue
                        seen_keys.add(key)
                        normalized = {field: row.get(field, "") for field in SCREENING_FIELDNAMES}
                        normalized["candidate_id"] = candidate_id
                        normalized["run_id"] = run_id
                        normalized["parent_id"] = row.get("parent_id") or row.get("parent") or ""
                        normalized["parent_smiles"] = row.get("parent_smiles", "")
                        normalized["status"] = "ok"
                        normalized.update(ensure_light_smiles_identity({**row, "smiles": can, "input_smiles": row.get("input_smiles") or smiles}))
                        normalized["affinity_kcal_mol"] = row.get("affinity_kcal_mol") or row.get("best_affinity_kcal_mol") or ""
                        normalized["live_min_pose_rmsd_to_lead"] = row.get("live_min_pose_rmsd_to_lead") or row.get("lead_difference_rmsd") or row.get("best_rmsd_to_lead") or ""
                        normalized["live_inner_rmsd"] = row.get("live_inner_rmsd") or row.get("inner_rmsd") or ""
                        normalized["live_inner_cluster_fraction"] = row.get("live_inner_cluster_fraction") or row.get("inner_cluster_fraction") or ""
                        normalized["live_whole_rmsd"] = row.get("live_whole_rmsd") or row.get("whole_rmsd") or ""
                        normalized["affinity_recheck_delta"] = row.get("affinity_recheck_delta", "")
                        normalized["affinity_recheck_status"] = row.get("affinity_recheck_status", "")
                        if "classification" in normalized:
                            normalized["classification"] = row.get("classification", "")
                        normalized["pose_path"] = str(resolved_pose_path) if resolved_pose_path else pose_path
                        loaded.append(normalized)
            except Exception:
                continue
    return loaded


def merge_screening_rows(history_rows: list[dict[str, str | float | int]], current_rows: list[dict[str, str | float | int]]) -> list[dict[str, str | float | int]]:
    merged: dict[tuple[str, str], dict[str, str | float | int]] = {}
    for row in history_rows + current_rows:
        if not is_successful_pose_screening_row(row):
            continue
        candidate_id = str(row.get("candidate_id", ""))
        run_id = str(row.get("run_id", ""))
        pose_path = str(row.get("pose_path", ""))
        smiles = str(row.get("smiles", ""))
        key = (candidate_id or smiles, run_id or pose_path or smiles)
        merged[key] = row
    return list(merged.values())


def is_successful_pose_screening_row(row: dict[str, str | float | int]) -> bool:
    status = str(row.get("status", "") or "ok")
    if status not in {"ok", "reference"}:
        return False
    if not row.get("smiles"):
        return False
    if safe_float(row.get("affinity_kcal_mol")) is None and safe_float(row.get("best_affinity_kcal_mol")) is None:
        return False
    pose_path_value = str(row.get("pose_path", "") or "")
    pose_path = host_path(pose_path_value) if pose_path_value else None
    return bool(pose_path and pose_path.exists())


def enrich_screening_rows_with_pose_analysis(
    rows: list[dict[str, str | float | int]],
    analysis_rows: list[dict[str, str | float | int]],
) -> list[dict[str, str | float | int]]:
    by_candidate: dict[str, dict[str, str | float | int]] = {
        str(row.get("candidate_id", "")): row for row in analysis_rows if row.get("candidate_id")
    }
    enriched_rows: list[dict[str, str | float | int]] = []
    for row in rows:
        if not is_successful_pose_screening_row(row):
            continue
        analysis = by_candidate.get(str(row.get("candidate_id", "")))
        if analysis is None:
            continue
        updated = dict(row)
        updated["status"] = "ok"
        updated["live_min_pose_rmsd_to_lead"] = analysis.get("min_pose_rmsd_to_lead", "")
        updated["live_near_lead_mode_count"] = analysis.get("near_lead_mode_count", "")
        updated["live_mode_count"] = analysis.get("mode_count", "")
        updated["live_inner_rmsd"] = analysis.get("best_affinity_mode_mean_rmsd", "")
        updated["live_inner_cluster_fraction"] = analysis.get("best_affinity_mode_cluster_fraction", "")
        updated["live_whole_rmsd"] = analysis.get("whole_rmsd", "")
        enriched_rows.append(updated)
    return enriched_rows


def load_parent_smiles_file(path: Path) -> list[tuple[str, str]]:
    parents: list[tuple[str, str]] = []
    if not path.exists():
        raise FileNotFoundError(f"parent smiles file not found: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        handle.seek(0)
        if "," in first and "smiles" in first.lower():
            for index, row in enumerate(csv.DictReader(handle), start=1):
                smiles = row.get("smiles") or row.get("mol_smiles") or ""
                can = canonical(smiles) if smiles else None
                if can:
                    parents.append((row.get("candidate_id") or row.get("parent_id") or f"P{index:04d}", can))
        else:
            for index, line in enumerate(handle, start=1):
                parts = line.strip().split()
                if not parts:
                    continue
                can = canonical(parts[0])
                if can:
                    parents.append((parts[1] if len(parts) > 1 else f"P{index:04d}", can))
    return parents


def resolve_active_parent_file(args: argparse.Namespace, outdir: Path, diagnostics: list[str]) -> Path | None:
    if args.parent_smiles_file:
        parent_file = args.parent_smiles_file.expanduser().resolve()
        if parent_file.exists():
            diagnostics.append(f"parent_smiles_file={parent_file}")
            return parent_file
        diagnostics.append(f"parent_smiles_file_missing={parent_file}; fallback_to_default_or_lead=true")
        print(f"parent_smiles_file_missing={parent_file}; falling back to default parents or LEAD", flush=True)
    if not args.no_outdir_history:
        default_parent_file = outdir / "next_round_parents.smi"
        if default_parent_file.exists():
            diagnostics.append(f"default_next_round_parents_loaded={default_parent_file}")
            return default_parent_file
    return None


def parent_source_label(parent_file: Path | None) -> str:
    if parent_file is None:
        return "parent_source=LEAD"
    if parent_file.name == "next_round_parents.smi":
        return f"default_next_round_parents_loaded={parent_file}"
    return f"parent_smiles_file={parent_file}"


def ensure_lead_parent(
    parents: list[tuple[str, str]],
    lead_smiles: str,
    diagnostics: list[str] | None = None,
    *,
    label: str = "LEAD_FORCED",
) -> list[tuple[str, str]]:
    lead_can = canonical(lead_smiles) or lead_smiles
    kept: list[tuple[str, str]] = [(parent_id, parent_smiles) for parent_id, parent_smiles in parents if (canonical(parent_smiles) or parent_smiles) != lead_can]
    if diagnostics is not None and len(kept) == len(parents):
        diagnostics.append("lead_smiles_forced_into_active_parents=true")
    return [(label, lead_can), *kept]


def parse_config(
    path: Path,
    cpu_override: int | None,
    exhaustiveness_override: int | None,
    center_x_override: float | None = None,
    center_y_override: float | None = None,
    center_z_override: float | None = None,
    size_x_override: float | None = None,
    size_y_override: float | None = None,
    size_z_override: float | None = None,
    energy_range_override: int | None = None,
) -> ConfigBox:
    values: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    required = ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]
    missing = [key for key in required if key not in values]
    if missing:
        raise ValueError(f"config missing required Vina box fields: {', '.join(missing)}")
    return ConfigBox(
        center_x=center_x_override if center_x_override is not None else float(values["center_x"]),
        center_y=center_y_override if center_y_override is not None else float(values["center_y"]),
        center_z=center_z_override if center_z_override is not None else float(values["center_z"]),
        size_x=size_x_override if size_x_override is not None else float(values["size_x"]),
        size_y=size_y_override if size_y_override is not None else float(values["size_y"]),
        size_z=size_z_override if size_z_override is not None else float(values["size_z"]),
        exhaustiveness=exhaustiveness_override or int(values.get("exhaustiveness", "8")),
        cpu=cpu_override or int(values.get("cpu", "1")),
        energy_range=energy_range_override or int(values.get("energy_range", "3")),
    )


def write_effective_config(path: Path, receptor: Path, box: ConfigBox) -> None:
    lines = [
        f"receptor = {receptor}",
        "",
        f"center_x = {box.center_x}",
        f"center_y = {box.center_y}",
        f"center_z = {box.center_z}",
        "",
        f"size_x = {box.size_x}",
        f"size_y = {box.size_y}",
        f"size_z = {box.size_z}",
        "",
        f"cpu = {box.cpu}",
        f"exhaustiveness = {box.exhaustiveness}",
        f"energy_range = {box.energy_range}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def properties(smiles: str) -> dict[str, float | int]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    return {
        "mw": round(float(Descriptors.MolWt(mol)), 3),
        "logp": round(float(Crippen.MolLogP(mol)), 3),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(mol)), 3),
        "rot_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "formal_charge": int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms())),
    }


def passes_filters(props: dict[str, float | int], parent_props: dict[str, float | int]) -> tuple[bool, str]:
    max_heavy_atoms = max(45, int(parent_props["heavy_atoms"]) + 4)
    max_mw = max(650.0, float(parent_props["mw"]) + 120.0)
    max_tpsa = max(150.0, float(parent_props["tpsa"]) + 40.0)
    if not 12 <= int(props["heavy_atoms"]) <= max_heavy_atoms:
        return False, "heavy_atom_limit"
    if not 180 <= float(props["mw"]) <= max_mw:
        return False, "mw_limit"
    if not -1.0 <= float(props["logp"]) <= 6.0:
        return False, "logp_limit"
    max_hbd = max(5, int(parent_props["hbd"]) + 1)
    max_hba = max(10, int(parent_props["hba"]) + 2)
    if int(props["hbd"]) > max_hbd or int(props["hba"]) > max_hba:
        return False, "hbond_count_limit"
    if float(props["tpsa"]) > max_tpsa:
        return False, "tpsa_limit"
    if int(props["rot_bonds"]) > min(10, int(parent_props["rot_bonds"]) + 2):
        return False, "rotatable_bond_limit"
    if abs(int(props["formal_charge"])) > 1:
        return False, "charge_limit"
    if int(props["heavy_atoms"]) > int(parent_props["heavy_atoms"]) + 4:
        return False, "too_much_growth"
    return True, ""


def aromatic_single_substitutions(parent: str) -> list[tuple[str, str]]:
    mol = Chem.MolFromSmiles(parent)
    if mol is None:
        return []
    replacements = [
        ("F", "F"),
        ("Cl", "Cl"),
        ("Br", "Br"),
        ("Me", "C"),
        ("Et", "CC"),
        ("OH", "O"),
        ("NH2", "N"),
        ("OMe", "OC"),
        ("SMe", "SC"),
        ("Vinyl", "C=C"),
        ("Ethynyl", "C#C"),
        ("CN", "C#N"),
        ("CHO", "C=O"),
        ("COMe", "C(C)=O"),
        ("CONH2", "C(N)=O"),
    ]
    analogs: list[tuple[str, str]] = []
    rw_parent = Chem.RWMol(mol)
    for atom in mol.GetAtoms():
        if not atom.GetIsAromatic() or atom.GetAtomicNum() != 6 or atom.GetTotalNumHs() < 1:
            continue
        atom_idx = atom.GetIdx()
        for label, frag_smiles in replacements:
            rw = Chem.RWMol(rw_parent)
            frag = Chem.MolFromSmiles(frag_smiles)
            if frag is None:
                continue
            combo = Chem.CombineMols(rw.GetMol(), frag)
            edit = Chem.RWMol(combo)
            new_atom_idx = mol.GetNumAtoms()
            edit.AddBond(atom_idx, new_atom_idx, Chem.BondType.SINGLE)
            smiles = canonical(Chem.MolToSmiles(edit.GetMol()))
            if smiles:
                analogs.append((smiles, f"aromatic_C{atom_idx}_{label}"))
    return analogs


def aromatic_n_swaps(parent: str) -> list[tuple[str, str]]:
    mol = Chem.MolFromSmiles(parent)
    if mol is None:
        return []
    analogs: list[tuple[str, str]] = []
    for atom in mol.GetAtoms():
        if not atom.GetIsAromatic() or atom.GetAtomicNum() != 6 or atom.GetTotalNumHs() < 1:
            continue
        atom_idx = atom.GetIdx()
        edit = Chem.RWMol(mol)
        edit_atom = edit.GetAtomWithIdx(atom_idx)
        edit_atom.SetAtomicNum(7)
        edit_atom.SetNoImplicit(True)
        edit_atom.SetNumExplicitHs(0)
        smiles = canonical(Chem.MolToSmiles(edit.GetMol()))
        if smiles:
            analogs.append((smiles, f"aromatic_C{atom_idx}_Nswap"))
    return analogs


def terminal_hetero_deletions(parent: str) -> list[tuple[str, str]]:
    mol = Chem.MolFromSmiles(parent)
    if mol is None:
        return []
    analogs: list[tuple[str, str]] = []
    for atom in mol.GetAtoms():
        if atom.GetDegree() != 1 or atom.GetAtomicNum() not in {9, 17, 35}:
            continue
        atom_idx = atom.GetIdx()
        label = atom.GetSymbol()
        edit = Chem.RWMol(mol)
        edit.RemoveAtom(atom_idx)
        smiles = canonical(Chem.MolToSmiles(edit.GetMol()))
        if smiles:
            analogs.append((smiles, f"delete_terminal_{atom_idx}_{label}"))
    return analogs


def terminal_aliphatic_substitutions(parent: str) -> list[tuple[str, str]]:
    mol = Chem.MolFromSmiles(parent)
    if mol is None:
        return []
    replacements = [
        ("OH", "O"),
        ("NH2", "N"),
        ("F", "F"),
    ]
    analogs: list[tuple[str, str]] = []
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic() or atom.GetAtomicNum() != 6 or atom.GetTotalNumHs() < 1:
            continue
        if sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1) > 1:
            continue
        atom_idx = atom.GetIdx()
        for label, frag_smiles in replacements:
            frag = Chem.MolFromSmiles(frag_smiles)
            if frag is None:
                continue
            combo = Chem.CombineMols(mol, frag)
            edit = Chem.RWMol(combo)
            edit.AddBond(atom_idx, mol.GetNumAtoms(), Chem.BondType.SINGLE)
            smiles = canonical(Chem.MolToSmiles(edit.GetMol()))
            if smiles:
                analogs.append((smiles, f"terminal_C{atom_idx}_{label}"))
    return analogs


def connected_fragment_without_atom(mol: Chem.Mol, start_idx: int, blocked_idx: int) -> set[int]:
    visited: set[int] = set()
    stack = [start_idx]
    while stack:
        idx = stack.pop()
        if idx in visited or idx == blocked_idx:
            continue
        visited.add(idx)
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in visited and nbr_idx != blocked_idx:
                stack.append(nbr_idx)
    return visited


def remove_atoms(parent_mol: Chem.Mol, atom_indices: set[int]) -> str | None:
    if not atom_indices or len(atom_indices) >= parent_mol.GetNumAtoms():
        return None
    edit = Chem.RWMol(parent_mol)
    for idx in sorted(atom_indices, reverse=True):
        edit.RemoveAtom(idx)
    return canonical(Chem.MolToSmiles(edit.GetMol()))


def heavy_atom_count_in_indices(mol: Chem.Mol, atom_indices: set[int]) -> int:
    return sum(1 for idx in atom_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1)


def replace_removed_atoms_with_linker(
    mol: Chem.Mol,
    remove_indices: set[int],
    left_idx: int,
    right_idx: int,
    linker_symbols: list[str],
    label: str,
) -> tuple[str, str] | None:
    if left_idx in remove_indices or right_idx in remove_indices or left_idx == right_idx:
        return None
    old_to_new: dict[int, int] = {}
    edit = Chem.RWMol()
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in remove_indices:
            continue
        new_atom = Chem.Atom(atom.GetAtomicNum())
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        new_atom.SetIsAromatic(atom.GetIsAromatic())
        old_to_new[idx] = edit.AddAtom(new_atom)
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin in remove_indices or end in remove_indices:
            continue
        edit.AddBond(old_to_new[begin], old_to_new[end], bond.GetBondType())
    left_new = old_to_new[left_idx]
    right_new = old_to_new[right_idx]
    if not linker_symbols:
        edit.AddBond(left_new, right_new, Chem.BondType.SINGLE)
    else:
        previous = left_new
        for symbol in linker_symbols:
            atom = Chem.Atom(symbol)
            new_idx = edit.AddAtom(atom)
            edit.AddBond(previous, new_idx, Chem.BondType.SINGLE)
            previous = new_idx
        edit.AddBond(previous, right_new, Chem.BondType.SINGLE)
    try:
        candidate = edit.GetMol()
        Chem.SanitizeMol(candidate)
    except Exception:
        return None
    smiles = canonical(Chem.MolToSmiles(candidate))
    return (smiles, label) if smiles else None


def aromatic_substituent_deletions(parent: str, max_removed_heavy_atoms: int = 8) -> list[tuple[str, str]]:
    mol = Chem.MolFromSmiles(parent)
    if mol is None:
        return []
    analogs: list[tuple[str, str]] = []
    for atom in mol.GetAtoms():
        if not atom.GetIsAromatic() or not atom.IsInRing():
            continue
        anchor_idx = atom.GetIdx()
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            bond = mol.GetBondBetweenAtoms(anchor_idx, nbr_idx)
            if nbr.IsInRing() or bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            fragment = connected_fragment_without_atom(mol, nbr_idx, anchor_idx)
            heavy_count = sum(1 for idx in fragment if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1)
            if heavy_count == 0 or heavy_count > max_removed_heavy_atoms:
                continue
            smiles = remove_atoms(mol, fragment)
            if smiles:
                analogs.append((smiles, f"delete_aromatic_C{anchor_idx}_substituent_{heavy_count}ha"))
    return analogs


def shrink_aromatic_substituents(parent: str) -> list[tuple[str, str]]:
    mol = Chem.MolFromSmiles(parent)
    if mol is None:
        return []
    analogs: list[tuple[str, str]] = []
    for atom in mol.GetAtoms():
        if not atom.GetIsAromatic() or not atom.IsInRing():
            continue
        anchor_idx = atom.GetIdx()
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            bond = mol.GetBondBetweenAtoms(anchor_idx, nbr_idx)
            if nbr.IsInRing() or bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                tail: set[int] = set()
                for tail_atom in nbr.GetNeighbors():
                    tail_idx = tail_atom.GetIdx()
                    if tail_idx == anchor_idx:
                        continue
                    tail |= connected_fragment_without_atom(mol, tail_idx, nbr_idx)
                if tail and heavy_atom_count_in_indices(mol, tail) <= 4:
                    smiles = remove_atoms(mol, tail)
                    if smiles:
                        analogs.append((smiles, f"shrink_aromatic_C{anchor_idx}_alkyl_to_Me"))
            if nbr.GetAtomicNum() in {8, 16}:
                terminal_carbons = [
                    x.GetIdx()
                    for x in nbr.GetNeighbors()
                    if x.GetIdx() != anchor_idx and x.GetAtomicNum() == 6 and x.GetDegree() == 1
                ]
                for terminal_idx in terminal_carbons:
                    smiles = remove_atoms(mol, {terminal_idx})
                    if smiles:
                        label = "alkoxy_to_OH" if nbr.GetAtomicNum() == 8 else "thioether_to_SH"
                        analogs.append((smiles, f"shrink_aromatic_C{anchor_idx}_{label}"))
                    edit = Chem.RWMol(mol)
                    edit.GetAtomWithIdx(nbr_idx).SetAtomicNum(6)
                    smiles = remove_atoms(edit.GetMol(), {terminal_idx})
                    if smiles:
                        analogs.append((smiles, f"shrink_aromatic_C{anchor_idx}_heteroMe_to_Me"))
            if nbr.GetAtomicNum() == 6:
                has_carbonyl_o = any(
                    mol.GetBondBetweenAtoms(nbr_idx, x.GetIdx()) is not None
                    and mol.GetBondBetweenAtoms(nbr_idx, x.GetIdx()).GetBondType() == Chem.BondType.DOUBLE
                    and x.GetAtomicNum() == 8
                    for x in nbr.GetNeighbors()
                )
                if has_carbonyl_o:
                    terminal_carbons = [
                        x.GetIdx()
                        for x in nbr.GetNeighbors()
                        if x.GetIdx() != anchor_idx and x.GetAtomicNum() == 6 and x.GetDegree() == 1
                    ]
                    for terminal_idx in terminal_carbons:
                        smiles = remove_atoms(mol, {terminal_idx})
                        if smiles:
                            analogs.append((smiles, f"shrink_aromatic_C{anchor_idx}_acyl_to_CHO"))
    return analogs


def amide_linker_replacements(parent: str) -> list[tuple[str, str]]:
    mol = Chem.MolFromSmiles(parent)
    if mol is None:
        return []
    analogs: list[tuple[str, str]] = []
    seen: set[str] = set()
    for carbonyl in mol.GetAtoms():
        if carbonyl.GetAtomicNum() != 6:
            continue
        carbonyl_idx = carbonyl.GetIdx()
        oxygen_neighbors = [
            nbr.GetIdx()
            for nbr in carbonyl.GetNeighbors()
            if nbr.GetAtomicNum() == 8
            and (bond := mol.GetBondBetweenAtoms(carbonyl_idx, nbr.GetIdx())) is not None
            and bond.GetBondType() == Chem.BondType.DOUBLE
        ]
        if not oxygen_neighbors:
            continue
        amide_nitrogens = [
            nbr
            for nbr in carbonyl.GetNeighbors()
            if nbr.GetAtomicNum() == 7
            and (bond := mol.GetBondBetweenAtoms(carbonyl_idx, nbr.GetIdx())) is not None
            and bond.GetBondType() == Chem.BondType.SINGLE
        ]
        if not amide_nitrogens:
            continue
        carbonyl_side = [
            nbr.GetIdx()
            for nbr in carbonyl.GetNeighbors()
            if nbr.GetAtomicNum() > 1 and nbr.GetIdx() not in set(oxygen_neighbors + [n.GetIdx() for n in amide_nitrogens])
        ]
        if len(carbonyl_side) != 1:
            continue
        left_idx = carbonyl_side[0]
        left_atom = mol.GetAtomWithIdx(left_idx)
        if left_atom.GetAtomicNum() == 1:
            continue
        for nitrogen in amide_nitrogens:
            nitrogen_idx = nitrogen.GetIdx()
            right_neighbors = [
                nbr.GetIdx()
                for nbr in nitrogen.GetNeighbors()
                if nbr.GetIdx() != carbonyl_idx and nbr.GetAtomicNum() > 1
            ]
            if len(right_neighbors) != 1:
                continue
            right_idx = right_neighbors[0]
            right_atom = mol.GetAtomWithIdx(right_idx)
            if right_atom.GetAtomicNum() == 1:
                continue
            if mol.GetBondBetweenAtoms(left_idx, right_idx) is not None:
                continue
            ketone_result = replace_removed_atoms_with_linker(
                mol,
                {nitrogen_idx},
                carbonyl_idx,
                right_idx,
                [],
                f"amide_to_ketone_CO_C{carbonyl_idx}_N{nitrogen_idx}",
            )
            if ketone_result:
                smiles, label = ketone_result
                if smiles != canonical(parent) and smiles not in seen:
                    seen.add(smiles)
                    analogs.append((smiles, label))
            variants = [
                ("amide_to_direct", {carbonyl_idx, nitrogen_idx, *oxygen_neighbors}, []),
                ("amide_to_NH", {carbonyl_idx, nitrogen_idx, *oxygen_neighbors}, ["N"]),
                ("amide_to_ether_O", {carbonyl_idx, nitrogen_idx, *oxygen_neighbors}, ["O"]),
                ("amide_to_thioether_S", {carbonyl_idx, nitrogen_idx, *oxygen_neighbors}, ["S"]),
                ("amide_to_methylene_CH2", {carbonyl_idx, nitrogen_idx, *oxygen_neighbors}, ["C"]),
            ]
            for variant_label, remove_indices, linker_symbols in variants:
                result = replace_removed_atoms_with_linker(
                    mol,
                    set(remove_indices),
                    left_idx,
                    right_idx,
                    linker_symbols,
                    f"{variant_label}_C{carbonyl_idx}_N{nitrogen_idx}",
                )
                if not result:
                    continue
                smiles, label = result
                if smiles == canonical(parent) or smiles in seen:
                    continue
                seen.add(smiles)
                analogs.append((smiles, label))
    return analogs


def ring_replacement_mol(
    mol: Chem.Mol,
    ring: tuple[int, ...],
    external: list[tuple[int, int, Chem.BondType]],
    atom_symbols: list[str],
    label: str,
) -> tuple[str, str] | None:
    ring_set = set(ring)
    old_to_new: dict[int, int] = {}
    edit = Chem.RWMol()
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in ring_set:
            continue
        old_to_new[idx] = edit.AddAtom(Chem.Atom(atom))
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin in ring_set or end in ring_set:
            continue
        edit.AddBond(old_to_new[begin], old_to_new[end], bond.GetBondType())
    ring_new_indices: list[int] = []
    for symbol in atom_symbols:
        new_atom = Chem.Atom(symbol)
        new_atom.SetIsAromatic(True)
        ring_new_indices.append(edit.AddAtom(new_atom))
    for index, new_idx in enumerate(ring_new_indices):
        edit.AddBond(new_idx, ring_new_indices[(index + 1) % len(ring_new_indices)], Chem.BondType.AROMATIC)
    used_positions: set[int] = set()
    for old_pos, outside_idx, bond_type in external:
        new_pos = round(old_pos * (len(atom_symbols) - 1) / (len(ring) - 1))
        if new_pos in used_positions or atom_symbols[new_pos] != "C":
            return None
        used_positions.add(new_pos)
        edit.AddBond(old_to_new[outside_idx], ring_new_indices[new_pos], bond_type)
    try:
        candidate = edit.GetMol()
        Chem.SanitizeMol(candidate)
    except Exception:
        return None
    smiles = canonical(Chem.MolToSmiles(candidate))
    return (smiles, label) if smiles else None


def aromatic_six_to_five_ring_replacements(parent: str) -> list[tuple[str, str]]:
    mol = Chem.MolFromSmiles(parent)
    if mol is None:
        return []
    templates = [
        ("pyrrole_like_5N", ["N", "C", "C", "C", "C"]),
        ("furan", ["O", "C", "C", "C", "C"]),
        ("furan_alt", ["C", "O", "C", "C", "C"]),
        ("isoxazole", ["O", "N", "C", "C", "C"]),
        ("isoxazole_alt", ["N", "O", "C", "C", "C"]),
        ("oxadiazole_1_2_4", ["O", "N", "C", "N", "C"]),
        ("oxadiazole_1_3_4", ["O", "C", "N", "N", "C"]),
        ("thiophene", ["S", "C", "C", "C", "C"]),
        ("oxazole", ["O", "C", "N", "C", "C"]),
        ("thiazole", ["S", "C", "N", "C", "C"]),
    ]
    analogs: list[tuple[str, str]] = []
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        ring_set = set(ring)
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        if any(
            nbr.IsInRing() and nbr.GetIdx() not in ring_set
            for idx in ring
            for nbr in mol.GetAtomWithIdx(idx).GetNeighbors()
        ):
            continue
        external: list[tuple[int, int, Chem.BondType]] = []
        for pos, idx in enumerate(ring):
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in ring_set:
                    continue
                bond = mol.GetBondBetweenAtoms(idx, nbr_idx)
                if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                    external = []
                    break
                external.append((pos, nbr_idx, bond.GetBondType()))
            if len(external) > 3:
                break
        if not external or len(external) > 3:
            continue
        for template_label, atom_symbols in templates:
            for offset in range(len(atom_symbols)):
                rotated = atom_symbols[offset:] + atom_symbols[:offset]
                result = ring_replacement_mol(
                    mol,
                    ring,
                    external,
                    rotated,
                    f"ring6to5_{template_label}_ring{ring[0]}_rot{offset}",
                )
                if result:
                    analogs.append(result)
    return analogs


def intramolecular_ring_closures(parent: str, max_candidates: int = 80) -> list[tuple[str, str]]:
    mol = Chem.MolFromSmiles(parent)
    if mol is None:
        return []
    ring_info = mol.GetRingInfo()
    atom_rings = [set(ring) for ring in ring_info.AtomRings()]
    ring_atoms = {idx for ring in atom_rings for idx in ring}
    candidate_atoms = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 6
        and atom.GetTotalNumHs() > 0
        and atom.GetIdx() in ring_atoms
        and atom.GetDegree() <= 3
    ]
    analogs: list[tuple[str, str]] = []
    seen: set[str] = set()
    for pos, begin in enumerate(candidate_atoms):
        begin_rings = [ring for ring in atom_rings if begin in ring]
        for end in candidate_atoms[pos + 1 :]:
            if mol.GetBondBetweenAtoms(begin, end) is not None:
                continue
            end_rings = [ring for ring in atom_rings if end in ring]
            if any(ring_a is ring_b for ring_a in begin_rings for ring_b in end_rings):
                continue
            try:
                path_len = len(Chem.rdmolops.GetShortestPath(mol, begin, end)) - 1
            except Exception:
                continue
            if path_len < 3 or path_len > 8:
                continue
            edit = Chem.RWMol(mol)
            edit.AddBond(begin, end, Chem.BondType.SINGLE)
            try:
                candidate = edit.GetMol()
                Chem.SanitizeMol(candidate)
            except Exception:
                continue
            smiles = canonical(Chem.MolToSmiles(candidate))
            if not smiles or smiles == canonical(parent) or smiles in seen:
                continue
            seen.add(smiles)
            analogs.append((smiles, f"ring_closure_C{begin}_C{end}_path{path_len}"))
            if len(analogs) >= max_candidates:
                return analogs
    return analogs


def partial_aromatic_saturations(parent: str, max_candidates: int = 80) -> list[tuple[str, str]]:
    mol = Chem.MolFromSmiles(parent)
    if mol is None:
        return []
    try:
        kek = Chem.Mol(mol)
        Chem.Kekulize(kek, clearAromaticFlags=True)
    except Exception:
        return []
    analogs: list[tuple[str, str]] = []
    seen: set[str] = set()
    for bond in kek.GetBonds():
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        if begin.GetAtomicNum() != 6 or end.GetAtomicNum() != 6:
            continue
        if not begin.IsInRing() or not end.IsInRing():
            continue
        if begin.GetDegree() > 3 or end.GetDegree() > 3:
            continue
        edit = Chem.RWMol(kek)
        edit_bond = edit.GetBondBetweenAtoms(begin.GetIdx(), end.GetIdx())
        if edit_bond is None:
            continue
        edit_bond.SetBondType(Chem.BondType.SINGLE)
        for atom_idx in (begin.GetIdx(), end.GetIdx()):
            atom = edit.GetAtomWithIdx(atom_idx)
            atom.SetIsAromatic(False)
            atom.SetNoImplicit(False)
        try:
            candidate = edit.GetMol()
            Chem.SanitizeMol(candidate)
        except Exception:
            continue
        smiles = canonical(Chem.MolToSmiles(candidate))
        if not smiles or smiles == canonical(parent) or smiles in seen:
            continue
        seen.add(smiles)
        analogs.append((smiles, f"dearomatize_saturate_bond_{begin.GetIdx()}_{end.GetIdx()}"))
        if len(analogs) >= max_candidates:
            return analogs
    return analogs


def atom_deletions(parent: str) -> list[tuple[str, str]]:
    mol = Chem.MolFromSmiles(parent)
    if mol is None:
        return []
    analogs: list[tuple[str, str]] = []
    for atom in mol.GetAtoms():
        if atom.GetDegree() != 1 or atom.GetAtomicNum() in {1, 6, 7, 8, 9, 17}:
            continue
        edit = Chem.RWMol(mol)
        idx = atom.GetIdx()
        label = atom.GetSymbol()
        edit.RemoveAtom(idx)
        smiles = canonical(Chem.MolToSmiles(edit.GetMol()))
        if smiles:
            analogs.append((smiles, f"delete_terminal_{idx}_{label}"))
    return analogs


EDIT_MODE_GROUPS = {
    "add": (
        aromatic_single_substitutions,
        terminal_aliphatic_substitutions,
    ),
    "delete": (
        aromatic_substituent_deletions,
        terminal_hetero_deletions,
        atom_deletions,
    ),
    "shrink": (
        shrink_aromatic_substituents,
    ),
    "drastic": (
        aromatic_n_swaps,
        amide_linker_replacements,
        aromatic_six_to_five_ring_replacements,
        intramolecular_ring_closures,
        partial_aromatic_saturations,
    ),
}
DEFAULT_EDIT_MODES = ("add", "delete", "shrink", "drastic")


def parse_edit_modes(edit_mode_text: str) -> tuple[str, ...]:
    normalized = [
        item.strip().lower()
        for item in re.split(r"[,/+\s]+", edit_mode_text or "")
        if item.strip()
    ]
    if not normalized or normalized == ["default"]:
        return DEFAULT_EDIT_MODES
    aliases = {
        "adds": "add",
        "addition": "add",
        "deletion": "delete",
        "del": "delete",
        "remove": "delete",
        "reduction": "shrink",
        "shirnk": "shrink",
        "small": "shrink",
        "ring": "drastic",
        "scaffold": "drastic",
    }
    modes: list[str] = []
    for item in normalized:
        mode = aliases.get(item, item)
        if mode not in EDIT_MODE_GROUPS:
            raise ValueError(f"unknown --edit-mode value: {item}; use add, delete, shrink, drastic, or comma combinations")
        if mode not in modes:
            modes.append(mode)
    return tuple(modes)


def selected_edit_mode_text(args: argparse.Namespace) -> str:
    shortcut_modes = [
        mode
        for flag_name, mode in (
            ("add_only", "add"),
            ("delete_only", "delete"),
            ("shrink_only", "shrink"),
            ("drastic_only", "drastic"),
        )
        if getattr(args, flag_name, False)
    ]
    if shortcut_modes:
        return ",".join(shortcut_modes)
    return str(args.edit_mode)


def drastic_two_step_edits(parent: str, max_first: int = 40, max_total: int = 160) -> list[tuple[str, str]]:
    first_steps = (
        shrink_aromatic_substituents(parent)
        + aromatic_substituent_deletions(parent, max_removed_heavy_atoms=5)
        + aromatic_n_swaps(parent)
        + amide_linker_replacements(parent)
        + aromatic_six_to_five_ring_replacements(parent)
        + intramolecular_ring_closures(parent, max_candidates=40)
        + partial_aromatic_saturations(parent, max_candidates=40)
    )
    analogs: list[tuple[str, str]] = []
    seen: set[str] = set()
    for first_smiles, first_label in first_steps[:max_first]:
        second_steps = (
            aromatic_single_substitutions(first_smiles)
            + terminal_aliphatic_substitutions(first_smiles)
            + shrink_aromatic_substituents(first_smiles)
            + aromatic_substituent_deletions(first_smiles, max_removed_heavy_atoms=4)
            + aromatic_n_swaps(first_smiles)
        )
        for second_smiles, second_label in second_steps:
            can = canonical(second_smiles)
            if not can or can in seen or can == canonical(parent):
                continue
            seen.add(can)
            analogs.append((can, f"drastic2_{first_label}__{second_label}"))
            if len(analogs) >= max_total:
                return analogs
    return analogs


def enumerate_raw_edits(parent: str, edit_modes: tuple[str, ...]) -> list[tuple[str, str]]:
    raw: list[tuple[str, str]] = []
    for mode in edit_modes:
        for generator in EDIT_MODE_GROUPS[mode]:
            raw.extend(generator(parent))
        if mode == "drastic":
            raw.extend(drastic_two_step_edits(parent))
    return raw


def enumerate_round(
    parent: str,
    parent_id: str,
    round_index: int,
    seen: set[str],
    history_seen_smiles: set[str],
    history_seen_edit_labels: set[str],
    limit: int,
    rng: random.Random,
    edit_modes: tuple[str, ...],
    deterministic_batch: bool = False,
) -> tuple[list[Candidate], list[dict[str, str | int | float]]]:
    raw = enumerate_raw_edits(parent, edit_modes)
    if not deterministic_batch:
        rng.shuffle(raw)
    parent_props = properties(parent)
    candidates: list[Candidate] = []
    rejected: list[dict[str, str | int | float]] = []
    for index, (smiles, label) in enumerate(raw, start=1):
        can = canonical(smiles)
        if can is None:
            rejected.append(
                {
                    "round": round_index,
                    "parent_id": parent_id,
                    "edit_label": label,
                    "input_smiles": smiles,
                    "canonical_smiles": "",
                    "reason": "invalid_after_edit",
                }
            )
            continue
        if can in seen:
            rejected.append(
                {
                    "round": round_index,
                    "parent_id": parent_id,
                    "edit_label": label,
                    "input_smiles": smiles,
                    "canonical_smiles": can,
                    "reason": "duplicate",
                }
            )
            continue
        if can in history_seen_smiles:
            rejected.append(
                {
                    "round": round_index,
                    "parent_id": parent_id,
                    "edit_label": label,
                    "input_smiles": smiles,
                    "canonical_smiles": can,
                    "reason": "history_duplicate_smiles",
                }
            )
            continue
        props = properties(can)
        ok, reason = passes_filters(props, parent_props)
        if not ok:
            rejected.append(
                {
                    "round": round_index,
                    "parent_id": parent_id,
                    "edit_label": label,
                    "input_smiles": smiles,
                    "canonical_smiles": can,
                    "reason": reason,
                    **props,
                }
            )
            continue
        seen.add(can)
        candidates.append(Candidate(f"R{round_index:02d}_{len(candidates) + 1:04d}", parent_id, can, label, round_index))
        if len(candidates) >= limit:
            break
    return candidates, rejected


def embed_sdf(smiles: str, name: str, sdf_path: Path, seed: int) -> None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    status = AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=True, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    if status != 0:
        raise RuntimeError(f"RDKit EmbedMolecule failed: {status}")
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=300)
    mol.SetProp("_Name", name)
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol)
    writer.close()


def meeko_commands(meeko: str, sdf_path: Path, pdbqt_path: Path) -> list[tuple[str, list[str]]]:
    commands: list[tuple[str, list[str]]] = []
    if command_exists(meeko):
        commands.append(("meeko_script", [meeko, "-i", str(sdf_path), "-o", str(pdbqt_path)]))
    commands.append(("meeko_module", [sys.executable, "-m", "meeko.cli.mk_prepare_ligand", "-i", str(sdf_path), "-o", str(pdbqt_path)]))
    commands.append(("meeko_module_alt", [sys.executable, "-m", "meeko.cli.prepare_ligand", "-i", str(sdf_path), "-o", str(pdbqt_path)]))
    return commands


def prepare_pdbqt(sdf_path: Path, pdbqt_path: Path, meeko: str, obabel: str) -> tuple[bool, str, str]:
    logs: list[str] = []
    for mode, cmd in meeko_commands(meeko, sdf_path, pdbqt_path):
        if pdbqt_path.exists():
            pdbqt_path.unlink()
        code, output = run(cmd)
        logs.append(f"## {mode}\n+ {' '.join(cmd)}\nexit_code={code}\n{output}")
        if code == 0 and pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
            return True, "\n".join(logs), mode

    fallback = [obabel, "-isdf", str(sdf_path), "-opdbqt", "-O", str(pdbqt_path), "--partialcharge", "gasteiger"]
    if pdbqt_path.exists():
        pdbqt_path.unlink()
    code2, output2 = run(fallback)
    logs.append(f"## obabel\n+ {' '.join(fallback)}\nexit_code={code2}\n{output2}")
    if code2 == 0 and pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
        return True, "\n".join(logs), "obabel"
    return False, "\n".join(logs), "failed"


def dock(vina: str, receptor: Path, ligand: Path, out_pose: Path, log_path: Path, box: ConfigBox, seed: int, num_modes: int = 9) -> tuple[bool, float | None, str]:
    cmd = [
        vina,
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
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
        "--seed",
        str(seed),
        "--num_modes",
        str(num_modes),
        "--energy_range",
        str(box.energy_range),
        "--out",
        str(out_pose),
    ]
    code, output = run(cmd, cwd=out_pose.parent)
    log_path.write_text(output, encoding="utf-8", errors="replace")
    match = AFFINITY_LOG_RE.search(output)
    affinity = float(match.group(1)) if match else None
    return code == 0 and affinity is not None, affinity, output


def read_best_affinity(path: Path) -> float:
    text = path.read_text(encoding="utf-8", errors="replace")
    pdbqt_match = AFFINITY_PDBQT_RE.search(text)
    if pdbqt_match:
        return float(pdbqt_match.group(1))
    log_match = AFFINITY_LOG_RE.search(text)
    if log_match:
        return float(log_match.group(1))
    raise ValueError(f"could not parse affinity from {path}")


def is_heavy_pdbqt_atom(line: str) -> bool:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return False
    atom_name = line[12:16].strip()
    atom_type = line[77:].strip().split()[0] if len(line) > 77 and line[77:].strip() else atom_name
    element_hint = re.sub(r"[^A-Za-z]", "", atom_type or atom_name).upper()
    return bool(element_hint) and not element_hint.startswith("H")


def pdbqt_atom_record(line: str) -> dict[str, str | float] | None:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None
    atom_name = line[12:16].strip()
    residue_name = line[17:20].strip()
    residue_id = line[22:26].strip()
    atom_type = line[77:].strip().split()[0] if len(line) > 77 and line[77:].strip() else atom_name
    raw_element = re.sub(r"[^A-Za-z]", "", atom_type or atom_name).upper()
    element = PDBQT_ELEMENT_FALLBACKS.get(raw_element, PDBQT_ELEMENT_FALLBACKS.get(raw_element[:1], raw_element[:1].title()))
    return {
        "x": x,
        "y": y,
        "z": z,
        "element": element,
        "atom_name": atom_name,
        "residue_name": residue_name,
        "residue_id": residue_id,
        "atom_type": atom_type,
    }


def parse_pdbqt_modes(path: Path) -> list[PdbqtMode]:
    modes: list[PdbqtMode] = []
    current: list[tuple[float, float, float]] | None = None
    current_atoms: list[dict[str, str | float]] | None = None
    current_index = 1
    current_affinity: float | None = None
    saw_model = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            if current is not None:
                modes.append(PdbqtMode(current_index, current_affinity, current, current_atoms or []))
            current = []
            current_atoms = []
            parts = line.split()
            current_index = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else len(modes) + 1
            current_affinity = None
            saw_model = True
            continue
        affinity_match = AFFINITY_PDBQT_RE.search(line)
        if affinity_match:
            current_affinity = float(affinity_match.group(1))
            continue
        if line.startswith("ENDMDL"):
            if current is not None:
                modes.append(PdbqtMode(current_index, current_affinity, current, current_atoms or []))
                current = None
                current_atoms = None
            continue
        atom_record = pdbqt_atom_record(line)
        if atom_record is not None and is_heavy_pdbqt_atom(line):
            target = current
            if target is None and not saw_model:
                current = []
                current_atoms = []
                target = current
            if target is not None:
                target.append((float(atom_record["x"]), float(atom_record["y"]), float(atom_record["z"])))
                if current_atoms is not None:
                    current_atoms.append(atom_record)
    if current is not None:
        modes.append(PdbqtMode(current_index, current_affinity, current, current_atoms or []))
    return [mode for mode in modes if mode.coords]


def coordinate_rmsd(reference: list[tuple[float, float, float]], pose: list[tuple[float, float, float]]) -> tuple[float, int]:
    matched = min(len(reference), len(pose))
    if matched == 0:
        return float("nan"), 0
    total = 0.0
    for ref_xyz, pose_xyz in zip(reference[:matched], pose[:matched]):
        total += sum((a - b) ** 2 for a, b in zip(ref_xyz, pose_xyz))
    return math.sqrt(total / matched), matched


def squared_distance(atom_a: dict[str, str | float], atom_b: dict[str, str | float]) -> float:
    return (
        (float(atom_a["x"]) - float(atom_b["x"])) ** 2
        + (float(atom_a["y"]) - float(atom_b["y"])) ** 2
        + (float(atom_a["z"]) - float(atom_b["z"])) ** 2
    )


def parse_receptor_atoms(path: Path) -> list[dict[str, str | float]]:
    atoms: list[dict[str, str | float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        atom = pdbqt_atom_record(line)
        if atom is not None and atom.get("element") != "H":
            atoms.append(atom)
    return atoms


def interaction_summary(
    ligand_atoms: list[dict[str, str | float]],
    receptor_atoms: list[dict[str, str | float]],
) -> dict[str, int]:
    hbonds: set[tuple[int, int]] = set()
    hydrophobic: set[tuple[int, int]] = set()
    vdw: set[tuple[int, int]] = set()
    pi_contacts: set[tuple[int, str]] = set()
    ch_pi_contacts: set[tuple[int, str]] = set()
    for ligand_index, ligand_atom in enumerate(ligand_atoms):
        ligand_element = str(ligand_atom.get("element", ""))
        ligand_type = str(ligand_atom.get("atom_type", ""))
        ligand_name = str(ligand_atom.get("atom_name", ""))
        for receptor_index, receptor_atom in enumerate(receptor_atoms):
            receptor_element = str(receptor_atom.get("element", ""))
            distance2 = squared_distance(ligand_atom, receptor_atom)
            if ligand_element in HBOND_ELEMENTS and receptor_element in HBOND_ELEMENTS and distance2 <= 3.5**2:
                hbonds.add((ligand_index, receptor_index))
            if ligand_element in HYDROPHOBIC_ELEMENTS and receptor_element in HYDROPHOBIC_ELEMENTS and distance2 <= 4.5**2:
                hydrophobic.add((ligand_index, receptor_index))
            if distance2 <= 4.0**2:
                vdw.add((ligand_index, receptor_index))
            residue_name = str(receptor_atom.get("residue_name", ""))
            if ligand_element == "C" and residue_name in AROMATIC_RESIDUES and distance2 <= 5.0**2:
                pi_contacts.add((ligand_index, f"{residue_name}{receptor_atom.get('residue_id', '')}"))
            if (
                ligand_element == "C"
                and residue_name in AROMATIC_RESIDUES
                and distance2 <= 4.2**2
                and not ligand_type.startswith("A")
                and not ligand_name.upper().startswith("A")
            ):
                ch_pi_contacts.add((ligand_index, f"{residue_name}{receptor_atom.get('residue_id', '')}"))
    return {
        "hbond_count": len(hbonds),
        "hydrophobic_count": len(hydrophobic),
        "vdw_contact_count": len(vdw),
        "pi_contact_count": len(pi_contacts),
        "ch_pi_count": len(ch_pi_contacts),
    }


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[int(position)]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def normal_scores(values: list[float | None], *, lower_is_better: bool) -> list[float | None]:
    finite = [value for value in values if value is not None and math.isfinite(value)]
    if not finite:
        return [None for _ in values]
    mean = sum(finite) / len(finite)
    variance = sum((value - mean) ** 2 for value in finite) / len(finite)
    stdev = math.sqrt(variance)
    if stdev < 1e-9:
        return [0.5 if value is not None and math.isfinite(value) else None for value in values]
    scores: list[float | None] = []
    for value in values:
        if value is None or not math.isfinite(value):
            scores.append(None)
            continue
        z = (value - mean) / stdev
        score = 1.0 - normal_cdf(z) if lower_is_better else normal_cdf(z)
        scores.append(round(score, 6))
    return scores


def best_affinity_mode_summary(
    modes: list[PdbqtMode],
    cluster_rmsd_cutoff: float,
) -> dict[str, str | float | int]:
    if not modes:
        return {
            "best_affinity_mode": "",
            "best_affinity_mode_cluster_size": 0,
            "best_affinity_mode_cluster_fraction": 0.0,
            "best_affinity_mode_mean_rmsd": "",
        }
    best_mode = min(
        modes,
        key=lambda mode: mode.affinity if mode.affinity is not None else float("inf"),
    )
    rmsds: list[float] = []
    near_count = 0
    for mode in modes:
        rmsd, _ = coordinate_rmsd(best_mode.coords, mode.coords)
        if not math.isfinite(rmsd):
            continue
        rmsds.append(rmsd)
        if rmsd <= cluster_rmsd_cutoff:
            near_count += 1
    mean_rmsd = sum(rmsds) / len(rmsds) if rmsds else float("nan")
    return {
        "best_affinity_mode": best_mode.mode_index,
        "best_affinity_mode_cluster_size": near_count,
        "best_affinity_mode_cluster_fraction": round(near_count / len(modes), 4),
        "best_affinity_mode_mean_rmsd": round(mean_rmsd, 4) if math.isfinite(mean_rmsd) else "",
    }


def whole_pose_rmsd_summary(modes: list[PdbqtMode]) -> dict[str, str | float]:
    if len(modes) == 1:
        return {"whole_rmsd": 0.0}
    pairwise: list[float] = []
    for index, mode in enumerate(modes):
        for other in modes[index + 1 :]:
            rmsd, _ = coordinate_rmsd(mode.coords, other.coords)
            if math.isfinite(rmsd):
                pairwise.append(rmsd)
    mean_rmsd = sum(pairwise) / len(pairwise) if pairwise else float("nan")
    return {"whole_rmsd": round(mean_rmsd, 4) if math.isfinite(mean_rmsd) else ""}


def internal_pose_cluster_summary(
    modes: list[PdbqtMode],
    cluster_rmsd_cutoff: float,
) -> dict[str, str | float | int]:
    if not modes:
        return {
            "best_internal_cluster_size": 0,
            "best_internal_cluster_fraction": 0.0,
            "best_internal_cluster_representative_mode": "",
            "best_internal_cluster_mean_rmsd": "",
        }
    best_mode = modes[0].mode_index
    best_neighbors: list[float] = []
    for mode in modes:
        neighbor_rmsds: list[float] = []
        for other in modes:
            rmsd, _ = coordinate_rmsd(mode.coords, other.coords)
            if math.isfinite(rmsd) and rmsd <= cluster_rmsd_cutoff:
                neighbor_rmsds.append(rmsd)
        if len(neighbor_rmsds) > len(best_neighbors):
            best_mode = mode.mode_index
            best_neighbors = neighbor_rmsds
        elif len(neighbor_rmsds) == len(best_neighbors) and neighbor_rmsds:
            current_mean = sum(neighbor_rmsds) / len(neighbor_rmsds)
            best_mean = sum(best_neighbors) / len(best_neighbors) if best_neighbors else float("inf")
            if current_mean < best_mean:
                best_mode = mode.mode_index
                best_neighbors = neighbor_rmsds
    mean_rmsd = sum(best_neighbors) / len(best_neighbors) if best_neighbors else float("nan")
    return {
        "best_internal_cluster_size": len(best_neighbors),
        "best_internal_cluster_fraction": round(len(best_neighbors) / len(modes), 4),
        "best_internal_cluster_representative_mode": best_mode,
        "best_internal_cluster_mean_rmsd": round(mean_rmsd, 4) if math.isfinite(mean_rmsd) else "",
    }


def analyze_pose_rmsd(
    lead_pose: Path,
    receptor: Path,
    rows: list[dict[str, str | float | int]],
    lead_smiles: str,
    quantile: float,
    near_rmsd_cutoff: float,
    internal_cluster_rmsd_cutoff: float,
    affinity_threshold: float,
    affinity_penalty_scale: float,
) -> tuple[
    list[dict[str, str | float | int]],
    list[dict[str, str | float | int]],
    list[dict[str, str | float | int]],
    list[dict[str, str | float | int]],
]:
    lead_modes = parse_pdbqt_modes(lead_pose)
    if not lead_modes:
        raise ValueError(f"no heavy-atom pose modes found in lead pose: {lead_pose}")
    lead_reference = lead_modes[0].coords
    receptor_atoms = parse_receptor_atoms(receptor)
    parent_reference_by_id: dict[str, list[tuple[float, float, float]]] = {}
    for history_row in rows:
        candidate_id = str(history_row.get("candidate_id", ""))
        pose_path_value = str(history_row.get("pose_path", ""))
        if not candidate_id or not pose_path_value:
            continue
        pose_path = host_path(pose_path_value)
        if not pose_path or not pose_path.exists():
            continue
        modes = parse_pdbqt_modes(pose_path)
        if not modes:
            continue
        best_mode = min(modes, key=lambda mode: mode.affinity if mode.affinity is not None else float("inf"))
        parent_reference_by_id.setdefault(candidate_id, best_mode.coords)
    analysis_rows: list[dict[str, str | float | int]] = []
    mode_rows: list[dict[str, str | float | int]] = []
    lead_rmsd_values: list[tuple[int, float, int]] = []
    for mode in lead_modes:
        rmsd, matched = coordinate_rmsd(lead_reference, mode.coords)
        if math.isfinite(rmsd):
            lead_rmsd_values.append((mode.mode_index, rmsd, matched))
            mode_rows.append(
                {
                    "candidate_id": "LEAD_REFERENCE",
                    "mode_index": mode.mode_index,
                    "mode_affinity_kcal_mol": mode.affinity if mode.affinity is not None else "",
                    "mode_rmsd_to_lead": round(rmsd, 4),
                    "within_rmsd_cutoff": "true" if rmsd <= near_rmsd_cutoff else "false",
                    "matched_heavy_atoms": matched,
                    "mode_heavy_atoms": len(mode.coords),
                    "lead_heavy_atoms": len(lead_reference),
                    "pose_path": str(lead_pose),
                    "smiles": lead_smiles,
                }
            )
    lead_near_count = sum(1 for _, rmsd, _ in lead_rmsd_values if rmsd <= near_rmsd_cutoff)
    lead_internal_cluster = internal_pose_cluster_summary(lead_modes, internal_cluster_rmsd_cutoff)
    lead_best_mode_summary = best_affinity_mode_summary(lead_modes, internal_cluster_rmsd_cutoff)
    lead_whole_rmsd = whole_pose_rmsd_summary(lead_modes)
    lead_interactions = interaction_summary(lead_modes[0].atoms, receptor_atoms)
    analysis_rows.append(
        {
            "candidate_id": "LEAD_REFERENCE",
            "status": "reference",
            "affinity_kcal_mol": lead_modes[0].affinity if lead_modes[0].affinity is not None else "",
            "delta_vs_baseline": 0.0,
            "pose_path": str(lead_pose),
            "mode_count": len(lead_modes),
            "near_lead_mode_count": lead_near_count,
            "near_lead_mode_fraction": round(lead_near_count / len(lead_modes), 4) if lead_modes else 0.0,
            **lead_internal_cluster,
            **lead_best_mode_summary,
            **lead_whole_rmsd,
            **lead_interactions,
            "matched_heavy_atoms": len(lead_reference),
            "lead_heavy_atoms": len(lead_reference),
            "candidate_mode_heavy_atoms_min": min(len(mode.coords) for mode in lead_modes),
            "min_pose_rmsd_to_lead": 0.0,
            "best_rmsd_mode": 1,
            "rmsd_q25_cutoff": "",
            "rmsd_q25_keep": "true",
            "smiles": lead_smiles,
        }
    )
    for row in rows:
        if row.get("status") != "ok":
            continue
        pose_path_value = str(row.get("pose_path", ""))
        if not pose_path_value:
            continue
        pose_path = host_path(pose_path_value)
        if not pose_path or not pose_path.exists():
            continue
        modes = parse_pdbqt_modes(pose_path)
        if not modes:
            continue
        rmsd_values: list[tuple[int, float, int]] = []
        for mode_index, mode in enumerate(modes, start=1):
            rmsd, matched = coordinate_rmsd(lead_reference, mode.coords)
            if math.isfinite(rmsd):
                rmsd_values.append((mode.mode_index, rmsd, matched))
                mode_rows.append(
                    {
                        "candidate_id": str(row.get("candidate_id", "")),
                        "mode_index": mode.mode_index,
                        "mode_affinity_kcal_mol": mode.affinity if mode.affinity is not None else "",
                        "mode_rmsd_to_lead": round(rmsd, 4),
                        "within_rmsd_cutoff": "true" if rmsd <= near_rmsd_cutoff else "false",
                        "matched_heavy_atoms": matched,
                        "mode_heavy_atoms": len(mode.coords),
                        "lead_heavy_atoms": len(lead_reference),
                        "pose_path": pose_path_value,
                        "smiles": str(row.get("smiles", "")),
                    }
                )
        if not rmsd_values:
            continue
        best_mode, best_rmsd, matched_atoms = min(rmsd_values, key=lambda item: item[1])
        family_references = [lead_reference]
        parent_reference = parent_reference_by_id.get(str(row.get("parent_id", "")))
        if parent_reference is not None:
            family_references.append(parent_reference)
        family_rmsd_values: list[tuple[int, float, int]] = []
        for mode in modes:
            best_family_rmsd = float("nan")
            best_family_matched = 0
            for reference in family_references:
                rmsd, matched = coordinate_rmsd(reference, mode.coords)
                if math.isfinite(rmsd) and (not math.isfinite(best_family_rmsd) or rmsd < best_family_rmsd):
                    best_family_rmsd = rmsd
                    best_family_matched = matched
            if math.isfinite(best_family_rmsd):
                family_rmsd_values.append((mode.mode_index, best_family_rmsd, best_family_matched))
        _, best_family_rmsd, _ = min(family_rmsd_values, key=lambda item: item[1]) if family_rmsd_values else (best_mode, best_rmsd, matched_atoms)
        candidate_heavy_min = min(len(mode.coords) for mode in modes)
        near_mode_count = sum(1 for _, rmsd, _ in rmsd_values if rmsd <= near_rmsd_cutoff)
        internal_cluster = internal_pose_cluster_summary(modes, internal_cluster_rmsd_cutoff)
        best_mode_cluster = best_affinity_mode_summary(modes, internal_cluster_rmsd_cutoff)
        whole_rmsd = whole_pose_rmsd_summary(modes)
        best_affinity_mode = min(modes, key=lambda mode: mode.affinity if mode.affinity is not None else float("inf"))
        interactions = interaction_summary(best_affinity_mode.atoms, receptor_atoms)
        analysis_rows.append(
            {
                "candidate_id": str(row.get("candidate_id", "")),
                "status": str(row.get("status", "")),
                "affinity_kcal_mol": row.get("affinity_kcal_mol", ""),
                "delta_vs_baseline": row.get("delta_vs_baseline", ""),
                "pose_path": pose_path_value,
                "mode_count": len(modes),
                "near_lead_mode_count": near_mode_count,
                "near_lead_mode_fraction": round(near_mode_count / len(modes), 4) if modes else 0.0,
                **internal_cluster,
                **best_mode_cluster,
                **whole_rmsd,
                **interactions,
                "matched_heavy_atoms": matched_atoms,
                "lead_heavy_atoms": len(lead_reference),
                "candidate_mode_heavy_atoms_min": candidate_heavy_min,
                "min_pose_rmsd_to_lead": round(best_family_rmsd, 4),
                "min_pose_rmsd_to_origin_lead": round(best_rmsd, 4),
                "best_rmsd_mode": best_mode,
                "rmsd_q25_cutoff": "",
                "rmsd_q25_keep": "",
                "smiles": str(row.get("smiles", "")),
            }
        )
    finite_rmsds = [float(row["min_pose_rmsd_to_lead"]) for row in analysis_rows]
    cutoff = percentile(finite_rmsds, quantile)
    for row in analysis_rows:
        row["rmsd_q25_cutoff"] = round(cutoff, 4) if math.isfinite(cutoff) else ""
        row["rmsd_q25_keep"] = "true" if math.isfinite(cutoff) and float(row["min_pose_rmsd_to_lead"]) <= cutoff else "false"
    analysis_rows.sort(key=lambda row: (row["rmsd_q25_keep"] != "true", float(row["min_pose_rmsd_to_lead"]), float(row["affinity_kcal_mol"])))
    ranked_rows: list[dict[str, str | float | int]] = []
    alternative_rows: list[dict[str, str | float | int]] = []
    for row in analysis_rows:
        if row["candidate_id"] == "LEAD_REFERENCE":
            ranked_rows.append(
                {
                    "rank": 0,
                    "candidate_id": row["candidate_id"],
                    "final_score": 1.0,
                    "affinity_score": 1.0,
                    "inner_score": 1.0,
                    "lead_family_score": 1.0,
                    "whole_score": 1.0,
                    "interaction_score": 1.0,
                    "affinity_kcal_mol": row["affinity_kcal_mol"],
                    "affinity_penalty": 0.0,
                    "min_pose_rmsd_to_lead": 0.0,
                    "rmsd_score": 1.0,
                    "near_lead_mode_count": row.get("near_lead_mode_count", 0),
                    "near_lead_mode_fraction": row.get("near_lead_mode_fraction", 0.0),
                    "best_rmsd_mode": 1,
                    "mode_count": row["mode_count"],
                    "best_affinity_mode": row.get("best_affinity_mode", ""),
                    "best_affinity_mode_mean_rmsd": row.get("best_affinity_mode_mean_rmsd", ""),
                    "whole_rmsd": row.get("whole_rmsd", ""),
                    "hbond_count": row.get("hbond_count", ""),
                    "hydrophobic_count": row.get("hydrophobic_count", ""),
                    "vdw_contact_count": row.get("vdw_contact_count", ""),
                    "pi_contact_count": row.get("pi_contact_count", ""),
                    "ch_pi_count": row.get("ch_pi_count", ""),
                    "rmsd_q25_keep": "true",
                    "pose_path": row["pose_path"],
                    "smiles": row["smiles"],
                }
            )
            continue
        ranked_rows.append(
            {
                "rank": 0,
                "candidate_id": row["candidate_id"],
                "final_score": "",
                "affinity_score": "",
                "inner_score": "",
                "lead_family_score": "",
                "whole_score": "",
                "interaction_score": "",
                "affinity_kcal_mol": row["affinity_kcal_mol"],
                "affinity_penalty": "",
                "min_pose_rmsd_to_lead": row["min_pose_rmsd_to_lead"],
                "rmsd_score": "",
                "near_lead_mode_count": row.get("near_lead_mode_count", 0),
                "near_lead_mode_fraction": row.get("near_lead_mode_fraction", 0.0),
                "best_rmsd_mode": row["best_rmsd_mode"],
                "mode_count": row["mode_count"],
                "best_affinity_mode": row.get("best_affinity_mode", ""),
                "best_affinity_mode_cluster_fraction": row.get("best_affinity_mode_cluster_fraction", ""),
                "best_affinity_mode_mean_rmsd": row.get("best_affinity_mode_mean_rmsd", ""),
                "whole_rmsd": row.get("whole_rmsd", ""),
                "hbond_count": row.get("hbond_count", ""),
                "hydrophobic_count": row.get("hydrophobic_count", ""),
                "vdw_contact_count": row.get("vdw_contact_count", ""),
                "pi_contact_count": row.get("pi_contact_count", ""),
                "ch_pi_count": row.get("ch_pi_count", ""),
                "rmsd_q25_keep": row["rmsd_q25_keep"],
                "pose_path": row["pose_path"],
                "smiles": row["smiles"],
            }
        )
        min_rmsd = float(row["min_pose_rmsd_to_lead"])
        affinity = float(row["affinity_kcal_mol"])
        cluster_fraction = float(row.get("best_internal_cluster_fraction", 0.0))
        cluster_size = int(row.get("best_internal_cluster_size", 0))
        lead_like_count = int(row.get("near_lead_mode_count", 0))
        alternative_score = 0.65 * cluster_fraction + 0.35 * (1.0 / (1.0 + max(0.0, min_rmsd - near_rmsd_cutoff)))
        if lead_like_count == 0 and cluster_size >= 2:
            alternative_rows.append(
                {
                    "rank": 0,
                    "candidate_id": row["candidate_id"],
                    "alternative_score": round(alternative_score, 6),
                    "affinity_kcal_mol": row["affinity_kcal_mol"],
                    "min_pose_rmsd_to_lead": row["min_pose_rmsd_to_lead"],
                    "near_lead_mode_count": row.get("near_lead_mode_count", 0),
                    "mode_count": row["mode_count"],
                    "best_internal_cluster_size": row.get("best_internal_cluster_size", 0),
                    "best_internal_cluster_fraction": row.get("best_internal_cluster_fraction", 0.0),
                    "best_internal_cluster_representative_mode": row.get("best_internal_cluster_representative_mode", ""),
                    "best_internal_cluster_mean_rmsd": row.get("best_affinity_mode_mean_rmsd", row.get("best_internal_cluster_mean_rmsd", "")),
                    "pose_path": row["pose_path"],
                    "smiles": row["smiles"],
                }
            )
    affinity_scores = normal_scores([safe_float(row.get("affinity_kcal_mol")) for row in ranked_rows], lower_is_better=True)
    inner_scores = normal_scores([safe_float(row.get("best_affinity_mode_mean_rmsd")) for row in ranked_rows], lower_is_better=True)
    lead_scores = normal_scores([safe_float(row.get("min_pose_rmsd_to_lead")) for row in ranked_rows], lower_is_better=True)
    whole_scores = normal_scores([safe_float(row.get("whole_rmsd")) for row in ranked_rows], lower_is_better=True)
    interaction_values = [
        (
            (safe_float(row.get("hbond_count")) or 0.0) * 2.0
            + (safe_float(row.get("ch_pi_count")) or 0.0) * 1.8
            + (safe_float(row.get("pi_contact_count")) or 0.0) * 1.2
            + (safe_float(row.get("hydrophobic_count")) or 0.0) * 0.35
            + (safe_float(row.get("vdw_contact_count")) or 0.0) * 0.05
        )
        for row in ranked_rows
    ]
    interaction_scores = normal_scores(interaction_values, lower_is_better=False)
    for row, affinity_score, inner_score, lead_score, whole_score, interaction_score in zip(
        ranked_rows, affinity_scores, inner_scores, lead_scores, whole_scores, interaction_scores
    ):
        if row.get("candidate_id") == "LEAD_REFERENCE":
            row["affinity_score"] = 1.0
            row["inner_score"] = 1.0
            row["lead_family_score"] = 1.0
            row["whole_score"] = 1.0
            row["interaction_score"] = 1.0
            row["final_score"] = 1.0
            continue
        row["affinity_score"] = affinity_score if affinity_score is not None else ""
        row["inner_score"] = inner_score if inner_score is not None else ""
        row["lead_family_score"] = lead_score if lead_score is not None else ""
        row["whole_score"] = whole_score if whole_score is not None else ""
        row["interaction_score"] = interaction_score if interaction_score is not None else ""
        weighted_scores = [
            (interaction_score, 0.38),
            (inner_score, 0.24),
            (affinity_score, 0.22),
            (whole_score, 0.16),
        ]
        present_weight = sum(weight for score, weight in weighted_scores if score is not None)
        row["final_score"] = (
            round(sum(float(score) * weight for score, weight in weighted_scores if score is not None) / present_weight, 6)
            if present_weight > 0
            else ""
        )
    ranked_rows.sort(
        key=lambda row: (
            -(safe_float(row.get("final_score")) if safe_float(row.get("final_score")) is not None else -1.0),
            -(safe_float(row.get("interaction_score")) if safe_float(row.get("interaction_score")) is not None else -1.0),
            safe_float(row.get("best_affinity_mode_mean_rmsd")) if safe_float(row.get("best_affinity_mode_mean_rmsd")) is not None else 999.0,
            float(row["affinity_kcal_mol"]),
        )
    )
    for rank, row in enumerate(ranked_rows, start=1):
        row["rank"] = rank
    alternative_rows.sort(
        key=lambda row: (
            -float(row["alternative_score"]),
            -float(row["best_internal_cluster_fraction"]),
            float(row["affinity_kcal_mol"]),
        )
    )
    for rank, row in enumerate(alternative_rows, start=1):
        row["rank"] = rank
    mode_rows.sort(key=lambda row: (str(row["candidate_id"]), float(row["mode_rmsd_to_lead"])))
    return analysis_rows, mode_rows, ranked_rows, alternative_rows


def summarize_pose_against_lead(
    lead_reference: list[tuple[float, float, float]],
    pose_path: Path,
    near_rmsd_cutoff: float,
    internal_cluster_rmsd_cutoff: float,
) -> dict[str, str | float | int]:
    if not pose_path.exists():
        return {
            "live_min_pose_rmsd_to_lead": "",
            "live_near_lead_mode_count": "",
            "live_mode_count": "",
            "live_inner_rmsd": "",
            "live_inner_cluster_fraction": "",
            "live_whole_rmsd": "",
        }
    modes = parse_pdbqt_modes(pose_path)
    if not modes:
        return {
            "live_min_pose_rmsd_to_lead": "",
            "live_near_lead_mode_count": "",
            "live_mode_count": 0,
            "live_inner_rmsd": "",
            "live_inner_cluster_fraction": "",
            "live_whole_rmsd": "",
        }
    rmsd_values: list[tuple[int, float]] = []
    for mode in modes:
        rmsd, _ = coordinate_rmsd(lead_reference, mode.coords)
        if math.isfinite(rmsd):
            rmsd_values.append((mode.mode_index, rmsd))
    if not rmsd_values:
        return {
            "live_min_pose_rmsd_to_lead": "",
            "live_near_lead_mode_count": "",
            "live_mode_count": len(modes),
            "live_inner_rmsd": "",
            "live_inner_cluster_fraction": "",
            "live_whole_rmsd": whole_pose_rmsd_summary(modes).get("whole_rmsd", ""),
        }
    _, best_rmsd = min(rmsd_values, key=lambda item: item[1])
    near_count = sum(1 for _, rmsd in rmsd_values if rmsd <= near_rmsd_cutoff)
    best_mode_cluster = best_affinity_mode_summary(modes, internal_cluster_rmsd_cutoff)
    whole = whole_pose_rmsd_summary(modes)
    return {
        "live_min_pose_rmsd_to_lead": round(best_rmsd, 4),
        "live_near_lead_mode_count": near_count,
        "live_mode_count": len(modes),
        "live_inner_rmsd": best_mode_cluster.get("best_affinity_mode_mean_rmsd", ""),
        "live_inner_cluster_fraction": best_mode_cluster.get("best_affinity_mode_cluster_fraction", ""),
        "live_whole_rmsd": whole.get("whole_rmsd", ""),
    }


def build_combined_pose_summary(
    analysis_rows: list[dict[str, str | float | int]],
    ranked_rows: list[dict[str, str | float | int]],
    alternative_rows: list[dict[str, str | float | int]],
) -> list[dict[str, str | float | int]]:
    lead_like_by_id = {str(row["candidate_id"]): row for row in ranked_rows}
    alternative_by_id = {str(row["candidate_id"]): row for row in alternative_rows}
    combined: list[dict[str, str | float | int]] = []
    for analysis in analysis_rows:
        candidate_id = str(analysis["candidate_id"])
        lead_like = lead_like_by_id.get(candidate_id, {})
        alternative = alternative_by_id.get(candidate_id, {})
        is_reference = candidate_id == "LEAD_REFERENCE"
        is_lead_like = str(analysis.get("rmsd_q25_keep", "")) == "true" or int(analysis.get("near_lead_mode_count", 0)) > 0
        is_alternative = bool(alternative)
        if is_reference:
            classification = "lead_reference"
        elif is_lead_like and is_alternative:
            classification = "lead_like_and_alternative_cluster"
        elif is_lead_like:
            classification = "lead_like"
        elif is_alternative:
            classification = "alternative_stable_pose"
        else:
            classification = "weak_or_unstable_pose"
        combined.append(
            {
                "rank": 0,
                "candidate_id": candidate_id,
                "classification": classification,
                "final_score": lead_like.get("final_score", ""),
                "affinity_score": lead_like.get("affinity_score", ""),
                "inner_score": lead_like.get("inner_score", ""),
                "lead_family_score": lead_like.get("lead_family_score", ""),
                "whole_score": lead_like.get("whole_score", ""),
                "interaction_score": lead_like.get("interaction_score", ""),
                "lead_like_score": lead_like.get("final_score", ""),
                "alternative_score": alternative.get("alternative_score", ""),
                "affinity_kcal_mol": analysis.get("affinity_kcal_mol", ""),
                "affinity_penalty": lead_like.get("affinity_penalty", ""),
                "lead_difference_rmsd": analysis.get("min_pose_rmsd_to_lead", ""),
                "best_lead_rmsd_mode": analysis.get("best_rmsd_mode", ""),
                "near_lead_mode_count": analysis.get("near_lead_mode_count", ""),
                "near_lead_mode_fraction": analysis.get("near_lead_mode_fraction", ""),
                "inner_rmsd": analysis.get("best_affinity_mode_mean_rmsd", ""),
                "whole_rmsd": analysis.get("whole_rmsd", ""),
                "inner_cluster_size": analysis.get("best_affinity_mode_cluster_size", ""),
                "inner_cluster_fraction": analysis.get("best_affinity_mode_cluster_fraction", ""),
                "hbond_count": analysis.get("hbond_count", ""),
                "hydrophobic_count": analysis.get("hydrophobic_count", ""),
                "vdw_contact_count": analysis.get("vdw_contact_count", ""),
                "pi_contact_count": analysis.get("pi_contact_count", ""),
                "ch_pi_count": analysis.get("ch_pi_count", ""),
                "inner_cluster_representative_mode": analysis.get("best_affinity_mode", ""),
                "mode_count": analysis.get("mode_count", ""),
                "best_affinity_mode": analysis.get("best_affinity_mode", ""),
                "rmsd_q25_keep": analysis.get("rmsd_q25_keep", ""),
                "pose_path": analysis.get("pose_path", ""),
                "smiles": analysis.get("smiles", ""),
            }
        )
    combined.sort(
        key=lambda row: (
            row["classification"] != "lead_reference",
            row["classification"] not in {"lead_like", "lead_like_and_alternative_cluster"},
            -safe_float(row.get("lead_like_score")) if safe_float(row.get("lead_like_score")) is not None else 999.0,
            safe_float(row.get("lead_difference_rmsd")) if safe_float(row.get("lead_difference_rmsd")) is not None else 999.0,
            -safe_float(row.get("alternative_score")) if safe_float(row.get("alternative_score")) is not None else 999.0,
        )
    )
    for rank, row in enumerate(combined, start=1):
        row["rank"] = rank
    return combined


def choose_next_round_parents(
    combined_rows: list[dict[str, str | float | int]],
    top_n: int,
) -> list[dict[str, str | float | int]]:
    if top_n <= 0:
        return []

    def lead_like_sort_key(row: dict[str, str | float | int]) -> tuple[float, float, float]:
        lead_score = safe_float(row.get("lead_like_score"))
        lead_rmsd = safe_float(row.get("lead_difference_rmsd"))
        affinity = safe_float(row.get("affinity_kcal_mol"))
        return (
            -(lead_score if lead_score is not None else -999.0),
            lead_rmsd if lead_rmsd is not None else 999.0,
            affinity if affinity is not None else 999.0,
        )

    def alternative_sort_key(row: dict[str, str | float | int]) -> tuple[float, float, float, float]:
        alternative_score = safe_float(row.get("alternative_score"))
        inner_fraction = safe_float(row.get("inner_cluster_fraction"))
        inner_rmsd = safe_float(row.get("inner_rmsd"))
        affinity = safe_float(row.get("affinity_kcal_mol"))
        return (
            -(alternative_score if alternative_score is not None else -999.0),
            -(inner_fraction if inner_fraction is not None else -999.0),
            inner_rmsd if inner_rmsd is not None else 999.0,
            affinity if affinity is not None else 999.0,
        )

    lead_like_rows = sorted(
        [
            row
            for row in combined_rows
            if row.get("candidate_id") != "LEAD_REFERENCE"
            and row.get("classification") in {"lead_like", "lead_like_and_alternative_cluster"}
            and row.get("smiles")
        ],
        key=lead_like_sort_key,
    )
    alternative_rows = sorted(
        [
            row
            for row in combined_rows
            if row.get("candidate_id") != "LEAD_REFERENCE"
            and row.get("classification") in {"alternative_stable_pose", "lead_like_and_alternative_cluster"}
            and row.get("smiles")
        ],
        key=alternative_sort_key,
    )

    lead_quota = max(1, math.ceil(top_n * 0.6))
    alternative_quota = max(1, top_n - lead_quota)
    selected: list[dict[str, str | float | int]] = []
    seen_smiles: set[str] = set()

    def add_row(row: dict[str, str | float | int], reason: str) -> None:
        smiles = str(row.get("smiles", ""))
        if not smiles:
            return
        can = canonical(smiles) or smiles
        if can in seen_smiles:
            return
        seen_smiles.add(can)
        selected.append(
            {
                "parent_rank": len(selected) + 1,
                "candidate_id": row.get("candidate_id", ""),
                "selection_reason": reason,
                "lead_difference_rmsd": row.get("lead_difference_rmsd", ""),
                "inner_rmsd": row.get("inner_rmsd", ""),
                "affinity_kcal_mol": row.get("affinity_kcal_mol", ""),
                "smiles": smiles,
            }
        )

    for row in lead_like_rows[:lead_quota]:
        add_row(row, f"lead_like_followup:{row.get('classification', '')}")
    for row in alternative_rows[:alternative_quota]:
        add_row(row, f"alternative_pose_followup:{row.get('classification', '')}")
    for row in lead_like_rows + alternative_rows:
        if len(selected) >= top_n:
            break
        add_row(row, f"fill_followup:{row.get('classification', '')}")
    return selected


def write_next_parent_files(outdir: Path, parent_rows: list[dict[str, str | float | int]]) -> None:
    valid_rows: list[dict[str, str | float | int]] = []
    lines: list[str] = []
    for row in parent_rows:
        smiles = canonical(str(row.get("smiles", "") or ""))
        if not smiles:
            print(f"warning: skipped next-round parent without valid SMILES: {row}", flush=True)
            continue
        enriched = enrich_display_row({**row, "smiles": smiles})
        label = str(enriched.get("nickname") or enriched.get("candidate_id") or enriched.get("event_id") or f"P{len(valid_rows) + 1:04d}")
        valid_rows.append(enriched)
        lines.append(f"{smiles}\t{label}")
    try_write_csv(outdir / "next_round_parents.csv", valid_rows, fieldnames=NEXT_PARENT_FIELDNAMES, label="next_round_parents.csv")
    (outdir / "next_round_parents.smi").write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def append_attempted_smiles(outdir: Path, row: dict[str, str | float | int]) -> None:
    smiles = str(row.get("smiles", "") or "")
    can = canonical(smiles) if smiles else None
    if not can:
        return
    attempted_path = outdir / ATTEMPTED_SMILES_FILENAME
    existing: set[str] = set()
    if attempted_path.exists():
        try:
            for line in attempted_path.read_text(encoding="utf-8", errors="replace").splitlines():
                parts = line.strip().split()
                if parts:
                    existing_can = canonical(parts[0])
                    if existing_can:
                        existing.add(existing_can)
        except Exception:
            existing = set()
    if can in existing:
        return
    fields = [
        can,
        str(row.get("candidate_id", "") or ""),
        str(row.get("run_id", "") or ""),
        str(row.get("status", "") or ""),
        str(row.get("edit_label", "") or ""),
    ]
    with attempted_path.open("a", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(fields) + "\n")


def append_successful_lead(outdir: Path, lead_smiles: str, row: dict[str, str | float | int]) -> bool:
    if str(row.get("status", "") or "") != "ok":
        return False
    can = canonical(lead_smiles)
    if not can:
        return False
    lead_path = outdir / "lead.csv"
    existing_rows = read_csv_rows(lead_path) if lead_path.exists() else []
    for existing in existing_rows:
        existing_can = canonical(existing.get("smiles", ""))
        if existing_can == can:
            return False
    new_row: dict[str, str | float | int] = {
        "lead_index": len(existing_rows) + 1,
        "first_seen_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "event_id": row_event_id(row),
        "candidate_id": row.get("candidate_id", ""),
        **smiles_identity_fields(can, lead_smiles),
        "smiles": can,
        "affinity_kcal_mol": row.get("affinity_kcal_mol", ""),
        "live_min_pose_rmsd_to_lead": row.get("live_min_pose_rmsd_to_lead", ""),
        "live_near_lead_mode_count": row.get("live_near_lead_mode_count", ""),
        "live_mode_count": row.get("live_mode_count", ""),
        "live_inner_rmsd": row.get("live_inner_rmsd", ""),
        "live_inner_cluster_fraction": row.get("live_inner_cluster_fraction", ""),
        "live_whole_rmsd": row.get("live_whole_rmsd", ""),
        "pose_path": row.get("pose_path", ""),
        "log_path": row.get("log_path", ""),
        "mw": row.get("mw", ""),
        "logp": row.get("logp", ""),
        "hbd": row.get("hbd", ""),
        "hba": row.get("hba", ""),
        "tpsa": row.get("tpsa", ""),
        "rot_bonds": row.get("rot_bonds", ""),
        "heavy_atoms": row.get("heavy_atoms", ""),
        "formal_charge": row.get("formal_charge", ""),
    }
    existing_rows.append({key: str(value) for key, value in new_row.items()})
    try_write_csv(lead_path, existing_rows, fieldnames=LEAD_FIELDNAMES, label="lead.csv")
    return True


def write_ledger_outputs(
    *,
    outdir: Path,
    lead_pose: Path,
    receptor: Path,
    all_screening_rows: list[dict[str, str | float | int]],
    lead_smiles: str,
    rmsd_quantile: float,
    near_rmsd_cutoff: float,
    internal_cluster_rmsd_cutoff: float,
    affinity_threshold: float,
    affinity_penalty_scale: float,
    next_parent_count: int,
    write_debug_csvs: bool,
) -> dict[str, int | str | float]:
    nickname_map = load_user_nickname_map(outdir)
    rmsd_rows, mode_rmsd_rows, final_rank_rows, alternative_rows = analyze_pose_rmsd(
        lead_pose.resolve(),
        receptor.resolve(),
        all_screening_rows,
        lead_smiles,
        rmsd_quantile,
        near_rmsd_cutoff,
        internal_cluster_rmsd_cutoff,
        affinity_threshold,
        affinity_penalty_scale,
    )
    all_screening_rows = enrich_screening_rows_with_pose_analysis(all_screening_rows, rmsd_rows)
    all_screening_rows = apply_user_nicknames(all_screening_rows, nickname_map)
    try_write_csv(
        outdir / "all_screening_results.csv",
        all_screening_rows,
        fieldnames=SCREENING_FIELDNAMES,
        label="all_screening_results.csv",
    )
    combined_rows = build_combined_pose_summary(rmsd_rows, final_rank_rows, alternative_rows)
    molecule_rows = compact_molecule_rows(combined_rows, all_screening_rows)
    molecule_rows = apply_user_nicknames(molecule_rows, nickname_map)
    mode_rows = compact_mode_rows(mode_rmsd_rows, molecule_rows)
    history_rows = molecule_rows_to_history_rows(molecule_rows)
    history_rows = apply_user_nicknames(history_rows, nickname_map)
    try_write_csv(
        outdir / "final_ranked_candidates.csv",
        molecule_rows,
        fieldnames=FINAL_RANK_FIELDNAMES,
        label="final_ranked_candidates.csv",
    )
    try_write_csv(
        outdir / "history_summary.csv",
        history_rows,
        fieldnames=HISTORY_SUMMARY_FIELDNAMES,
        label="history_summary.csv",
    )
    try_write_csv(
        outdir / "mode_pose_details.csv",
        mode_rows,
        fieldnames=COMPACT_MODE_FIELDNAMES,
        label="mode_pose_details.csv",
    )
    next_parent_rows = choose_next_round_parents(combined_rows, next_parent_count)
    write_next_parent_files(outdir, next_parent_rows)
    if write_debug_csvs:
        try_write_csv(outdir / "pose_rmsd_analysis.csv", rmsd_rows, fieldnames=POSE_RMSD_FIELDNAMES, label="pose_rmsd_analysis.csv")
        try_write_csv(outdir / "pose_mode_rmsd.csv", mode_rmsd_rows, fieldnames=POSE_MODE_RMSD_FIELDNAMES, label="pose_mode_rmsd.csv")
        try_write_csv(outdir / "lead_like_ranked_candidates.csv", final_rank_rows, fieldnames=FINAL_RANK_FIELDNAMES, label="lead_like_ranked_candidates.csv")
        try_write_csv(outdir / "alternative_mode_candidates.csv", alternative_rows, fieldnames=ALTERNATIVE_MODE_FIELDNAMES, label="alternative_mode_candidates.csv")
        try_write_csv(
            outdir / "pose_rmsd_filtered.csv",
            [rmsd_row for rmsd_row in rmsd_rows if rmsd_row.get("rmsd_q25_keep") == "true"],
            fieldnames=POSE_RMSD_FIELDNAMES,
            label="pose_rmsd_filtered.csv",
        )
    return {
        "all_screening_rows": len(all_screening_rows),
        "all_screening_rows_with_pose_metrics": sum(1 for row in all_screening_rows if row.get("live_whole_rmsd") != ""),
        "molecule_rows": len(molecule_rows),
        "mode_rows": len(mode_rows),
        "pose_rmsd_rows": len(rmsd_rows),
        "raw_mode_rows": len(mode_rmsd_rows),
        "next_parent_rows": len(next_parent_rows),
        "rmsd_cutoff": rmsd_rows[0].get("rmsd_q25_cutoff", "") if rmsd_rows else "",
    }


def refresh_ledger_from_history_only(
    *,
    outdir: Path,
    history_dirs: list[Path],
    lead_pose: Path | None,
    receptor: Path | None,
    lead_smiles: str,
    rmsd_quantile: float,
    near_rmsd_cutoff: float,
    internal_cluster_rmsd_cutoff: float,
    affinity_threshold: float,
    affinity_penalty_scale: float,
    next_parent_count: int,
    write_debug_csvs: bool,
) -> dict[str, int | str | float]:
    history_rows = load_history_screening_rows(history_dirs)
    nickname_map = load_user_nickname_map(outdir)
    history_rows = apply_user_nicknames(history_rows, nickname_map)
    if lead_pose is not None:
        return write_ledger_outputs(
            outdir=outdir,
            lead_pose=lead_pose,
            receptor=receptor if receptor is not None else Path(""),
            all_screening_rows=history_rows,
            lead_smiles=lead_smiles,
            rmsd_quantile=rmsd_quantile,
            near_rmsd_cutoff=near_rmsd_cutoff,
            internal_cluster_rmsd_cutoff=internal_cluster_rmsd_cutoff,
            affinity_threshold=affinity_threshold,
            affinity_penalty_scale=affinity_penalty_scale,
            next_parent_count=next_parent_count,
            write_debug_csvs=write_debug_csvs,
        )
    try_write_csv(outdir / "all_screening_results.csv", history_rows, fieldnames=SCREENING_FIELDNAMES, label="all_screening_results.csv")
    molecule_rows = compact_molecule_rows([], history_rows)
    molecule_rows = apply_user_nicknames(molecule_rows, nickname_map)
    try_write_csv(outdir / "final_ranked_candidates.csv", molecule_rows, fieldnames=FINAL_RANK_FIELDNAMES, label="final_ranked_candidates.csv")
    try_write_csv(outdir / "history_summary.csv", apply_user_nicknames(molecule_rows_to_history_rows(molecule_rows), nickname_map), fieldnames=HISTORY_SUMMARY_FIELDNAMES, label="history_summary.csv")
    return {"history_rows": len(history_rows), "molecule_rows": len(molecule_rows)}


def write_csv(path: Path, rows: list[dict[str, str | float | int]], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = sorted({key for row in rows for key in row.keys()})
    if not fieldnames:
        raise ValueError(f"cannot write CSV without fieldnames: {path}")
    materialized_rows = [enrich_display_row(row) for row in rows]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows([{field: row.get(field, "") for field in fieldnames} for row in materialized_rows])


def try_write_csv(path: Path, rows: list[dict[str, str | float | int]], fieldnames: list[str], label: str) -> bool:
    try:
        write_csv(path, rows, fieldnames=fieldnames)
        return True
    except PermissionError:
        print(f"\nwarning: cannot update {label}; close the file if it is open: {path}", flush=True)
        return False


def write_progress_status(
    path: Path,
    *,
    started_at: float,
    round_index: int,
    current_index: int,
    total_current_round: int,
    candidate: Candidate | None,
    status: str,
    affinity: float | None,
    best_affinity: float,
    best_smiles: str,
    baseline: float,
    screened_count: int,
    rejected_count: int,
    rmsd_summary: dict[str, str | float | int] | None = None,
) -> None:
    elapsed = time.time() - started_at
    lines = [
        f"elapsed_seconds={elapsed:.1f}",
        f"round={round_index}",
        f"current={current_index}/{total_current_round}",
        f"screened_count={screened_count}",
        f"rejected_count={rejected_count}",
        f"baseline_affinity={baseline:.3f}",
        f"best_affinity={best_affinity:.3f}",
        f"best_delta_vs_baseline={best_affinity - baseline:.3f}",
        f"best_smiles={best_smiles}",
        f"current_status={status}",
    ]
    if candidate is not None:
        lines.extend(
            [
                f"current_candidate_id={candidate.candidate_id}",
                f"current_parent_id={candidate.parent_id}",
                f"current_edit_label={candidate.edit_label}",
                f"current_smiles={candidate.smiles}",
            ]
        )
    if affinity is not None:
        lines.append(f"current_affinity={affinity:.3f}")
        lines.append(f"current_delta_vs_baseline={affinity - baseline:.3f}")
    if rmsd_summary:
        for key in (
            "live_min_pose_rmsd_to_lead",
            "live_near_lead_mode_count",
            "live_mode_count",
            "live_inner_rmsd",
            "live_inner_cluster_fraction",
            "live_whole_rmsd",
        ):
            lines.append(f"{key}={rmsd_summary.get(key, '')}")
    try:
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    except PermissionError:
        pass


def print_progress_line(
    *,
    round_index: int,
    current_index: int,
    total_current_round: int,
    candidate_id: str,
    edit_label: str,
    status: str,
    affinity: float | None,
    best_affinity: float,
    baseline: float,
    rmsd_summary: dict[str, str | float | int] | None = None,
) -> None:
    affinity_text = "NA" if affinity is None else f"{affinity:.2f}"
    delta_text = "NA" if affinity is None else f"{affinity - baseline:+.2f}"
    rmsd_text = ""
    if rmsd_summary and rmsd_summary.get("live_min_pose_rmsd_to_lead") != "":
        rmsd_text = (
            f" rmsd={float(rmsd_summary['live_min_pose_rmsd_to_lead']):.2f}"
            f" near={rmsd_summary.get('live_near_lead_mode_count', '')}/{rmsd_summary.get('live_mode_count', '')}"
        )
        if rmsd_summary.get("live_inner_rmsd") != "":
            rmsd_text += (
                f" inner={float(rmsd_summary['live_inner_rmsd']):.2f}"
                f" inner_frac={rmsd_summary.get('live_inner_cluster_fraction', '')}"
            )
        if rmsd_summary.get("live_whole_rmsd") != "":
            rmsd_text += f" whole={float(rmsd_summary['live_whole_rmsd']):.2f}"
    message = (
        f"[round {round_index} {current_index}/{total_current_round}] "
        f"{candidate_id} {edit_label} status={status} affinity={affinity_text}{rmsd_text}"
    )
    print(message, flush=True)


def assert_writable_outdir(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    probe = outdir / ".lead_optimizer_write_test"
    try:
        probe.write_text("ok\n", encoding="utf-8")
        probe.unlink()
    except PermissionError as exc:
        raise PermissionError(
            f"Cannot write to output directory: {outdir}. "
            "Close any Windows program that has CSV/log files open there, or choose a new --outdir. "
            f"Original error: {exc}"
        ) from exc


def snapshot_existing_history(outdir: Path, run_id: str) -> Path | None:
    csv_paths = [path for path in outdir.glob("*.csv") if path.is_file()]
    if not csv_paths:
        return None
    snapshot_dir = outdir / "history_snapshots" / run_id
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    for csv_path in csv_paths:
        try:
            shutil.copy2(csv_path, snapshot_dir / csv_path.name)
        except PermissionError:
            print(f"warning: could not snapshot locked history file: {csv_path}", flush=True)
    return snapshot_dir


def migrate_history_csvs(history_dirs: list[Path], outdir: Path, run_id: str) -> Path | None:
    sources: list[Path] = []
    for history_dir in history_dirs:
        directory = history_dir.expanduser().resolve()
        if directory == outdir or not directory.exists():
            continue
        sources.extend(path for path in directory.rglob("*.csv") if path.is_file())
    if not sources:
        return None
    migrated_dir = outdir / "migrated_history" / run_id
    migrated_dir.mkdir(parents=True, exist_ok=True)
    for index, source in enumerate(sources, start=1):
        safe_parent = re.sub(r"[^A-Za-z0-9_.-]+", "_", source.parent.name)
        target = migrated_dir / f"{index:04d}_{safe_parent}_{source.name}"
        try:
            shutil.copy2(source, target)
        except PermissionError:
            print(f"warning: could not migrate locked history file: {source}", flush=True)
    return migrated_dir


def existing_csv_fieldnames(path: Path) -> list[str]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            return next(reader)
        except StopIteration:
            return []


def merged_fieldnames(existing: list[str], required: list[str]) -> list[str]:
    fieldnames = list(existing)
    for field in required:
        if field not in fieldnames:
            fieldnames.append(field)
    return fieldnames


def dock_single_smiles_to_ledger(
    *,
    smiles: str,
    outdir: Path,
    run_id: str,
    candidate_id: str,
    edit_label: str,
    receptor: Path,
    box: ConfigBox,
    vina: str,
    meeko: str,
    obabel: str,
    seed: int,
    baseline: float | None,
    num_modes: int = 9,
    input_smiles: str | None = None,
    update_all_screening: bool = True,
) -> dict[str, str | float | int]:
    can = canonical(smiles)
    if can is None:
        raise ValueError(f"single dock SMILES is invalid: {smiles}")
    sdf_dir = outdir / "sdf"
    pdbqt_dir = outdir / "pdbqt"
    pose_dir = outdir / "poses"
    log_dir = outdir / "logs"
    prep_log_dir = outdir / "prep_logs"
    for directory in (sdf_dir, pdbqt_dir, pose_dir, log_dir, prep_log_dir):
        directory.mkdir(parents=True, exist_ok=True)
    props = properties(can)
    safe_candidate_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", candidate_id) or "SINGLE"
    file_stem = f"{run_id}_{safe_candidate_id}"
    sdf_path = sdf_dir / f"{file_stem}.sdf"
    pdbqt_path = pdbqt_dir / f"{file_stem}.pdbqt"
    pose_path = pose_dir / f"{file_stem}.pdbqt"
    log_path = log_dir / f"{file_stem}.log"
    prep_log_path = prep_log_dir / f"{file_stem}.log"
    status = "ok"
    affinity: float | str = ""
    reason = ""
    prep_mode = ""
    try:
        embed_sdf(can, safe_candidate_id, sdf_path, seed)
        prep_ok, prep_output, prep_mode = prepare_pdbqt(sdf_path, pdbqt_path, meeko, obabel)
        prep_log_path.write_text(prep_output, encoding="utf-8", errors="replace")
        if not prep_ok:
            status = "prep_failed"
            reason = "pdbqt_preparation_failed"
        else:
            dock_ok, dock_affinity, _ = dock(vina, receptor, pdbqt_path, pose_path, log_path, box, seed, num_modes)
            if dock_ok and dock_affinity is not None:
                affinity = dock_affinity
            else:
                status = "dock_failed"
                reason = "vina_failed_or_no_affinity"
    except Exception as exc:
        status = "failed"
        reason = str(exc)
    row: dict[str, str | float | int] = {
        "candidate_id": safe_candidate_id,
        "run_id": run_id,
        "round": 0,
        "parent_id": "MANUAL",
        "parent_smiles": "",
        "edit_label": edit_label,
        **smiles_identity_fields(can, input_smiles or smiles),
        "smiles": can,
        "status": status,
        "affinity_kcal_mol": affinity,
        "delta_vs_baseline": round(float(affinity) - baseline, 3) if isinstance(affinity, float) and baseline is not None else "",
        "live_min_pose_rmsd_to_lead": "",
        "live_near_lead_mode_count": "",
        "live_mode_count": "",
        "live_inner_rmsd": "",
        "live_inner_cluster_fraction": "",
        "live_whole_rmsd": "",
        "reason": reason,
        "prep_mode": prep_mode,
        "sdf_path": str(sdf_path),
        "pdbqt_path": str(pdbqt_path),
        "pose_path": str(pose_path),
        "log_path": str(log_path),
        **props,
    }
    if update_all_screening:
        all_screening_path = outdir / "all_screening_results.csv"
        existing_rows = read_csv_rows(all_screening_path) if all_screening_path.exists() else []
        merged_rows: dict[tuple[str, str], dict[str, str | float | int]] = {}
        for existing in existing_rows:
            key = (existing.get("candidate_id", ""), existing.get("smiles", ""))
            merged_rows[key] = existing
        merged_rows[(safe_candidate_id, can)] = row
        write_csv(all_screening_path, list(merged_rows.values()), fieldnames=SCREENING_FIELDNAMES)
    return row


def historical_rows_needing_redock(rows: list[dict[str, str | float | int]]) -> list[dict[str, str | float | int]]:
    needed: list[dict[str, str | float | int]] = []
    seen_smiles: set[str] = set()
    for row in rows:
        smiles = str(row.get("smiles", ""))
        can = canonical(smiles) if smiles else None
        if not can or can in seen_smiles:
            continue
        pose_path_value = str(row.get("pose_path", ""))
        pose_path = host_path(pose_path_value) if pose_path_value else None
        pose_ok = bool(pose_path and pose_path.exists() and parse_pdbqt_modes(pose_path))
        missing_metrics = not row.get("live_inner_rmsd") and not row.get("inner_rmsd")
        history_unscored = str(row.get("classification", "")) == "history_unscored_pose"
        if history_unscored or not pose_ok or missing_metrics:
            needed.append(row)
            seen_smiles.add(can)
    return needed


def redock_history_rows(
    *,
    rows: list[dict[str, str | float | int]],
    outdir: Path,
    run_id: str,
    max_count: int,
    receptor: Path,
    box: ConfigBox,
    vina: str,
    meeko: str,
    obabel: str,
    seed: int,
    baseline: float | None,
    num_modes: int = 9,
) -> list[dict[str, str | float | int]]:
    if max_count <= 0:
        return []
    redocked: list[dict[str, str | float | int]] = []
    for index, source_row in enumerate(historical_rows_needing_redock(rows)[:max_count], start=1):
        smiles = row_smiles_text(source_row)
        can = canonical(smiles)
        if not can:
            continue
        old_affinity = safe_float(source_row.get("affinity_kcal_mol") or source_row.get("best_affinity_kcal_mol"))
        row = dock_single_smiles_to_ledger(
            smiles=can,
            outdir=outdir,
            run_id=run_id,
            candidate_id=f"HRE_{index:04d}",
            edit_label=f"history_redock:{source_row.get('candidate_id', '')}",
            receptor=receptor,
            box=box,
            vina=vina,
            meeko=meeko,
            obabel=obabel,
            seed=seed + index,
            baseline=baseline,
            num_modes=num_modes,
            input_smiles=str(source_row.get("input_smiles", "") or smiles),
            update_all_screening=False,
        )
        new_affinity = safe_float(row.get("affinity_kcal_mol"))
        if old_affinity is not None and new_affinity is not None:
            delta = round(new_affinity - old_affinity, 3)
            row["affinity_recheck_delta"] = delta
            row["affinity_recheck_status"] = "consistent" if abs(delta) <= 0.2 else "changed_gt_0.2"
        else:
            row["affinity_recheck_delta"] = ""
            row["affinity_recheck_status"] = "missing_old_or_new_affinity"
        redocked.append(row)
    return redocked


def history_pose_for_smiles(
    smiles: str,
    rows: list[dict[str, str | float | int]],
) -> Path | None:
    can = canonical(smiles)
    if can is None:
        return None
    best_pose: Path | None = None
    best_affinity = float("inf")
    for row in rows:
        row_smiles = row_smiles_text(row)
        row_can = canonical(row_smiles) if row_smiles else None
        if row_can != can:
            continue
        pose_path_value = str(row.get("pose_path", ""))
        pose_path = host_path(pose_path_value) if pose_path_value else None
        if not pose_path or not pose_path.exists():
            continue
        try:
            if not parse_pdbqt_modes(pose_path):
                continue
        except Exception:
            continue
        affinity = safe_float(row.get("affinity_kcal_mol") or row.get("best_affinity_kcal_mol"))
        affinity_key = affinity if affinity is not None else 999.0
        if best_pose is None or affinity_key < best_affinity:
            best_pose = pose_path
            best_affinity = affinity_key
    return best_pose


def resolve_lead_pose_reference(
    *,
    lead_smiles: str,
    explicit_lead_pose: Path | None,
    history_screening_rows: list[dict[str, str | float | int]],
    outdir: Path,
    run_id: str,
    receptor: Path,
    box: ConfigBox,
    vina: str,
    meeko: str,
    obabel: str,
    seed: int,
    baseline: float | None,
    num_modes: int,
    allow_auto_dock: bool,
    lead_input_smiles: str | None = None,
) -> tuple[Path | None, list[dict[str, str | float | int]], list[str]]:
    notes: list[str] = []
    if explicit_lead_pose is not None:
        lead_pose = explicit_lead_pose.resolve()
        notes.append(f"lead_pose_source=explicit:{lead_pose}")
        return lead_pose, [], notes

    history_pose = history_pose_for_smiles(lead_smiles, history_screening_rows)
    if history_pose is not None:
        notes.append(f"lead_pose_source=history:{history_pose}")
        print(f"lead_pose_source=history {history_pose}", flush=True)
        return history_pose, [], notes

    if not allow_auto_dock:
        notes.append("lead_pose_source=none; refresh_only_did_not_auto_dock=true")
        return None, [], notes

    print("lead_pose_source=auto_dock; docking --lead-smiles once as the RMSD reference", flush=True)
    row = dock_single_smiles_to_ledger(
        smiles=lead_smiles,
        outdir=outdir,
        run_id=f"{run_id}_lead_reference",
        candidate_id="LEAD_AUTO",
        edit_label="auto_lead_reference",
        receptor=receptor,
        box=box,
        vina=vina,
        meeko=meeko,
        obabel=obabel,
        seed=seed,
        baseline=baseline,
        num_modes=num_modes,
        input_smiles=lead_input_smiles or lead_smiles,
        update_all_screening=False,
    )
    if row.get("status") != "ok":
        notes.append(f"lead_pose_source=auto_dock_failed:{row.get('reason', '')}")
        print(f"warning: auto lead pose docking failed: {row.get('reason', '')}", flush=True)
        return None, [row], notes
    lead_pose = host_path(str(row.get("pose_path", "")))
    if lead_pose is None or not lead_pose.exists():
        notes.append("lead_pose_source=auto_dock_failed:missing_pose_path")
        return None, [row], notes
    notes.append(f"lead_pose_source=auto_dock:{lead_pose}")
    return lead_pose, [row], notes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lead-smiles")
    parser.add_argument("--parent-smiles-file", type=Path, help="Optional SMI/CSV file of parent SMILES for iterative expansion")
    parser.add_argument("--baseline-affinity", type=float)
    parser.add_argument("--baseline-output", type=Path, help="Existing Vina PDBQT/log to parse when --baseline-affinity is omitted")
    parser.add_argument("--stop-affinity", type=float, help="Optional early stop once a candidate reaches this affinity or better; disabled by default")
    parser.add_argument("--receptor", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--vina", default="vina")
    parser.add_argument("--meeko", default="mk_prepare_ligand.py")
    parser.add_argument("--obabel", default="obabel")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--max-rounds", type=int, default=2)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--exhaustiveness", type=int)
    parser.add_argument("--num-modes", type=int, default=9, help="Maximum Vina poses to request per ligand; Vina default is 9")
    parser.add_argument("--energy-range", type=int, help="Maximum kcal/mol range above the best mode for saved Vina poses; larger values can retain more diverse poses")
    parser.add_argument("--center-x", type=float)
    parser.add_argument("--center-y", type=float)
    parser.add_argument("--center-z", type=float)
    parser.add_argument("--size-x", type=float)
    parser.add_argument("--size-y", type=float)
    parser.add_argument("--size-z", type=float)
    parser.add_argument("--csv-update-interval", type=int, default=10, help="Refresh compact ledger CSVs every N docked ligands; 0 means only at the end")
    parser.add_argument("--lead-pose", type=Path, help="Lead _out.pdbqt file; mode 1 is used as the RMSD reference for post-run pose triage")
    parser.add_argument("--rmsd-quantile", type=float, default=0.25, help="Pose RMSD quantile to keep during post-run triage")
    parser.add_argument("--near-rmsd-cutoff", type=float, default=4.0, help="Mode RMSD cutoff in Angstrom for counting lead-like modes")
    parser.add_argument("--internal-cluster-rmsd-cutoff", type=float, default=3.0, help="Mode-to-mode RMSD cutoff in Angstrom for detecting stable alternative pose clusters")
    parser.add_argument("--affinity-threshold", type=float, default=-10.0, help="Deprecated compatibility option; final ranking now uses normal-distribution affinity scoring instead of hard affinity penalties")
    parser.add_argument("--affinity-penalty-scale", type=float, default=3.0, help="Deprecated compatibility option")
    parser.add_argument("--run-id", default="", help="Optional run id prefix for output pose and ligand filenames; defaults to timestamp")
    parser.add_argument("--history-dir", type=Path, action="append", default=[], help="Previous run directory to avoid rerunning known SMILES/edit labels; may be repeated")
    parser.add_argument("--no-outdir-history", action="store_true", help="Do not automatically use --outdir as a continuation history source")
    parser.add_argument("--snapshot-history", action="store_true", help="Optional manual CSV backup; disabled by default to keep ledger directories clean")
    parser.add_argument("--next-parent-count", type=int, default=3, help="Write the top N candidates as next-round parents after final ranking")
    parser.add_argument("--deterministic-batch", action="store_true", help="Use fixed enumeration order instead of the default randomized candidate batch")
    parser.add_argument("--edit-mode", default="default", help="Limit generated edits: add, delete, shrink, drastic, or comma combinations; default uses all")
    parser.add_argument("--add", dest="add_only", action="store_true", help="Shortcut for --edit-mode add")
    parser.add_argument("--delete", dest="delete_only", action="store_true", help="Shortcut for --edit-mode delete")
    parser.add_argument("--shrink", dest="shrink_only", action="store_true", help="Shortcut for --edit-mode shrink")
    parser.add_argument("--drastic", dest="drastic_only", action="store_true", help="Shortcut for --edit-mode drastic")
    parser.add_argument("--write-debug-csvs", action="store_true", help="Write verbose debug CSVs: screening_results.csv, rejected_candidates.csv, pose_rmsd_filtered.csv")
    parser.add_argument("--refresh-ledger-only", action="store_true", help="Merge and rewrite compact ledger CSVs from history without docking new molecules")
    parser.add_argument("--rescore-history", action="store_true", help="For --max-rounds 0, also recompute RMSD/ranking for the full historical ledger; slower on large ledgers")
    parser.add_argument("--redock-history-unscored", action="store_true", help="Redock historical rows with missing/unscored poses before ranking")
    parser.add_argument("--history-redock-limit", type=int, default=30, help="Maximum history_unscored/missing-pose SMILES to redock per run; 0 disables")
    args = parser.parse_args()

    if not args.lead_smiles:
        raise ValueError("provide --lead-smiles")
    if not args.receptor or not args.config:
        raise ValueError("provide --receptor and --config for docking mode")

    lead = canonical(args.lead_smiles)
    if lead is None:
        raise ValueError("lead SMILES is invalid")
    lead_input_smiles = clean_smiles_text(args.lead_smiles)
    baseline = args.baseline_affinity
    if baseline is None:
        if not args.baseline_output:
            raise ValueError("provide --baseline-affinity or --baseline-output")
        baseline = read_best_affinity(args.baseline_output)
    stop_affinity = args.stop_affinity

    outdir = args.outdir.resolve()
    assert_writable_outdir(outdir)
    user_nickname_map = load_user_nickname_map(outdir)
    run_id = args.run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
    continuation_history_dirs = list(args.history_dir)
    if not args.no_outdir_history and not continuation_history_dirs:
        continuation_history_dirs.append(outdir)
    snapshot_dir = snapshot_existing_history(outdir, run_id) if args.snapshot_history and continuation_history_dirs else None
    migrated_dir = None
    sdf_dir = outdir / "sdf"
    pdbqt_dir = outdir / "pdbqt"
    pose_dir = outdir / "poses"
    log_dir = outdir / "logs"
    prep_log_dir = outdir / "prep_logs"
    for directory in (sdf_dir, pdbqt_dir, pose_dir, log_dir, prep_log_dir):
        directory.mkdir(parents=True, exist_ok=True)

    box = parse_config(
        args.config.resolve(),
        args.cpu,
        args.exhaustiveness,
        center_x_override=args.center_x,
        center_y_override=args.center_y,
        center_z_override=args.center_z,
        size_x_override=args.size_x,
        size_y_override=args.size_y,
        size_z_override=args.size_z,
        energy_range_override=args.energy_range,
    )
    write_effective_config(outdir / "effective_vina_config.txt", args.receptor.resolve(), box)
    history_seen_smiles, history_seen_edit_labels, history_notes = load_history(continuation_history_dirs)
    if snapshot_dir is not None:
        history_notes.append(f"history_snapshot_dir={snapshot_dir}")
    if migrated_dir is not None:
        history_notes.append(f"migrated_history_dir={migrated_dir}")
    if any(history_dir.expanduser().resolve() == outdir for history_dir in continuation_history_dirs):
        history_notes.append("continuation_history_uses_outdir=true")
    history_screening_rows = load_history_screening_rows(continuation_history_dirs)
    history_screening_rows = apply_user_nicknames(history_screening_rows, user_nickname_map)
    redocked_history_rows: list[dict[str, str | float | int]] = []
    if args.redock_history_unscored and args.history_redock_limit > 0:
        redocked_history_rows = redock_history_rows(
            rows=history_screening_rows,
            outdir=outdir,
            run_id=f"{run_id}_history_redock",
            max_count=args.history_redock_limit,
            receptor=args.receptor.resolve(),
            box=box,
            vina=args.vina,
            meeko=args.meeko,
            obabel=args.obabel,
            seed=args.seed,
            baseline=baseline,
            num_modes=args.num_modes,
            input_smiles=lead_input_smiles,
            update_all_screening=False,
        )
        if redocked_history_rows:
            history_screening_rows = merge_screening_rows(history_screening_rows, redocked_history_rows)
    active_lead_pose, lead_reference_rows, lead_pose_notes = resolve_lead_pose_reference(
        lead_smiles=lead,
        explicit_lead_pose=args.lead_pose,
        history_screening_rows=history_screening_rows,
        outdir=outdir,
        run_id=run_id,
        receptor=args.receptor.resolve(),
        box=box,
        vina=args.vina,
        meeko=args.meeko,
        obabel=args.obabel,
        seed=args.seed,
        baseline=baseline,
        num_modes=args.num_modes,
        allow_auto_dock=not args.refresh_ledger_only and args.max_rounds != 0,
        lead_input_smiles=lead_input_smiles,
    )
    if lead_reference_rows:
        history_screening_rows = merge_screening_rows(history_screening_rows, lead_reference_rows)
        history_seen_smiles.add(lead)
    rows: list[dict[str, str | float | int]] = []
    rejected_rows: list[dict[str, str | float | int]] = []
    started_at = time.time()
    progress_path = outdir / "progress_status.txt"
    diagnostics: list[str] = [
        f"lead_smiles={lead}",
        f"run_id={run_id}",
        f"baseline_affinity={baseline:.3f}",
        f"stop_affinity={'disabled' if stop_affinity is None else f'{stop_affinity:.3f}'}",
        f"lead_properties={properties(lead)}",
        f"python_executable={sys.executable}",
        f"meeko_script_exists={command_exists(args.meeko)} ({args.meeko})",
        f"obabel_exists={command_exists(args.obabel)} ({args.obabel})",
        f"vina_exists={command_exists(args.vina)} ({args.vina})",
        f"history_unique_smiles={len(history_seen_smiles)}",
        f"history_unique_edit_labels={len(history_seen_edit_labels)}",
        f"history_screening_rows={len(history_screening_rows)}",
        f"history_redocked_rows={len(redocked_history_rows)}",
        *lead_pose_notes,
        *history_notes,
    ]
    best_smiles = lead
    best_affinity = baseline
    if continuation_history_dirs:
        print(
            f"loaded history: {len(history_seen_smiles)} smiles, {len(history_seen_edit_labels)} edit labels",
            flush=True,
        )
    live_lead_reference: list[tuple[float, float, float]] | None = None
    if active_lead_pose:
        live_lead_modes = parse_pdbqt_modes(active_lead_pose.resolve())
        if not live_lead_modes:
            raise ValueError(f"no heavy-atom pose modes found in lead pose: {active_lead_pose}")
        live_lead_reference = live_lead_modes[0].coords
    if args.refresh_ledger_only:
        ledger_counts = refresh_ledger_from_history_only(
            outdir=outdir,
            history_dirs=continuation_history_dirs,
            lead_pose=active_lead_pose,
            receptor=args.receptor.resolve(),
            lead_smiles=lead,
            rmsd_quantile=args.rmsd_quantile,
            near_rmsd_cutoff=args.near_rmsd_cutoff,
            internal_cluster_rmsd_cutoff=args.internal_cluster_rmsd_cutoff,
            affinity_threshold=args.affinity_threshold,
            affinity_penalty_scale=args.affinity_penalty_scale,
            next_parent_count=args.next_parent_count,
            write_debug_csvs=args.write_debug_csvs,
        )
        diagnostics.extend(f"{key}={value}" for key, value in ledger_counts.items())
        print(f"refreshed compact ledger: {ledger_counts}")
        return

    seen = {lead}
    parent_file = resolve_active_parent_file(args, outdir, diagnostics)
    if parent_file:
        parents = load_parent_smiles_file(parent_file)[: max(1, args.next_parent_count)]
        if not parents:
            diagnostics.append(f"parent_smiles_file_empty={parent_file}; fallback_to_lead=true")
            parents = [("LEAD", lead)]
    else:
        parents = [("LEAD", lead)]
    original_parent_count = len(parents)
    parents = ensure_lead_parent(parents, lead, diagnostics)
    if len(parents) != original_parent_count or parents[0][0] == "LEAD_FORCED":
        print("lead_smiles_forced_into_active_parents=true", flush=True)
    rng = random.Random(args.seed + sum(ord(ch) for ch in run_id))
    edit_modes = parse_edit_modes(selected_edit_mode_text(args))
    diagnostics.append(f"initial_parent_count={len(parents)}")
    diagnostics.append(f"batch_selection={'deterministic' if args.deterministic_batch else 'randomized'}")
    diagnostics.append(f"edit_modes={','.join(edit_modes)}")
    print(f"{parent_source_label(parent_file)} parents={len(parents)}", flush=True)
    print(
        f"iteration plan: max_rounds={args.max_rounds}, batch_size={args.batch_size}, "
        f"initial_parents={len(parents)}, edit_modes={','.join(edit_modes)}",
        flush=True,
    )
    write_progress_status(
        progress_path,
        started_at=started_at,
        round_index=0,
        current_index=0,
        total_current_round=0,
        candidate=None,
        status="initializing",
        affinity=None,
        best_affinity=best_affinity,
        best_smiles=best_smiles,
        baseline=baseline,
        screened_count=len(rows),
        rejected_count=len(rejected_rows),
        rmsd_summary=None,
    )

    if args.max_rounds == 0:
        print("lead_only_dock=true; docking --lead-smiles without generating analogs", flush=True)
        lead_row = dock_single_smiles_to_ledger(
            smiles=lead,
            outdir=outdir,
            run_id=run_id,
            candidate_id="LEAD_INPUT",
            edit_label="input_lead_only",
            receptor=args.receptor.resolve(),
            box=box,
            vina=args.vina,
            meeko=args.meeko,
            obabel=args.obabel,
            seed=args.seed,
            baseline=baseline,
            num_modes=args.num_modes,
        )
        if lead_row.get("status") == "ok":
            lead_pose_path = host_path(str(lead_row.get("pose_path", "")))
            if lead_pose_path and lead_pose_path.exists():
                active_lead_pose = lead_pose_path
                diagnostics.append(f"lead_pose_source=lead_only_dock:{active_lead_pose}")
                lead_modes = parse_pdbqt_modes(active_lead_pose)
                if lead_modes:
                    live_lead_reference = lead_modes[0].coords
                    lead_row.update(
                        summarize_pose_against_lead(
                            live_lead_reference,
                            active_lead_pose,
                            args.near_rmsd_cutoff,
                            args.internal_cluster_rmsd_cutoff,
                        )
                    )
        rows.append(lead_row)
        append_attempted_smiles(outdir, lead_row)
        affinity = safe_float(lead_row.get("affinity_kcal_mol"))
        if affinity is not None and affinity < best_affinity:
            best_affinity = affinity
            best_smiles = lead
        print_progress_line(
            round_index=0,
            current_index=1,
            total_current_round=1,
            candidate_id=str(lead_row.get("candidate_id", "")),
            edit_label=str(lead_row.get("edit_label", "")),
            status=str(lead_row.get("status", "")),
            affinity=affinity,
            best_affinity=best_affinity,
            baseline=baseline,
            rmsd_summary=lead_row,
        )

    for round_index in range(1, args.max_rounds + 1):
        round_candidates: list[Candidate] = []
        for parent_id, parent_smiles in parents:
            per_parent_limit = max(1, math.ceil(args.batch_size / len(parents)))
            candidates, rejected = enumerate_round(
                parent_smiles,
                parent_id,
                round_index,
                seen,
                history_seen_smiles,
                history_seen_edit_labels,
                per_parent_limit,
                rng,
                edit_modes,
                args.deterministic_batch,
            )
            round_candidates.extend(candidates)
            rejected_rows.extend(rejected)
        round_candidates = round_candidates[: args.batch_size]
        diagnostics.append(
            f"round_{round_index}: parents={len(parents)} accepted_for_docking={len(round_candidates)} rejected_before_docking={len(rejected_rows)}"
        )
        if len(round_candidates) < args.batch_size:
            diagnostics.append(f"round_{round_index}: underfilled_batch requested={args.batch_size} accepted={len(round_candidates)}")
            print(f"round {round_index}: only found {len(round_candidates)}/{args.batch_size} new candidates after history/filtering", flush=True)
        if args.write_debug_csvs:
            try_write_csv(
                outdir / "rejected_candidates.csv",
                rejected_rows,
                fieldnames=REJECTED_FIELDNAMES,
                label="rejected_candidates.csv",
            )
        if not round_candidates:
            diagnostics.append(f"round_{round_index}: no candidates survived enumeration and filters")
            counts = rejection_reason_counts(rejected_rows)
            diagnostics.append(f"round_{round_index}_rejection_counts={counts}")
            print(f"round {round_index}: no new candidates to dock; rejection counts: {counts}", flush=True)
            break

        next_parents: list[tuple[str, str]] = []
        parent_smiles_by_id = dict(parents)
        for current_index, candidate in enumerate(round_candidates, start=1):
            write_progress_status(
                progress_path,
                started_at=started_at,
                round_index=round_index,
                current_index=current_index,
                total_current_round=len(round_candidates),
                candidate=candidate,
                status="preparing",
                affinity=None,
                best_affinity=best_affinity,
                best_smiles=best_smiles,
                baseline=baseline,
                screened_count=len(rows),
                rejected_count=len(rejected_rows),
                rmsd_summary=None,
            )
            props = properties(candidate.smiles)
            file_stem = f"{run_id}_{candidate.candidate_id}"
            sdf_path = sdf_dir / f"{file_stem}.sdf"
            pdbqt_path = pdbqt_dir / f"{file_stem}.pdbqt"
            pose_path = pose_dir / f"{file_stem}.pdbqt"
            log_path = log_dir / f"{file_stem}.log"
            prep_log_path = prep_log_dir / f"{file_stem}.log"
            status = "ok"
            affinity: float | str = ""
            reason = ""
            prep_mode = ""
            try:
                embed_sdf(candidate.smiles, candidate.candidate_id, sdf_path, args.seed + round_index)
                prep_ok, prep_output, prep_mode = prepare_pdbqt(sdf_path, pdbqt_path, args.meeko, args.obabel)
                prep_log_path.write_text(prep_output, encoding="utf-8", errors="replace")
                if not prep_ok:
                    status = "prep_failed"
                    reason = "pdbqt_preparation_failed"
                else:
                    write_progress_status(
                        progress_path,
                        started_at=started_at,
                        round_index=round_index,
                        current_index=current_index,
                        total_current_round=len(round_candidates),
                        candidate=candidate,
                        status="docking",
                        affinity=None,
                        best_affinity=best_affinity,
                        best_smiles=best_smiles,
                        baseline=baseline,
                        screened_count=len(rows),
                        rejected_count=len(rejected_rows),
                        rmsd_summary=None,
                    )
                    dock_ok, dock_affinity, _ = dock(args.vina, args.receptor.resolve(), pdbqt_path, pose_path, log_path, box, args.seed, args.num_modes)
                    if dock_ok and dock_affinity is not None:
                        affinity = dock_affinity
                        if dock_affinity < best_affinity:
                            best_affinity = dock_affinity
                            best_smiles = candidate.smiles
                        if dock_affinity <= baseline - 0.15:
                            next_parents.append((candidate.candidate_id, candidate.smiles))
                    else:
                        status = "dock_failed"
                        reason = "vina_failed_or_no_affinity"
            except Exception as exc:
                status = "failed"
                reason = str(exc)

            rmsd_summary: dict[str, str | float | int] = {
                "live_min_pose_rmsd_to_lead": "",
                "live_near_lead_mode_count": "",
                "live_mode_count": "",
                "live_inner_rmsd": "",
                "live_inner_cluster_fraction": "",
                "live_whole_rmsd": "",
            }
            if status == "ok" and live_lead_reference is not None:
                rmsd_summary = summarize_pose_against_lead(
                    live_lead_reference,
                    pose_path,
                    args.near_rmsd_cutoff,
                    args.internal_cluster_rmsd_cutoff,
                )

            row: dict[str, str | float | int] = {
                "candidate_id": candidate.candidate_id,
                "run_id": run_id,
                "round": candidate.round_index,
                "parent_id": candidate.parent_id,
                "parent_smiles": parent_smiles_by_id.get(candidate.parent_id, ""),
                "edit_label": candidate.edit_label,
                **smiles_identity_fields(candidate.smiles, candidate.smiles),
                "smiles": candidate.smiles,
                "status": status,
                "affinity_kcal_mol": affinity,
                "delta_vs_baseline": round(float(affinity) - baseline, 3) if isinstance(affinity, float) else "",
                **rmsd_summary,
                "reason": reason,
                "prep_mode": prep_mode,
                "sdf_path": str(sdf_path),
                "pdbqt_path": str(pdbqt_path),
                "pose_path": str(pose_path),
                "log_path": str(log_path),
                **props,
            }
            rows.append(row)
            append_attempted_smiles(outdir, row)
            current_affinity = affinity if isinstance(affinity, float) else None
            print_progress_line(
                round_index=round_index,
                current_index=current_index,
                total_current_round=len(round_candidates),
                candidate_id=candidate.candidate_id,
                edit_label=candidate.edit_label,
                status=status,
                affinity=current_affinity,
                best_affinity=best_affinity,
                baseline=baseline,
                rmsd_summary=rmsd_summary,
            )
            write_progress_status(
                progress_path,
                started_at=started_at,
                round_index=round_index,
                current_index=current_index,
                total_current_round=len(round_candidates),
                candidate=candidate,
                status=status,
                affinity=current_affinity,
                best_affinity=best_affinity,
                best_smiles=best_smiles,
                baseline=baseline,
                screened_count=len(rows),
                rejected_count=len(rejected_rows),
                rmsd_summary=rmsd_summary,
            )
            if args.csv_update_interval > 0 and len(rows) % args.csv_update_interval == 0:
                all_screening_rows = merge_screening_rows(history_screening_rows, rows)
                if active_lead_pose:
                    try:
                        write_ledger_outputs(
                            outdir=outdir,
                            lead_pose=active_lead_pose,
                            receptor=args.receptor.resolve(),
                            all_screening_rows=all_screening_rows,
                            lead_smiles=lead,
                            rmsd_quantile=args.rmsd_quantile,
                            near_rmsd_cutoff=args.near_rmsd_cutoff,
                            internal_cluster_rmsd_cutoff=args.internal_cluster_rmsd_cutoff,
                            affinity_threshold=args.affinity_threshold,
                            affinity_penalty_scale=args.affinity_penalty_scale,
                            next_parent_count=args.next_parent_count,
                            write_debug_csvs=args.write_debug_csvs,
                        )
                    except Exception as exc:
                        print(f"warning: periodic ledger refresh failed: {exc}", flush=True)
                else:
                    try_write_csv(outdir / "all_screening_results.csv", all_screening_rows, fieldnames=SCREENING_FIELDNAMES, label="all_screening_results.csv")
                    compact_rows = compact_molecule_rows([], all_screening_rows)
                    try_write_csv(
                        outdir / "final_ranked_candidates.csv",
                        compact_rows,
                        fieldnames=FINAL_RANK_FIELDNAMES,
                        label="final_ranked_candidates.csv",
                    )
                    try_write_csv(
                        outdir / "history_summary.csv",
                        molecule_rows_to_history_rows(compact_rows),
                        fieldnames=HISTORY_SUMMARY_FIELDNAMES,
                        label="history_summary.csv",
                    )
            if stop_affinity is not None and isinstance(affinity, float) and affinity <= stop_affinity:
                (outdir / "best_hit.smi").write_text(f"{candidate.smiles}\t{candidate.candidate_id}\t{affinity}\n", encoding="utf-8")
                if args.write_debug_csvs:
                    try_write_csv(outdir / "screening_results.csv", rows, fieldnames=SCREENING_FIELDNAMES, label="screening_results.csv")
                if active_lead_pose:
                    try:
                        all_screening_rows = merge_screening_rows(history_screening_rows, rows)
                        ledger_counts = write_ledger_outputs(
                            outdir=outdir,
                            lead_pose=active_lead_pose,
                            receptor=args.receptor.resolve(),
                            all_screening_rows=all_screening_rows,
                            lead_smiles=lead,
                            rmsd_quantile=args.rmsd_quantile,
                            near_rmsd_cutoff=args.near_rmsd_cutoff,
                            internal_cluster_rmsd_cutoff=args.internal_cluster_rmsd_cutoff,
                            affinity_threshold=args.affinity_threshold,
                            affinity_penalty_scale=args.affinity_penalty_scale,
                            next_parent_count=args.next_parent_count,
                            write_debug_csvs=args.write_debug_csvs,
                        )
                        diagnostics.extend(f"{key}={value}" for key, value in ledger_counts.items())
                    except Exception as exc:
                        print(f"warning: pose RMSD analysis failed: {exc}", flush=True)
                print(f"STOP: {candidate.candidate_id} reached {affinity:.3f} kcal/mol")
                return

        parents = sorted(next_parents, key=lambda item: next(float(row["affinity_kcal_mol"]) for row in rows if row["candidate_id"] == item[0]))[: max(1, args.next_parent_count)]
        parents = ensure_lead_parent(parents, lead, diagnostics, label="LEAD_FORCED")
        if not parents:
            parents = [("BEST_OR_LEAD", best_smiles)]

    if args.write_debug_csvs:
        try_write_csv(
            outdir / "screening_results.csv",
            rows,
            fieldnames=SCREENING_FIELDNAMES,
            label="screening_results.csv",
        )
    if args.max_rounds == 0 and not args.rescore_history:
        all_screening_rows = apply_user_nicknames(merge_screening_rows(history_screening_rows, rows), user_nickname_map)
        try_write_csv(outdir / "all_screening_results.csv", all_screening_rows, fieldnames=SCREENING_FIELDNAMES, label="all_screening_results.csv")
        diagnostics.append("lead_only_fast_update=true")
        diagnostics.append("full_history_rescore=skipped")
    elif active_lead_pose:
        try:
            all_screening_rows = merge_screening_rows(history_screening_rows, rows)
            ledger_counts = write_ledger_outputs(
                outdir=outdir,
                lead_pose=active_lead_pose,
                receptor=args.receptor.resolve(),
                all_screening_rows=all_screening_rows,
                lead_smiles=lead,
                rmsd_quantile=args.rmsd_quantile,
                near_rmsd_cutoff=args.near_rmsd_cutoff,
                internal_cluster_rmsd_cutoff=args.internal_cluster_rmsd_cutoff,
                affinity_threshold=args.affinity_threshold,
                affinity_penalty_scale=args.affinity_penalty_scale,
                next_parent_count=args.next_parent_count,
                write_debug_csvs=args.write_debug_csvs,
            )
            diagnostics.extend(f"{key}={value}" for key, value in ledger_counts.items())
        except Exception as exc:
            diagnostics.append(f"pose_rmsd_analysis_failed={exc}")
            print(f"warning: pose RMSD analysis failed: {exc}", flush=True)
    else:
        all_screening_rows = merge_screening_rows(history_screening_rows, rows)
        try_write_csv(outdir / "all_screening_results.csv", all_screening_rows, fieldnames=SCREENING_FIELDNAMES, label="all_screening_results.csv")
        compact_rows = compact_molecule_rows([], all_screening_rows)
        try_write_csv(outdir / "final_ranked_candidates.csv", compact_rows, fieldnames=FINAL_RANK_FIELDNAMES, label="final_ranked_candidates.csv")
        try_write_csv(outdir / "history_summary.csv", molecule_rows_to_history_rows(compact_rows), fieldnames=HISTORY_SUMMARY_FIELDNAMES, label="history_summary.csv")
    if args.max_rounds == 0:
        lead_affinity = rows[0].get("affinity_kcal_mol", "") if rows else ""
        write_next_parent_files(
            outdir,
            [
                {
                    "smiles": lead,
                    "nickname": "LEAD_INPUT",
                    "event_id": row_event_id(rows[0]) if rows else "LEAD_INPUT",
                    "parent_rank": 1,
                    "candidate_id": "LEAD_INPUT",
                    "selection_reason": "current_input_lead",
                    "lead_difference_rmsd": 0,
                    "inner_rmsd": rows[0].get("live_inner_rmsd", "") if rows else "",
                    "whole_rmsd": rows[0].get("live_whole_rmsd", "") if rows else "",
                    "affinity_kcal_mol": lead_affinity,
                }
            ],
        )
        print(f"current lead saved as next_round parent: {lead}", flush=True)
    diagnostics.append(f"screened_rows={len(rows)}")
    diagnostics.append(f"rejected_rows={len(rejected_rows)}")
    diagnostics.append(f"best_affinity={best_affinity:.3f}")
    (outdir / "best_hit.smi").write_text(f"{best_smiles}\tBEST\t{best_affinity}\n", encoding="utf-8")
    if args.max_rounds == 0 and rows:
        (outdir / "current_lead.smi").write_text(f"{lead}\tLEAD_INPUT\t{rows[0].get('affinity_kcal_mol', '')}\n", encoding="utf-8")
    lead_row_for_stock = next(
        (
            row
            for row in rows
            if canonical(str(row.get("smiles", "") or "")) == lead and str(row.get("status", "") or "") == "ok"
        ),
        None,
    )
    if lead_row_for_stock is not None:
        if append_successful_lead(outdir, lead, lead_row_for_stock):
            diagnostics.append("lead_csv_appended=true")
        else:
            diagnostics.append("lead_csv_duplicate_or_unwritten=true")
    write_progress_status(
        progress_path,
        started_at=started_at,
        round_index=args.max_rounds,
        current_index=len(rows),
        total_current_round=len(rows),
        candidate=None,
        status="finished",
        affinity=None,
        best_affinity=best_affinity,
        best_smiles=best_smiles,
        baseline=baseline,
        screened_count=len(rows),
        rejected_count=len(rejected_rows),
    )
    if args.max_rounds == 0 and rows:
        lead_affinity = rows[0].get("affinity_kcal_mol", "")
        print(f"current lead affinity: {lead_affinity} kcal/mol")
        print(f"current lead smiles: {lead}")
        print(f"baseline/global reference affinity: {best_affinity:.3f} kcal/mol")
    else:
        print(f"best affinity: {best_affinity:.3f} kcal/mol")
        print(f"best smiles: {best_smiles}")
    print(f"wrote {outdir / 'all_screening_results.csv'}")


if __name__ == "__main__":
    main()
