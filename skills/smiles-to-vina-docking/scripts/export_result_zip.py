#!/usr/bin/env python
"""Export a Task2-style result.zip from a docking history and one SMILES."""

from __future__ import annotations

import argparse
import csv
import os
import re
import zipfile
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors

try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError:  # pragma: no cover - depends on RDKit build
    rdMolStandardize = None  # type: ignore[assignment]


LEDGER_SOURCES = [
    "final_ranked_candidates.csv",
    "history_summary.csv",
    "all_screening_results.csv",
    "candidate_pose_summary.csv",
    "screening_results.csv",
]
SMILES_COLUMNS = ["smiles", "mol_smiles", "canonical_smiles", "input_smiles"]


@dataclass(frozen=True)
class SmilesProfile:
    smiles: str
    canonical: str
    standard_canonical: str
    inchikey: str
    connectivity_key: str
    formula: str
    fingerprint: object


@dataclass
class LedgerSelection:
    row: dict[str, str] | None
    method: str
    similarity: float | None
    target_canonical: str | None
    nearest: list[dict[str, str]]


@dataclass
class RouteGuess:
    route: str
    source: str
    check: dict[str, str | bool | int]


def clean_smiles_text(smiles: str | None) -> str:
    """Remove shell/newline whitespace that often sneaks into copied SMILES."""
    if smiles is None:
        return ""
    return re.sub(r"\s+", "", str(smiles).strip())


def mol_from_smiles(smiles: str | None) -> Chem.Mol | None:
    text = clean_smiles_text(smiles)
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


def mol_inchikey(mol: Chem.Mol) -> str:
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:
        return ""


def profile_from_smiles(smiles: str | None) -> SmilesProfile | None:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
    standard = standardized_mol(mol)
    standard_can = Chem.MolToSmiles(standard, canonical=True)
    inchikey = mol_inchikey(standard)
    return SmilesProfile(
        smiles=canonical_smiles,
        canonical=canonical_smiles,
        standard_canonical=standard_can,
        inchikey=inchikey,
        connectivity_key=inchikey.split("-", 1)[0] if inchikey else "",
        formula=rdMolDescriptors.CalcMolFormula(standard),
        fingerprint=Chem.RDKFingerprint(standard),
    )


def safe_float(value: str | int | float | None) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def candidate_paths(path_text: str | None) -> list[Path]:
    if not path_text:
        return []
    text = str(path_text).strip().strip('"')
    paths: list[Path] = []

    def add(path: Path) -> None:
        if path not in paths:
            paths.append(path)

    add(Path(text))
    if text.startswith("/mnt/") and len(text) > 6 and text[6] == "/":
        drive = text[5].upper()
        add(Path(f"{drive}:/{text[7:]}"))
    drive_match = re.match(r"^([A-Za-z]):[\\/](.*)$", text)
    if drive_match:
        drive = drive_match.group(1).lower()
        rest = drive_match.group(2).replace("\\", "/")
        add(Path(f"/mnt/{drive}/{rest}"))
    return paths


def resolve_path(path_text: str | None) -> Path | None:
    paths = candidate_paths(path_text)
    for path in paths:
        if path.exists():
            return path
    return paths[0] if paths else None


def host_path(path_text: str | None) -> Path | None:
    """Backward-compatible path resolver used by notebook callers."""
    return resolve_path(path_text)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows([{field: row.get(field, "") for field in fieldnames} for row in rows])


def row_smiles(row: dict[str, str]) -> str:
    for column in SMILES_COLUMNS:
        value = clean_smiles_text(row.get(column))
        if value:
            return value
    return ""


def split_event_id(event_id: str | None) -> tuple[str, str]:
    text = str(event_id or "")
    if ":" in text:
        run_id, candidate_id = text.rsplit(":", 1)
        return run_id, candidate_id
    return "", text


def row_run_candidate(row: dict[str, str]) -> tuple[str, str]:
    event_run, event_candidate = split_event_id(row.get("event_id"))
    return row.get("run_id") or event_run, row.get("candidate_id") or event_candidate


def row_nickname(row: dict[str, str]) -> str:
    existing = str(row.get("nickname", "") or "")
    if not existing:
        return ""
    run_id, candidate_id = row_run_candidate(row)
    machine_aliases = {
        str(row.get("mol_id", "") or ""),
        str(row.get("candidate_id", "") or ""),
        str(row.get("event_id", "") or ""),
        candidate_id,
        f"{run_id}:{candidate_id}" if run_id and candidate_id else "",
    }
    parent = str(row.get("parent", "") or row.get("parent_id", "") or "")
    edit_label = str(row.get("edit_label", "") or "")
    if parent and edit_label:
        machine_aliases.add(f"{parent}+{edit_label}")
    return "" if existing in {alias for alias in machine_aliases if alias} else existing


def row_matches_requested_id(row: dict[str, str], candidate_id: str | None, run_id: str | None) -> bool:
    if not candidate_id and not run_id:
        return True
    event_run, event_candidate = row_run_candidate(row)
    possible_ids = {
        event_candidate,
        row.get("candidate_id", ""),
        row.get("mol_id", ""),
        row.get("event_id", ""),
        row.get("nickname", ""),
    }
    if candidate_id and candidate_id not in possible_ids:
        return False
    if run_id and event_run and event_run != run_id:
        return False
    return True


def source_rank(row: dict[str, str]) -> int:
    name = Path(row.get("_source_file", "")).name
    try:
        return LEDGER_SOURCES.index(name)
    except ValueError:
        return len(LEDGER_SOURCES)


def row_sort_key(row: dict[str, str]) -> tuple[int, int, float, float]:
    score = safe_float(row.get("final_score"))
    affinity = safe_float(row.get("affinity_kcal_mol") or row.get("best_affinity_kcal_mol"))
    has_pose_or_log = 0 if row.get("pose_path") or row.get("log_path") else 1
    return (
        source_rank(row),
        has_pose_or_log,
        -(score if score is not None else -1.0),
        affinity if affinity is not None else 999.0,
    )


def ledger_candidate_rows(ledger_dir: Path) -> list[dict[str, str]]:
    if not ledger_dir.exists():
        raise FileNotFoundError(f"history directory not found: {ledger_dir}")
    rows: list[dict[str, str]] = []
    for source_name in LEDGER_SOURCES:
        csv_path = ledger_dir / source_name
        if not csv_path.exists():
            continue
        for row in read_csv_rows(csv_path):
            row["_source_file"] = str(csv_path)
            rows.append(row)
    return rows


def nearest_rows(
    target: SmilesProfile,
    profiled_rows: list[tuple[dict[str, str], SmilesProfile]],
    limit: int = 8,
) -> list[dict[str, str]]:
    nearest: list[tuple[float, dict[str, str]]] = []
    for row, profile in profiled_rows:
        similarity = DataStructs.TanimotoSimilarity(target.fingerprint, profile.fingerprint)
        nearest.append((similarity, row))
    nearest.sort(key=lambda item: (-item[0], row_sort_key(item[1])))
    out: list[dict[str, str]] = []
    for similarity, row in nearest[:limit]:
        run_id, candidate_id = row_run_candidate(row)
        out.append(
            {
                "similarity": f"{similarity:.4f}",
                "source_file": Path(row.get("_source_file", "")).name,
                "event_id": row.get("event_id", ""),
                "candidate_id": candidate_id,
                "run_id": run_id,
                "affinity_kcal_mol": row.get("affinity_kcal_mol") or row.get("best_affinity_kcal_mol") or "",
                "smiles": row_smiles(row),
            }
        )
    return out


def nearest_rows_fast(
    target_mol: Chem.Mol,
    rows: list[dict[str, str]],
    limit: int = 8,
) -> list[dict[str, str]]:
    target_fp = Chem.RDKFingerprint(target_mol)
    nearest: list[tuple[float, dict[str, str]]] = []
    for row in rows:
        mol = mol_from_smiles(row_smiles(row))
        if mol is None:
            continue
        similarity = DataStructs.TanimotoSimilarity(target_fp, Chem.RDKFingerprint(mol))
        nearest.append((similarity, row))
    nearest.sort(key=lambda item: (-item[0], row_sort_key(item[1])))
    out: list[dict[str, str]] = []
    for similarity, row in nearest[:limit]:
        run_id, candidate_id = row_run_candidate(row)
        out.append(
            {
                "similarity": f"{similarity:.4f}",
                "source_file": Path(row.get("_source_file", "")).name,
                "event_id": row.get("event_id", ""),
                "candidate_id": candidate_id,
                "run_id": run_id,
                "affinity_kcal_mol": row.get("affinity_kcal_mol") or row.get("best_affinity_kcal_mol") or "",
                "smiles": row_smiles(row),
            }
        )
    return out


def find_ledger_candidate(
    ledger_dir: Path,
    target_smiles: str | None,
    candidate_id: str | None,
    run_id: str | None,
) -> LedgerSelection:
    rows = ledger_candidate_rows(ledger_dir)
    if not rows:
        raise FileNotFoundError(f"no readable history CSVs found in {ledger_dir}")

    rows = [row for row in rows if row_matches_requested_id(row, candidate_id, run_id)]
    if candidate_id and not target_smiles:
        selected = min(rows, key=row_sort_key, default=None)
        return LedgerSelection(selected, "candidate_id", None, None, [])

    target_mol = mol_from_smiles(target_smiles) if target_smiles else None
    if target_smiles and target_mol is None:
        raise ValueError(f"invalid SMILES: {target_smiles}")
    if target_mol is None:
        return LedgerSelection(None, "no_target_smiles", None, None, [])

    target_can = Chem.MolToSmiles(target_mol, canonical=True)
    exact_matches: list[dict[str, str]] = []
    row_mol_cache: dict[str, Chem.Mol | None] = {}
    for row in rows:
        smiles = row_smiles(row)
        if not smiles:
            continue
        if smiles in row_mol_cache:
            mol = row_mol_cache[smiles]
        else:
            mol = mol_from_smiles(smiles)
            row_mol_cache[smiles] = mol
        if mol is not None and Chem.MolToSmiles(mol, canonical=True) == target_can:
            exact_matches.append(row)

    if exact_matches:
        selected = min(exact_matches, key=row_sort_key)
        return LedgerSelection(
            row=selected,
            method="exact_canonical_smiles",
            similarity=None,
            target_canonical=target_can,
            nearest=[],
        )

    target_profile = profile_from_smiles(target_can)
    if target_profile is None:
        raise ValueError(f"invalid SMILES after canonicalization: {target_smiles}")

    profile_cache: dict[str, SmilesProfile | None] = {}
    profiled_rows: list[tuple[dict[str, str], SmilesProfile]] = []
    matches: list[tuple[int, str, float | None, dict[str, str]]] = []
    for row in rows:
        smiles = row_smiles(row)
        if smiles not in profile_cache:
            profile_cache[smiles] = profile_from_smiles(smiles)
        profile = profile_cache[smiles]
        if profile is None:
            continue
        profiled_rows.append((row, profile))
        if profile.standard_canonical == target_profile.standard_canonical:
            matches.append((1, "standardized_canonical_smiles", None, row))
        elif (
            profile.connectivity_key
            and profile.connectivity_key == target_profile.connectivity_key
            and profile.formula == target_profile.formula
        ):
            matches.append((2, "same_inchikey_connectivity", None, row))
        else:
            similarity = DataStructs.TanimotoSimilarity(target_profile.fingerprint, profile.fingerprint)
            if similarity >= 0.995 and profile.formula == target_profile.formula:
                matches.append((3, "high_fingerprint_similarity", similarity, row))

    if matches:
        matches.sort(key=lambda item: (item[0], row_sort_key(item[3])))
        _, method, similarity, row = matches[0]
        return LedgerSelection(
            row=row,
            method=method,
            similarity=similarity,
            target_canonical=target_profile.canonical,
            nearest=nearest_rows(target_profile, profiled_rows),
        )

    return LedgerSelection(
        row=None,
        method="no_match",
        similarity=None,
        target_canonical=target_profile.canonical,
        nearest=nearest_rows_fast(target_mol, rows),
    )


def infer_log_path_from_pose(pose_path_text: str | None) -> Path | None:
    pose_path = resolve_path(pose_path_text)
    if not pose_path:
        return None
    candidates = [
        pose_path.with_suffix(".log"),
        pose_path.parent.parent / "logs" / f"{pose_path.stem}.log" if pose_path.parent.name == "poses" else pose_path.with_suffix(".log"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def find_log_path(ledger_dir: Path, row: dict[str, str] | None) -> Path | None:
    if not row:
        return None
    direct = resolve_path(row.get("log_path"))
    if direct and direct.exists():
        return direct
    from_pose = infer_log_path_from_pose(row.get("pose_path"))
    if from_pose and from_pose.exists():
        return from_pose

    run_id, candidate_id = row_run_candidate(row)
    stems: list[str] = []
    pose_path = resolve_path(row.get("pose_path"))
    if pose_path:
        stems.append(pose_path.stem)
    if run_id and candidate_id:
        stems.append(f"{run_id}_{candidate_id}")
    if row.get("event_id"):
        stems.append(str(row["event_id"]).replace(":", "_"))
    if candidate_id:
        stems.append(candidate_id)

    log_dirs = [ledger_dir / "logs"]
    log_dirs.extend(path for path in ledger_dir.glob("*/logs") if path.is_dir())
    for log_dir in log_dirs:
        if not log_dir.exists():
            continue
        logs = list(log_dir.glob("*.log"))
        for stem in stems:
            exact = log_dir / f"{stem}.log"
            if exact.exists():
                return exact
            for log in logs:
                if stem and stem in log.stem:
                    return log
    return None


def split_mols(side: str) -> list[str]:
    return [item.strip() for item in side.split(".") if item.strip()]


def element_counts(smiles_list: list[str]) -> Counter[str] | None:
    counts: Counter[str] = Counter()
    for smiles in smiles_list:
        mol = mol_from_smiles(smiles)
        if mol is None:
            return None
        for atom in mol.GetAtoms():
            counts[atom.GetSymbol()] += 1
    return counts


def products_covered_by_reactants(reactants: list[str], products: list[str]) -> bool:
    reactant_counts = element_counts(reactants)
    product_counts = element_counts(products)
    if reactant_counts is None or product_counts is None:
        return False
    return all(reactant_counts[element] >= count for element, count in product_counts.items())


def check_route(mol_smiles: str, route: str) -> dict[str, str | bool | int]:
    reasons: list[str] = []
    steps = [step.strip() for step in route.split(",") if step.strip()]
    if not steps:
        return {
            "route_valid": False,
            "reason": "empty_route",
            "step_count": 0,
            "all_molecules_valid": False,
            "final_product_match": False,
            "balance_ok": False,
            "no_self_loop": False,
        }

    target_can = canonical(mol_smiles)
    all_molecules_valid = True
    balance_ok = True
    no_self_loop = True
    final_product_can = ""

    for step in steps:
        if ">>" not in step:
            reasons.append("missing_reaction_arrow")
            all_molecules_valid = False
            balance_ok = False
            continue
        left, right = step.split(">>", 1)
        reactants = split_mols(left)
        products = split_mols(right)
        reactant_cans = [canonical(item) for item in reactants]
        product_cans = [canonical(item) for item in products]
        if not reactants or not products or any(item is None for item in reactant_cans + product_cans):
            reasons.append("invalid_route_molecule")
            all_molecules_valid = False
            balance_ok = False
            continue
        if any(product in reactant_cans for product in product_cans):
            reasons.append("self_loop_product_in_reactants")
            no_self_loop = False
        if not products_covered_by_reactants(reactants, products):
            reasons.append("product_elements_not_covered_by_reactants")
            balance_ok = False
        final_product_can = ".".join(product_cans)

    final_product_match = bool(target_can and final_product_can == target_can)
    if not final_product_match:
        reasons.append("final_product_mismatch")
    route_valid = all_molecules_valid and final_product_match and balance_ok and no_self_loop
    return {
        "route_valid": route_valid,
        "reason": ";".join(dict.fromkeys(reasons)) if reasons else "ok",
        "step_count": len(steps),
        "all_molecules_valid": all_molecules_valid,
        "final_product_match": final_product_match,
        "balance_ok": balance_ok,
        "no_self_loop": no_self_loop,
    }


def heavy_fragment_sizes(mol: Chem.Mol) -> list[int]:
    sizes: list[int] = []
    for fragment in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True):
        sizes.append(sum(1 for atom in fragment.GetAtoms() if atom.GetAtomicNum() > 1))
    return sizes


def route_from_reactant_mol(target_can: str, reactant_mol: Chem.Mol, source: str) -> RouteGuess | None:
    try:
        Chem.SanitizeMol(reactant_mol)
        fragments = Chem.GetMolFrags(reactant_mol, asMols=True, sanitizeFrags=True)
    except Exception:
        return None
    reactants = [Chem.MolToSmiles(fragment, canonical=True) for fragment in fragments if fragment.GetNumAtoms() > 0]
    if len(reactants) < 2:
        return None
    route = f"{'.'.join(reactants)}>>{target_can}"
    check = check_route(target_can, route)
    if check.get("route_valid") is True:
        return RouteGuess(route=route, source=source, check=check)
    return None


def infer_amide_route(mol: Chem.Mol, target_can: str) -> RouteGuess | None:
    pattern = Chem.MolFromSmarts("[C:1](=[O:2])-[N:3]")
    if pattern is None:
        return None
    guesses: list[tuple[tuple[int, int], RouteGuess]] = []
    for carbon_idx, _oxygen_idx, nitrogen_idx in mol.GetSubstructMatches(pattern):
        carbon = mol.GetAtomWithIdx(carbon_idx)
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        n_extra = [
            atom.GetIdx()
            for atom in nitrogen.GetNeighbors()
            if atom.GetIdx() != carbon_idx and atom.GetAtomicNum() > 1
        ]
        c_extra = [
            atom.GetIdx()
            for atom in carbon.GetNeighbors()
            if atom.GetIdx() not in {nitrogen_idx} and atom.GetAtomicNum() > 1 and atom.GetSymbol() != "O"
        ]
        if not n_extra or not c_extra:
            continue
        rw = Chem.RWMol(mol)
        if rw.GetBondBetweenAtoms(carbon_idx, nitrogen_idx) is None:
            continue
        rw.RemoveBond(carbon_idx, nitrogen_idx)
        hydroxyl_idx = rw.AddAtom(Chem.Atom("O"))
        rw.AddBond(carbon_idx, hydroxyl_idx, Chem.BondType.SINGLE)
        reactant_mol = rw.GetMol()
        guess = route_from_reactant_mol(target_can, reactant_mol, "retrosynthetic_amide_coupling")
        if guess is None:
            continue
        sizes = sorted(heavy_fragment_sizes(reactant_mol), reverse=True)
        if len(sizes) < 2:
            continue
        guesses.append(((min(sizes[:2]), -abs(sizes[0] - sizes[1])), guess))
    if not guesses:
        return None
    guesses.sort(key=lambda item: item[0], reverse=True)
    return guesses[0][1]


def infer_sulfonamide_route(mol: Chem.Mol, target_can: str) -> RouteGuess | None:
    pattern = Chem.MolFromSmarts("[S:1](=[O:2])(=[O:3])-[N:4]")
    if pattern is None:
        return None
    for sulfur_idx, _o1_idx, _o2_idx, nitrogen_idx in mol.GetSubstructMatches(pattern):
        rw = Chem.RWMol(mol)
        if rw.GetBondBetweenAtoms(sulfur_idx, nitrogen_idx) is None:
            continue
        rw.RemoveBond(sulfur_idx, nitrogen_idx)
        chlorine_idx = rw.AddAtom(Chem.Atom("Cl"))
        rw.AddBond(sulfur_idx, chlorine_idx, Chem.BondType.SINGLE)
        guess = route_from_reactant_mol(target_can, rw.GetMol(), "retrosynthetic_sulfonamide")
        if guess is not None:
            return guess
    return None


def infer_route(smiles: str) -> RouteGuess:
    target_can = canonical(smiles)
    mol = mol_from_smiles(smiles)
    if target_can is None or mol is None:
        route = default_route(smiles)
        return RouteGuess(route, "placeholder_invalid_smiles", check_route(smiles, route))

    for infer in (infer_amide_route, infer_sulfonamide_route):
        guess = infer(mol, target_can)
        if guess is not None:
            return guess

    route = default_route(target_can)
    return RouteGuess(route, "placeholder_no_supported_disconnection", check_route(target_can, route))


def default_route(smiles: str) -> str:
    return f"PROPOSE_ROUTE_MANUALLY>>{smiles}"


def export_result_zip(
    *,
    ledger_dir: Path,
    smiles: str | None,
    candidate_id: str | None,
    run_id: str | None,
    route: str | None,
    result_dir: Path,
    zip_name: str,
) -> Path:
    selection = find_ledger_candidate(ledger_dir, smiles, candidate_id, run_id)
    selected = selection.row
    selected_smiles = row_smiles(selected) if selected else clean_smiles_text(smiles)
    if not selected_smiles:
        raise ValueError("provide --smiles or a --candidate-id found in the history")
    mol_smiles = canonical(selected_smiles)
    if mol_smiles is None:
        raise ValueError(f"invalid selected SMILES: {selected_smiles}")

    if route and route.strip():
        route_text = route.strip()
        route_guess = RouteGuess(route_text, "user_supplied", check_route(mol_smiles, route_text))
    else:
        route_guess = infer_route(mol_smiles)
        route_text = route_guess.route

    result_dir.mkdir(parents=True, exist_ok=True)
    result_csv = result_dir / "result.csv"
    result_log = result_dir / "result.log"
    write_csv(result_csv, [{"mol_smiles": mol_smiles, "route": route_text}], ["mol_smiles", "route"])

    log_text = ""
    log_path = find_log_path(ledger_dir, selected)
    pose_path = resolve_path(selected.get("pose_path")) if selected else None
    if log_path and log_path.exists():
        log_text = log_path.read_text(encoding="utf-8", errors="replace")

    lines = [
        "AI4S Task2 result export log",
        f"timestamp_local: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"history_dir: {ledger_dir}",
        "",
        "[Selected Molecule]",
        f"input_smiles: {clean_smiles_text(smiles)}",
        f"mol_smiles: {mol_smiles}",
        f"route: {route_text}",
        f"route_source: {route_guess.source}",
        f"route_valid: {route_guess.check.get('route_valid')}",
        f"route_check_reason: {route_guess.check.get('reason')}",
    ]
    if route_text.startswith("PROPOSE_ROUTE_MANUALLY>>"):
        lines.append("route_warning: placeholder route; provide or verify a chemically realistic route before official submission.")

    if selected:
        event_run, event_candidate = row_run_candidate(selected)
        lines.extend(
            [
                "",
                "[History Match]",
                f"match_method: {selection.method}",
                f"match_similarity: {selection.similarity if selection.similarity is not None else ''}",
                f"target_canonical_smiles: {selection.target_canonical or ''}",
                f"source_file: {selected.get('_source_file', '')}",
                f"event_id: {selected.get('event_id', '')}",
                f"candidate_id: {event_candidate}",
                f"run_id: {event_run}",
                f"mol_id: {selected.get('mol_id', '')}",
                f"nickname: {row_nickname(selected)}",
                f"edit_label: {selected.get('edit_label', '')}",
                f"affinity_kcal_mol: {selected.get('affinity_kcal_mol') or selected.get('best_affinity_kcal_mol') or ''}",
                f"final_score: {selected.get('final_score', '')}",
                f"lead_difference_rmsd: {selected.get('lead_difference_rmsd') or selected.get('best_rmsd_to_lead') or ''}",
                f"inner_rmsd: {selected.get('inner_rmsd') or selected.get('live_inner_rmsd') or ''}",
                f"pose_path: {pose_path or selected.get('pose_path', '')}",
                f"log_path: {log_path or selected.get('log_path', '')}",
            ]
        )
    else:
        lines.extend(
            [
                "",
                "[History Match]",
                f"match_method: {selection.method}",
                f"target_canonical_smiles: {selection.target_canonical or ''}",
                "status: no exact/standardized/InChIKey match found in history",
            ]
        )

    if selection.nearest:
        lines.extend(["", "[Nearest History Candidates]"])
        for idx, near in enumerate(selection.nearest, start=1):
            lines.append(
                " | ".join(
                    [
                        f"{idx}",
                        f"similarity={near['similarity']}",
                        f"source={near['source_file']}",
                        f"event_id={near['event_id']}",
                        f"candidate_id={near['candidate_id']}",
                        f"affinity={near['affinity_kcal_mol']}",
                        f"smiles={near['smiles']}",
                    ]
                )
            )

    lines.extend(["", "[Vina Log]", log_text.strip() if log_text.strip() else "No Vina log found for this molecule."])
    result_log.write_text("\n".join(lines) + "\n", encoding="utf-8")

    zip_path = result_dir / zip_name
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.write(result_csv, "result.csv")
        archive.write(result_log, "result.log")
    return zip_path


def export_smiles_result_zip(
    smiles: str,
    ledger_dir: str | Path,
    outdir: str | Path | None = None,
    route: str | None = None,
    zip_name: str = "result.zip",
) -> Path:
    """Create result.zip for one SMILES from an existing docking history."""
    ledger_path = Path(ledger_dir).resolve()
    result_dir = Path(outdir).resolve() if outdir else ledger_path / "result_export"
    return export_result_zip(
        ledger_dir=ledger_path,
        smiles=smiles,
        candidate_id=None,
        run_id=None,
        route=route,
        result_dir=result_dir,
        zip_name=zip_name,
    )


def export_candidate_result_zip(
    candidate_id: str,
    ledger_dir: str | Path,
    outdir: str | Path | None = None,
    run_id: str | None = None,
    route: str | None = None,
    zip_name: str = "result.zip",
) -> Path:
    """Create result.zip for one candidate_id from an existing docking history."""
    ledger_path = Path(ledger_dir).resolve()
    result_dir = Path(outdir).resolve() if outdir else ledger_path / "result_export"
    return export_result_zip(
        ledger_dir=ledger_path,
        smiles=None,
        candidate_id=candidate_id,
        run_id=run_id,
        route=route,
        result_dir=result_dir,
        zip_name=zip_name,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--history-dir", "--ledger-dir", dest="ledger_dir", metavar="HISTORY_DIR", type=Path, required=True, help="Docking history directory")
    parser.add_argument("--smiles", help="SMILES to export; matched by canonical/standardized identity when present in the history")
    parser.add_argument("--candidate-id", help="Candidate id, mol_id, event_id, or nickname to export instead of matching only by SMILES")
    parser.add_argument("--run-id", help="Optional run id disambiguator for --candidate-id")
    parser.add_argument("--route", help="Reaction route SMILES for result.csv; if omitted, a simple RDKit retrosynthetic route is attempted")
    parser.add_argument("--outdir", type=Path, help="Output directory; defaults to <history-dir>/result_export")
    parser.add_argument("--zip-name", default="result.zip")
    args = parser.parse_args()

    ledger_dir = args.ledger_dir.resolve()
    outdir = args.outdir.resolve() if args.outdir else ledger_dir / "result_export"
    zip_path = export_result_zip(
        ledger_dir=ledger_dir,
        smiles=args.smiles,
        candidate_id=args.candidate_id,
        run_id=args.run_id,
        route=args.route,
        result_dir=outdir,
        zip_name=args.zip_name,
    )
    print(f"wrote {zip_path}")
    print(f"wrote {outdir / 'result.csv'}")
    print(f"wrote {outdir / 'result.log'}")


if __name__ == "__main__":
    main()
