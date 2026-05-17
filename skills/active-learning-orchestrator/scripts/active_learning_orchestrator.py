#!/usr/bin/env python
"""Orchestrate an active-learning loop across the docking/refinement skills.

This is intentionally a glue layer, not a fourth mega-skill. It calls existing
skill scripts, then builds frontier/gap CSVs that a human or the refiner can use
for the next iteration.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import importlib.util
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_PYTHON = "/mnt/l/WSL/softwares/conda_envs/md/bin/python"
DEFAULT_PROJECT = Path("~/vina_task2")
DEFAULT_POSE_ANALYZER = Path("/mnt/e/Codex/skills/pose-analyzer")
DEFAULT_REFINER = Path("/mnt/e/Codex/skills/druglike-pocket-refiner")
SMILES_DOCK_UTILS = Path("/mnt/e/Codex/skills/smiles-to-vina-docking/scripts/dock_utils.py")
_DOCK_UTILS_MODULE = None

ID_FIELDS = {
    "seq_id",
    "nickname",
    "score_set",
    "input_smiles",
    "canonical_smiles",
    "smiles",
    "murcko_scaffold",
}
EXCLUDED_MODEL_FIELDS = {
    "official_binding_score",
    "predicted_binding_score",
    "prediction_residual",
    "pure_smiles_pc1",
    "pure_smiles_pc2",
    "structural_interaction_pc1",
    "structural_interaction_pc2",
}


def read_csv(path: Path) -> list[dict[str, str]]:
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


def is_finite_number(value: object) -> bool:
    parsed = safe_float(value)
    return not math.isnan(parsed) and math.isfinite(parsed)


def clamp01(value: float) -> float:
    if math.isnan(value):
        return 0.0
    return max(0.0, min(1.0, value))


def percentile(values: list[float], p: float) -> float:
    finite = sorted(v for v in values if not math.isnan(v))
    if not finite:
        return math.nan
    idx = max(0, min(len(finite) - 1, round((len(finite) - 1) * p)))
    return finite[idx]


def robust_scale(value: float, lo: float, hi: float) -> float:
    if math.isnan(value) or math.isnan(lo) or math.isnan(hi) or abs(hi - lo) < 1e-12:
        return 0.5
    return clamp01((value - lo) / (hi - lo))


def euclidean(a: list[float], b: list[float]) -> float:
    total = 0.0
    for x, y in zip(a, b):
        total += (x - y) ** 2
    return math.sqrt(total)


def cosine_distance(a: list[float], b: list[float]) -> float:
    dot = 0.0
    norm_a = 0.0
    norm_b = 0.0
    for x, y in zip(a, b):
        dot += x * y
        norm_a += x * x
        norm_b += y * y
    if norm_a <= 1e-12 or norm_b <= 1e-12:
        return 1.0
    return 1.0 - (dot / (math.sqrt(norm_a) * math.sqrt(norm_b)))


def vector_distance(a: list[float], b: list[float], metric: str) -> float:
    if metric == "cosine":
        return cosine_distance(a, b)
    return euclidean(a, b)


def run_command(command: list[str], dry_run: bool) -> None:
    print("+ " + " ".join(command), flush=True)
    if dry_run:
        return
    subprocess.run(command, check=True)


def run_command_capture(command: list[str], dry_run: bool) -> subprocess.CompletedProcess[str] | None:
    print("+ " + " ".join(command), flush=True)
    if dry_run:
        return None
    return subprocess.run(command, check=True, text=True, capture_output=True)


def shell_path(path: Path) -> str:
    return str(path.expanduser())


def load_dock_utils():
    global _DOCK_UTILS_MODULE
    if _DOCK_UTILS_MODULE is not None:
        return _DOCK_UTILS_MODULE
    spec = importlib.util.spec_from_file_location("active_learning_dock_utils", SMILES_DOCK_UTILS)
    if spec is None or spec.loader is None:
        return None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    _DOCK_UTILS_MODULE = module
    return module


def canonical_smiles(smiles: str) -> str:
    smiles = str(smiles or "").strip()
    if not smiles:
        return ""
    du = load_dock_utils()
    if du is None:
        return smiles
    try:
        return du.canonical(smiles) or ""
    except Exception:
        return ""


def numeric_feature_fields(rows: list[dict[str, str]]) -> list[str]:
    fields: list[str] = []
    if not rows:
        return fields
    for field in rows[0].keys():
        if field in ID_FIELDS or field in EXCLUDED_MODEL_FIELDS:
            continue
        values = [safe_float(row.get(field)) for row in rows]
        finite = [v for v in values if not math.isnan(v)]
        if len(finite) >= max(5, len(rows) // 20):
            fields.append(field)
    return fields


def scaled_vectors(rows: list[dict[str, str]], fields: list[str]) -> list[list[float]]:
    ranges: dict[str, tuple[float, float]] = {}
    for field in fields:
        values = [safe_float(row.get(field)) for row in rows]
        lo = percentile(values, 0.05)
        hi = percentile(values, 0.95)
        ranges[field] = (lo, hi)
    out: list[list[float]] = []
    for row in rows:
        out.append([robust_scale(safe_float(row.get(field)), *ranges[field]) for field in fields])
    return out


def mean_k_distance(index: int, vectors: list[list[float]], k: int, metric: str = "euclidean") -> float:
    if len(vectors) <= 1:
        return 0.0
    dists = [vector_distance(vectors[index], other, metric) for j, other in enumerate(vectors) if j != index]
    dists.sort()
    use = dists[: max(1, min(k, len(dists)))]
    return sum(use) / len(use)


def nearest_scored_distance(index: int, vectors: list[list[float]], scored_indexes: set[int], metric: str = "euclidean") -> float:
    if not scored_indexes:
        return 0.0
    dists = [vector_distance(vectors[index], vectors[j], metric) for j in scored_indexes if j != index]
    if not dists:
        return 0.0
    return min(dists)


def druglike_score(row: dict[str, str]) -> float:
    qed = safe_float(row.get("qed"))
    logp = safe_float(row.get("logp"))
    tpsa = safe_float(row.get("tpsa"))
    rot = safe_float(row.get("rot_bonds"))
    score = 0.0
    score += 0.45 * clamp01(qed)
    score += 0.20 * (1.0 - clamp01(abs((logp if not math.isnan(logp) else 2.5) - 2.5) / 4.0))
    score += 0.20 * (1.0 - clamp01(abs((tpsa if not math.isnan(tpsa) else 90.0) - 90.0) / 120.0))
    score += 0.15 * (1.0 - clamp01((rot if not math.isnan(rot) else 8.0) / 16.0))
    return clamp01(score)


def quality_score(row: dict[str, str], official_min: float, official_max: float) -> float:
    official = safe_float(row.get("official_binding_score"))
    pred = safe_float(row.get("predicted_binding_score"))
    affinity = safe_float(row.get("affinity_kcal_mol"))
    if not math.isnan(official):
        score = robust_scale(official, official_min, official_max)
    elif not math.isnan(pred):
        score = 0.55 * robust_scale(pred, official_min, official_max)
    else:
        score = 0.0
    if not math.isnan(affinity):
        # Vina is better when more negative. Map -7..-13 roughly to 0..1.
        score = 0.75 * score + 0.25 * clamp01((-affinity - 7.0) / 6.0)
    return clamp01(score)


def selection_origin(row: dict[str, str]) -> str:
    edit_label = str(row.get("edit_label", "") or "").lower()
    nickname = str(row.get("nickname", "") or "").lower()
    refinement = str(row.get("druglike_refinement_score", "") or "").strip()
    if "refine" in edit_label or "druglike_refine" in edit_label or nickname.startswith("refine") or refinement:
        return "iteration_generated"
    return "history_existing"


def selection_reason(
    row: dict[str, str],
    novelty: float,
    uncertainty: float,
    quality: float,
    druglike: float,
    chembl: float,
    contact: float,
    metric: str,
) -> str:
    parts = [
        f"origin={selection_origin(row)}",
        f"distance_metric={metric}",
        f"novelty={novelty:.3f}",
        f"uncertainty={uncertainty:.3f}",
        f"quality={quality:.3f}",
        f"druglike={druglike:.3f}",
        f"chembl={chembl:.3f}",
        f"contact={contact:.3f}",
    ]
    if is_finite_number(row.get("affinity_kcal_mol")):
        parts.append(f"affinity={safe_float(row.get('affinity_kcal_mol')):.2f}")
    return ";".join(parts)


def count_successful_docked(history_csv: Path) -> int:
    if not history_csv.exists():
        return 0
    rows = read_csv(history_csv)
    return sum(1 for row in rows if is_finite_number(row.get("affinity_kcal_mol")))


def adaptive_weights(iteration: int, docked_count: int, args: argparse.Namespace) -> dict[str, float]:
    if not args.auto_weights:
        return {
            "novelty": args.weight_novelty,
            "uncertainty": args.weight_uncertainty,
            "quality": args.weight_quality,
            "druglike": args.weight_druglike,
            "chembl": args.weight_chembl,
            "contact": args.weight_contact,
            "phase": "manual",
        }
    if iteration <= 3:
        weights = {
            "novelty": 0.38,
            "uncertainty": 0.24,
            "quality": 0.10,
            "druglike": 0.12,
            "chembl": 0.08,
            "contact": 0.08,
            "phase": "early_exploration",
        }
    elif docked_count > 200:
        weights = {
            "novelty": 0.18,
            "uncertainty": 0.30,
            "quality": 0.26,
            "druglike": 0.10,
            "chembl": 0.08,
            "contact": 0.08,
            "phase": "late_exploitation",
        }
    else:
        weights = {
            "novelty": 0.28,
            "uncertainty": 0.24,
            "quality": 0.18,
            "druglike": 0.12,
            "chembl": 0.09,
            "contact": 0.09,
            "phase": "balanced",
        }
    return weights


def build_frontier(args: argparse.Namespace) -> None:
    feature_matrix = args.feature_matrix.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    rows = read_csv(feature_matrix)
    if not rows:
        raise SystemExit(f"No rows found in {feature_matrix}")
    fields = numeric_feature_fields(rows)
    vectors = scaled_vectors(rows, fields)
    metric = getattr(args, "distance_metric", "cosine")
    scored_indexes = {i for i, row in enumerate(rows) if not math.isnan(safe_float(row.get("official_binding_score")))}
    official_values = [safe_float(row.get("official_binding_score")) for row in rows]
    official_min = percentile(official_values, 0.02)
    official_max = percentile(official_values, 0.98)
    all_novelty = [mean_k_distance(i, vectors, args.neighbor_k, metric) for i in range(len(rows))]
    nov_lo = percentile(all_novelty, 0.05)
    nov_hi = percentile(all_novelty, 0.95)
    scored_dists = [nearest_scored_distance(i, vectors, scored_indexes, metric) for i in range(len(rows))]
    sd_lo = percentile(scored_dists, 0.05)
    sd_hi = percentile(scored_dists, 0.95)
    ranked: list[dict[str, object]] = []
    for i, row in enumerate(rows):
        novelty = robust_scale(all_novelty[i], nov_lo, nov_hi)
        uncertainty = robust_scale(scored_dists[i], sd_lo, sd_hi)
        q = quality_score(row, official_min, official_max)
        d = druglike_score(row)
        chembl = clamp01(safe_float(row.get("chembl_scaffold_similarity")))
        contact = clamp01(safe_float(row.get("surface_contact_fraction_4a")))
        if contact == 0.0:
            contact = clamp01(safe_float(row.get("ligand_atom_contact_fraction_4a")))
        frontier = (
            args.weight_novelty * novelty
            + args.weight_uncertainty * uncertainty
            + args.weight_quality * q
            + args.weight_druglike * d
            + args.weight_chembl * chembl
            + args.weight_contact * contact
        )
        ranked.append(
            {
                "seq_id": row.get("seq_id", ""),
                "nickname": row.get("nickname", ""),
                "canonical_smiles": row.get("canonical_smiles", ""),
                "official_binding_score": row.get("official_binding_score", ""),
                "predicted_binding_score": row.get("predicted_binding_score", ""),
                "affinity_kcal_mol": row.get("affinity_kcal_mol", ""),
                "pure_smiles_pc1": row.get("pure_smiles_pc1", ""),
                "pure_smiles_pc2": row.get("pure_smiles_pc2", ""),
                "structural_interaction_pc1": row.get("structural_interaction_pc1", ""),
                "structural_interaction_pc2": row.get("structural_interaction_pc2", ""),
                "qed": row.get("qed", ""),
                "chembl_scaffold_similarity": row.get("chembl_scaffold_similarity", ""),
                "surface_contact_fraction_4a": row.get("surface_contact_fraction_4a", ""),
                "novelty_score": round(novelty, 6),
                "uncertainty_proxy": round(uncertainty, 6),
                "quality_score": round(q, 6),
                "druglike_score": round(d, 6),
                "contact_score": round(contact, 6),
                "frontier_score": round(frontier, 6),
                "selection_origin": selection_origin(row),
                "selection_reason": selection_reason(row, novelty, uncertainty, q, d, chembl, contact, metric),
            }
        )
    ranked.sort(key=lambda item: safe_float(item.get("frontier_score")), reverse=True)
    for rank, item in enumerate(ranked, start=1):
        item["active_learning_rank"] = rank
        item["needs_official_binding_score"] = "yes" if str(item.get("official_binding_score", "")).strip() in {"", "nan"} else "no"
    frontier_rows = ranked[: args.top_frontier_seeds]
    official_candidates = [
        row for row in ranked
        if row.get("needs_official_binding_score") == "yes" and is_finite_number(row.get("affinity_kcal_mol"))
    ]
    generated_candidates = [row for row in official_candidates if row.get("selection_origin") == "iteration_generated"]
    if len(generated_candidates) >= args.top_official_recommendations:
        official_rows = generated_candidates[: args.top_official_recommendations]
    else:
        official_rows = [*generated_candidates]
        seen_ids = {str(row.get("seq_id", "")) for row in official_rows}
        for row in official_candidates:
            if str(row.get("seq_id", "")) in seen_ids:
                continue
            official_rows.append(row)
            if len(official_rows) >= args.top_official_recommendations:
                break
    fields_out = [
        "active_learning_rank",
        "seq_id",
        "nickname",
        "canonical_smiles",
        "needs_official_binding_score",
        "official_binding_score",
        "predicted_binding_score",
        "affinity_kcal_mol",
        "pure_smiles_pc1",
        "pure_smiles_pc2",
        "structural_interaction_pc1",
        "structural_interaction_pc2",
        "qed",
        "chembl_scaffold_similarity",
        "surface_contact_fraction_4a",
        "novelty_score",
        "uncertainty_proxy",
        "quality_score",
        "druglike_score",
        "contact_score",
        "frontier_score",
        "selection_origin",
        "selection_reason",
    ]
    write_csv(outdir / "frontier_seeds.csv", frontier_rows, fields_out)
    write_csv(outdir / "official_score_recommendations.csv", official_rows, fields_out)
    write_csv(outdir / "active_learning_shortlist.csv", ranked[: args.top_shortlist], fields_out)
    with (outdir / "active_learning_shortlist.smi").open("w", encoding="utf-8") as handle:
        for row in ranked[: args.top_shortlist]:
            smiles = str(row.get("canonical_smiles", "")).strip()
            if smiles:
                handle.write(f"{smiles}\t{row.get('seq_id', '')}\n")
    build_gap_maps(rows, outdir, args.grid_bins)
    dashboard_dir = getattr(args, "dashboard_dir", None)
    if dashboard_dir:
        export_dashboard_recommendations(Path(dashboard_dir), frontier_rows, ranked[: args.top_shortlist], official_rows)
    print(f"frontier_seeds={outdir / 'frontier_seeds.csv'} rows={len(frontier_rows)}")
    print(f"official_score_recommendations={outdir / 'official_score_recommendations.csv'} rows={len(official_rows)}")
    print(f"shortlist={outdir / 'active_learning_shortlist.csv'} rows={min(args.top_shortlist, len(ranked))}")


def export_dashboard_recommendations(dashboard_dir: Path, frontier_rows: list[dict[str, object]], shortlist_rows: list[dict[str, object]], official_rows: list[dict[str, object]]) -> None:
    data_dir = dashboard_dir.expanduser().resolve() / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    frontier_ids = {str(row.get("seq_id", "")) for row in frontier_rows if row.get("seq_id")}
    official_ids = {str(row.get("seq_id", "")) for row in official_rows if row.get("seq_id")}
    payload = {
        "frontier_seeds": frontier_rows,
        "official_score_recommendations": official_rows,
        "shortlist": shortlist_rows,
        "frontier_seq_ids": sorted(frontier_ids),
        "official_score_seq_ids": sorted(official_ids),
        "purpose": "Recommended molecules/seeds for the next official binding score cycle.",
    }
    with (data_dir / "active_learning_recommendations.json").open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False)


def deduplicate_refiner_outputs(outdir: Path) -> dict[str, int]:
    outdir = outdir.expanduser().resolve()
    files = [
        "druglike_refinement_candidates.csv",
        "druglike_refinement_ranked.csv",
        "druglike_refinement_docked_ranked.csv",
    ]
    summary = {"files_checked": 0, "rows_before": 0, "rows_after": 0, "duplicates_removed": 0}
    for name in files:
        path = outdir / name
        if not path.exists():
            continue
        rows = read_csv(path)
        if not rows:
            continue
        fields = list(rows[0].keys())
        seen: set[str] = set()
        kept: list[dict[str, object]] = []
        removed = 0
        for row in rows:
            smiles = row.get("candidate_smiles") or row.get("canonical_smiles") or row.get("smiles")
            canonical = canonical_smiles(str(smiles or ""))
            key = canonical or str(smiles or "").strip()
            if not key:
                kept.append(row)
                continue
            if key in seen:
                removed += 1
                continue
            seen.add(key)
            kept.append(row)
        if removed:
            write_csv(path, kept, fields)
        summary["files_checked"] += 1
        summary["rows_before"] += len(rows)
        summary["rows_after"] += len(kept)
        summary["duplicates_removed"] += removed
    summary_path = outdir / "deduplication_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(
        f"deduplicated_refiner_outputs={outdir} duplicates_removed={summary['duplicates_removed']} "
        f"rows_before={summary['rows_before']} rows_after={summary['rows_after']}",
        flush=True,
    )
    return summary


def build_gap_maps(rows: list[dict[str, str]], outdir: Path, bins: int) -> None:
    for prefix, xfield, yfield in (
        ("pure_smiles", "pure_smiles_pc1", "pure_smiles_pc2"),
        ("structural_interaction", "structural_interaction_pc1", "structural_interaction_pc2"),
    ):
        xs = [safe_float(row.get(xfield)) for row in rows]
        ys = [safe_float(row.get(yfield)) for row in rows]
        finite = [(x, y) for x, y in zip(xs, ys) if not math.isnan(x) and not math.isnan(y)]
        if not finite:
            continue
        min_x, max_x = min(x for x, _ in finite), max(x for x, _ in finite)
        min_y, max_y = min(y for _, y in finite), max(y for _, y in finite)
        cells: dict[tuple[int, int], list[dict[str, str]]] = {}
        for row, x, y in zip(rows, xs, ys):
            if math.isnan(x) or math.isnan(y):
                continue
            bx = min(bins - 1, max(0, int((x - min_x) / (max_x - min_x or 1.0) * bins)))
            by = min(bins - 1, max(0, int((y - min_y) / (max_y - min_y or 1.0) * bins)))
            cells.setdefault((bx, by), []).append(row)
        out: list[dict[str, object]] = []
        for bx in range(bins):
            for by in range(bins):
                members = cells.get((bx, by), [])
                scored = [safe_float(m.get("official_binding_score")) for m in members]
                scored = [v for v in scored if not math.isnan(v)]
                preds = [safe_float(m.get("predicted_binding_score")) for m in members]
                preds = [v for v in preds if not math.isnan(v)]
                count = len(members)
                priority = (1.0 / (1.0 + count)) + (0.25 if scored else 0.0)
                out.append(
                    {
                        "x_bin": bx,
                        "y_bin": by,
                        "molecule_count": count,
                        "scored_count": len(scored),
                        "mean_official_binding_score": round(sum(scored) / len(scored), 6) if scored else "",
                        "mean_predicted_binding_score": round(sum(preds) / len(preds), 6) if preds else "",
                        "gap_priority": round(priority, 6),
                    }
                )
        out.sort(key=lambda item: safe_float(item.get("gap_priority")), reverse=True)
        write_csv(
            outdir / f"gap_regions_{prefix}.csv",
            out,
            [
                "x_bin",
                "y_bin",
                "molecule_count",
                "scored_count",
                "mean_official_binding_score",
                "mean_predicted_binding_score",
                "gap_priority",
            ],
        )


def frontier_weight_namespace(iteration: int, docked_count: int, args: argparse.Namespace) -> dict[str, object]:
    weights = adaptive_weights(iteration, docked_count, args)
    return {
        "weight_novelty": weights["novelty"],
        "weight_uncertainty": weights["uncertainty"],
        "weight_quality": weights["quality"],
        "weight_druglike": weights["druglike"],
        "weight_chembl": weights["chembl"],
        "weight_contact": weights["contact"],
        "weight_phase": weights["phase"],
        "docked_count": docked_count,
    }


def read_metrics(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


@dataclass
class Paths:
    project: Path
    python: str
    pose_analyzer: Path
    refiner: Path

    @property
    def history_csv(self) -> Path:
        return self.project / "dock_history" / "dock_history.csv"

    @property
    def history_dir(self) -> Path:
        return self.project / "dock_history"

    @property
    def analysis_csv(self) -> Path:
        return self.project / "dock_history" / "binding_score_analysis.csv"

    @property
    def reference_csv(self) -> Path:
        return self.project / "references" / "chembl_kinase_smiles.csv"

    @property
    def score_dir(self) -> Path:
        return self.project / "score_space_model"

    @property
    def dashboard_dir(self) -> Path:
        return self.project / "score_space_dashboard"

    @property
    def active_dir(self) -> Path:
        return self.project / "active_learning"

    @property
    def receptor(self) -> Path:
        return self.project / "vina_bin" / "target.pdbqt"

    @property
    def config(self) -> Path:
        return self.project / "vina_bin" / "idz6F_config.txt"

    @property
    def vina(self) -> Path:
        return self.project / "vina_bin" / "vina"


def select_seed_smiles(frontier_csv: Path, max_count: int) -> list[str]:
    rows = read_csv(frontier_csv)
    smiles: list[str] = []
    seen: set[str] = set()
    for row in rows:
        value = (row.get("canonical_smiles") or "").strip()
        if not value or value in seen:
            continue
        seen.add(value)
        smiles.append(value)
        if len(smiles) >= max_count:
            break
    return smiles


def run_refiner_generation(paths: Paths, args: argparse.Namespace, iter_dir: Path) -> None:
    frontier_csv = iter_dir / "frontier_seeds.csv"
    seed_smiles = select_seed_smiles(frontier_csv, args.frontier_seeds_to_refine)
    if not seed_smiles:
        print(f"no frontier seed SMILES found in {frontier_csv}; skipping generation")
        return
    outdir = iter_dir / "generated_refinement"
    command = [
        paths.python,
        shell_path(paths.refiner / "scripts" / "refine_druglike_candidates.py"),
        "--history-dir",
        shell_path(paths.history_dir),
        "--reference-smiles-csv",
        shell_path(paths.reference_csv),
        "--outdir",
        shell_path(outdir),
        "--receptor",
        shell_path(args.receptor.expanduser() if args.receptor else paths.receptor),
        "--config",
        shell_path(args.config.expanduser() if args.config else paths.config),
        "--vina",
        shell_path(args.vina.expanduser() if args.vina else paths.vina),
        "--max-rounds",
        str(args.refiner_rounds),
        "--batch-size",
        str(args.refiner_batch_size),
        "--beam-size",
        str(args.refiner_beam_size),
        "--target-refinement-score",
        str(args.target_refinement_score),
        "--min-generated-before-target-stop",
        str(args.min_generated_before_target_stop),
        "--dock-workers",
        str(args.dock_workers),
    ]
    for smiles in seed_smiles:
        command.extend(["--expert-smiles", smiles])
    if args.dock_generated:
        command.extend(["--dock-top-candidates", str(args.dock_top_candidates)])
    if args.allow_add:
        command.append("--allow-add")
    if args.allow_qve_loss:
        command.append("--allow-qve-loss")
    if args.redock_refined:
        command.append("--redock-refined")
    print(f"generation_seeds={len(seed_smiles)} outdir={outdir}")
    run_command(command, args.dry_run)
    if not args.dry_run:
        deduplicate_refiner_outputs(outdir)


def run_score_model(paths: Paths, args: argparse.Namespace) -> dict[str, object]:
    model_result = run_command_capture(
        [
            paths.python,
            shell_path(paths.pose_analyzer / "scripts" / "score_space_model.py"),
            "--history-csv",
            shell_path(paths.history_csv),
            "--analysis-csv",
            shell_path(paths.analysis_csv),
            "--reference-smiles-csv",
            shell_path(paths.reference_csv),
            "--outdir",
            shell_path(paths.score_dir),
            "--train-size",
            str(args.train_size),
            "--test-size",
            str(args.test_size),
            *(["--require-shap"] if args.require_shap else []),
        ],
        args.dry_run,
    )
    if model_result and model_result.stdout:
        print(model_result.stdout, end="" if model_result.stdout.endswith("\n") else "\n", flush=True)
    if model_result and model_result.stderr:
        print(model_result.stderr, end="" if model_result.stderr.endswith("\n") else "\n", flush=True)
    return read_metrics(paths.score_dir / "binding_score_model_metrics.json")


def rebuild_dashboard(paths: Paths, args: argparse.Namespace) -> None:
    run_command(
        [
            paths.python,
            shell_path(paths.pose_analyzer / "scripts" / "build_score_space_dashboard.py"),
            "--feature-matrix",
            shell_path(paths.score_dir / "binding_score_feature_matrix.csv"),
            "--outdir",
            shell_path(paths.dashboard_dir),
        ],
        args.dry_run,
    )


def build_iteration_frontier(
    paths: Paths,
    args: argparse.Namespace,
    iter_dir: Path,
    iteration: int,
    metrics: dict[str, object],
) -> argparse.Namespace | None:
    if args.dry_run:
        print("dry-run: skipped frontier calculation")
        return None
    gap_args = argparse.Namespace(
        feature_matrix=paths.score_dir / "binding_score_feature_matrix.csv",
        outdir=iter_dir,
        top_frontier_seeds=args.top_frontier_seeds,
        top_shortlist=args.top_shortlist,
        top_official_recommendations=args.top_official_recommendations,
        neighbor_k=args.neighbor_k,
        distance_metric=args.distance_metric,
        grid_bins=args.grid_bins,
        **frontier_weight_namespace(iteration, count_successful_docked(paths.history_csv), args),
        dashboard_dir=paths.dashboard_dir,
    )
    build_frontier(gap_args)
    write_iteration_manifest(iter_dir, iteration, args, gap_args, metrics)
    return gap_args


def run_pipeline(args: argparse.Namespace) -> None:
    paths = Paths(
        project=args.project.expanduser(),
        python=args.python,
        pose_analyzer=args.pose_analyzer.expanduser(),
        refiner=args.refiner.expanduser(),
    )
    previous_mae: float | None = None
    low_improvement_streak = 0
    for iteration in range(args.start_iteration, args.start_iteration + args.iterations):
        iter_dir = paths.active_dir / f"iteration_{iteration:02d}"
        print(f"\n=== active learning iteration {iteration} ===", flush=True)
        metrics = run_score_model(paths, args)
        current_mae = safe_float(metrics.get("test_mae")) if metrics else math.nan
        if previous_mae is not None and not math.isnan(current_mae):
            relative_improvement = (previous_mae - current_mae) / max(abs(previous_mae), 1e-12)
            if relative_improvement < args.convergence_min_relative_mae_improvement:
                low_improvement_streak += 1
            else:
                low_improvement_streak = 0
            print(
                f"mae_convergence_check previous={previous_mae:.6f} current={current_mae:.6f} "
                f"relative_improvement={relative_improvement:.4f} streak={low_improvement_streak}/{args.convergence_patience}",
                flush=True,
            )
        if not math.isnan(current_mae):
            previous_mae = current_mae
        if low_improvement_streak >= args.convergence_patience:
            write_convergence_report(iter_dir, iteration, previous_mae, low_improvement_streak, args)
            print(
                f"convergence_reached=true iteration={iteration} reason=mae_relative_improvement_below_"
                f"{args.convergence_min_relative_mae_improvement} for {low_improvement_streak} checks",
                flush=True,
            )
            break
        rebuild_dashboard(paths, args)
        build_iteration_frontier(paths, args, iter_dir, iteration, metrics)
        if args.generate_new_smiles:
            run_refiner_generation(paths, args, iter_dir)
            if args.dock_generated:
                print("[active-learning] dock_generated=true; refreshing score model, dashboard, and frontier after new dock_history rows", flush=True)
                refreshed_metrics = run_score_model(paths, args)
                rebuild_dashboard(paths, args)
                build_iteration_frontier(paths, args, iter_dir, iteration, refreshed_metrics)
        if args.stop_after_frontier:
            print("stop_after_frontier=true; review frontier seeds before generation/docking.")
            break
    print("dashboard command:")
    print(f"  cd {paths.dashboard_dir} && {paths.python} -m http.server 8765 --bind 127.0.0.1")
    print("open: http://127.0.0.1:8765")


def write_iteration_manifest(
    iter_dir: Path,
    iteration: int,
    args: argparse.Namespace,
    gap_args: argparse.Namespace | None = None,
    metrics: dict[str, object] | None = None,
) -> None:
    iter_dir.mkdir(parents=True, exist_ok=True)
    weight_source = gap_args or args
    metrics = metrics or {}
    manifest = {
        "iteration": iteration,
        "project": str(args.project.expanduser()),
        "top_frontier_seeds": args.top_frontier_seeds,
        "top_shortlist": args.top_shortlist,
        "top_official_recommendations": args.top_official_recommendations,
        "distance_metric": getattr(weight_source, "distance_metric", getattr(args, "distance_metric", "cosine")),
        "frontier_weight_phase": getattr(weight_source, "weight_phase", "manual"),
        "successful_docked_count": getattr(weight_source, "docked_count", ""),
        "weights": {
            "novelty": getattr(weight_source, "weight_novelty", args.weight_novelty),
            "uncertainty": getattr(weight_source, "weight_uncertainty", args.weight_uncertainty),
            "quality": getattr(weight_source, "weight_quality", args.weight_quality),
            "druglike": getattr(weight_source, "weight_druglike", args.weight_druglike),
            "chembl": getattr(weight_source, "weight_chembl", args.weight_chembl),
            "contact": getattr(weight_source, "weight_contact", args.weight_contact),
        },
        "model_metrics": {
            "test_mae": metrics.get("test_mae", ""),
            "test_r2": metrics.get("test_r2", ""),
            "cv_mae_mean": metrics.get("cv_mae_mean", ""),
            "cv_r2_mean": metrics.get("cv_r2_mean", ""),
            "model_interpretation": metrics.get("model_interpretation", ""),
        },
        "human_next_step": "Review frontier_seeds.csv. Pick seeds for refiner or official testing before running the next iteration.",
    }
    with (iter_dir / "iteration_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)


def write_convergence_report(
    iter_dir: Path,
    iteration: int,
    current_mae: float | None,
    streak: int,
    args: argparse.Namespace,
) -> None:
    iter_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "convergence_reached": True,
        "iteration": iteration,
        "criterion": "test_mae_relative_improvement",
        "min_relative_improvement": args.convergence_min_relative_mae_improvement,
        "patience": args.convergence_patience,
        "low_improvement_streak": streak,
        "current_test_mae": current_mae if current_mae is not None else "",
    }
    with (iter_dir / "convergence_report.json").open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    def add_common(p: argparse.ArgumentParser) -> None:
        p.add_argument("--project", type=Path, default=DEFAULT_PROJECT)
        p.add_argument("--python", default=DEFAULT_PYTHON)
        p.add_argument("--pose-analyzer", type=Path, default=DEFAULT_POSE_ANALYZER)
        p.add_argument("--refiner", type=Path, default=DEFAULT_REFINER)
        p.add_argument("--top-frontier-seeds", type=int, default=30)
        p.add_argument("--top-shortlist", type=int, default=100)
        p.add_argument("--top-official-recommendations", type=int, default=3)
        p.add_argument("--neighbor-k", type=int, default=8)
        p.add_argument("--distance-metric", choices=["cosine", "euclidean"], default="cosine", help="Distance metric for high-dimensional frontier/uncertainty sampling")
        p.add_argument("--grid-bins", type=int, default=16)
        p.add_argument("--weight-novelty", type=float, default=0.28)
        p.add_argument("--weight-uncertainty", type=float, default=0.22)
        p.add_argument("--weight-quality", type=float, default=0.18)
        p.add_argument("--weight-druglike", type=float, default=0.14)
        p.add_argument("--weight-chembl", type=float, default=0.10)
        p.add_argument("--weight-contact", type=float, default=0.08)
        p.add_argument("--auto-weights", dest="auto_weights", action="store_true", default=True, help="Adapt frontier weights by iteration phase and successful docked count")
        p.add_argument("--manual-weights", dest="auto_weights", action="store_false", help="Use the explicit --weight-* values without phase adaptation")

    run = sub.add_parser("run", help="Run score model, dashboard build, and frontier analysis")
    add_common(run)
    run.add_argument("--train-size", type=int, default=60)
    run.add_argument("--test-size", type=int, default=20)
    run.add_argument("--iterations", type=int, default=1, help="Number of project iterations to plan/run. Use 10 for a full active-learning campaign; default 1 is safer for review.")
    run.add_argument("--start-iteration", type=int, default=1)
    run.add_argument("--require-shap", action="store_true")
    run.add_argument("--generate-new-smiles", action="store_true", help="After frontier selection, call druglike-pocket-refiner on frontier seeds to create new SMILES")
    run.add_argument("--dock-generated", action="store_true", help="When generating new SMILES, dock top refined candidates and write them back to dock_history.csv")
    run.add_argument("--frontier-seeds-to-refine", type=int, default=5)
    run.add_argument("--refiner-rounds", type=int, default=2)
    run.add_argument("--refiner-batch-size", type=int, default=200)
    run.add_argument("--refiner-beam-size", type=int, default=3)
    run.add_argument("--target-refinement-score", type=float, default=5.0)
    run.add_argument("--min-generated-before-target-stop", type=int, default=0, help="Pass through to refiner; 0 means at least max(batch-size, dock-top-candidates) new candidates before target-score early stop")
    run.add_argument("--dock-top-candidates", type=int, default=20)
    run.add_argument("--dock-workers", type=int, default=1)
    run.add_argument("--allow-add", action="store_true")
    run.add_argument("--allow-qve-loss", action="store_true")
    run.add_argument("--redock-refined", action="store_true")
    run.add_argument("--receptor", type=Path)
    run.add_argument("--config", type=Path)
    run.add_argument("--vina", type=Path)
    run.add_argument("--convergence-patience", type=int, default=2, help="Stop after this many consecutive MAE checks with too little improvement")
    run.add_argument("--convergence-min-relative-mae-improvement", type=float, default=0.05, help="Minimum relative test-MAE improvement required to reset convergence streak")
    run.add_argument("--stop-after-frontier", action="store_true", help="Stop after the first frontier seed set so a human can review before later iterations")
    run.add_argument("--dry-run", action="store_true")

    frontier = sub.add_parser("frontier", help="Only build frontier/gap outputs from an existing feature matrix")
    frontier.add_argument("--feature-matrix", type=Path, default=DEFAULT_PROJECT / "score_space_model" / "binding_score_feature_matrix.csv")
    frontier.add_argument("--outdir", type=Path, default=DEFAULT_PROJECT / "active_learning")
    frontier.add_argument("--dashboard-dir", type=Path, help="Optional dashboard directory; writes data/active_learning_recommendations.json for highlighting")
    frontier.add_argument("--top-frontier-seeds", type=int, default=30)
    frontier.add_argument("--top-shortlist", type=int, default=100)
    frontier.add_argument("--top-official-recommendations", type=int, default=3)
    frontier.add_argument("--neighbor-k", type=int, default=8)
    frontier.add_argument("--distance-metric", choices=["cosine", "euclidean"], default="cosine")
    frontier.add_argument("--grid-bins", type=int, default=16)
    frontier.add_argument("--weight-novelty", type=float, default=0.28)
    frontier.add_argument("--weight-uncertainty", type=float, default=0.22)
    frontier.add_argument("--weight-quality", type=float, default=0.18)
    frontier.add_argument("--weight-druglike", type=float, default=0.14)
    frontier.add_argument("--weight-chembl", type=float, default=0.10)
    frontier.add_argument("--weight-contact", type=float, default=0.08)
    return parser


def main(argv: Iterable[str] | None = None) -> None:
    parser = build_parser()
    raw_args = list(argv) if argv is not None else sys.argv[1:]
    if raw_args and raw_args[0] not in {"run", "frontier", "-h", "--help"}:
        raw_args = ["run", *raw_args]
    args = parser.parse_args(raw_args)
    if args.command == "run":
        run_pipeline(args)
    elif args.command == "frontier":
        build_frontier(args)
    else:
        parser.error(f"unknown command: {args.command}")


if __name__ == "__main__":
    main()
