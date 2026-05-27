#!/usr/bin/env python
"""Bridge: import pose-analyzer structural-interaction-space frontier seeds into teammate navigation memory.

Reads pose-analyzer feature matrix (and optionally frontier_seeds.csv), converts
high-frontier-score molecules into teammate-compatible navigation candidates, and
writes them to teammate's internal navigation_memory directory so the CommitteeJudge
can assign navigation_model_support_score to them.

Also writes a pending-runtime-evidence record describing the structural gap analysis.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def safe_float(v: str) -> float:
    try:
        if v in ("", None):
            return math.nan
        return float(v)
    except Exception:
        return math.nan


def fmt(v: float | str) -> str:
    if isinstance(v, float):
        if math.isnan(v):
            return ""
        return f"{v:.4f}"
    return str(v) if v is not None else ""


def build_nav_candidates(feature_rows: list[dict[str, str]], frontier_rows: list[dict[str, str]] | None, top_n: int) -> list[dict[str, str]]:
    """Build teammate navigation-memory rows from pose-analyzer frontier seeds."""
    # If frontier_rows provided, use them; otherwise fall back to feature_rows sorted by predicted novelty
    sources = frontier_rows if frontier_rows else feature_rows
    scored: list[tuple[float, dict[str, str]]] = []
    for row in sources:
        fs = safe_float(row.get("frontier_score"))
        if not math.isnan(fs):
            scored.append((fs, row))
            continue
        # fallback: use structural-interaction-space distance from origin as a proxy for novelty
        s1 = safe_float(row.get("structural_interaction_pc1"))
        s2 = safe_float(row.get("structural_interaction_pc2"))
        if not math.isnan(s1) and not math.isnan(s2):
            scored.append((math.sqrt(s1 * s1 + s2 * s2), row))
    scored.sort(key=lambda x: x[0], reverse=True)

    nav: list[dict[str, str]] = []
    for i, (_, row) in enumerate(scored[:top_n], start=1):
        smi = row.get("canonical_smiles", "") or row.get("mol_smiles", "")
        if not smi:
            continue
        cid = row.get("seq_id", "") or row.get("candidate_id", "") or f"PA_{i:04d}"
        ob = row.get("official_binding_score", "")
        pred = row.get("predicted_binding_score", "")
        qed = row.get("qed", "")
        aff = row.get("affinity_kcal_mol", "")
        chembl = row.get("chembl_scaffold_similarity", "")
        tags = ";".join(t for t in [
            "embedding_navigation",
            "pose_analyzer_frontier",
            "structural_gap",
            "diversity_quota",
        ] if t)
        nav.append({
            "candidate_id": f"PA_{cid}",
            "source_candidate_id": cid,
            "mol_smiles": smi,
            "route": row.get("route", ""),
            "navigation_family": "pose_analyzer_structural_frontier",
            "transformer_frontier_type": "structural_interaction_space",
            "nearest_positive_id": row.get("nearest_positive_id", ""),
            "nearest_positive_similarity": row.get("nearest_positive_similarity", ""),
            "nearest_negative_id": row.get("nearest_negative_id", ""),
            "nearest_negative_similarity": row.get("nearest_negative_similarity", ""),
            "committee_expected_binding": ob if ob else pred,
            "internal_priority_score": fmt(row.get("frontier_score", "")),
            "qed": qed,
            "affinity_kcal_mol": aff,
            "chembl_scaffold_similarity": chembl,
            "route_valid": "unknown",
            "score_allowed": "False",
            "allowed_use": "candidate_navigation_and_diversity_quota_only",
            "forbidden_use": "final_judge_or_binding_proof",
            "source_path": str(Path("~/vina_task2/score_space_model/binding_score_feature_matrix.csv").expanduser()),
            "tags": tags,
        })
    return nav


def write_nav_csv(nav_dir: Path, rows: list[dict[str, str]]) -> None:
    nav_dir.mkdir(parents=True, exist_ok=True)
    path = nav_dir / "pose_analyzer_frontier_candidates.csv"
    fields = [
        "candidate_id", "source_candidate_id", "mol_smiles", "route",
        "navigation_family", "transformer_frontier_type",
        "nearest_positive_id", "nearest_positive_similarity",
        "nearest_negative_id", "nearest_negative_similarity",
        "committee_expected_binding", "internal_priority_score",
        "qed", "affinity_kcal_mol", "chembl_scaffold_similarity",
        "route_valid", "score_allowed", "allowed_use", "forbidden_use",
        "source_path", "tags",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})
    print(f"wrote={path} rows={len(rows)}", flush=True)


def write_pending_evidence(pending_dir: Path, nav_count: int, gap_count: int) -> None:
    pending_dir.mkdir(parents=True, exist_ok=True)
    evidence = {
        "pending_rule_id": "PEND_POSE_ANALYZER_001",
        "rule_text": (
            "Pose-analyzer structural-interaction-space frontier analysis identified "
            f"{nav_count} frontier seeds and {gap_count} under-explored grid regions. "
            "These candidates are positioned at the boundary between known dense clusters "
            "and empty space in the docking-contact-descriptor embedding."
        ),
        "source_type": "internal_runtime_observation",
        "source_title": "pose-analyzer explore_space_gaps.py",
        "source_url_or_citation": "~/vina_task2/space_exploration/gap_regions.csv",
        "evidence_quote_or_summary": json.dumps({
            "method": "grid_density + knn_mean_distance on structural_interaction_pc1/pc2",
            "policy": "frontier seeds are navigation-only; they may seed diversity quotas but must not override hard gates or analog-pair evidence",
            "space_priority": "structural_interaction_space over pure_smiles_space",
        }, ensure_ascii=False),
        "applicable_skill": "CandidateGenerator;CommitteeJudge;DynamicKBExpander",
        "confidence": "pending",
        "limitations": (
            "Space coordinates are PCA projections of docking/contact/descriptors; "
            "they compress information and may miss subtle pose differences. "
            "Frontier seeds require docking validation before binding claims."
        ),
        "pending_status": "review_only_no_score",
        "score_allowed": False,
    }
    path = pending_dir / "pending_runtime_evidence.jsonl"
    # Append mode: teammate may already have other pending evidence
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(evidence, ensure_ascii=False) + "\n")
    print(f"appended={path}", flush=True)


def write_bridge_report(output_dir: Path, nav_rows: list[dict[str, str]], gap_path: Path | None) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Pose-Analyzer → Teammate Navigation Bridge Report",
        "",
        f"Imported {len(nav_rows)} frontier candidates into teammate navigation memory.",
        "",
        "## What was imported",
        "",
        "- Source: pose-analyzer `binding_score_feature_matrix.csv` (structural-interaction-space)",
        "- Destination: `agent/kb/internal/navigation_memory/pose_analyzer_frontier_candidates.csv`",
        "- Candidate tags: `embedding_navigation;pose_analyzer_frontier;structural_gap;diversity_quota`",
        "",
        "## Why this matters",
        "",
        "CommitteeJudge assigns `navigation_model_support_score` to candidates tagged with `embedding_navigation`.",
        "By feeding pose-analyzer frontier seeds into this channel, we tell the committee:",
        "'These molecules live at the edge of known space; they are worth exploring for diversity.'",
        "",
        "## Limitations",
        "",
        "- `score_allowed` is explicitly `False` for all imported candidates.",
        "- They are navigation signals only; they cannot override hard gates, route validity, or analog-pair evidence.",
        "- Official binding scores from these candidates must be obtained experimentally before they can be promoted.",
        "",
    ]
    if gap_path and gap_path.exists():
        lines.append(f"- Gap analysis source: `{gap_path}`")
    lines.append("")
    (output_dir / "POSE_ANALYZER_BRIDGE_REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote={output_dir / 'POSE_ANALYZER_BRIDGE_REPORT.md'}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-matrix", type=Path, required=True,
                        help="pose-analyzer binding_score_feature_matrix.csv")
    parser.add_argument("--frontier-csv", type=Path,
                        help="Optional explore_space_gaps.py frontier_seeds.csv (fallback to feature matrix)")
    parser.add_argument("--teammate-root", type=Path, required=True,
                        help="Path to teammate_portable_package root (contains Code/agent/kb/)")
    parser.add_argument("--top-candidates", type=int, default=30,
                        help="Number of frontier seeds to import")
    parser.add_argument("--gap-csv", type=Path,
                        help="Optional explore_space_gaps.py gap_regions.csv for audit report")
    args = parser.parse_args()

    feature_rows = read_csv(args.feature_matrix.expanduser())
    frontier_rows = read_csv(args.frontier_csv.expanduser()) if args.frontier_csv else None

    nav_rows = build_nav_candidates(feature_rows, frontier_rows, args.top_candidates)
    if not nav_rows:
        raise ValueError("no valid frontier candidates found")

    teammate_root = args.teammate_root.expanduser().resolve()
    nav_dir = teammate_root / "Code" / "agent" / "kb" / "internal" / "navigation_memory"
    pending_dir = teammate_root / "Code" / "agent" / "kb" / "pending"
    output_dir = teammate_root / "Output"

    write_nav_csv(nav_dir, nav_rows)
    write_pending_evidence(pending_dir, len(nav_rows), len(read_csv(args.gap_csv.expanduser())) if args.gap_csv else 0)
    write_bridge_report(output_dir, nav_rows, args.gap_csv.expanduser() if args.gap_csv else None)

    print("bridge_complete", flush=True)


if __name__ == "__main__":
    main()
