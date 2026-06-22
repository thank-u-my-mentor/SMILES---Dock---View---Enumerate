#!/usr/bin/env python
"""Build a small static HTML report for docking and pose-analysis outputs."""

from __future__ import annotations

import argparse
import csv
import html
import os
from pathlib import Path


DEFAULT_TASK_ROOT = Path("~/vina_task2")
DEFAULT_HISTORY_CSV = DEFAULT_TASK_ROOT / "dock_history" / "dock_history.csv"
DEFAULT_POSE_ANALYSIS_DIR = DEFAULT_TASK_ROOT / "pose_analysis"
DEFAULT_REPORT_DIR = DEFAULT_TASK_ROOT / "pose_report"


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


def safe_float(value: object) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except Exception:
        return None


def existing_success_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for row in rows:
        affinity = safe_float(row.get("affinity_kcal_mol"))
        if affinity is None:
            continue
        if row.get("reason"):
            continue
        if not row.get("pose_path"):
            continue
        out.append(row)
    return out


def rel_link(value: object, report_dir: Path, fallback_label: str = "open") -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    path = Path(text).expanduser()
    if not path.is_absolute():
        href = text.replace("\\", "/")
    else:
        try:
            href = os.path.relpath(path, report_dir).replace("\\", "/")
        except ValueError:
            href = path.as_posix()
    label = fallback_label if path.is_dir() else (Path(text).name or fallback_label)
    return f'<a href="{html.escape(href)}">{html.escape(label)}</a>'


def cell(row: dict[str, str], field: str, report_dir: Path) -> str:
    if field == "smiles":
        value = row.get("smiles") or row.get("canonical_smiles") or row.get("input_smiles") or ""
    else:
        value = row.get(field, "")
    if field in {"pose_path", "log_path", "prep_log_path", "view_pml", "view_pse", "pocket_dir", "standard_pose_path", "ligand_pose_path"}:
        return rel_link(value, report_dir, field)
    return html.escape(str(value))


def table_html(
    rows: list[dict[str, str]],
    fields: list[str],
    report_dir: Path,
    empty_text: str,
    limit: int = 50,
) -> str:
    if not rows:
        return f'<p class="empty">{html.escape(empty_text)}</p>'
    header = "".join(f"<th>{html.escape(field)}</th>" for field in fields)
    body_rows = []
    for row in rows[:limit]:
        body_rows.append("<tr>" + "".join(f"<td>{cell(row, field, report_dir)}</td>" for field in fields) + "</tr>")
    more = ""
    if len(rows) > limit:
        more = f'<p class="empty">Showing {limit} of {len(rows)} rows.</p>'
    return f"<div class=\"table-wrap\"><table><thead><tr>{header}</tr></thead><tbody>{''.join(body_rows)}</tbody></table></div>{more}"


def score_key(row: dict[str, str]) -> tuple[float, float]:
    official = safe_float(row.get("official_binding_score"))
    affinity = safe_float(row.get("affinity_kcal_mol"))
    if official is not None:
        return (-official, affinity if affinity is not None else 999.0)
    return (0.0, affinity if affinity is not None else 999.0)


def build_html(
    *,
    history_rows: list[dict[str, str]],
    family_rows: list[dict[str, str]],
    residue_rows: list[dict[str, str]],
    surface_rows: list[dict[str, str]],
    report_dir: Path,
    pose_analysis_dir: Path,
) -> str:
    successful = sorted(existing_success_rows(history_rows), key=score_key)
    failed = [row for row in history_rows if row.get("reason")]
    best_affinity = min((safe_float(row.get("affinity_kcal_mol")) for row in successful), default=None)
    pse_count = sum(1 for row in family_rows if row.get("view_pse"))
    pml_count = sum(1 for row in family_rows if row.get("view_pml"))
    mode_fields = [
        "standard_seq_id",
        "standard_analysis_score",
        "analysis_score_source",
        "standard_mode",
        "standard_mode_affinity",
        "support_count",
        "weighted_support",
        "standard_binding_surface_coverage_4a",
        "standard_surface_contact_fraction_4a",
        "view_pml",
        "view_pse",
    ]
    dock_fields = [
        "seq_id",
        "nickname",
        "smiles",
        "affinity_kcal_mol",
        "inner_rmsd",
        "cnn_pose_score",
        "hbond_count",
        "hydrophobic_count",
        "pi_contact_count",
        "pose_path",
    ]
    residue_fields = ["residue", "score_weighted_support"]
    surface_fields = ["residue", "residue_name", "chain", "residue_id", "surface_like", "source"]
    failed_fields = ["seq_id", "smiles", "reason", "prep_log_path", "log_path"]
    mode_table = table_html(family_rows, mode_fields, report_dir, "No mode-family summary was found.")
    dock_table = table_html(successful, dock_fields, report_dir, "No successful docked rows were found.", limit=100)
    residue_table = table_html(residue_rows, residue_fields, report_dir, "No consensus residue table was found.", limit=30)
    surface_table = table_html(surface_rows, surface_fields, report_dir, "No predicted binding-surface table was found.", limit=40)
    failed_table = table_html(failed, failed_fields, report_dir, "No failed rows recorded.", limit=30)
    pose_dir_link = rel_link(pose_analysis_dir, report_dir, "pose_analysis")
    best_text = "" if best_affinity is None else f"{best_affinity:.3f} kcal/mol"
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Docking Pose Report</title>
  <style>
    :root {{
      color-scheme: light;
      --ink: #182026;
      --muted: #5e6a70;
      --line: #d9e0df;
      --panel: #f7f9f8;
      --accent: #0b6f6a;
      --accent-2: #7a4d12;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: var(--ink);
      background: #ffffff;
      line-height: 1.45;
    }}
    header {{
      padding: 28px clamp(18px, 4vw, 48px) 18px;
      border-bottom: 1px solid var(--line);
      background: linear-gradient(180deg, #f8fbfa 0%, #ffffff 100%);
    }}
    h1 {{
      margin: 0 0 12px;
      font-size: clamp(26px, 4vw, 42px);
      font-weight: 720;
      letter-spacing: 0;
    }}
    h2 {{
      margin: 0 0 12px;
      font-size: 18px;
      letter-spacing: 0;
    }}
    a {{ color: var(--accent); text-decoration-thickness: 1px; text-underline-offset: 2px; }}
    main {{ padding: 20px clamp(18px, 4vw, 48px) 48px; }}
    .summary {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
      gap: 10px;
      max-width: 980px;
    }}
    .metric {{
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 12px;
      background: var(--panel);
      min-height: 74px;
    }}
    .metric span {{
      display: block;
      color: var(--muted);
      font-size: 12px;
      text-transform: uppercase;
    }}
    .metric strong {{
      display: block;
      margin-top: 6px;
      font-size: 20px;
      font-weight: 700;
      color: var(--ink);
    }}
    section {{ margin-top: 28px; }}
    .empty {{ color: var(--muted); margin: 8px 0 0; }}
    .table-wrap {{
      width: 100%;
      overflow-x: auto;
      border: 1px solid var(--line);
      border-radius: 6px;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 13px;
      min-width: 860px;
    }}
    th, td {{
      padding: 8px 10px;
      border-bottom: 1px solid var(--line);
      text-align: left;
      vertical-align: top;
      max-width: 320px;
      overflow-wrap: anywhere;
    }}
    th {{
      position: sticky;
      top: 0;
      background: #edf4f3;
      font-size: 12px;
      color: #263338;
    }}
    tbody tr:nth-child(even) {{ background: #fbfcfc; }}
    .note {{
      color: var(--muted);
      max-width: 920px;
      margin: 8px 0 0;
    }}
  </style>
</head>
<body>
  <header>
    <h1>Docking Pose Report</h1>
    <div class="summary">
      <div class="metric"><span>Docking rows</span><strong>{len(history_rows)}</strong></div>
      <div class="metric"><span>Successful docks</span><strong>{len(successful)}</strong></div>
      <div class="metric"><span>Best affinity</span><strong>{html.escape(best_text or "n/a")}</strong></div>
      <div class="metric"><span>Mode families</span><strong>{len(family_rows)}</strong></div>
      <div class="metric"><span>PyMOL PML</span><strong>{pml_count}</strong></div>
      <div class="metric"><span>PyMOL PSE</span><strong>{pse_count}</strong></div>
    </div>
    <p class="note">Pose-analysis folder: {pose_dir_link}</p>
  </header>
  <main>
    <section>
      <h2>Reference Mode Families</h2>
      {mode_table}
    </section>
    <section>
      <h2>Docked Molecules</h2>
      {dock_table}
    </section>
    <section>
      <h2>Consensus Contact Residues</h2>
      {residue_table}
    </section>
    <section>
      <h2>Predicted Binding Surface</h2>
      {surface_table}
    </section>
    <section>
      <h2>Failed Rows</h2>
      {failed_table}
    </section>
  </main>
</body>
</html>
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--history-csv", type=Path, default=DEFAULT_HISTORY_CSV)
    parser.add_argument("--pose-analysis-dir", type=Path, default=DEFAULT_POSE_ANALYSIS_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_REPORT_DIR)
    args = parser.parse_args()

    history_csv = args.history_csv.expanduser().resolve()
    pose_analysis_dir = args.pose_analysis_dir.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    history_rows = read_csv(history_csv)
    family_rows = read_csv(pose_analysis_dir / "mode_family_summary.csv")
    residue_rows = read_csv(pose_analysis_dir / "consensus_residue_contacts.csv")
    surface_rows = read_csv(pose_analysis_dir / "predicted_binding_surface_residues.csv")
    html_text = build_html(
        history_rows=history_rows,
        family_rows=family_rows,
        residue_rows=residue_rows,
        surface_rows=surface_rows,
        report_dir=outdir,
        pose_analysis_dir=pose_analysis_dir,
    )
    report_path = outdir / "index.html"
    report_path.write_text(html_text, encoding="utf-8")
    print(f"wrote={report_path}")


if __name__ == "__main__":
    main()
