#!/usr/bin/env python
"""Build an HTML report focused on residue-ligand velocity distributions."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rdc-dir", required=True, type=Path, help="Directory containing the original RDC report")
    parser.add_argument("--velocity-dir", required=True, type=Path, help="Directory containing velocity distribution outputs")
    parser.add_argument("--title", default="Distribution Density of Distance Fluctuations Between Residues and Substrate")
    parser.add_argument(
        "--subtitle",
        default="Residue-ligand distance-fluctuation velocity distribution; ligand selection: <code>resname LIG</code>",
    )
    parser.add_argument("--top-n", type=int, default=12)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def copy_outputs(velocity_dir: Path, rdc_dir: Path) -> None:
    for name in [
        "velocity_distribution_density.png",
        "residue_ligand_velocity_distribution.csv",
        "residue_ligand_velocity_distribution_summary.csv",
    ]:
        src = velocity_dir / name
        if src.exists():
            shutil.copy2(src, rdc_dir / name)


def table(rows: list[dict[str, str]], top_n: int) -> str:
    columns = [
        "residue",
        "mean_velocity_x1e_minus2_A_per_ps",
        "std_velocity_x1e_minus2_A_per_ps",
        "fraction_approaching",
        "fraction_departing",
        "p05_x1e_minus2",
        "p95_x1e_minus2",
    ]
    head = "".join(f"<th>{column}</th>" for column in columns)
    body = []
    for row in rows[:top_n]:
        cls = "approach" if float(row["fraction_approaching"]) > float(row["fraction_departing"]) else "depart"
        body.append("<tr class='" + cls + "'>" + "".join(f"<td>{row.get(column, '')}</td>" for column in columns) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def main() -> None:
    args = parse_args()
    args.rdc_dir.mkdir(parents=True, exist_ok=True)
    copy_outputs(args.velocity_dir, args.rdc_dir)
    summary_path = args.rdc_dir / "residue_ligand_velocity_distribution_summary.csv"
    rows = read_rows(summary_path) if summary_path.exists() else []
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{args.title}</title>
  <style>
    :root {{
      --ink: #1f2328;
      --muted: #6b7280;
      --red: #c1121f;
      --line: #d6d6d6;
      --green: #2f7d59;
      --blue: #355c7d;
      --paper: #ffffff;
    }}
    body {{
      margin: 0;
      background: var(--paper);
      color: var(--ink);
      font-family: Arial, Helvetica, sans-serif;
    }}
    main {{
      max-width: 1120px;
      margin: 0 auto;
      padding: 26px 28px 38px;
    }}
    .panel-label {{
      font-size: 42px;
      font-weight: 800;
      line-height: 1;
      margin-bottom: 8px;
    }}
    h1 {{
      color: var(--red);
      font-size: 23px;
      line-height: 1.25;
      text-align: center;
      margin: 0 0 8px;
      font-weight: 800;
    }}
    .subtitle {{
      text-align: center;
      color: var(--muted);
      font-size: 13px;
      margin: 0 0 18px;
    }}
    .figure-wrap {{
      position: relative;
      border-top: 1px solid transparent;
    }}
    .figure-wrap img {{
      display: block;
      width: min(100%, 980px);
      margin: 0 auto;
    }}
    .axis-note {{
      display: grid;
      grid-template-columns: repeat(3, 1fr);
      gap: 12px;
      margin: 14px auto 20px;
      max-width: 980px;
      font-size: 14px;
      color: #333;
    }}
    .axis-note div {{
      border-top: 1px solid var(--line);
      padding-top: 8px;
      text-align: center;
    }}
    .axis-note strong {{
      display: block;
      margin-bottom: 3px;
      font-size: 15px;
    }}
    .definition {{
      max-width: 980px;
      margin: 0 auto 20px;
      color: #333;
      font-size: 14px;
      line-height: 1.55;
    }}
    code {{
      background: #f3f4f6;
      border-radius: 4px;
      padding: 1px 5px;
    }}
    table {{
      border-collapse: collapse;
      width: 100%;
      margin-top: 18px;
      font-size: 13px;
    }}
    th, td {{
      border-bottom: 1px solid #e5e7eb;
      padding: 7px 8px;
      text-align: right;
      white-space: nowrap;
    }}
    th:first-child, td:first-child {{
      text-align: left;
    }}
    th {{
      color: #374151;
      font-weight: 700;
      background: #fafafa;
    }}
    tr.approach td:first-child {{
      color: var(--green);
      font-weight: 700;
    }}
    tr.depart td:first-child {{
      color: var(--blue);
      font-weight: 700;
    }}
    .files {{
      margin-top: 18px;
      color: var(--muted);
      font-size: 13px;
    }}
  </style>
</head>
<body>
<main>
  <div class="panel-label">C</div>
  <h1>{args.title}</h1>
  <p class="subtitle">{args.subtitle}</p>
  <section class="figure-wrap">
    <img src="velocity_distribution_density.png" alt="Velocity distribution density">
  </section>
  <section class="axis-note">
    <div><strong>Negative V</strong>Residue-substrate distance is shortening.</div>
    <div><strong>Near zero</strong>Low distance fluctuation over adjacent frames.</div>
    <div><strong>Positive V</strong>Residue-substrate distance is increasing.</div>
  </section>
  <section class="definition">
    <p><code>V = delta(distance) / delta(time)</code>, plotted as <code>V x 10^-2 (A/ps)</code>. A wider distribution indicates stronger local distance fluctuation. This velocity-distribution view captures short-timescale dynamic changes; the older linear RDC table remains available in the CSV files for long-timescale drift.</p>
  </section>
  <h2>Residue Velocity Summary</h2>
  {table(rows, args.top_n)}
  <p class="files">Files: <code>velocity_distribution_density.png</code>, <code>residue_ligand_velocity_distribution.csv</code>, <code>residue_ligand_velocity_distribution_summary.csv</code>, plus the original RDC CSV files in this directory.</p>
</main>
</body>
</html>
"""
    (args.rdc_dir / "index.html").write_text(html, encoding="utf-8")
    print(args.rdc_dir / "index.html")


if __name__ == "__main__":
    main()
