#!/usr/bin/env python
"""Build a small HTML report for residue-ligand RDC outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rdc-dir", required=True, type=Path)
    parser.add_argument("--title", default="Residue-Ligand RDC Report")
    parser.add_argument("--top-n", type=int, default=30)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def td(text: object) -> str:
    return f"<td>{str(text)}</td>"


def main() -> None:
    args = parse_args()
    rows = read_rows(args.rdc_dir / "residue_ligand_rdc_within_cutoff.csv")
    rows = rows[: args.top_n]
    columns = [
        "residue",
        "min_distance_A",
        "initial_distance_A",
        "final_distance_A",
        "slope_A_per_ns",
        "rdc_per_ns",
        "fit_r2",
    ]
    table_rows = []
    for row in rows:
        rdc = float(row["rdc_per_ns"])
        cls = "approach" if rdc < 0 else "depart"
        table_rows.append(
            f"<tr class='{cls}'>" + "".join(td(row.get(column, "")) for column in columns) + "</tr>"
        )
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{args.title}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; color: #222; }}
    h1 {{ font-size: 22px; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(360px, 1fr)); gap: 18px; }}
    img {{ max-width: 100%; border: 1px solid #ddd; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 18px; font-size: 13px; }}
    th, td {{ border-bottom: 1px solid #ddd; padding: 6px 8px; text-align: right; }}
    th:first-child, td:first-child {{ text-align: left; }}
    tr.approach td:first-child {{ color: #087f5b; font-weight: 700; }}
    tr.depart td:first-child {{ color: #a61e4d; font-weight: 700; }}
    code {{ background: #f3f3f3; padding: 1px 4px; }}
  </style>
</head>
<body>
  <h1>{args.title}</h1>
  <p><code>rdc_per_ns = slope_A_per_ns / reference_distance_A</code>. Positive values mean moving away from the ligand; negative values mean approaching.</p>
  <div class="grid">
    <div><h2>Top Absolute RDC</h2><img src="top_abs_rdc_per_ns.png"></div>
    <div><h2>Top Absolute Distance Slope</h2><img src="top_abs_slope_A_per_ns.png"></div>
    <div><h2>Distance Traces</h2><img src="top_rdc_distance_traces.png"></div>
  </div>
  <h2>Top Residues</h2>
  <table>
    <thead><tr>{''.join(f'<th>{column}</th>' for column in columns)}</tr></thead>
    <tbody>
      {''.join(table_rows)}
    </tbody>
  </table>
  <p>Files: <code>residue_ligand_rdc_all.csv</code>, <code>residue_ligand_rdc_within_cutoff.csv</code>, <code>residue_ligand_distance_timeseries.csv</code>.</p>
</body>
</html>
"""
    (args.rdc_dir / "index.html").write_text(html, encoding="utf-8")
    print(args.rdc_dir / "index.html")


if __name__ == "__main__":
    main()
