#!/usr/bin/env python3
"""Render every scene in a PyMOL PSE to per-PDB PNG files."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from pymol import cmd


DEFAULT_REVIEW_DIR = Path("/mnt/e/QZHAO-LAB/Metal-F_plasmid_candidate_review")


def safe_name(text: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")
    return value or "scene"


def read_summary(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle))


def read_review(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with open(path, newline="", encoding="utf-8", errors="replace") as handle:
        return {row.get("pdb_id", "").upper(): row for row in csv.DictReader(handle) if row.get("pdb_id")}


def rel_for_html(path: Path, base: Path) -> str:
    try:
        return path.relative_to(base).as_posix()
    except ValueError:
        return path.as_posix()


def render(
    pse: Path,
    summary_csv: Path,
    outdir: Path,
    index_csv: Path,
    review_csv: Path | None,
    width: int,
    height: int,
    ray: bool,
    max_scenes: int,
) -> int:
    rows = read_summary(summary_csv)
    if max_scenes:
        rows = rows[:max_scenes]
    review = read_review(review_csv) if review_csv else {}

    cmd.delete("all")
    cmd.load(str(pse))
    scenes = set(cmd.get_scene_list())
    outdir.mkdir(parents=True, exist_ok=True)
    index_rows: list[dict[str, str]] = []

    cmd.bg_color("white")
    cmd.set("ray_shadow", "off")
    cmd.set("antialias", 2)
    cmd.set("orthoscopic", 1)
    cmd.set("ray_opaque_background", 1)
    cmd.set("label_color", "yellow")

    for i, row in enumerate(rows, 1):
        pdb_id = str(row.get("pdb_id", "")).upper()
        scene = str(row.get("scene", ""))
        if not pdb_id or not scene:
            continue
        if scene not in scenes:
            print(f"[render_pse_scene_pngs] skip missing scene: {scene}", flush=True)
            continue

        target_dir = outdir / pdb_id
        target_dir.mkdir(parents=True, exist_ok=True)
        png = target_dir / f"{safe_name(scene)}.png"

        print(f"[render_pse_scene_pngs] {i}/{len(rows)} {pdb_id} {scene}", flush=True)
        cmd.scene(scene, "recall")
        cmd.refresh()
        if ray:
            cmd.ray(width, height)
        cmd.png(str(png), width=width, height=height, dpi=180, ray=0)

        meta = review.get(pdb_id, {})
        index_rows.append(
            {
                "pdb_id": pdb_id,
                "scene": scene,
                "png_path": str(png),
                "png_relpath": rel_for_html(png, index_csv.parent),
                "nearest_F_metal_A": row.get("nearest_F_metal_A", ""),
                "fluorine_comp_ids": row.get("fluorine_comp_ids", ""),
                "metal_comp_ids": row.get("metal_comp_ids", ""),
                "final_decision": meta.get("final_decision", ""),
                "source_organism_corrected": meta.get("source_organism_corrected", ""),
                "expression_curation_status": meta.get("expression_curation_status", ""),
            }
        )

    index_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "pdb_id",
        "scene",
        "png_path",
        "png_relpath",
        "nearest_F_metal_A",
        "fluorine_comp_ids",
        "metal_comp_ids",
        "final_decision",
        "source_organism_corrected",
        "expression_curation_status",
    ]
    with open(index_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(index_rows)
    print(f"[render_pse_scene_pngs] rendered={len(index_rows)} index={index_csv}", flush=True)
    return len(index_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pse", default=str(DEFAULT_REVIEW_DIR / "Metal-F_plasmid_all_nonduplicate_candidates.pse"))
    parser.add_argument("--summary-csv", default=str(DEFAULT_REVIEW_DIR / "Metal-F_plasmid_all_nonduplicate_candidates_pse_summary.csv"))
    parser.add_argument("--review-csv", default=str(DEFAULT_REVIEW_DIR / "Metal-F_plasmid_candidate_review.csv"))
    parser.add_argument("--outdir", default=str(DEFAULT_REVIEW_DIR / "scene_png"))
    parser.add_argument("--index-csv", default=str(DEFAULT_REVIEW_DIR / "scene_image_index.csv"))
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=1200)
    parser.add_argument("--max-scenes", type=int, default=0)
    parser.add_argument("--no-ray", action="store_true", help="Use OpenGL png export without ray tracing.")
    args = parser.parse_args()
    render(
        Path(args.pse),
        Path(args.summary_csv),
        Path(args.outdir),
        Path(args.index_csv),
        Path(args.review_csv) if args.review_csv else None,
        args.width,
        args.height,
        not args.no_ray,
        args.max_scenes,
    )


if __name__ == "__main__":
    main()
