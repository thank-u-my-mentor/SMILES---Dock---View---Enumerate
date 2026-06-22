#!/usr/bin/env python
"""Run a simple receptor + starting-SMILES docking and pose-analysis project."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
SKILL_DIR = SCRIPT_DIR.parent
CODEX_SKILLS_DIR = SKILL_DIR.parent
POSE_ANALYZER_DIR = CODEX_SKILLS_DIR / "pose-analyzer"
POSE_ANALYZER_SCRIPTS = POSE_ANALYZER_DIR / "scripts"
DEFAULT_PROJECT_DIR = Path("~/vina_task2")


def run(cmd: list[str], *, dry_run: bool = False) -> None:
    print("$ " + " ".join(cmd), flush=True)
    if dry_run:
        return
    subprocess.run(cmd, check=True)


def resolve_engine(engine: str, explicit: str | None) -> str:
    if explicit:
        return explicit
    if shutil.which(engine):
        return engine
    if engine == "glina" and shutil.which("gnina"):
        return "gnina"
    if engine == "gnina" and shutil.which("glina"):
        return "glina"
    if engine != "vina" and shutil.which("vina"):
        return "vina"
    return engine


def default_rounds(args: argparse.Namespace) -> int:
    if args.max_rounds is not None:
        return args.max_rounds
    return 1 if args.analogs else 0


def ensure_pose_analyzer() -> None:
    needed = POSE_ANALYZER_SCRIPTS / "infer_binding_mode_families.py"
    if not needed.exists():
        raise FileNotFoundError(f"pose-analyzer sibling skill not found at {needed}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smiles", required=True, help="Starting ligand SMILES")
    parser.add_argument("--receptor", type=Path, required=True, help="Prepared receptor; use PDBQT for Vina-style docking")
    parser.add_argument("--config", type=Path, required=True, help="Vina-style config with center_x/y/z and size_x/y/z")
    parser.add_argument("--project-dir", type=Path, default=DEFAULT_PROJECT_DIR, help="Project root; defaults to ~/vina_task2")
    parser.add_argument("--nickname", default="lead_001")
    parser.add_argument("--engine", choices=["vina", "gnina", "glina"], default="vina", help="Docking executable family")
    parser.add_argument("--vina", "--engine-exe", dest="vina", help="Explicit docking executable path; also used for gnina/glina if supplied")
    parser.add_argument("--meeko", help="Optional mk_prepare_ligand.py path")
    parser.add_argument("--obabel", help="Optional obabel path")
    parser.add_argument("--analogs", action="store_true", help="Generate and dock one conservative analog batch")
    parser.add_argument("--max-rounds", type=int, help="Analog generation rounds; overrides --analogs")
    parser.add_argument("--batch-size", type=int, default=24, help="Maximum generated analogs per round")
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--save-pse", action="store_true", help="Ask PyMOL to save .pse sessions in addition to .pml scripts")
    parser.add_argument("--tune-grid-box", action="store_true", help="Run pose-analyzer local grid-box tuning")
    parser.add_argument("--score-column", default="auto", help="Pose-analysis score column; auto uses official score if present, otherwise affinity")
    parser.add_argument("--skip-env-check", action="store_true")
    parser.add_argument("--skip-pose-analysis", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    ensure_pose_analyzer()
    project_dir = args.project_dir.expanduser().resolve()
    history_dir = project_dir / "dock_history"
    pose_analysis_dir = project_dir / "pose_analysis"
    report_dir = project_dir / "pose_report"
    docking_executable = resolve_engine(args.engine, args.vina)

    if not args.skip_env_check:
        env_cmd = [
            sys.executable,
            str(SCRIPT_DIR / "check_local_env.py"),
            "--docking-engine",
            docking_executable,
            "--strict",
        ]
        run(env_cmd, dry_run=args.dry_run)

    dock_cmd = [
        sys.executable,
        str(SCRIPT_DIR / "dock_smiles.py"),
        "--smiles",
        args.smiles,
        "--nickname",
        args.nickname,
        "--history-dir",
        str(history_dir),
        "--receptor",
        str(args.receptor.expanduser().resolve()),
        "--config",
        str(args.config.expanduser().resolve()),
        "--vina",
        docking_executable,
        "--max-rounds",
        str(default_rounds(args)),
        "--batch-size",
        str(args.batch_size),
        "--exhaustiveness",
        str(args.exhaustiveness),
        "--num-modes",
        str(args.num_modes),
    ]
    if args.cpu is not None:
        dock_cmd.extend(["--cpu", str(args.cpu)])
    if args.meeko:
        dock_cmd.extend(["--meeko", args.meeko])
    if args.obabel:
        dock_cmd.extend(["--obabel", args.obabel])
    run(dock_cmd, dry_run=args.dry_run)

    if args.skip_pose_analysis:
        print(f"history_dir={history_dir}", flush=True)
        return

    infer_cmd = [
        sys.executable,
        str(POSE_ANALYZER_SCRIPTS / "infer_binding_mode_families.py"),
        "--analysis-csv",
        str(history_dir / "dock_history.csv"),
        "--history-csv",
        str(history_dir / "dock_history.csv"),
        "--receptor",
        str(args.receptor.expanduser().resolve()),
        "--config",
        str(args.config.expanduser().resolve()),
        "--outdir",
        str(pose_analysis_dir),
        "--score-column",
        args.score_column,
    ]
    if args.save_pse:
        infer_cmd.append("--save-pse")
    if args.tune_grid_box:
        infer_cmd.append("--tune-grid-box")
    run(infer_cmd, dry_run=args.dry_run)

    report_cmd = [
        sys.executable,
        str(POSE_ANALYZER_SCRIPTS / "build_pose_report.py"),
        "--history-csv",
        str(history_dir / "dock_history.csv"),
        "--pose-analysis-dir",
        str(pose_analysis_dir),
        "--outdir",
        str(report_dir),
    ]
    run(report_cmd, dry_run=args.dry_run)
    print(f"report={report_dir / 'index.html'}", flush=True)


if __name__ == "__main__":
    main()
