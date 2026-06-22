#!/usr/bin/env python3
"""Write GAMESS runner scripts for an MCPB.py directory."""

from __future__ import annotations

import argparse
from pathlib import Path


COORD_MARKER = "COORDINATES OF ALL ATOMS ARE (ANGS)"


def extract_last_gamess_coords(log: Path) -> list[str]:
    lines = log.read_text(encoding="utf-8", errors="replace").splitlines()
    indices = [idx for idx, line in enumerate(lines) if COORD_MARKER in line]
    if not indices:
        raise RuntimeError(f"no coordinate block found in {log}")
    coords: list[str] = []
    for line in lines[indices[-1] + 1 :]:
        stripped = line.strip()
        if not stripped:
            if coords:
                break
            continue
        parts = stripped.split()
        if len(parts) == 5 and parts[0][0].isalpha():
            try:
                coords.append(
                    f"{parts[0]:<2s} {float(parts[1]):8.1f} "
                    f"{float(parts[2]):14.6f} {float(parts[3]):14.6f} {float(parts[4]):14.6f}"
                )
            except ValueError:
                if coords:
                    break
        elif coords:
            break
    if not coords:
        raise RuntimeError(f"failed to parse coordinates from {log}")
    return coords


def fill_data(template: Path, out: Path, coords: list[str]) -> None:
    lines = template.read_text(encoding="utf-8", errors="replace").splitlines()
    data_idx = next(i for i, line in enumerate(lines) if line.strip().upper() == "$DATA")
    end_idx = next(i for i in range(data_idx + 1, len(lines)) if lines[i].strip().upper() == "$END")
    new_lines = lines[: data_idx + 3] + coords + [" $END"] + lines[end_idx + 1 :]
    out.write_text("\n".join(new_lines) + "\n", encoding="utf-8", newline="\n")


def normalize_open_shell_gamess_input(path: Path) -> None:
    """Make MCPB-generated GAMESS input explicit for open-shell DFT."""
    text = path.read_text(encoding="utf-8", errors="replace")
    if not path.with_suffix(path.suffix + ".bak_mcpb_default").exists():
        path.with_suffix(path.suffix + ".bak_mcpb_default").write_text(text, encoding="utf-8", newline="\n")
    lines = text.splitlines()
    out: list[str] = []
    inserted_scf = any(line.strip().upper().startswith("$SCF") for line in lines)
    for line in lines:
        stripped = line.strip()
        if stripped.upper().startswith("$CONTRL"):
            body = stripped
            if "SCFTYP=" not in body.upper():
                body = body.replace("$CONTRL", "$CONTRL SCFTYP=UHF", 1)
            if "MAXIT=" not in body.upper():
                body = body.replace("$END", "MAXIT=200 $END")
            out.append(" " + " ".join(body.split()))
            if not inserted_scf:
                out.append(" $SCF DIRSCF=.T. DIIS=.T. DAMP=.T. ETHRSH=2.0 MAXDII=20 $END")
                inserted_scf = True
            continue
        out.append(line)
    path.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Write MCPB GAMESS opt/fc/large-mk runner scripts.")
    parser.add_argument("--mcpb-dir", type=Path, required=True)
    parser.add_argument("--group", default="TJ_M7_FE")
    parser.add_argument("--gamess", type=Path, default=Path("/home/qin/softwares/gamess/rungms"))
    parser.add_argument("--cores", type=int, default=16)
    args = parser.parse_args()

    root = args.mcpb_dir
    group = args.group
    opt = f"{group}_small_opt"
    fc = f"{group}_small_fc"
    large = f"{group}_large_mk"
    for stem in (opt, fc, large):
        if not (root / f"{stem}.inp").exists():
            raise FileNotFoundError(root / f"{stem}.inp")
        normalize_open_shell_gamess_input(root / f"{stem}.inp")

    (root / "gamess_logs").mkdir(exist_ok=True)
    (root / "02_run_mcpb_gamess_qm.sh").write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
mkdir -p gamess_logs
GAMESS={args.gamess}
VERSION=${{VERSION:-00}}
NCORES=${{NCORES:-{args.cores}}}

run_job() {{
  local stem="$1"
  echo "[GAMESS] $stem start $(date) NCORES=$NCORES" | tee "gamess_logs/${{stem}}.status"
  "$GAMESS" "$stem" "$VERSION" "$NCORES" > "gamess_logs/${{stem}}.log" 2>&1
  local rc=$?
  echo "[GAMESS] $stem exit_code=$rc end $(date)" | tee -a "gamess_logs/${{stem}}.status"
  grep -q 'EXECUTION OF GAMESS TERMINATED NORMALLY' "gamess_logs/${{stem}}.log"
}}

run_job {opt}

/mnt/l/WSL/softwares/conda_envs/md/bin/python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/write_mcpb_gamess_pipeline.py \\
  --mcpb-dir "$PWD" --group {group} --fill-only

run_job {fc}
run_job {large}

echo all_mcpb_gamess_qm_finished
echo next='bash 03_after_qm_run_mcpb_steps.sh'
""",
        encoding="utf-8",
        newline="\n",
    )
    (root / "02_run_mcpb_gamess_qm.sh").chmod(0o755)

    (root / "monitor_mcpb_gamess_qm.sh").write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
for stem in {opt} {fc} {large}; do
  echo "===== $stem ====="
  tail -n 30 "gamess_logs/${{stem}}.log" 2>/dev/null || true
  grep -E 'NSERCH:|TOTAL ENERGY|S-SQUARED|TERMINATED|ERROR|FAILURE|HESSIAN|PDC|ELPOT' "gamess_logs/${{stem}}.log" 2>/dev/null | tail -n 30 || true
done
""",
        encoding="utf-8",
        newline="\n",
    )
    (root / "monitor_mcpb_gamess_qm.sh").chmod(0o755)
    (root / "start_mcpb_qm_tmux.sh").write_text(
        """#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
SESSION=${SESSION:-tj_mcpb_m7_qm}
if tmux has-session -t "$SESSION" 2>/dev/null; then
  echo "tmux_session_exists=$SESSION"
  echo "attach=tmux attach -t $SESSION"
  exit 0
fi
tmux new-session -d -s "$SESSION" "cd '$PWD'; env NCORES=${NCORES:-16} bash 02_run_mcpb_gamess_qm.sh 2>&1 | tee logs/tj_mcpb_m7_qm.tmux.log; echo; echo '[tmux] finished'; exec bash"
echo "tmux_session=$SESSION"
echo "attach=tmux attach -t $SESSION"
echo "monitor=cd $PWD && bash monitor_mcpb_gamess_qm.sh"
""",
        encoding="utf-8",
        newline="\n",
    )
    (root / "start_mcpb_qm_tmux.sh").chmod(0o755)

    # Optional immediate fill mode.
    if getattr(args, "fill_only", False):
        coords = extract_last_gamess_coords(root / "gamess_logs" / f"{opt}.log")
        fill_data(root / f"{fc}.inp", root / f"{fc}.inp", coords)
        print(f"filled={root / (fc + '.inp')} atoms={len(coords)}")
    else:
        print(f"runner={root / '02_run_mcpb_gamess_qm.sh'}")
        print(f"monitor={root / 'monitor_mcpb_gamess_qm.sh'}")
        print(f"tmux_starter={root / 'start_mcpb_qm_tmux.sh'}")
    return 0


if __name__ == "__main__":
    # Add a hidden flag without complicating the visible CLI help too much.
    import sys

    if "--fill-only" in sys.argv:
        sys.argv.remove("--fill-only")
        namespace_marker = True
    else:
        namespace_marker = False
    original_parse_args = argparse.ArgumentParser.parse_args

    def parse_args_with_marker(self, *args, **kwargs):
        ns = original_parse_args(self, *args, **kwargs)
        ns.fill_only = namespace_marker
        return ns

    argparse.ArgumentParser.parse_args = parse_args_with_marker
    raise SystemExit(main())
