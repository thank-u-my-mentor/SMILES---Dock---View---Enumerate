#!/usr/bin/env python3
"""Write a GAMESS SCF runner that tries one job, then a fallback job.

The generated script is intentionally plain bash because it is meant to run
inside tmux on WSL where long nested quoting has caused false failures.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def q(text: str) -> str:
    return "'" + text.replace("'", "'\"'\"'") + "'"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mcpb-dir", required=True)
    parser.add_argument("--first-job", required=True)
    parser.add_argument("--fallback-job", required=True)
    parser.add_argument("--session", required=True)
    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument(
        "--prefix",
        default="run_m5_then_m7_scf",
        help="Output script prefix written under --mcpb-dir.",
    )
    args = parser.parse_args()

    mcpb_dir = Path(args.mcpb_dir).resolve()
    first_fg = mcpb_dir / f"run_{args.first_job}_foreground.sh"
    fallback_fg = mcpb_dir / f"run_{args.fallback_job}_foreground.sh"
    if not first_fg.exists():
        raise SystemExit(f"missing foreground script: {first_fg}")
    if not fallback_fg.exists():
        raise SystemExit(f"missing fallback foreground script: {fallback_fg}")

    first_log = f"gamess_logs/{args.first_job}.log"
    fallback_log = f"gamess_logs/{args.fallback_job}.log"
    chain_status = f"gamess_logs/{args.first_job}_to_{args.fallback_job}.chain_status"

    foreground = mcpb_dir / f"{args.prefix}_foreground.sh"
    tmux_script = mcpb_dir / f"{args.prefix}_tmux.sh"
    monitor = mcpb_dir / f"monitor_{args.prefix}.sh"

    foreground.write_text(
        f"""#!/usr/bin/env bash
set -uo pipefail
cd {q(str(mcpb_dir))}
mkdir -p gamess_logs
STATUS={q(chain_status)}
: > "$STATUS"
echo "[chain] first={args.first_job} fallback={args.fallback_job} cores=${{NCORES:-{args.cores}}}" | tee -a "$STATUS"
echo "[chain] starting first job at $(date)" | tee -a "$STATUS"
env NCORES="${{NCORES:-{args.cores}}}" bash {q(first_fg.name)}
FIRST_RC=$?
echo "[chain] first exit=$FIRST_RC at $(date)" | tee -a "$STATUS"
if grep -q "DENSITY CONVERGED" {q(first_log)} && ! grep -q "SCF IS UNCONVERGED\\|SCF HAS NOT CONVERGED" {q(first_log)}; then
  echo "[chain] first job reached usable DENSITY CONVERGED; fallback skipped" | tee -a "$STATUS"
  exit 0
fi
echo "[chain] first job not usable; launching fallback at $(date)" | tee -a "$STATUS"
env NCORES="${{NCORES:-{args.cores}}}" bash {q(fallback_fg.name)}
FALLBACK_RC=$?
echo "[chain] fallback exit=$FALLBACK_RC at $(date)" | tee -a "$STATUS"
if grep -q "DENSITY CONVERGED" {q(fallback_log)} && ! grep -q "SCF IS UNCONVERGED\\|SCF HAS NOT CONVERGED" {q(fallback_log)}; then
  echo "[chain] fallback job reached usable DENSITY CONVERGED" | tee -a "$STATUS"
  exit 0
fi
echo "[chain] fallback also did not reach usable DENSITY CONVERGED" | tee -a "$STATUS"
exit 1
""",
        encoding="utf-8",
    )
    tmux_script.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd {q(str(mcpb_dir))}
tmux has-session -t {q(args.session)} 2>/dev/null && tmux kill-session -t {q(args.session)}
tmux new-session -d -s {q(args.session)} "cd {q(str(mcpb_dir))}; env NCORES={args.cores} bash {q(foreground.name)}; echo; echo '[chain] finished'; exec bash"
echo "session={args.session}"
echo "status={chain_status}"
""",
        encoding="utf-8",
    )
    monitor.write_text(
        f"""#!/usr/bin/env bash
set -uo pipefail
cd {q(str(mcpb_dir))}
echo "== chain status =="
test -f {q(chain_status)} && tail -n 40 {q(chain_status)} || true
echo
for job in {q(args.first_job)} {q(args.fallback_job)}; do
  log="gamess_logs/${{job}}.log"
  echo "== $job =="
  if [ -f "$log" ]; then
    /mnt/l/WSL/softwares/conda_envs/md/bin/python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/summarize_gamess_scf.py "$log" --last 20 || tail -n 40 "$log"
  else
    echo "log not created yet: $log"
  fi
  echo
done
echo "== live GAMESS/DDI processes =="
pgrep -af "{args.first_job}|{args.fallback_job}|ddikick.x|gamess.00.x" || true
""",
        encoding="utf-8",
    )
    for path in (foreground, tmux_script, monitor):
        path.chmod(0o755)
    print(f"foreground={foreground}")
    print(f"tmux={tmux_script}")
    print(f"monitor={monitor}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
