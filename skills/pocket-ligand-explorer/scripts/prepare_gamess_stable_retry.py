#!/usr/bin/env python3
"""Prepare a more conservative GAMESS retry input from an existing MCPB input."""

from __future__ import annotations

import argparse
from pathlib import Path


def replace_card(lines: list[str], prefix: str, replacement: str) -> list[str]:
    out: list[str] = []
    replaced = False
    for line in lines:
        if line.strip().upper().startswith(prefix):
            out.extend(replacement.splitlines())
            replaced = True
        else:
            out.append(line)
    if not replaced:
        out = replacement.splitlines() + out
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a stable GAMESS retry and tmux runner.")
    parser.add_argument("--mcpb-dir", type=Path, required=True)
    parser.add_argument("--source-job", required=True)
    parser.add_argument("--retry-job", required=True)
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--mult", type=int, required=True)
    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--session", default="zfy_mcpb_m7_stable")
    parser.add_argument("--gamess", default="/home/qin/softwares/gamess/rungms")
    parser.add_argument("--nstep", type=int, default=300)
    args = parser.parse_args()

    root = args.mcpb_dir
    source = root / f"{args.source_job}.inp"
    retry = root / f"{args.retry_job}.inp"
    if not source.exists():
        raise FileNotFoundError(source)

    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    lines = replace_card(lines, "$SYSTEM", " $SYSTEM MEMDDI=200 MWORDS=120 $END")
    lines = replace_card(
        lines,
        "$CONTRL",
        (
            " $CONTRL\n"
            "  SCFTYP=UHF DFTTYP=B3LYP RUNTYP=OPTIMIZE\n"
            f"  ICHARG={args.charge} MULT={args.mult} MAXIT=200\n"
            " $END"
        ),
    )
    lines = replace_card(
        lines,
        "$SCF",
        (
            " $SCF\n"
            "  DIRSCF=.T. DIIS=.T. DAMP=.T. SHIFT=.T.\n"
            "  ETHRSH=10.0 MAXDII=30\n"
            " $END"
        ),
    )
    lines = replace_card(lines, "$STATPT", f" $STATPT NSTEP={args.nstep} OPTTOL=0.0002 $END")
    retry.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")

    (root / "gamess_logs").mkdir(exist_ok=True)
    foreground = root / f"run_{args.retry_job}_foreground.sh"
    foreground.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={args.retry_job}
NCORES=${{NCORES:-{args.cores}}}
GAMESS=${{GAMESS:-{args.gamess}}}
VERSION=${{VERSION:-00}}
test -s "${{JOB}}.inp" || {{ echo "missing=${{JOB}}.inp"; exit 2; }}
mkdir -p gamess_logs logs
echo "[GAMESS] $JOB start $(date) NCORES=$NCORES" | tee "gamess_logs/${{JOB}}.status"
set +e
"$GAMESS" "$JOB" "$VERSION" "$NCORES" > "gamess_logs/${{JOB}}.log" 2>&1
rc=$?
set -e
echo "[GAMESS] $JOB exit_code=$rc end $(date)" | tee -a "gamess_logs/${{JOB}}.status"
if grep -q 'EXECUTION OF GAMESS TERMINATED NORMALLY' "gamess_logs/${{JOB}}.log"; then
  echo "[GAMESS] $JOB status=normal" | tee -a "gamess_logs/${{JOB}}.status"
else
  echo "[GAMESS] $JOB status=not_normal" | tee -a "gamess_logs/${{JOB}}.status"
fi
exit "$rc"
""",
        encoding="utf-8",
        newline="\n",
    )
    foreground.chmod(0o755)

    tmux_runner = root / f"run_{args.retry_job}_tmux.sh"
    tmux_runner.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
SESSION=${{SESSION:-{args.session}}}
NCORES=${{NCORES:-{args.cores}}}
if tmux has-session -t "$SESSION" 2>/dev/null; then
  echo "tmux_session_exists=$SESSION"
  echo "monitor=cd $PWD && bash monitor_{args.retry_job}.sh"
  exit 0
fi
tmux new-session -d -s "$SESSION" "cd '$PWD'; env NCORES='$NCORES' bash {foreground.name}; echo; echo '[tmux] finished'; exec bash"
echo "tmux_session=$SESSION"
echo "monitor=cd $PWD && bash monitor_{args.retry_job}.sh"
echo "tail=tail -f $PWD/gamess_logs/{args.retry_job}.log"
""",
        encoding="utf-8",
        newline="\n",
    )
    tmux_runner.chmod(0o755)

    monitor = root / f"monitor_{args.retry_job}.sh"
    monitor.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={args.retry_job}
echo "===== status ====="
cat "gamess_logs/${{JOB}}.status" 2>/dev/null || true
echo
echo "===== recent log ====="
tail -n 80 "gamess_logs/${{JOB}}.log" 2>/dev/null || true
echo
echo "===== markers ====="
grep -E 'RUNTYP=|ITER EX|NSERCH:|TOTAL ENERGY|S-SQUARED|END OF GEOMETRY SEARCH|FAILURE|ERROR|TERMINATED' "gamess_logs/${{JOB}}.log" 2>/dev/null | tail -n 80 || true
""",
        encoding="utf-8",
        newline="\n",
    )
    monitor.chmod(0o755)

    print(f"retry_input={retry}")
    print(f"foreground={foreground}")
    print(f"tmux_runner={tmux_runner}")
    print(f"monitor={monitor}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
