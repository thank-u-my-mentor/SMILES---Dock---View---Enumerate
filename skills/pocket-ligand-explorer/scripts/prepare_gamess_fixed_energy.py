#!/usr/bin/env python3
"""Create a fresh fixed-geometry GAMESS ENERGY/SCF job from an MCPB input.

This is for Fe/radical MCPB workflows where an earlier $VEC is either missing
or not trustworthy. It intentionally starts from the input coordinates and a
fresh HUCKEL guess, while forcing the fine DFT grid from the first iteration.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def _starts_card(line: str, card: str) -> bool:
    return line.strip().upper().startswith(card.upper())


def _skip_card(lines: list[str], i: int) -> int:
    """Return the first index after a GAMESS card starting at i."""
    if "$END" in lines[i].upper():
        return i + 1
    i += 1
    while i < len(lines) and "$END" not in lines[i].upper():
        i += 1
    return min(i + 1, len(lines))


def remove_cards(lines: list[str], cards: list[str]) -> list[str]:
    out: list[str] = []
    i = 0
    upper_cards = [c.upper() for c in cards]
    while i < len(lines):
        if any(_starts_card(lines[i], card) for card in upper_cards):
            i = _skip_card(lines, i)
            continue
        out.append(lines[i])
        i += 1
    return out


def upsert_card(lines: list[str], card: str, replacement: str) -> list[str]:
    out: list[str] = []
    i = 0
    replaced = False
    while i < len(lines):
        if _starts_card(lines[i], card):
            out.extend(replacement.splitlines())
            replaced = True
            i = _skip_card(lines, i)
            continue
        out.append(lines[i])
        i += 1
    if replaced:
        return out

    insert_at = next(
        (idx for idx, line in enumerate(out) if _starts_card(line, "$DATA")),
        len(out),
    )
    return out[:insert_at] + replacement.splitlines() + out[insert_at:]


def scf_card(mode: str) -> str:
    if mode == "soscf":
        return (
            " $SCF\n"
            "  DIRSCF=.T. DIIS=.F. SOSCF=.T. DAMP=.T. SHIFT=.T.\n"
            "  ETHRSH=10.0 MAXDII=30\n"
            " $END"
        )
    if mode == "damp":
        return (
            " $SCF\n"
            "  DIRSCF=.T. DIIS=.F. SOSCF=.F. DAMP=.T. SHIFT=.T.\n"
            " $END"
        )
    if mode == "very-damp":
        return (
            " $SCF\n"
            "  DIRSCF=.T. DIIS=.F. SOSCF=.F. DAMP=.T. SHIFT=.T.\n"
            "  ETHRSH=0.0 MAXDII=1\n"
            " $END"
        )
    return (
        " $SCF\n"
        "  DIRSCF=.T. DIIS=.T. SOSCF=.F. DAMP=.T. SHIFT=.T.\n"
        "  ETHRSH=10.0 MAXDII=30\n"
        " $END"
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Create a fresh fixed-geometry GAMESS ENERGY input and tmux/tail "
            "helpers from an MCPB-generated small_opt input."
        )
    )
    parser.add_argument("--mcpb-dir", type=Path, required=True)
    parser.add_argument("--source-job", required=True)
    parser.add_argument("--out-job", required=True)
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--mult", type=int, required=True)
    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--session", default=None)
    parser.add_argument("--maxit", type=int, default=200)
    parser.add_argument("--memddi", type=int, default=400)
    parser.add_argument("--mwords", type=int, default=200)
    parser.add_argument(
        "--scf-mode",
        choices=["diis", "soscf", "damp", "very-damp"],
        default="diis",
    )
    parser.add_argument("--gamess", default="/home/qin/softwares/gamess/rungms")
    args = parser.parse_args()
    if not (0 <= args.maxit <= 200):
        raise SystemExit(
            f"GAMESS $CONTRL MAXIT must be between 0 and 200 for this build; got {args.maxit}."
        )

    root = args.mcpb_dir
    source = root / f"{args.source_job}.inp"
    if not source.exists():
        raise FileNotFoundError(source)

    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    lines = remove_cards(lines, ["$STATPT", "$GUESS", "$VEC", "$ELPOT", "$PDC"])
    lines = upsert_card(lines, "$SYSTEM", f" $SYSTEM MEMDDI={args.memddi} MWORDS={args.mwords} $END")
    lines = upsert_card(
        lines,
        "$CONTRL",
        (
            " $CONTRL\n"
            "  SCFTYP=UHF DFTTYP=B3LYP RUNTYP=ENERGY\n"
            f"  ICHARG={args.charge} MULT={args.mult} MAXIT={args.maxit}\n"
            "  COORD=UNIQUE UNITS=ANGS\n"
            " $END"
        ),
    )
    lines = upsert_card(lines, "$SCF", scf_card(args.scf_mode))
    lines = upsert_card(
        lines,
        "$DFT",
        (
            " $DFT\n"
            "  NRAD=96 NLEB=302 NRAD0=96 NLEB0=302 SWOFF=0.0\n"
            " $END"
        ),
    )
    lines = upsert_card(lines, "$GUESS", " $GUESS GUESS=HUCKEL $END")

    out_inp = root / f"{args.out_job}.inp"
    out_inp.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")

    (root / "gamess_logs").mkdir(exist_ok=True)
    (root / "logs").mkdir(exist_ok=True)

    foreground = root / f"run_{args.out_job}_foreground.sh"
    foreground.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={args.out_job}
NCORES=${{NCORES:-{args.cores}}}
GAMESS=${{GAMESS:-{args.gamess}}}
VERSION=${{VERSION:-00}}
mkdir -p gamess_logs logs
echo "[GAMESS] $JOB start $(date) NCORES=$NCORES fresh_guess=HUCKEL dft_grid=fine_from_start scf_mode={args.scf_mode}" | tee "gamess_logs/${{JOB}}.status"
set +e
"$GAMESS" "$JOB" "$VERSION" "$NCORES" > "gamess_logs/${{JOB}}.log" 2>&1
rc=$?
set -e
echo "[GAMESS] $JOB exit_code=$rc end $(date)" | tee -a "gamess_logs/${{JOB}}.status"
if grep -q 'DENSITY CONVERGED' "gamess_logs/${{JOB}}.log"; then
  echo "[GAMESS] $JOB scf=density_converged" | tee -a "gamess_logs/${{JOB}}.status"
fi
if grep -q 'SCF IS UNCONVERGED\\|SCF HAS NOT CONVERGED' "gamess_logs/${{JOB}}.log"; then
  echo "[GAMESS] $JOB scf=unconverged" | tee -a "gamess_logs/${{JOB}}.status"
fi
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

    session = args.session or args.out_job.lower()
    tmux_runner = root / f"run_{args.out_job}_tmux.sh"
    tmux_runner.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
SESSION=${{SESSION:-{session}}}
NCORES=${{NCORES:-{args.cores}}}
if tmux has-session -t "$SESSION" 2>/dev/null; then
  echo "tmux_session_exists=$SESSION"
  echo "monitor=cd $PWD && bash monitor_{args.out_job}.sh"
  echo "tail=tail -f $PWD/gamess_logs/{args.out_job}.log"
  exit 0
fi
tmux new-session -d -s "$SESSION" "cd '$PWD'; env NCORES='$NCORES' bash {foreground.name}; echo; echo '[tmux] finished'; exec bash"
echo "tmux_session=$SESSION"
echo "monitor=cd $PWD && bash monitor_{args.out_job}.sh"
echo "tail=tail -f $PWD/gamess_logs/{args.out_job}.log"
""",
        encoding="utf-8",
        newline="\n",
    )
    tmux_runner.chmod(0o755)

    monitor = root / f"monitor_{args.out_job}.sh"
    monitor.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={args.out_job}
echo "===== status ====="
cat "gamess_logs/${{JOB}}.status" 2>/dev/null || true
echo
echo "===== running ====="
pgrep -af "${{JOB}}|ddikick|gamess.00.x" || true
echo
echo "===== markers ====="
grep -E 'GUESS =|RUNTYP=|NRAD0|NLEB0|DFT CODE IS SWITCHING|ITER EX|DENSITY CONVERGED|SCF IS UNCONVERGED|SCF HAS NOT CONVERGED|FINAL U-B3LYP|S-SQUARED|TERMINATED' "gamess_logs/${{JOB}}.log" 2>/dev/null | tail -n 140 || true
echo
echo "===== recent log ====="
tail -n 90 "gamess_logs/${{JOB}}.log" 2>/dev/null || true
""",
        encoding="utf-8",
        newline="\n",
    )
    monitor.chmod(0o755)

    guarded = root / f"run_{args.out_job}_guarded.sh"
    guarded.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={args.out_job}
NCORES=${{NCORES:-{args.cores}}}
PYTHON=${{PYTHON:-/mnt/l/WSL/softwares/conda_envs/md/bin/python}}
GUARD=${{GUARD:-/mnt/e/Codex/skills/pocket-ligand-explorer/scripts/guard_gamess_scf.py}}
LOG="gamess_logs/${{JOB}}.log"
STATUS="gamess_logs/${{JOB}}.guard.status"

echo "[guarded] launching $JOB with NCORES=$NCORES"
env NCORES="$NCORES" bash {foreground.name} &
run_pid=$!

set +e
"$PYTHON" "$GUARD" \\
  --log "$LOG" \\
  --kill-pattern "$JOB" \\
  --status "$STATUS" \\
  --poll 30 \\
  --kill
guard_rc=$?
wait "$run_pid"
run_rc=$?
set -e

echo "[guarded] guard_exit=$guard_rc run_exit=$run_rc"
if [[ "$guard_rc" -ne 0 ]]; then
  exit "$guard_rc"
fi
exit "$run_rc"
""",
        encoding="utf-8",
        newline="\n",
    )
    guarded.chmod(0o755)

    print(f"input={out_inp}")
    print(f"foreground={foreground}")
    print(f"guarded={guarded}")
    print(f"tmux={tmux_runner}")
    print(f"monitor={monitor}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
