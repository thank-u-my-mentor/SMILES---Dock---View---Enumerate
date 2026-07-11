#!/usr/bin/env python3
"""Build a GAMESS GRADIENT diagnostic job from a converged MOREAD vector."""

from __future__ import annotations

import argparse
from pathlib import Path


def replace_card(lines: list[str], card: str, replacement: str) -> list[str]:
    out: list[str] = []
    i = 0
    replaced = False
    target = card.upper()
    while i < len(lines):
        line = lines[i]
        if line.strip().upper().startswith(target):
            out.extend(replacement.splitlines())
            replaced = True
            if "$END" not in line.upper():
                i += 1
                while i < len(lines) and "$END" not in lines[i].upper():
                    i += 1
            i += 1
            continue
        out.append(line)
        i += 1
    if not replaced:
        out = replacement.splitlines() + out
    return out


def remove_vec(lines: list[str]) -> list[str]:
    out: list[str] = []
    i = 0
    while i < len(lines):
        if lines[i].strip().upper() == "$VEC":
            i += 1
            while i < len(lines) and lines[i].strip().upper() != "$END":
                i += 1
            i += 1
            continue
        out.append(lines[i])
        i += 1
    return out


def extract_vec(dat_path: Path) -> list[str]:
    lines = dat_path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip().upper() == "$VEC":
            start = i
            break
    if start is None:
        raise SystemExit(f"No $VEC block found in {dat_path}")
    end = None
    for i in range(start + 1, len(lines)):
        if lines[i].strip().upper() == "$END":
            end = i
            break
    if end is None:
        raise SystemExit(f"No $END after $VEC in {dat_path}")
    return lines[start : end + 1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mcpb-dir", type=Path, required=True)
    parser.add_argument("--source-job", required=True)
    parser.add_argument("--out-job", required=True)
    parser.add_argument("--dat", type=Path, required=True)
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--mult", type=int, required=True)
    parser.add_argument("--norb", type=int, required=True)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--session", default=None)
    parser.add_argument("--gamess", default="/home/qin/softwares/gamess/rungms")
    parser.add_argument("--scf-mode", choices=["soscf", "diis"], default="soscf")
    args = parser.parse_args()

    root = args.mcpb_dir
    source = root / f"{args.source_job}.inp"
    if not source.exists():
        raise FileNotFoundError(source)
    vec = extract_vec(args.dat)
    lines = remove_vec(source.read_text(encoding="utf-8", errors="replace").splitlines())

    lines = replace_card(
        lines,
        "$CONTRL",
        (
            " $CONTRL\n"
            "  SCFTYP=UHF DFTTYP=B3LYP RUNTYP=GRADIENT\n"
            f"  ICHARG={args.charge} MULT={args.mult} MAXIT=200 COORD=UNIQUE UNITS=ANGS\n"
            " $END"
        ),
    )
    if args.scf_mode == "soscf":
        scf = (
            " $SCF\n"
            "  DIRSCF=.T. DIIS=.F. SOSCF=.T. DAMP=.T. SHIFT=.T.\n"
            "  ETHRSH=10.0 MAXDII=30\n"
            " $END"
        )
    else:
        scf = (
            " $SCF\n"
            "  DIRSCF=.T. DIIS=.T. DAMP=.T. SHIFT=.T.\n"
            "  ETHRSH=10.0 MAXDII=20\n"
            " $END"
        )
    lines = replace_card(lines, "$SCF", scf)
    lines = replace_card(
        lines,
        "$DFT",
        (
            " $DFT\n"
            "  NRAD=96 NLEB=302 NRAD0=96 NLEB0=302 SWOFF=0.0\n"
            " $END"
        ),
    )
    lines = replace_card(lines, "$GUESS", f" $GUESS GUESS=MOREAD NORB={args.norb} $END")
    out_inp = root / f"{args.out_job}.inp"
    out_inp.write_text("\n".join(lines + vec) + "\n", encoding="utf-8", newline="\n")

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
echo "[GAMESS] $JOB start $(date) NCORES=$NCORES MOREAD={args.dat} runtype=GRADIENT scf_mode={args.scf_mode} dft_grid=fine_from_start" | tee "gamess_logs/${{JOB}}.status"
set +e
"$GAMESS" "$JOB" "$VERSION" "$NCORES" > "gamess_logs/${{JOB}}.log" 2>&1
rc=$?
set -e
echo "[GAMESS] $JOB exit_code=$rc end $(date)" | tee -a "gamess_logs/${{JOB}}.status"
if grep -q 'SCF IS UNCONVERGED\\|SCF HAS NOT CONVERGED' "gamess_logs/${{JOB}}.log"; then
  echo "[GAMESS] $JOB scf=unconverged" | tee -a "gamess_logs/${{JOB}}.status"
fi
if grep -q 'DENSITY CONVERGED' "gamess_logs/${{JOB}}.log"; then
  echo "[GAMESS] $JOB scf=density_converged" | tee -a "gamess_logs/${{JOB}}.status"
fi
if grep -q 'GRADIENT OF THE ENERGY' "gamess_logs/${{JOB}}.log"; then
  echo "[GAMESS] $JOB gradient=printed" | tee -a "gamess_logs/${{JOB}}.status"
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
grep -E 'RUNTYP=|GUESS =|NRAD0|NLEB0|SWOFF|DFT CODE IS SWITCHING|DENSITY CONVERGED|SCF IS UNCONVERGED|SCF HAS NOT CONVERGED|FINAL U-B3LYP|S-SQUARED|GRADIENT|MAXIMUM GRADIENT|TERMINATED' "gamess_logs/${{JOB}}.log" 2>/dev/null | tail -n 120 || true
echo
echo "===== recent log ====="
tail -n 80 "gamess_logs/${{JOB}}.log" 2>/dev/null || true
""",
        encoding="utf-8",
        newline="\n",
    )
    monitor.chmod(0o755)

    print(f"input={out_inp}")
    print(f"vec_lines={len(vec)}")
    print(f"foreground={foreground}")
    print(f"tmux={tmux_runner}")
    print(f"monitor={monitor}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
