#!/usr/bin/env python3
"""Build a GAMESS HESSIAN job from an evaluated OPT geometry plus MOREAD.

The input geometry is taken from an evaluated ``NSERCH`` point in a previous
GAMESS optimization log. The orbital vector is read from the matching restart
``.dat`` file. This avoids hand-copying coordinates and keeps the Fe-radical
fine-grid/MOREAD convention consistent across MCPB projects.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from prepare_gamess_moread_optimize_from_log import (
    extract_vec,
    parse_data_block,
    remove_vec,
    replace_card,
    replace_data_block,
    write_text,
)
from export_gamess_optimized_model import parse_coords_for_nserch


def remove_card(lines: list[str], card: str) -> list[str]:
    out: list[str] = []
    i = 0
    target = card.upper()
    while i < len(lines):
        line = lines[i]
        if line.strip().upper().startswith(target):
            if "$END" not in line.upper():
                i += 1
                while i < len(lines) and "$END" not in lines[i].upper():
                    i += 1
            i += 1
            continue
        out.append(line)
        i += 1
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mcpb-dir", type=Path, required=True)
    parser.add_argument("--source-job", required=True)
    parser.add_argument("--source-log", type=Path, required=True)
    parser.add_argument("--out-job", required=True)
    parser.add_argument("--dat", type=Path, required=True)
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--mult", type=int, required=True)
    parser.add_argument("--norb", type=int, required=True)
    parser.add_argument("--from-nserch", type=int, required=True)
    parser.add_argument("--cores", type=int, default=16)
    parser.add_argument("--session", default=None)
    parser.add_argument("--gamess", default="/home/qin/softwares/gamess/rungms")
    parser.add_argument("--scf-mode", choices=["soscf", "diis"], default="soscf")
    args = parser.parse_args()

    root = args.mcpb_dir
    source = root / f"{args.source_job}.inp"
    if not source.exists():
        raise FileNotFoundError(source)
    lines = remove_vec(source.read_text(encoding="utf-8", errors="replace").splitlines())
    _, _, _, atom_lines = parse_data_block(lines)
    nserch, coords, _metadata = parse_coords_for_nserch(
        args.source_log, len(atom_lines), args.from_nserch
    )
    vec = extract_vec(args.dat)

    lines = replace_data_block(lines, coords, f"hessian from evaluated NSERCH={nserch}")
    lines = remove_card(lines, "$STATPT")
    lines = replace_card(lines, "$SYSTEM", " $SYSTEM MEMDDI=400 MWORDS=200 $END")
    lines = replace_card(
        lines,
        "$CONTRL",
        (
            " $CONTRL\n"
            "  SCFTYP=UHF DFTTYP=B3LYP RUNTYP=HESSIAN\n"
            f"  ICHARG={args.charge} MULT={args.mult} MAXIT=200 COORD=UNIQUE UNITS=ANGS\n"
            " $END"
        ),
    )
    if args.scf_mode == "soscf":
        scf_card = (
            " $SCF\n"
            "  DIRSCF=.T. DIIS=.F. SOSCF=.T. DAMP=.T. SHIFT=.T.\n"
            "  ETHRSH=10.0 MAXDII=30\n"
            " $END"
        )
    else:
        scf_card = (
            " $SCF\n"
            "  DIRSCF=.T. DIIS=.T. DAMP=.T. SHIFT=.T.\n"
            "  ETHRSH=10.0 MAXDII=20\n"
            " $END"
        )
    lines = replace_card(lines, "$SCF", scf_card)
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
    write_text(out_inp, "\n".join(lines + vec) + "\n")

    (root / "gamess_logs").mkdir(exist_ok=True)
    (root / "logs").mkdir(exist_ok=True)
    foreground = root / f"run_{args.out_job}_foreground.sh"
    write_text(
        foreground,
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={args.out_job}
NCORES=${{NCORES:-{args.cores}}}
GAMESS=${{GAMESS:-{args.gamess}}}
VERSION=${{VERSION:-00}}
mkdir -p gamess_logs logs
echo "[GAMESS] $JOB start $(date) NCORES=$NCORES from_job={args.source_job} from_nserch={nserch} MOREAD={args.dat} runtype=HESSIAN scf_mode={args.scf_mode} dft_grid=fine_from_start no_harmonic_constraints" | tee "gamess_logs/${{JOB}}.status"
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
if grep -q 'FREQUENCIES IN CM' "gamess_logs/${{JOB}}.log"; then
  echo "[GAMESS] $JOB hessian=frequencies_printed" | tee -a "gamess_logs/${{JOB}}.status"
fi
if grep -q 'EXECUTION OF GAMESS TERMINATED NORMALLY' "gamess_logs/${{JOB}}.log"; then
  echo "[GAMESS] $JOB status=normal" | tee -a "gamess_logs/${{JOB}}.status"
else
  echo "[GAMESS] $JOB status=not_normal" | tee -a "gamess_logs/${{JOB}}.status"
fi
exit "$rc"
""",
    )
    foreground.chmod(0o755)

    session = args.session or args.out_job.lower()
    tmux_runner = root / f"run_{args.out_job}_tmux.sh"
    write_text(
        tmux_runner,
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
    )
    tmux_runner.chmod(0o755)

    monitor = root / f"monitor_{args.out_job}.sh"
    write_text(
        monitor,
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
grep -E 'RUNTYP=|GUESS =|NRAD0|NLEB0|DENSITY CONVERGED|SCF IS UNCONVERGED|SCF HAS NOT CONVERGED|FINAL U-B3LYP|S-SQUARED|HESSIAN|FORCE CONSTANT|FREQUENCIES|IMAGINARY|TERMINATED' "gamess_logs/${{JOB}}.log" 2>/dev/null | tail -n 160 || true
echo
echo "===== recent log ====="
tail -n 100 "gamess_logs/${{JOB}}.log" 2>/dev/null || true
""",
    )
    monitor.chmod(0o755)

    print(f"input={out_inp}")
    print(f"from_nserch={nserch}")
    print(f"vec_lines={len(vec)}")
    print(f"foreground={foreground}")
    print(f"tmux={tmux_runner}")
    print(f"monitor={monitor}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
