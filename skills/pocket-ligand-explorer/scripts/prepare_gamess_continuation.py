#!/usr/bin/env python3
"""Prepare a GAMESS geometry-optimization continuation from an evaluated NSERCH."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from export_gamess_optimized_model import parse_coords_for_nserch, parse_xyz  # noqa: E402

ATOMIC_NUMBERS = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "S": 16,
    "FE": 26,
    "Fe": 26,
}


def write_text_lf(path: Path, text: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def write_input(
    path: Path,
    atoms,
    coords,
    *,
    charge: int,
    mult: int,
    nstep: int,
    title: str,
    mwords: int,
    memddi: int,
) -> None:
    scftyp = "RHF" if mult == 1 else "UHF"
    lines = [
        f" $SYSTEM MWORDS={mwords} MEMDDI={memddi} $END",
        " $CONTRL",
        f"  SCFTYP={scftyp} DFTTYP=B3LYP RUNTYP=OPTIMIZE",
        f"  ICHARG={charge} MULT={mult} MAXIT=200 COORD=CART UNITS=ANGS",
        " $END",
        f" $STATPT NSTEP={nstep} OPTTOL=0.001 $END",
        " $SCF DIRSCF=.T. DIIS=.T. DAMP=.T. ETHRSH=2.0 MAXDII=20 $END",
        " $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=1 $END",
        " $DATA",
        title,
        "C1",
    ]
    for atom, coord_atom in zip(atoms, coords):
        element = "Fe" if atom.element.upper() == "FE" else atom.element.capitalize()
        atomic_number = ATOMIC_NUMBERS.get(atom.element.upper())
        if atomic_number is None:
            raise SystemExit(f"Unsupported element {atom.element}")
        x, y, z = coord_atom.xyz
        lines.append(f"{element:<2s} {float(atomic_number):8.1f} {x:14.6f} {y:14.6f} {z:14.6f}")
    lines += [" $END", ""]
    write_text_lf(path, "\n".join(lines))


def write_scripts(outdir: Path, job_name: str) -> None:
    write_text_lf(
        outdir / "run_gamess_job.sh",
        """#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")"
mkdir -p logs
job="${1:?usage: run_gamess_job.sh JOB [NCORES]}"
GAMESS=${GAMESS:-/home/qin/softwares/gamess/rungms}
VERSION=${VERSION:-00}
NCORES=${2:-${NCORES:-16}}
log="logs/${job}.log"
status="logs/${job}.status"
echo "[GAMESS] ${job} start $(date) NCORES=${NCORES}" | tee "$status"
"$GAMESS" "$job" "$VERSION" "$NCORES" > "$log" 2>&1
rc=$?
echo "[GAMESS] ${job} exit_code=${rc} end $(date)" | tee -a "$status"
if grep -q 'EXECUTION OF GAMESS TERMINATED NORMALLY' "$log"; then
  echo "[GAMESS] ${job} status=normal" | tee -a "$status"
elif grep -q 'EXECUTION OF GAMESS TERMINATED -ABNORMALLY-' "$log"; then
  echo "[GAMESS] ${job} status=abnormal" | tee -a "$status"
else
  echo "[GAMESS] ${job} status=no_termination_marker" | tee -a "$status"
fi
exit "$rc"
""",
    )
    write_text_lf(
        outdir / "run_continuation.sh",
        f"""#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")"
NCORES=${{NCORES:-16}}
bash run_gamess_job.sh {job_name} "$NCORES"
""",
    )
    write_text_lf(
        outdir / "start_continuation_tmux.sh",
        """#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
SESSION=${SESSION:-feiii_m7_cont60}
NCORES=${NCORES:-16}
tmux has-session -t "$SESSION" 2>/dev/null && {
  echo "tmux_session_exists=$SESSION"
  echo "attach=tmux attach -t $SESSION"
  exit 0
}
tmux new-session -d -s "$SESSION" "cd '$PWD'; env NCORES='$NCORES' bash run_continuation.sh; echo; echo '[tmux] finished'; exec bash"
echo "tmux_session=$SESSION"
echo "attach=tmux attach -t $SESSION"
echo "monitor=bash monitor_continuation.sh"
""",
    )
    write_text_lf(
        outdir / "monitor_continuation.sh",
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
log="logs/{job_name}.log"
status="logs/{job_name}.status"
test -f "$status" && cat "$status" || true
test -f "$log" && grep -E 'RUNTYP=|BEGINNING GEOMETRY SEARCH|NSERCH:|TOTAL ENERGY|S-SQUARED|END OF GEOMETRY SEARCH|FAILURE|TERMINATED' "$log" | tail -n 60 || true
""",
    )
    for script in (
        "run_gamess_job.sh",
        "run_continuation.sh",
        "start_continuation_tmux.sh",
        "monitor_continuation.sh",
    ):
        (outdir / script).chmod(0o755)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--source-log", default="logs/FEIII_radical_M7_fastopt.log")
    parser.add_argument("--initial-xyz", default="feiii_radical_truncated_model.xyz")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--from-nserch", type=int, default=None)
    parser.add_argument("--nstep", type=int, default=60)
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--mult", type=int, default=7)
    parser.add_argument("--job-name", default="FEIII_radical_M7_cont60")
    parser.add_argument("--mwords", type=int, default=50)
    parser.add_argument("--memddi", type=int, default=10)
    args = parser.parse_args()

    workdir = Path(args.workdir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(exist_ok=True)

    atoms = parse_xyz(workdir / args.initial_xyz)
    nserch, coords, metadata = parse_coords_for_nserch(
        workdir / args.source_log,
        len(atoms),
        args.from_nserch,
    )
    title = f"FeIII radical M7 continuation from evaluated NSERCH={nserch}"
    write_input(
        outdir / f"{args.job_name}.inp",
        atoms,
        coords,
        charge=args.charge,
        mult=args.mult,
        nstep=args.nstep,
        title=title,
        mwords=args.mwords,
        memddi=args.memddi,
    )
    write_scripts(outdir, args.job_name)
    write_text_lf(
        outdir / "README_continuation.md",
        "\n".join(
            [
                "# M7 Continuation",
                "",
                f"- Source log: `{workdir / args.source_log}`",
                f"- Starting evaluated NSERCH: `{nserch}`",
                f"- Previous energy: `{metadata.get('energy_hartree', 'NA')}` Hartree",
                f"- Previous grad max/RMS: `{metadata.get('grad_max', 'NA')}` / `{metadata.get('grad_rms', 'NA')}`",
                f"- New job: `{args.job_name}`",
                f"- New `NSTEP`: `{args.nstep}`",
                f"- Charge/multiplicity: `{args.charge}` / `{args.mult}`",
                "",
                "Run:",
                "",
                "```bash",
                f"cd {outdir.as_posix()}",
                "env NCORES=16 SESSION=feiii_m7_cont60 bash start_continuation_tmux.sh",
                "bash monitor_continuation.sh",
                "```",
                "",
                "This continuation uses the last evaluated geometry, not the unknown next predicted geometry printed after an NSTEP-limit warning.",
            ]
        )
        + "\n",
    )

    print(f"outdir={outdir}")
    print(f"job={args.job_name}")
    print(f"input={outdir / (args.job_name + '.inp')}")
    print(f"from_nserch={nserch}")
    print(f"previous_energy={metadata.get('energy_hartree', 'NA')}")
    print(f"previous_grad_max={metadata.get('grad_max', 'NA')}")
    print(f"previous_grad_rms={metadata.get('grad_rms', 'NA')}")
    print(f"run={outdir / 'start_continuation_tmux.sh'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
