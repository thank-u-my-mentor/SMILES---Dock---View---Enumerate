#!/usr/bin/env python3
"""Convert a guarded MOREAD GAMESS OPT input into a HUCKEL-start input.

This is useful when the harmonic OPT restraints are accepted but GAMESS fails
while reading a previous UHF $VEC block. It keeps the same geometry, STATPT
restraints, charge, multiplicity, and basis, but removes $VEC and changes
$GUESS to HUCKEL.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def write_text_lf(path: Path, text: str) -> None:
    path.write_bytes(text.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8"))


def remove_vec_and_set_huckel(lines: list[str]) -> list[str]:
    out: list[str] = []
    i = 0
    dollar = "$"
    replaced_guess = False
    while i < len(lines):
        stripped = lines[i].strip().upper()
        if stripped == dollar + "VEC":
            i += 1
            while i < len(lines) and lines[i].strip().upper() != dollar + "END":
                i += 1
            i += 1
            continue
        if stripped.startswith(dollar + "GUESS"):
            out.append(f" {dollar}GUESS GUESS=HUCKEL {dollar}END")
            replaced_guess = True
            i += 1
            continue
        out.append(lines[i])
        i += 1
    if not replaced_guess:
        out.insert(0, f" {dollar}GUESS GUESS=HUCKEL {dollar}END")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mcpb-dir", type=Path, required=True)
    parser.add_argument("--source-job", required=True)
    parser.add_argument("--out-job", required=True)
    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--session", default=None)
    args = parser.parse_args()

    root = args.mcpb_dir
    source = root / f"{args.source_job}.inp"
    if not source.exists():
        raise FileNotFoundError(source)
    out_inp = root / f"{args.out_job}.inp"
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    write_text_lf(out_inp, "\n".join(remove_vec_and_set_huckel(lines)) + "\n")

    source_constraints = root / f"{args.source_job}_harmonic_constraints.tsv"
    out_constraints = root / f"{args.out_job}_harmonic_constraints.tsv"
    if source_constraints.exists():
        shutil.copyfile(source_constraints, out_constraints)

    (root / "gamess_logs").mkdir(exist_ok=True)
    (root / "logs").mkdir(exist_ok=True)

    foreground = root / f"run_{args.out_job}_foreground.sh"
    write_text_lf(
        foreground,
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={args.out_job}
NCORES=${{NCORES:-{args.cores}}}
GAMESS=${{GAMESS:-/home/qin/softwares/gamess/rungms}}
VERSION=${{VERSION:-00}}
mkdir -p gamess_logs logs
echo "[GAMESS] $JOB start $(date) NCORES=$NCORES GUESS=HUCKEL guarded" | tee "gamess_logs/${{JOB}}.status"
set +e
"$GAMESS" "$JOB" "$VERSION" "$NCORES" > "gamess_logs/${{JOB}}.log" 2>&1
rc=$?
set -e
echo "[GAMESS] $JOB exit_code=$rc end $(date)" | tee -a "gamess_logs/${{JOB}}.status"
if grep -q 'SCF IS UNCONVERGED\\|SCF HAS NOT CONVERGED' "gamess_logs/${{JOB}}.log"; then
  echo "[GAMESS] $JOB scf=unconverged" | tee -a "gamess_logs/${{JOB}}.status"
fi
if grep -q 'END OF GEOMETRY SEARCH' "gamess_logs/${{JOB}}.log"; then
  echo "[GAMESS] $JOB geometry=end_of_search" | tee -a "gamess_logs/${{JOB}}.status"
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
    write_text_lf(
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
    write_text_lf(
        monitor,
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={args.out_job}
echo "===== constraints ====="
cat "${{JOB}}_harmonic_constraints.tsv" 2>/dev/null || true
echo
echo "===== status ====="
cat "gamess_logs/${{JOB}}.status" 2>/dev/null || true
echo
echo "===== running ====="
pgrep -af "${{JOB}}|ddikick|gamess.00.x" || true
echo
echo "===== markers ====="
grep -E 'GUESS =|RUNTYP=|HARMONIC|BEGINNING GEOMETRY|NSERCH|GRAD\\. MAX|CONVERGED|SCF IS UNCONVERGED|SCF HAS NOT CONVERGED|FINAL U-B3LYP|S-SQUARED|END OF GEOMETRY SEARCH|TERMINATED' "gamess_logs/${{JOB}}.log" 2>/dev/null | tail -n 160 || true
echo
echo "===== recent log ====="
tail -n 80 "gamess_logs/${{JOB}}.log" 2>/dev/null || true
""",
    )
    monitor.chmod(0o755)

    print(f"input={out_inp}")
    print(f"foreground={foreground}")
    print(f"tmux={tmux_runner}")
    print(f"monitor={monitor}")
    print(f"constraints={out_constraints if out_constraints.exists() else 'none'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
