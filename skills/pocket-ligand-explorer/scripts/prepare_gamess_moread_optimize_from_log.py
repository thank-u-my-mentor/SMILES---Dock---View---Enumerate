#!/usr/bin/env python3
"""Continue a GAMESS OPTIMIZE job from an evaluated NSERCH geometry plus MOREAD.

This is safer than reusing an old input file directly: the new job starts from
the last evaluated geometry in a previous GAMESS log and reads the converged
orbital vector from a matching restart ``.dat`` file.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from export_gamess_optimized_model import parse_coords_for_nserch  # noqa: E402


ATOMIC_NUMBERS = {
    "H": 1.0,
    "C": 6.0,
    "N": 7.0,
    "O": 8.0,
    "S": 16.0,
    "FE": 26.0,
}


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


def parse_data_block(lines: list[str]) -> tuple[int, int, list[str], list[str]]:
    start = None
    for i, line in enumerate(lines):
        if line.strip().upper() == "$DATA":
            start = i
            break
    if start is None:
        raise SystemExit("No $DATA block found in source input")
    end = None
    for i in range(start + 1, len(lines)):
        if lines[i].strip().upper() == "$END":
            end = i
            break
    if end is None:
        raise SystemExit("No $END after $DATA in source input")
    header = lines[start + 1 : min(start + 3, end)]
    atom_lines = lines[start + 3 : end]
    atoms = [line for line in atom_lines if len(line.split()) >= 5]
    if not atoms:
        raise SystemExit("No atom lines found in $DATA block")
    return start, end, header, atoms


def element_token(element: str) -> str:
    return "Fe" if element.upper() == "FE" else element.capitalize()


def replace_data_block(lines: list[str], coords, title_suffix: str) -> list[str]:
    start, end, header, atom_lines = parse_data_block(lines)
    if len(atom_lines) != len(coords):
        raise SystemExit(
            f"$DATA atom count {len(atom_lines)} does not match log coordinate count {len(coords)}"
        )
    title = header[0] if header else "GAMESS continuation"
    symmetry = header[1] if len(header) > 1 else "C1"
    new_block = [" $DATA", f"{title} {title_suffix}".strip(), symmetry]
    for source_line, coord_atom in zip(atom_lines, coords):
        parts = source_line.split()
        source_element = parts[0]
        element = element_token(coord_atom.element or source_element)
        atomic_number = None
        if len(parts) >= 2:
            try:
                atomic_number = float(parts[1])
            except ValueError:
                atomic_number = None
        if atomic_number is None:
            atomic_number = ATOMIC_NUMBERS.get(coord_atom.element.upper())
        if atomic_number is None:
            raise SystemExit(f"Unsupported element in $DATA continuation: {coord_atom.element}")
        x, y, z = coord_atom.xyz
        new_block.append(f"{element:<2s} {atomic_number:8.1f} {x:14.6f} {y:14.6f} {z:14.6f}")
    new_block.append(" $END")
    return lines[:start] + new_block + lines[end + 1 :]


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8", newline="\n")


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
    parser.add_argument("--from-nserch", type=int, default=None)
    parser.add_argument("--nstep", type=int, default=30)
    parser.add_argument("--opttol", default="0.0002")
    parser.add_argument("--scf-mode", choices=["diis", "soscf"], default="soscf")
    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--session", default=None)
    parser.add_argument(
        "--harmonic-distance",
        nargs=4,
        action="append",
        metavar=("ATOM_I", "ATOM_J", "TARGET_A", "FORCE_KCAL_MOL"),
        help=(
            "Add one GAMESS harmonic distance guard to $STATPT. "
            "Example: --harmonic-distance 1 4 1.54 250"
        ),
    )
    args = parser.parse_args()

    root = args.mcpb_dir
    source = root / f"{args.source_job}.inp"
    if not source.exists():
        raise FileNotFoundError(source)
    lines = remove_vec(source.read_text(encoding="utf-8", errors="replace").splitlines())
    _, _, _, atom_lines = parse_data_block(lines)
    nserch, coords, metadata = parse_coords_for_nserch(args.source_log, len(atom_lines), args.from_nserch)
    vec = extract_vec(args.dat)

    lines = replace_data_block(lines, coords, f"continuation from NSERCH={nserch}")
    lines = replace_card(
        lines,
        "$CONTRL",
        (
            " $CONTRL\n"
            "  SCFTYP=UHF DFTTYP=B3LYP RUNTYP=OPTIMIZE\n"
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
    statpt_lines = [f" $STATPT NSTEP={args.nstep} OPTTOL={args.opttol}"]
    guard_descriptions: list[str] = []
    if args.harmonic_distance:
        ihmcon: list[str] = []
        shmcon: list[str] = []
        fhmcon: list[str] = []
        for idx, (atom_i, atom_j, target, force) in enumerate(args.harmonic_distance, start=1):
            try:
                ai = int(atom_i)
                aj = int(atom_j)
                target_f = float(target)
                force_f = float(force)
            except ValueError as exc:
                raise SystemExit(f"Invalid --harmonic-distance values: {(atom_i, atom_j, target, force)}") from exc
            ihmcon.extend(["1", str(ai), str(aj)])
            shmcon.append(f"{target_f:.4f}")
            fhmcon.append(f"{force_f:.1f}")
            guard_descriptions.append(f"{idx}: atoms {ai}-{aj} target {target_f:.4f} A force {force_f:.1f} kcal/mol")
        statpt_lines.append(f"  IHMCON(1)={','.join(ihmcon)}")
        statpt_lines.append(f"  SHMCON(1)={','.join(shmcon)}")
        statpt_lines.append(f"  FHMCON(1)={','.join(fhmcon)}")
    statpt_lines.append(" $END")
    lines = replace_card(lines, "$STATPT", "\n".join(statpt_lines))
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
GAMESS=${{GAMESS:-/home/qin/softwares/gamess/rungms}}
VERSION=${{VERSION:-00}}
mkdir -p gamess_logs logs
echo "[GAMESS] $JOB start $(date) NCORES=$NCORES from_nserch={nserch} MOREAD={args.dat} scf_mode={args.scf_mode} dft_grid=fine_from_start" | tee "gamess_logs/${{JOB}}.status"
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
grep -E 'GUESS =|RUNTYP=|NRAD0|NLEB0|SWOFF|DFT CODE IS SWITCHING|BEGINNING GEOMETRY|NSERCH|GRAD\\. MAX|CONVERGED|SCF IS UNCONVERGED|SCF HAS NOT CONVERGED|FINAL U-B3LYP|S-SQUARED|END OF GEOMETRY SEARCH|THE GEOMETRY SEARCH IS NOT CONVERGED|TERMINATED' "gamess_logs/${{JOB}}.log" 2>/dev/null | tail -n 140 || true
echo
echo "===== recent log ====="
tail -n 80 "gamess_logs/${{JOB}}.log" 2>/dev/null || true
""",
    )
    monitor.chmod(0o755)

    readme = root / f"README_{args.out_job}.md"
    write_text(
        readme,
        "\n".join(
            [
                f"# {args.out_job}",
                "",
                f"- Starts from evaluated NSERCH: `{nserch}`",
                f"- Previous energy: `{metadata.get('energy_hartree', 'NA')}` Hartree",
                f"- Previous GRAD.MAX/RMS: `{metadata.get('grad_max', 'NA')}` / `{metadata.get('grad_rms', 'NA')}`",
                f"- Uses MOREAD vector: `{args.dat}`",
                f"- New NSTEP: `{args.nstep}`",
                f"- Harmonic distance guards: `{'; '.join(guard_descriptions) if guard_descriptions else 'none'}`",
                "",
                "Run/monitor:",
                "",
                "```bash",
                f"cd {root.as_posix()}",
                f"bash run_{args.out_job}_tmux.sh",
                f"bash monitor_{args.out_job}.sh",
                "```",
            ]
        )
        + "\n",
    )

    print(f"input={out_inp}")
    print(f"from_nserch={nserch}")
    print(f"previous_energy={metadata.get('energy_hartree', 'NA')}")
    print(f"previous_grad_max={metadata.get('grad_max', 'NA')}")
    print(f"previous_grad_rms={metadata.get('grad_rms', 'NA')}")
    print(f"vec_lines={len(vec)}")
    print(f"foreground={foreground}")
    print(f"tmux={tmux_runner}")
    print(f"monitor={monitor}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
