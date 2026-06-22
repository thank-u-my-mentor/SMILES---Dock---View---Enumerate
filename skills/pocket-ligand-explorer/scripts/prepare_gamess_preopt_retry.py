#!/usr/bin/env python3
"""Prepare a cheap GAMESS pre-optimization input from an MCPB small-opt input."""

from __future__ import annotations

import argparse
from pathlib import Path


def _card_start(line: str, card: str) -> bool:
    return line.strip().upper().startswith(card.upper())


def replace_card(lines: list[str], card: str, replacement: str) -> list[str]:
    """Replace a GAMESS $CARD, including simple multi-line cards ending in $END."""
    out: list[str] = []
    i = 0
    replaced = False
    while i < len(lines):
        line = lines[i]
        if _card_start(line, card):
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


def method_cards(method: str, charge: int, mult: int, runtype: str) -> tuple[str, str]:
    method = method.lower()
    runtype = runtype.upper()
    if method == "uhf":
        control = (
            " $CONTRL\n"
            f"  SCFTYP=UHF RUNTYP={runtype}\n"
            f"  ICHARG={charge} MULT={mult} MAXIT=200\n"
            " $END"
        )
    elif method == "rohf":
        control = (
            " $CONTRL\n"
            f"  SCFTYP=ROHF RUNTYP={runtype}\n"
            f"  ICHARG={charge} MULT={mult} MAXIT=200\n"
            " $END"
        )
    elif method == "svwn":
        control = (
            " $CONTRL\n"
            f"  SCFTYP=UHF DFTTYP=SVWN RUNTYP={runtype}\n"
            f"  ICHARG={charge} MULT={mult} MAXIT=200\n"
            " $END"
        )
    elif method == "b3lyp":
        control = (
            " $CONTRL\n"
            f"  SCFTYP=UHF DFTTYP=B3LYP RUNTYP={runtype}\n"
            f"  ICHARG={charge} MULT={mult} MAXIT=200\n"
            " $END"
        )
    else:
        raise ValueError(f"Unsupported method: {method}")
    return control, method


def basis_card(basis: str) -> str:
    basis = basis.lower()
    if basis in {"sto3g", "sto-3g"}:
        return " $BASIS GBASIS=STO NGAUSS=3 $END"
    if basis in {"321g", "3-21g", "n21"}:
        return " $BASIS GBASIS=N21 NGAUSS=3 $END"
    if basis in {"631g", "6-31g", "n31"}:
        return " $BASIS GBASIS=N31 NGAUSS=6 $END"
    if basis in {"631gd", "6-31g(d)", "n31d"}:
        return " $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END"
    raise ValueError(f"Unsupported basis: {basis}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a cheap GAMESS pre-optimization job and monitor scripts."
    )
    parser.add_argument("--mcpb-dir", type=Path, required=True)
    parser.add_argument("--source-job", required=True)
    parser.add_argument("--preopt-job", required=True)
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--mult", type=int, required=True)
    parser.add_argument("--method", choices=["uhf", "rohf", "svwn", "b3lyp"], default="uhf")
    parser.add_argument("--basis", default="sto3g")
    parser.add_argument("--runtype", choices=["optimize", "energy"], default="optimize")
    parser.add_argument("--scf-mode", choices=["diis", "soscf"], default="diis")
    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--session", default="zfy_mcpb_m7_preopt")
    parser.add_argument("--gamess", default="/home/qin/softwares/gamess/rungms")
    parser.add_argument("--nstep", type=int, default=160)
    parser.add_argument("--opttol", default="0.001")
    args = parser.parse_args()

    root = args.mcpb_dir
    source = root / f"{args.source_job}.inp"
    preopt = root / f"{args.preopt_job}.inp"
    if not source.exists():
        raise FileNotFoundError(source)

    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    control, _ = method_cards(args.method, args.charge, args.mult, args.runtype)
    lines = replace_card(lines, "$SYSTEM", " $SYSTEM MEMDDI=160 MWORDS=80 $END")
    lines = replace_card(lines, "$CONTRL", control)
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
            "  ETHRSH=10.0 MAXDII=30\n"
            " $END"
        )
    lines = replace_card(
        lines,
        "$SCF",
        scf,
    )
    lines = replace_card(lines, "$STATPT", f" $STATPT NSTEP={args.nstep} OPTTOL={args.opttol} $END")
    lines = replace_card(lines, "$BASIS", basis_card(args.basis))
    preopt.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")

    (root / "gamess_logs").mkdir(exist_ok=True)
    foreground = root / f"run_{args.preopt_job}_foreground.sh"
    foreground.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={args.preopt_job}
NCORES=${{NCORES:-{args.cores}}}
GAMESS=${{GAMESS:-{args.gamess}}}
VERSION=${{VERSION:-00}}
test -s "${{JOB}}.inp" || {{ echo "missing=${{JOB}}.inp"; exit 2; }}
mkdir -p gamess_logs logs
echo "[GAMESS] $JOB start $(date) NCORES=$NCORES method={args.method} basis={args.basis}" | tee "gamess_logs/${{JOB}}.status"
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

    tmux_runner = root / f"run_{args.preopt_job}_tmux.sh"
    tmux_runner.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
SESSION=${{SESSION:-{args.session}}}
NCORES=${{NCORES:-{args.cores}}}
if tmux has-session -t "$SESSION" 2>/dev/null; then
  echo "tmux_session_exists=$SESSION"
  echo "monitor=cd $PWD && bash monitor_{args.preopt_job}.sh"
  exit 0
fi
tmux new-session -d -s "$SESSION" "cd '$PWD'; env NCORES='$NCORES' bash {foreground.name}; echo; echo '[tmux] finished'; exec bash"
echo "tmux_session=$SESSION"
echo "monitor=cd $PWD && bash monitor_{args.preopt_job}.sh"
echo "tail=tail -f $PWD/gamess_logs/{args.preopt_job}.log"
""",
        encoding="utf-8",
        newline="\n",
    )
    tmux_runner.chmod(0o755)

    monitor = root / f"monitor_{args.preopt_job}.sh"
    monitor.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={args.preopt_job}
echo "===== status ====="
cat "gamess_logs/${{JOB}}.status" 2>/dev/null || true
echo
echo "===== running ====="
pgrep -af "${{JOB}}|ddikick|gamess.00.x" || true
echo
echo "===== recent log ====="
tail -n 80 "gamess_logs/${{JOB}}.log" 2>/dev/null || true
echo
echo "===== markers ====="
grep -E 'RUNTYP=|GBASIS|ITER EX|NSERCH:|TOTAL ENERGY|S-SQUARED|END OF GEOMETRY SEARCH|FAILURE|ERROR|TERMINATED' "gamess_logs/${{JOB}}.log" 2>/dev/null | tail -n 100 || true
""",
        encoding="utf-8",
        newline="\n",
    )
    monitor.chmod(0o755)

    print(f"preopt_input={preopt}")
    print(f"foreground={foreground}")
    print(f"tmux_runner={tmux_runner}")
    print(f"monitor={monitor}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
