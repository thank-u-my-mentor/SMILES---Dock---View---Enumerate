#!/usr/bin/env python3
"""Build a MOREAD GAMESS geometry optimization with Fe-donor harmonic restraints."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Atom:
    index: int
    name: str
    resname: str
    chain: str
    resseq: str
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class Constraint:
    kind: str
    atom_a: Atom
    atom_b: Atom
    target: float
    force: float


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


def parse_pdb_atoms(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            atom = Atom(
                index=len(atoms) + 1,
                name=line[12:16].strip(),
                resname=line[17:20].strip(),
                chain=line[21:22].strip(),
                resseq=line[22:26].strip(),
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
            )
        except ValueError as exc:
            raise SystemExit(f"Could not parse PDB atom line in {path}: {line}") from exc
        atoms.append(atom)
    if not atoms:
        raise SystemExit(f"No ATOM/HETATM records found in {path}")
    return atoms


def parse_selector(text: str) -> tuple[str, str, str, str]:
    parts = text.split(":")
    if len(parts) != 4:
        raise argparse.ArgumentTypeError(
            f"selector must be RESNAME:CHAIN:RESSEQ:ATOM, got {text!r}"
        )
    return tuple(part.strip() for part in parts)  # type: ignore[return-value]


def parse_pair(text: str) -> tuple[tuple[str, str, str, str], tuple[str, str, str, str]]:
    if "," not in text:
        raise argparse.ArgumentTypeError(
            "pair must be SEL1,SEL2 where each selector is RESNAME:CHAIN:RESSEQ:ATOM"
        )
    left, right = text.split(",", 1)
    return parse_selector(left.strip()), parse_selector(right.strip())


def select_atom(atoms: list[Atom], selector: tuple[str, str, str, str]) -> Atom:
    resname, chain, resseq, atom_name = selector
    matches = [
        atom
        for atom in atoms
        if (resname in {"", "*"} or atom.resname == resname)
        and (chain in {"", "*"} or atom.chain == chain)
        and (resseq in {"", "*"} or atom.resseq == resseq)
        and (atom_name in {"", "*"} or atom.name == atom_name)
    ]
    if len(matches) != 1:
        pretty = ":".join(selector)
        detail = ", ".join(
            f"{a.index}:{a.resname}:{a.chain}:{a.resseq}:{a.name}" for a in matches[:20]
        )
        raise SystemExit(
            f"Selector {pretty!r} matched {len(matches)} atom(s), expected 1. {detail}"
        )
    return matches[0]


def distance(a: Atom, b: Atom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def statpt_card(
    *,
    nstep: int,
    opttol: str,
    pairs: list[tuple[int, int]],
    targets: list[float],
    forces: list[float],
) -> str:
    ih_lines: list[str] = []
    ih_values: list[int] = []
    for atom_a, atom_b in pairs:
        ih_values.extend([1, atom_a, atom_b])
    for start in range(0, len(ih_values), 9):
        array_index = start + 1
        chunk = ",".join(str(v) for v in ih_values[start : start + 9])
        ih_lines.append(f"         IHMCON({array_index})={chunk}")
    sh_lines: list[str] = []
    fh_lines: list[str] = []
    for start in range(0, len(targets), 6):
        array_index = start + 1
        sh_values = ",".join(f"{target:.4f}" for target in targets[start : start + 6])
        fh_values = ",".join(f"{force:.1f}" for force in forces[start : start + 6])
        sh_lines.append(f"         SHMCON({array_index})={sh_values}")
        fh_lines.append(f"         FHMCON({array_index})={fh_values}")
    return "\n".join(
        [
            f" $STATPT NSTEP={nstep} OPTTOL={opttol}",
            *ih_lines,
            *sh_lines,
            *fh_lines,
            " $END",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mcpb-dir", type=Path, required=True)
    parser.add_argument("--source-job", required=True)
    parser.add_argument("--out-job", required=True)
    parser.add_argument("--dat", type=Path, required=True)
    parser.add_argument("--small-pdb", type=Path, required=True)
    parser.add_argument("--metal", type=parse_selector, required=True)
    parser.add_argument("--donor", type=parse_selector, action="append", required=True)
    parser.add_argument(
        "--protect-bond",
        type=parse_pair,
        action="append",
        default=[],
        help=(
            "Additional harmonic distance constraint as "
            "RES:CHAIN:RESSEQ:ATOM,RES:CHAIN:RESSEQ:ATOM. "
            "Use this to protect nonreactive ligand covalent bonds during short OPT."
        ),
    )
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--mult", type=int, required=True)
    parser.add_argument("--norb", type=int, required=True)
    parser.add_argument("--nstep", type=int, default=60)
    parser.add_argument("--opttol", default="0.0002")
    parser.add_argument("--force", type=float, default=500.0)
    parser.add_argument("--metal-force", type=float, default=None)
    parser.add_argument("--protect-force", type=float, default=None)
    parser.add_argument("--target-distance", type=float, default=None)
    parser.add_argument("--scf-mode", choices=["diis", "soscf"], default="soscf")
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--session", default=None)
    args = parser.parse_args()

    root = args.mcpb_dir
    source = root / f"{args.source_job}.inp"
    if not source.exists():
        raise FileNotFoundError(source)
    vec = extract_vec(args.dat)
    lines = remove_vec(source.read_text(encoding="utf-8", errors="replace").splitlines())

    atoms = parse_pdb_atoms(args.small_pdb)
    metal = select_atom(atoms, args.metal)
    donors = [select_atom(atoms, selector) for selector in args.donor]
    metal_force = args.force if args.metal_force is None else args.metal_force
    protect_force = args.force if args.protect_force is None else args.protect_force
    constraints: list[Constraint] = []
    for donor in donors:
        target = args.target_distance if args.target_distance is not None else distance(metal, donor)
        constraints.append(Constraint("metal", metal, donor, target, metal_force))
    for selector_a, selector_b in args.protect_bond:
        atom_a = select_atom(atoms, selector_a)
        atom_b = select_atom(atoms, selector_b)
        constraints.append(
            Constraint("protected_bond", atom_a, atom_b, distance(atom_a, atom_b), protect_force)
        )

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
    lines = replace_card(
        lines,
        "$STATPT",
        statpt_card(
            nstep=args.nstep,
            opttol=args.opttol,
            pairs=[(item.atom_a.index, item.atom_b.index) for item in constraints],
            targets=[item.target for item in constraints],
            forces=[item.force for item in constraints],
        ),
    )
    lines = replace_card(
        lines,
        "$GUESS",
        f" $GUESS GUESS=MOREAD NORB={args.norb} $END",
    )

    out_inp = root / f"{args.out_job}.inp"
    out_inp.write_text("\n".join(lines + vec) + "\n", encoding="utf-8", newline="\n")

    report = root / f"{args.out_job}_harmonic_constraints.tsv"
    report.write_text(
        "constraint\tkind\tatom_a_index\tatom_a\tatom_b_index\tatom_b\ttarget_A\tforce_kcal_mol_A2\tcurrent_A\n"
        + "\n".join(
            (
                f"{i}\t{kind}"
                f"\t{atom_a.index}\t{atom_a.resname}:{atom_a.chain}:{atom_a.resseq}:{atom_a.name}"
                f"\t{atom_b.index}\t{atom_b.resname}:{atom_b.chain}:{atom_b.resseq}:{atom_b.name}"
                f"\t{target:.4f}\t{force:.1f}\t{distance(atom_a, atom_b):.4f}"
            )
            for i, (kind, atom_a, atom_b, target, force) in enumerate(
                (
                    (
                        item.kind,
                        item.atom_a,
                        item.atom_b,
                        item.target,
                        item.force,
                    )
                    for item in constraints
                ),
                start=1,
            )
        )
        + "\n",
        encoding="utf-8",
        newline="\n",
    )

    (root / "gamess_logs").mkdir(exist_ok=True)
    foreground = root / f"run_{args.out_job}_foreground.sh"
    foreground.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={args.out_job}
NCORES=${{NCORES:-{args.cores}}}
GAMESS=${{GAMESS:-/home/qin/softwares/gamess/rungms}}
VERSION=${{VERSION:-00}}
mkdir -p gamess_logs logs
echo "[GAMESS] $JOB start $(date) NCORES=$NCORES MOREAD={args.dat} metal_force={metal_force} protect_force={protect_force} dft_grid=fine_from_start" | tee "gamess_logs/${{JOB}}.status"
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
        encoding="utf-8",
        newline="\n",
    )
    foreground.chmod(0o755)

    background = root / f"run_{args.out_job}_background.sh"
    background.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
NCORES=${{NCORES:-{args.cores}}}
nohup env NCORES="$NCORES" bash {foreground.name} > logs/{args.out_job}.nohup.log 2>&1 &
pid=$!
echo "$pid" > logs/{args.out_job}.pid
echo "pid=$pid"
echo "tail=tail -f $PWD/gamess_logs/{args.out_job}.log"
echo "monitor=bash $PWD/monitor_{args.out_job}.sh"
""",
        encoding="utf-8",
        newline="\n",
    )
    background.chmod(0o755)

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
grep -E 'GUESS =|RUNTYP=|NRAD0|NLEB0|SWOFF|DFT CODE IS SWITCHING|HARMONIC|BEGINNING GEOMETRY|NSERCH|GRAD\\. MAX|CONVERGED|SCF IS UNCONVERGED|SCF HAS NOT CONVERGED|FINAL U-B3LYP|S-SQUARED|END OF GEOMETRY SEARCH|TERMINATED' "gamess_logs/${{JOB}}.log" 2>/dev/null | tail -n 140 || true
echo
echo "===== recent log ====="
tail -n 80 "gamess_logs/${{JOB}}.log" 2>/dev/null || true
""",
        encoding="utf-8",
        newline="\n",
    )
    monitor.chmod(0o755)

    print(f"input={out_inp}")
    print(f"constraints={report}")
    print(f"vec_lines={len(vec)}")
    print(f"foreground={foreground}")
    print(f"background={background}")
    print(f"tmux={tmux_runner}")
    print(f"monitor={monitor}")
    for item in constraints:
        print(
            "constraint "
            f"{item.kind} "
            f"{item.atom_a.resname}:{item.atom_a.chain}:{item.atom_a.resseq}:{item.atom_a.name}#{item.atom_a.index}-"
            f"{item.atom_b.resname}:{item.atom_b.chain}:{item.atom_b.resseq}:{item.atom_b.name}#{item.atom_b.index} "
            f"target={item.target:.4f} force={item.force:.1f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
