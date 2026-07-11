#!/usr/bin/env python3
"""Write a GROMACS EM/NVT/NPT/production runner for MCPB/ACPYPE outputs."""

from __future__ import annotations

import argparse
import math
from pathlib import Path


EM_MDP = """integrator              = steep
emtol                   = 1000.0
emstep                  = 0.01
nsteps                  = 50000
nstenergy               = 100
nstlog                  = 100
constraints             = h-bonds
constraint_algorithm    = lincs
lincs_iter              = 2
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
DispCorr                = EnerPres
pbc                     = xyz
"""

NVT_MDP = """integrator              = md
dt                      = 0.002
nsteps                  = 250000
nstxout-compressed      = 5000
nstenergy               = 1000
nstlog                  = 1000
continuation            = no
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 2
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
DispCorr                = EnerPres
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = no
pbc                     = xyz
gen_vel                 = yes
gen_temp                = 300
gen_seed                = -1
"""

NPT_MDP = """integrator              = md
dt                      = 0.002
nsteps                  = 250000
nstxout-compressed      = 5000
nstenergy               = 1000
nstlog                  = 1000
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 2
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
DispCorr                = EnerPres
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
pbc                     = xyz
gen_vel                 = no
"""


def production_mdp(ns: float, xtc_ps: float) -> str:
    nsteps = int(round(ns * 500_000))
    nst = max(1, int(round(xtc_ps / 0.002)))
    return f"""integrator              = md
dt                      = 0.002
nsteps                  = {nsteps}
nstxout-compressed      = {nst}
nstenergy               = {nst}
nstlog                  = {nst}
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 2
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
DispCorr                = EnerPres
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
pbc                     = xyz
gen_vel                 = no
"""


def write_runner(
    outdir: Path,
    basename: str,
    production_name: str,
    gmx: str,
    gmxrc: str,
    extra_ld: str,
    ntomp: int,
    maxwarn: int,
) -> None:
    runner = f"""#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

export LD_LIBRARY_PATH="{extra_ld}:${{LD_LIBRARY_PATH:-}}"
if [ -n "{gmxrc}" ] && [ -f "{gmxrc}" ]; then
  # shellcheck disable=SC1090
  set +u
  source "{gmxrc}" >/dev/null 2>&1 || true
  set -u
fi

GMX="{gmx}"
TOP="{basename}_GMX.top"
START="{basename}_GMX.gro"

run_mdrun_gpu() {{
  local deffnm="$1"
  "$GMX" mdrun -deffnm "$deffnm" -ntomp {ntomp} -nb gpu -pme gpu -bonded gpu -update gpu -v
}}

if [ ! -s em.gro ]; then
  "$GMX" grompp -f em.mdp -c "$START" -p "$TOP" -o em.tpr -maxwarn {maxwarn} > logs/grompp_em.log 2>&1
  "$GMX" mdrun -deffnm em -ntomp {ntomp} -v > logs/mdrun_em.log 2>&1
fi

if [ ! -s nvt.gro ]; then
  "$GMX" grompp -f nvt.mdp -c em.gro -r em.gro -p "$TOP" -o nvt.tpr -maxwarn {maxwarn} > logs/grompp_nvt.log 2>&1
  run_mdrun_gpu nvt > logs/mdrun_nvt.log 2>&1
fi

if [ ! -s npt.gro ]; then
  "$GMX" grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p "$TOP" -o npt.tpr -maxwarn {maxwarn} > logs/grompp_npt.log 2>&1
  run_mdrun_gpu npt > logs/mdrun_npt.log 2>&1
fi

if [ ! -s {production_name}.gro ]; then
  "$GMX" grompp -f {production_name}.mdp -c npt.gro -t npt.cpt -p "$TOP" -o {production_name}.tpr -maxwarn {maxwarn} > logs/grompp_{production_name}.log 2>&1
  run_mdrun_gpu {production_name} > logs/mdrun_{production_name}.log 2>&1
fi

if [ -s {production_name}.xtc ] && [ -s {production_name}.tpr ]; then
  echo "[PLE] fixing trajectory PBC for visualization"
  printf "System\\n" | "$GMX" trjconv -s {production_name}.tpr -f {production_name}.xtc -o {production_name}_nojump.xtc -pbc nojump > logs/trjconv_{production_name}_nojump.log 2>&1
  printf "Protein\\nSystem\\n" | "$GMX" trjconv -s {production_name}.tpr -f {production_name}_nojump.xtc -o {production_name}_centered_protein_compact.xtc -pbc mol -center -ur compact > logs/trjconv_{production_name}_centered.log 2>&1
  printf "System\\n" | "$GMX" trjconv -s {production_name}.tpr -f {production_name}_centered_protein_compact.xtc -o {production_name}_centered_protein_compact.gro -dump 0 > logs/trjconv_{production_name}_centered_gro.log 2>&1
  echo "visualization_structure=$(pwd)/{production_name}_centered_protein_compact.gro"
  echo "visualization_trajectory=$(pwd)/{production_name}_centered_protein_compact.xtc"
fi
"""
    path = outdir / f"run_{production_name}.sh"
    path.write_text(runner)
    path.chmod(0o755)

    starter = f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
SESSION="${{SESSION:-tj_mcpb_{production_name}}}"
tmux has-session -t "$SESSION" 2>/dev/null && {{
  echo "tmux session already exists: $SESSION"
  exit 0
}}
tmux new-session -d -s "$SESSION" "bash run_{production_name}.sh"
echo "started tmux session: $SESSION"
echo "monitor: tail -f $(pwd)/logs/mdrun_{production_name}.log"
"""
    starter_path = outdir / f"start_{production_name}_tmux.sh"
    starter_path.write_text(starter)
    starter_path.chmod(0o755)

    monitor = f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
echo "== active GROMACS processes =="
pgrep -af "{production_name}|gmx_mpi|mdrun" || true
echo
for f in logs/grompp_em.log logs/mdrun_em.log logs/grompp_nvt.log logs/mdrun_nvt.log logs/grompp_npt.log logs/mdrun_npt.log logs/grompp_{production_name}.log logs/mdrun_{production_name}.log; do
  [ -f "$f" ] || continue
  echo "===== $f ====="
  tail -n 30 "$f"
done
"""
    monitor_path = outdir / f"monitor_{production_name}.sh"
    monitor_path.write_text(monitor)
    monitor_path.chmod(0o755)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmxdir", required=True, help="ACPYPE .amb2gmx directory.")
    parser.add_argument("--basename", required=True, help="ACPYPE basename before _GMX.gro/top.")
    parser.add_argument("--production-ns", type=float, default=200.0)
    parser.add_argument("--xtc-ps", type=float, default=5.0, help="Compressed trajectory interval in ps.")
    parser.add_argument("--gmx", default="/home/qin/softwares/gromacs-2026.0-gpu/bin/gmx_mpi")
    parser.add_argument("--gmxrc", default="/home/qin/softwares/gromacs-2026.0-gpu/bin/GMXRC")
    parser.add_argument(
        "--extra-ld-library-path",
        default="/home/qin/softwares/cp2k-2026.1/build_minlib/src:/home/qin/softwares/gromacs-2026.0-gpu/lib",
    )
    parser.add_argument("--ntomp", type=int, default=16)
    parser.add_argument("--maxwarn", type=int, default=2)
    args = parser.parse_args()

    outdir = Path(args.gmxdir)
    gro = outdir / f"{args.basename}_GMX.gro"
    top = outdir / f"{args.basename}_GMX.top"
    if not gro.exists() or not top.exists():
        raise SystemExit(f"Missing ACPYPE GROMACS files: {gro} / {top}")

    ns_label = str(int(args.production_ns)) if math.isclose(args.production_ns, round(args.production_ns)) else str(args.production_ns).replace(".", "p")
    production_name = f"production_{ns_label}ns"

    (outdir / "logs").mkdir(exist_ok=True)
    (outdir / "em.mdp").write_text(EM_MDP)
    (outdir / "nvt.mdp").write_text(NVT_MDP)
    (outdir / "npt.mdp").write_text(NPT_MDP)
    (outdir / f"{production_name}.mdp").write_text(production_mdp(args.production_ns, args.xtc_ps))
    write_runner(
        outdir=outdir,
        basename=args.basename,
        production_name=production_name,
        gmx=args.gmx,
        gmxrc=args.gmxrc,
        extra_ld=args.extra_ld_library_path,
        ntomp=args.ntomp,
        maxwarn=args.maxwarn,
    )
    print(f"wrote={outdir / 'em.mdp'}")
    print(f"wrote={outdir / 'nvt.mdp'}")
    print(f"wrote={outdir / 'npt.mdp'}")
    print(f"wrote={outdir / f'{production_name}.mdp'}")
    print(f"wrote={outdir / f'run_{production_name}.sh'}")
    print(f"wrote={outdir / f'start_{production_name}_tmux.sh'}")
    print(f"wrote={outdir / f'monitor_{production_name}.sh'}")


if __name__ == "__main__":
    main()
