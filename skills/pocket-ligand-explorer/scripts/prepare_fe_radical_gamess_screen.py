from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path


ATOMIC_NUMBERS = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "S": 16,
    "Fe": 26,
}


@dataclass(frozen=True)
class Atom:
    serial: int
    name: str
    resn: str
    chain: str
    resi: str
    elem: str
    x: float
    y: float
    z: float


@dataclass
class ModelAtom:
    name: str
    resn: str
    chain: str
    resi: str
    elem: str
    x: float
    y: float
    z: float
    note: str = ""


Vec = tuple[float, float, float]


def v_add(a: Vec, b: Vec) -> Vec:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def v_sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def v_scale(a: Vec, s: float) -> Vec:
    return (a[0] * s, a[1] * s, a[2] * s)


def v_dot(a: Vec, b: Vec) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def v_cross(a: Vec, b: Vec) -> Vec:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def v_norm(a: Vec) -> float:
    return math.sqrt(v_dot(a, a))


def v_unit(a: Vec, fallback: Vec = (1.0, 0.0, 0.0)) -> Vec:
    n = v_norm(a)
    if n < 1e-8:
        return fallback
    return (a[0] / n, a[1] / n, a[2] / n)


def coord(atom: Atom | ModelAtom) -> Vec:
    return (atom.x, atom.y, atom.z)


def distance(a: Atom | ModelAtom, b: Atom | ModelAtom) -> float:
    return v_norm(v_sub(coord(a), coord(b)))


def perpendicular_basis(axis: Vec) -> tuple[Vec, Vec]:
    axis = v_unit(axis)
    trial = (1.0, 0.0, 0.0)
    if abs(v_dot(axis, trial)) > 0.85:
        trial = (0.0, 1.0, 0.0)
    p = v_unit(v_cross(axis, trial))
    q = v_unit(v_cross(axis, p))
    return p, q


def placed_atom(base: Atom | ModelAtom, name: str, direction: Vec, length: float, *, resn: str, chain: str, resi: str, note: str) -> ModelAtom:
    direction = v_unit(direction)
    pos = v_add(coord(base), v_scale(direction, length))
    return ModelAtom(name=name, resn=resn, chain=chain, resi=resi, elem="H", x=pos[0], y=pos[1], z=pos[2], note=note)


def methyl_hydrogens(cap: Atom | ModelAtom, anchor: Atom | ModelAtom, *, prefix: str, resn: str, chain: str, resi: str) -> list[ModelAtom]:
    to_anchor = v_unit(v_sub(coord(anchor), coord(cap)))
    p, q = perpendicular_basis(to_anchor)
    hydrogens: list[ModelAtom] = []
    for i, angle in enumerate((0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0), start=1):
        radial = v_add(v_scale(p, math.cos(angle)), v_scale(q, math.sin(angle)))
        direction = v_add(v_scale(to_anchor, -1.0 / 3.0), v_scale(radial, math.sqrt(8.0 / 9.0)))
        hydrogens.append(placed_atom(cap, f"{prefix}{i}", direction, 1.09, resn=resn, chain=chain, resi=resi, note="added methyl cap H"))
    return hydrogens


def single_hydrogen(atom: Atom | ModelAtom, neighbors: list[Atom | ModelAtom], *, name: str, resn: str, chain: str, resi: str, length: float = 1.01) -> ModelAtom:
    direction = (0.0, 0.0, 0.0)
    for neighbor in neighbors:
        direction = v_add(direction, v_sub(coord(atom), coord(neighbor)))
    return placed_atom(atom, name, direction, length, resn=resn, chain=chain, resi=resi, note="added valence H")


def water_hydrogens(oxygen: Atom, metal: Atom) -> list[ModelAtom]:
    away_from_fe = v_unit(v_sub(coord(oxygen), coord(metal)))
    p, _ = perpendicular_basis(away_from_fe)
    half_angle = math.radians(104.5 / 2.0)
    hydrogens: list[ModelAtom] = []
    for name, sign in (("H1", 1.0), ("H2", -1.0)):
        direction = v_add(v_scale(away_from_fe, math.cos(half_angle)), v_scale(p, sign * math.sin(half_angle)))
        hydrogens.append(
            placed_atom(oxygen, name, direction, 0.96, resn="HOH", chain=oxygen.chain, resi=oxygen.resi, note="added water H")
        )
    return hydrogens


def parse_pdb(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        raw_elem = line[76:80].strip() if len(line) >= 78 else ""
        elem = "".join(ch for ch in raw_elem if ch.isalpha())
        if not elem:
            elem = "".join(ch for ch in line[12:16].strip() if ch.isalpha())[:1]
        elem = "Fe" if elem.upper() == "FE" else elem.capitalize()
        atoms.append(
            Atom(
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                resn=line[17:20].strip(),
                chain=line[21].strip(),
                resi=line[22:26].strip(),
                elem=elem,
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
            )
        )
    return atoms


def find_atom(atoms: list[Atom], *, resn: str, chain: str, resi: str, name: str | None = None, elem: str | None = None) -> Atom:
    matches = [
        atom
        for atom in atoms
        if atom.resn == resn
        and atom.chain == chain
        and atom.resi == resi
        and (name is None or atom.name == name)
        and (elem is None or atom.elem.lower() == elem.lower())
    ]
    if len(matches) != 1:
        label = f"{resn}:{chain}:{resi}:{name or '*'}:{elem or '*'}"
        raise SystemExit(f"Expected one atom for {label}, found {len(matches)}")
    return matches[0]


def residue_atoms(atoms: list[Atom], *, resn: str, chain: str, resi: str) -> list[Atom]:
    found = [atom for atom in atoms if atom.resn == resn and atom.chain == chain and atom.resi == resi]
    if not found:
        raise SystemExit(f"Missing residue {resn}:{chain}:{resi}")
    return found


def copy_atom(atom: Atom, *, resn: str | None = None, chain: str | None = None, resi: str | None = None, note: str = "") -> ModelAtom:
    return ModelAtom(
        name=atom.name,
        resn=resn or atom.resn,
        chain=chain if chain is not None else atom.chain,
        resi=resi if resi is not None else atom.resi,
        elem=atom.elem,
        x=atom.x,
        y=atom.y,
        z=atom.z,
        note=note,
    )


def build_his_fragment(atoms: list[Atom], *, chain: str, resi: str, tag: str) -> list[ModelAtom]:
    source = {atom.name: atom for atom in residue_atoms(atoms, resn="HIS", chain=chain, resi=resi)}
    needed = ["CB", "CG", "ND1", "CD2", "CE1", "NE2"]
    missing = [name for name in needed if name not in source]
    if missing:
        raise SystemExit(f"HIS {chain}:{resi} missing atoms: {', '.join(missing)}")
    resn = f"H{tag}"[:3]
    model = [copy_atom(source[name], resn=resn, resi=resi, note="His->methyl-imidazole") for name in needed]
    by_name = {atom.name: atom for atom in model}
    model += methyl_hydrogens(by_name["CB"], by_name["CG"], prefix="HB", resn=resn, chain=chain, resi=resi)
    model.append(single_hydrogen(by_name["ND1"], [by_name["CG"], by_name["CE1"]], name="HD1", resn=resn, chain=chain, resi=resi))
    model.append(single_hydrogen(by_name["CD2"], [by_name["CG"], by_name["NE2"]], name="HD2", resn=resn, chain=chain, resi=resi))
    model.append(single_hydrogen(by_name["CE1"], [by_name["ND1"], by_name["NE2"]], name="HE1", resn=resn, chain=chain, resi=resi))
    return model


def build_glu_acetate(atoms: list[Atom], *, chain: str, resi: str) -> list[ModelAtom]:
    source = {atom.name: atom for atom in residue_atoms(atoms, resn="GLU", chain=chain, resi=resi)}
    needed = ["CG", "CD", "OE1", "OE2"]
    missing = [name for name in needed if name not in source]
    if missing:
        raise SystemExit(f"GLU {chain}:{resi} missing atoms: {', '.join(missing)}")
    resn = "EAC"
    model = [copy_atom(source[name], resn=resn, resi=resi, note="Glu->acetate") for name in needed]
    by_name = {atom.name: atom for atom in model}
    model += methyl_hydrogens(by_name["CG"], by_name["CD"], prefix="HG", resn=resn, chain=chain, resi=resi)
    return model


def build_model(args: argparse.Namespace) -> list[ModelAtom]:
    atoms = parse_pdb(Path(args.pdb))
    metal = find_atom(atoms, resn=args.metal_resn, chain=args.metal_chain, resi=args.metal_resi, name=args.metal_atom, elem="Fe")
    model: list[ModelAtom] = [copy_atom(metal, resn="FE", note="Fe(III)")]
    model += build_his_fragment(atoms, chain=args.chain, resi="187", tag="1")
    model += build_his_fragment(atoms, chain=args.chain, resi="270", tag="2")
    model += build_glu_acetate(atoms, chain=args.chain, resi="349")
    model += [copy_atom(atom, note="ACT acetate") for atom in residue_atoms(atoms, resn="ACT", chain=args.chain, resi="501")]
    model += [copy_atom(atom, note="full UNL radical ligand") for atom in residue_atoms(atoms, resn="UNL", chain=args.chain, resi="1")]
    water_o = find_atom(atoms, resn="HOH", chain="B", resi="875", name="O", elem="O")
    model.append(copy_atom(water_o, note="coordinating water O"))
    model += water_hydrogens(water_o, metal)
    return model


def write_model_pdb(path: Path, atoms: list[ModelAtom]) -> None:
    lines: list[str] = []
    for serial, atom in enumerate(atoms, start=1):
        lines.append(
            f"HETATM{serial:5d} {atom.name[:4]:>4s} {atom.resn[:3]:>3s} {atom.chain[:1] or 'A'}"
            f"{int(atom.resi):4d}    {atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
            f"  1.00  0.00          {atom.elem:>2s}"
        )
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")


def write_model_xyz(path: Path, atoms: list[ModelAtom]) -> None:
    lines = [str(len(atoms)), "Fe(III)-radical truncated screening model"]
    for atom in atoms:
        lines.append(f"{atom.elem:<2s} {atom.x:14.6f} {atom.y:14.6f} {atom.z:14.6f}  # {atom.resn}:{atom.resi}:{atom.name} {atom.note}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_gamess_input(path: Path, atoms: list[ModelAtom], *, charge: int, mult: int, nstep: int) -> None:
    scftyp = "RHF" if mult == 1 else "UHF"
    lines = [
        " $SYSTEM MWORDS=50 MEMDDI=10 $END",
        " $CONTRL",
        f"  SCFTYP={scftyp} DFTTYP=B3LYP RUNTYP=OPTIMIZE",
        f"  ICHARG={charge} MULT={mult} MAXIT=200 COORD=CART UNITS=ANGS",
        " $END",
        f" $STATPT NSTEP={nstep} OPTTOL=0.001 $END",
        " $SCF DIRSCF=.T. DIIS=.T. DAMP=.T. ETHRSH=2.0 MAXDII=20 $END",
        " $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=1 $END",
        " $DATA",
        f"FeIII radical six-coordinate screen MULT={mult}",
        "C1",
    ]
    for atom in atoms:
        atomic_number = ATOMIC_NUMBERS.get(atom.elem)
        if atomic_number is None:
            raise SystemExit(f"Unsupported element {atom.elem!r} in {atom}")
        lines.append(f"{atom.elem:<2s} {float(atomic_number):8.1f} {atom.x:14.6f} {atom.y:14.6f} {atom.z:14.6f}")
    lines += [" $END", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def write_run_script(outdir: Path, job_names: list[str]) -> None:
    job_script = outdir / "run_gamess_job.sh"
    job_script.write_text(
        "#!/usr/bin/env bash\n"
        "set -uo pipefail\n"
        "cd \"$(dirname \"$0\")\"\n"
        "mkdir -p logs\n"
        "job=\"${1:?usage: run_gamess_job.sh JOB [NCORES]}\"\n"
        "GAMESS=${GAMESS:-/home/qin/softwares/gamess/rungms}\n"
        "VERSION=${VERSION:-00}\n"
        "NCORES=${2:-${NCORES:-4}}\n"
        "log=\"logs/${job}.log\"\n"
        "status=\"logs/${job}.status\"\n"
        "echo \"[GAMESS] ${job} start $(date) NCORES=${NCORES}\" | tee \"$status\"\n"
        "\"$GAMESS\" \"$job\" \"$VERSION\" \"$NCORES\" > \"$log\" 2>&1\n"
        "rc=$?\n"
        "echo \"[GAMESS] ${job} exit_code=${rc} end $(date)\" | tee -a \"$status\"\n"
        "if grep -q 'EXECUTION OF GAMESS TERMINATED NORMALLY' \"$log\"; then\n"
        "  echo \"[GAMESS] ${job} status=normal\" | tee -a \"$status\"\n"
        "elif grep -q 'EXECUTION OF GAMESS TERMINATED -ABNORMALLY-' \"$log\"; then\n"
        "  echo \"[GAMESS] ${job} status=abnormal\" | tee -a \"$status\"\n"
        "else\n"
        "  echo \"[GAMESS] ${job} status=no_termination_marker\" | tee -a \"$status\"\n"
        "fi\n"
        "exit \"$rc\"\n",
        encoding="utf-8",
    )
    job_script.chmod(0o755)

    script = outdir / "run_fast_spin_screen.sh"
    script.write_text(
        "#!/usr/bin/env bash\n"
        "set -uo pipefail\n"
        "cd \"$(dirname \"$0\")\"\n"
        "mkdir -p logs\n"
        "NCORES=${NCORES:-8}\n"
        "echo \"[spin-screen] jobs: " + " ".join(job_names) + "\"\n"
        "echo \"[spin-screen] NCORES=$NCORES\"\n"
        "overall=0\n"
        + "".join(f"bash run_gamess_job.sh {job} \"$NCORES\" || overall=$?\n" for job in job_names)
        + "echo \"[spin-screen] overall_exit=${overall}\"\n"
        + "echo '[spin-screen] done'\n",
        encoding="utf-8",
    )
    script.chmod(0o755)

    tmux_script = outdir / "start_fast_spin_screen_tmux.sh"
    tmux_script.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "cd \"$(dirname \"$0\")\"\n"
        "SESSION=${SESSION:-feiii_spin_screen}\n"
        "NCORES=${NCORES:-4}\n"
        "tmux has-session -t \"$SESSION\" 2>/dev/null && {\n"
        "  echo \"tmux_session_exists=$SESSION\"\n"
        "  echo \"attach=tmux attach -t $SESSION\"\n"
        "  exit 0\n"
        "}\n"
        "tmux new-session -d -s \"$SESSION\" \"cd '$PWD'; env NCORES='$NCORES' bash run_fast_spin_screen.sh; echo; echo '[tmux] finished, press Ctrl-b d to detach or exit to close'; exec bash\"\n"
        "echo \"tmux_session=$SESSION\"\n"
        "echo \"attach=tmux attach -t $SESSION\"\n"
        "echo \"monitor=bash monitor_fast_spin_screen.sh\"\n",
        encoding="utf-8",
    )
    tmux_script.chmod(0o755)
    monitor = outdir / "monitor_fast_spin_screen.sh"
    monitor.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "cd \"$(dirname \"$0\")\"\n"
        "for log in logs/*.log; do\n"
        "  echo \"===== $log =====\"\n"
        "  grep -E 'INPUT CARD>|RUNTYP=|NSERCH=|STATIONARY POINT|END OF GEOMETRY SEARCH|TOTAL ENERGY|TERMINATED|FAILURE|ERROR|S\\*\\*2|S-SQUARED' \"$log\" | tail -n 40 || true\n"
        "done\n",
        encoding="utf-8",
    )
    monitor.chmod(0o755)


def write_workflow_markdown(outdir: Path, atoms: list[ModelAtom], *, charge: int) -> None:
    counts: dict[str, int] = {}
    for atom in atoms:
        counts[atom.elem] = counts.get(atom.elem, 0) + 1
    count_text = ", ".join(f"{elem}={counts[elem]}" for elem in sorted(counts))
    (outdir / "WORKFLOW_FE_RADICAL_MCPB.md").write_text(
        f"""# TJ Fe(III)-Radical MCPB Working Model

This document records the current idea before the expensive MCPB Hessian step.

## Current Chemical Model

- Metal: Fe(III), chain A active-site Fe.
- Coordination model: six-coordinate non-heme Fe.
- Donors: His187 NE2, His270 NE2, Glu349 OE1, ACT501 O2, UNL N1, HOH875 O.
- Radical: UNL C4/benzyl-side radical candidate, represented by deleting one C4 hydrogen.
- Net charge used for the screening QM model: `{charge}`.
- Spin states to screen first: `MULT=5` and `MULT=7`.
- Atom count in truncated QM model: `{len(atoms)}` ({count_text}).

Important parity check: `MULT=5` and `MULT=7` require an even total electron
count. If GAMESS prints `CHECK YOUR INPUT CHARGE AND MULTIPLICITY`, revise the
model charge, ligand protonation, or selected spin states before running a long
job.

## Workflow Diagram

```mermaid
flowchart TD
  A[Curated active-site PDB\\nH++ protonation + Fe + ACT + UNL] --> B[Radical preflight\\nremove one C4 H from UNL]
  B --> C[AmberTools ligand-only check\\nantechamber + parmchk2 + tleap]
  C --> D{{Ligand format OK?}}
  D -- no --> C1[Fix names, charges, atom graph, radical atom]
  C1 --> C
  D -- yes --> E[Build truncated Fe QM model\\nHis->methyl-imidazole\\nGlu->acetate\\nUNL kept full]
  E --> F[Fast GAMESS B3LYP/6-31G(d,p) screen\\nMULT=5 and MULT=7]
  F --> G{{Stable coordination?\\nreasonable spin density?}}
  G -- no --> E1[Revise oxidation state, ligand set, protonation, or spin]
  E1 --> E
  G -- yes --> H[MCPB.py production parameterization\\nsmall-model Hessian + large-model charges]
  H --> I[Amber bonded Fe force field\\nfrcmod/mol2/lib/prep]
  I --> J[Amber tleap -> ACPYPE -> GROMACS]
  J --> K[Production MD + PBC fix + RDC dashboard]
```

## What MCPB.py Is Doing Here

Classical MD knows only atoms, charges, bonds, angles, torsions, and nonbonded
parameters. A catalytic Fe center is not an ordinary ion in water: His/Glu/ACT/UNL/HOH
create a ligand field, and the Fe-donor distances should be maintained by a
bonded metal-site force field rather than by a naked Fe nonbonded model.

MCPB.py bridges QM and classical MD:

1. Define the metal center and its coordinating atoms.
2. Build a small QM model of Fe plus truncated donor ligands.
3. Optimize that model and compute force constants with QM.
4. Use the Seminario-style analysis to turn the QM Hessian into classical
   Fe-donor bond and angle parameters.
5. Fit or transfer charges for the larger metal site.
6. Write Amber-readable parameter files that can later be converted to GROMACS.

The expensive part is not "running MD"; it is asking QM to describe the metal
coordination surface well enough that the later Newtonian MD has meaningful
bond/angle parameters around Fe.

## How To Monitor

```bash
cd {outdir.as_posix()}
bash monitor_fast_spin_screen.sh
tail -f logs/FEIII_radical_M5_fastopt.log
tail -f logs/FEIII_radical_M7_fastopt.log
```

## CPU Use And Faster Runs

This screen is launched through `run_fast_spin_screen.sh` or, more robustly, a
tmux session:

```bash
cd {outdir.as_posix()}
env NCORES=16 SESSION=feiii_spin_screen bash start_fast_spin_screen_tmux.sh
tmux attach -t feiii_spin_screen
```

On a 32-thread CPU, `NCORES=4` naturally looks like about 12-15% total CPU even
when GAMESS is working. Use `NCORES=12-16` when the input has passed the early
charge/multiplicity and SCF-start checks. Do not expect perfect scaling; DFT has
serial and communication-heavy phases. If a higher-core run becomes unstable,
fall back to `NCORES=8`.

Good signs:

- `EXECUTION OF GAMESS TERMINATED NORMALLY`
- `END OF GEOMETRY SEARCH`
- Fe-donor bonds do not break or drift into absurd distances

Bad signs:

- SCF non-convergence
- optimized structure loses an Fe-His/Glu/ACT/UNL/HOH contact
- spin density ends up on the wrong chemical fragment
- one multiplicity is much less stable or chemically distorted
""",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Build Fe(III)-radical MULT=5/7 GAMESS screening inputs.")
    parser.add_argument("--pdb", default="/mnt/e/TJ/feiii_radical_preflight/complex_AFeIII_UNL_C4radical_dropH4A.pdb")
    parser.add_argument("--outdir", default="/mnt/e/TJ/feiii_radical_preflight/gamess_spin_screen")
    parser.add_argument("--chain", default="A")
    parser.add_argument("--metal-resn", default="FE2")
    parser.add_argument("--metal-chain", default="A")
    parser.add_argument("--metal-resi", default="431")
    parser.add_argument("--metal-atom", default="FE")
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--nstep", type=int, default=80)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    model = build_model(args)
    write_model_pdb(outdir / "feiii_radical_truncated_model.pdb", model)
    write_model_xyz(outdir / "feiii_radical_truncated_model.xyz", model)

    job_names: list[str] = []
    for mult in (5, 7):
        job = f"FEIII_radical_M{mult}_fastopt"
        write_gamess_input(outdir / f"{job}.inp", model, charge=args.charge, mult=mult, nstep=args.nstep)
        job_names.append(job)
    write_run_script(outdir, job_names)
    write_workflow_markdown(outdir, model, charge=args.charge)

    print(f"outdir={outdir}")
    print(f"model_pdb={outdir / 'feiii_radical_truncated_model.pdb'}")
    print(f"model_xyz={outdir / 'feiii_radical_truncated_model.xyz'}")
    for job in job_names:
        print(f"gamess_input={outdir / (job + '.inp')}")
    print(f"run_script={outdir / 'run_fast_spin_screen.sh'}")
    print(f"workflow_md={outdir / 'WORKFLOW_FE_RADICAL_MCPB.md'}")
    print(f"atoms={len(model)} charge={args.charge} mults=5,7")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
