#!/usr/bin/env python3
"""Guarded GAMESS OPT loop for Fe-radical MCPB small models.

The loop is deliberately conservative:

1. Wait for a GAMESS OPT job to finish.
2. Audit every evaluated NSERCH geometry for ligand bond identity and required
   Fe-donor distances.
3. Pick the chemically safe geometry with the best gradient.
4. If it is good enough for Hessian, stop and write HESSIAN_CANDIDATE.txt.
5. Otherwise continue OPT from that safe NSERCH with MOREAD and adjusted
   harmonic distance guards.

This is not a black-box optimizer.  It writes a round-by-round Markdown report
so the chemical assumptions remain visible.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from export_gamess_optimized_model import (  # noqa: E402
    apply_transform,
    best_rigid_transform,
    parse_coords_for_nserch,
)
from export_mcpb_opt_identity_pdb import (  # noqa: E402
    audit_mol2_bonds,
    coords_by_atom_name,
    distance,
    parse_pdb,
    to_model_atoms,
    write_displacements,
    write_fe_distances,
    write_pdb,
)


RESTART_DIR = Path("/home/qin/softwares/gamess/restart")


@dataclass
class HarmonicGuard:
    atom_i: int
    atom_j: int
    target: float
    force: float

    def key(self) -> tuple[int, int]:
        return tuple(sorted((self.atom_i, self.atom_j)))

    def as_args(self) -> list[str]:
        return [
            "--harmonic-distance",
            str(self.atom_i),
            str(self.atom_j),
            f"{self.target:.4f}",
            f"{self.force:.1f}",
        ]


@dataclass
class FeGuard:
    label: str
    min_a: float
    max_a: float
    matching_harmonic_key: tuple[int, int] | None = None


@dataclass
class Candidate:
    nserch: int
    energy: float | None
    grad_max: float | None
    grad_rms: float | None
    s2: float | None
    ligand_problems: int
    fe_failures: list[str]
    outdir: Path

    @property
    def safe(self) -> bool:
        return self.ligand_problems == 0 and not self.fe_failures and self.grad_max is not None


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def append(path: Path, line: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(f"- `{stamp}` {line}\n")


def evaluated_nserches(log: Path) -> list[int]:
    values = [
        int(match.group(1))
        for match in re.finditer(r"NSERCH:\s*(\d+)\s+E=", read_text(log))
    ]
    return sorted(set(values))


def finished(log: Path, status: Path) -> bool:
    text = read_text(log) + "\n" + read_text(status)
    return (
        "EXECUTION OF GAMESS TERMINATED NORMALLY" in text
        or "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-" in text
        or "exit_code=" in text
    )


def normal(log: Path) -> bool:
    return "EXECUTION OF GAMESS TERMINATED NORMALLY" in read_text(log)


def wait_for_finish(root: Path, job: str, poll_seconds: int, timeout_hours: float, report: Path) -> None:
    log = root / "gamess_logs" / f"{job}.log"
    status = root / "gamess_logs" / f"{job}.status"
    append(report, f"WAIT `{job}`")
    deadline = time.time() + timeout_hours * 3600.0
    while time.time() < deadline:
        if finished(log, status):
            append(report, f"FINISHED `{job}` normal={normal(log)} evaluated={evaluated_nserches(log)[-5:]}")
            return
        time.sleep(poll_seconds)
    raise SystemExit(f"Timeout waiting for {job}")


def parse_ligands(items: list[list[str]]) -> list[tuple[str, Path, str]]:
    ligands = []
    for resname, mol2, label in items:
        ligands.append((resname, Path(mol2), label))
    return ligands


def parse_guard(values: list[str]) -> HarmonicGuard:
    if len(values) != 4:
        raise argparse.ArgumentTypeError("guard must be ATOM_I ATOM_J TARGET_A FORCE")
    return HarmonicGuard(int(values[0]), int(values[1]), float(values[2]), float(values[3]))


def parse_fe_guard(values: list[str]) -> FeGuard:
    if len(values) < 3:
        raise argparse.ArgumentTypeError("fe guard must be LABEL MIN_A MAX_A [ATOM_I ATOM_J]")
    key = None
    if len(values) >= 5:
        key = tuple(sorted((int(values[3]), int(values[4]))))
    return FeGuard(values[0], float(values[1]), float(values[2]), key)


def audit_candidate(
    *,
    root: Path,
    template_pdb: Path,
    log: Path,
    nserch: int,
    ligands: list[tuple[str, Path, str]],
    fe_guards: list[FeGuard],
    outdir: Path,
) -> Candidate:
    atoms = parse_pdb(template_pdb)
    initial_model = to_model_atoms(atoms)
    _, final_model, metadata = parse_coords_for_nserch(log, len(atoms), nserch)
    rot, trans, rmsd = best_rigid_transform(
        [atom.xyz for atom in final_model],
        [atom.xyz for atom in initial_model],
    )
    initial_coords = [atom.xyz for atom in initial_model]
    optimized_coords = [apply_transform(atom.xyz, rot, trans) for atom in final_model]

    outdir.mkdir(parents=True, exist_ok=True)
    write_pdb(outdir / "initial_small_model.pdb", atoms, initial_coords, "Initial MCPB small model")
    write_pdb(
        outdir / "optimized_small_model_aligned.pdb",
        atoms,
        optimized_coords,
        f"Guarded GAMESS OPT aligned, NSERCH={nserch}",
    )
    write_displacements(outdir / "atom_displacement.tsv", atoms, initial_coords, optimized_coords)
    write_fe_distances(outdir / "fe_donor_distance_compare.tsv", atoms, initial_coords, optimized_coords, 3.5)

    ligand_problems = 0
    problem_lines: list[str] = []
    for resname, mol2, label in ligands:
        mol2_path = mol2 if mol2.is_absolute() else root / mol2
        rows, problems = audit_mol2_bonds(
            mol2=mol2_path,
            resname=resname,
            atoms=atoms,
            initial_coords=initial_coords,
            optimized_coords=optimized_coords,
            long_cutoff=2.2,
            short_cutoff=0.85,
            large_change_cutoff=0.35,
        )
        ligand_problems += problems
        (outdir / f"{resname}_mol2_bond_identity_audit.tsv").write_text(
            "\n".join(rows) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        if problems:
            problem_lines.extend([f"{label}: {row}" for row in rows[1:] if row.endswith("BROKEN_LONG_BOND,LARGE_CHANGE") or "BROKEN" in row or "LARGE_CHANGE" in row])

    fe_failures: list[str] = []
    fe_indices = [idx for idx, atom in enumerate(atoms) if atom.element.upper() == "FE" or atom.name.upper() == "FE"]
    if fe_indices:
        fe_idx = fe_indices[0]
        for guard in fe_guards:
            matched = [idx for idx, atom in enumerate(atoms) if atom.label == guard.label]
            if not matched:
                fe_failures.append(f"{guard.label}:missing")
                continue
            d = distance(optimized_coords[fe_idx], optimized_coords[matched[0]])
            if d < guard.min_a or d > guard.max_a:
                fe_failures.append(f"{guard.label}:{d:.4f}_outside_{guard.min_a:.2f}-{guard.max_a:.2f}")

    summary = [
        "# Guarded OPT Candidate",
        "",
        f"- NSERCH: `{nserch}`",
        f"- Energy: `{metadata.get('energy_hartree', 'NA')}`",
        f"- GRAD.MAX/RMS: `{metadata.get('grad_max', 'NA')}` / `{metadata.get('grad_rms', 'NA')}`",
        f"- S-squared: `{metadata.get('s2', 'NA')}`",
        f"- Alignment RMSD: `{rmsd:.6f} A`",
        f"- Ligand problem count: `{ligand_problems}`",
        f"- Fe guard failures: `{'; '.join(fe_failures) if fe_failures else 'none'}`",
        "",
        "## Ligand Problems",
        "",
        *(problem_lines or ["none"]),
    ]
    (outdir / "README_guard_audit.md").write_text("\n".join(summary) + "\n", encoding="utf-8", newline="\n")

    return Candidate(
        nserch=nserch,
        energy=metadata.get("energy_hartree") if isinstance(metadata.get("energy_hartree"), float) else None,
        grad_max=metadata.get("grad_max") if isinstance(metadata.get("grad_max"), float) else None,
        grad_rms=metadata.get("grad_rms") if isinstance(metadata.get("grad_rms"), float) else None,
        s2=metadata.get("s2") if isinstance(metadata.get("s2"), float) else None,
        ligand_problems=ligand_problems,
        fe_failures=fe_failures,
        outdir=outdir,
    )


def run(cmd: list[str], cwd: Path, report: Path) -> None:
    append(report, "RUN `" + " ".join(cmd) + "`")
    subprocess.run(cmd, cwd=str(cwd), check=True)


def start_next_job(
    *,
    root: Path,
    source_job: str,
    source_log: Path,
    source_dat: Path,
    out_job: str,
    charge: int,
    mult: int,
    norb: int,
    nserch: int,
    nstep: int,
    opttol: str,
    cores: int,
    session: str,
    guards: list[HarmonicGuard],
    report: Path,
) -> None:
    helper = SCRIPT_DIR / "prepare_gamess_moread_optimize_from_log.py"
    cmd = [
        sys.executable,
        str(helper),
        "--mcpb-dir",
        str(root),
        "--source-job",
        source_job,
        "--source-log",
        str(source_log),
        "--out-job",
        out_job,
        "--dat",
        str(source_dat),
        "--charge",
        str(charge),
        "--mult",
        str(mult),
        "--norb",
        str(norb),
        "--from-nserch",
        str(nserch),
        "--nstep",
        str(nstep),
        "--opttol",
        opttol,
        "--scf-mode",
        "soscf",
        "--cores",
        str(cores),
        "--session",
        session,
    ]
    for guard in guards:
        cmd.extend(guard.as_args())
    run(cmd, root, report)
    run(["bash", f"run_{out_job}_foreground.sh"], root, report)


def adjust_guards(candidates: list[Candidate], guards: list[HarmonicGuard], fe_guards: list[FeGuard], max_force: float) -> list[HarmonicGuard]:
    new = [HarmonicGuard(g.atom_i, g.atom_j, g.target, g.force) for g in guards]
    any_ligand_problem = any(candidate.ligand_problems for candidate in candidates)
    any_fe_problem_labels = {
        failure.split(":", 1)[0]
        for candidate in candidates
        for failure in candidate.fe_failures
    }
    for idx, guard in enumerate(new):
        increase = False
        # Current TJ nitrene critical covalent guard: UNL C1-C2.
        if any_ligand_problem and guard.key() == tuple(sorted((1, 4))):
            increase = True
        for fe_guard in fe_guards:
            if fe_guard.label in any_fe_problem_labels and fe_guard.matching_harmonic_key == guard.key():
                increase = True
        if increase:
            new[idx] = HarmonicGuard(guard.atom_i, guard.atom_j, guard.target, min(max_force, guard.force * 1.5))
    return new


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mcpb-dir", type=Path, required=True)
    parser.add_argument("--template-pdb", type=Path, required=True)
    parser.add_argument("--current-job", required=True)
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--mult", type=int, required=True)
    parser.add_argument("--norb", type=int, required=True)
    parser.add_argument("--ligand", nargs=3, action="append", required=True, metavar=("RESNAME", "MOL2", "LABEL"))
    parser.add_argument(
        "--fe-guard",
        nargs="+",
        action="append",
        default=[],
        metavar="FE_GUARD",
        help="Required Fe-donor distance window: LABEL MIN_A MAX_A [ATOM_I ATOM_J].",
    )
    parser.add_argument("--harmonic-distance", nargs=4, action="append", default=[], metavar=("ATOM_I", "ATOM_J", "TARGET", "FORCE"))
    parser.add_argument("--nstep", type=int, default=20)
    parser.add_argument("--opttol", default="0.0005")
    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--max-cycles", type=int, default=4)
    parser.add_argument("--poll-seconds", type=int, default=120)
    parser.add_argument("--timeout-hours", type=float, default=24.0)
    parser.add_argument("--hessian-grad-max", type=float, default=0.005)
    parser.add_argument("--hessian-grad-rms", type=float, default=0.002)
    parser.add_argument("--max-force", type=float, default=1000.0)
    parser.add_argument("--out-prefix", default="guard_loop")
    args = parser.parse_args()

    root = args.mcpb_dir
    template_pdb = args.template_pdb if args.template_pdb.is_absolute() else root / args.template_pdb
    ligands = parse_ligands(args.ligand)
    fe_guards = [parse_fe_guard(values) for values in args.fe_guard]
    guards = [parse_guard(values) for values in args.harmonic_distance]
    report = root / f"{args.out_prefix}_report.md"
    append(report, f"START current_job={args.current_job} max_cycles={args.max_cycles}")

    current_job = args.current_job
    for cycle in range(1, args.max_cycles + 1):
        log = root / "gamess_logs" / f"{current_job}.log"
        status = root / "gamess_logs" / f"{current_job}.status"
        wait_for_finish(root, current_job, args.poll_seconds, args.timeout_hours, report)
        nserches = evaluated_nserches(log)
        if not nserches:
            raise SystemExit(f"No evaluated NSERCH values in {log}")
        cycle_dir = root / f"{args.out_prefix}_cycle{cycle:02d}_{current_job}"
        candidates = [
            audit_candidate(
                root=root,
                template_pdb=template_pdb,
                log=log,
                nserch=nserch,
                ligands=ligands,
                fe_guards=fe_guards,
                outdir=cycle_dir / f"nserch_{nserch:03d}",
            )
            for nserch in nserches
        ]
        safe = [candidate for candidate in candidates if candidate.safe]
        rows = ["nserch\tenergy\tgrad_max\tgrad_rms\ts2\tligand_problems\tfe_failures\tsafe"]
        for candidate in candidates:
            rows.append(
                f"{candidate.nserch}\t{candidate.energy}\t{candidate.grad_max}\t{candidate.grad_rms}\t"
                f"{candidate.s2}\t{candidate.ligand_problems}\t{';'.join(candidate.fe_failures)}\t{candidate.safe}"
            )
        (cycle_dir / "candidate_summary.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
        if not safe:
            guards = adjust_guards(candidates, guards, fe_guards, args.max_force)
            append(report, f"CYCLE {cycle}: no safe candidate. adjusted guards={guards}")
            raise SystemExit("No chemically safe evaluated geometry; manual model review is needed before restarting.")
        best = sorted(safe, key=lambda item: (item.grad_max if item.grad_max is not None else 999.0, item.energy or 0.0))[0]
        append(
            report,
            f"CYCLE {cycle}: best_safe NSERCH={best.nserch} grad={best.grad_max} rms={best.grad_rms} outdir={best.outdir}",
        )
        if (
            best.grad_max is not None
            and best.grad_rms is not None
            and best.grad_max <= args.hessian_grad_max
            and best.grad_rms <= args.hessian_grad_rms
        ):
            marker = root / "HESSIAN_CANDIDATE.txt"
            marker.write_text(
                "\n".join(
                    [
                        f"job={current_job}",
                        f"nserch={best.nserch}",
                        f"energy={best.energy}",
                        f"grad_max={best.grad_max}",
                        f"grad_rms={best.grad_rms}",
                        f"s2={best.s2}",
                        f"audit_dir={best.outdir}",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            append(report, f"HESSIAN_CANDIDATE written `{marker}`")
            return 0

        guards = adjust_guards(candidates, guards, fe_guards, args.max_force)
        source_dat = RESTART_DIR / f"{current_job}.dat"
        if not source_dat.exists():
            raise SystemExit(f"Missing MOREAD dat for next cycle: {source_dat}")
        next_job = f"{args.out_prefix}_cycle{cycle:02d}_cont_from_{current_job}_n{best.nserch:03d}"
        append(report, f"CYCLE {cycle}: continuing as `{next_job}` with guards={guards}")
        start_next_job(
            root=root,
            source_job=current_job,
            source_log=log,
            source_dat=source_dat,
            out_job=next_job,
            charge=args.charge,
            mult=args.mult,
            norb=args.norb,
            nserch=best.nserch,
            nstep=args.nstep,
            opttol=args.opttol,
            cores=args.cores,
            session=next_job[:40],
            guards=guards,
            report=report,
        )
        current_job = next_job

    append(report, "STOP max cycles reached without Hessian candidate")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
