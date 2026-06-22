"""Rosetta-oriented loop repair helpers for PLE receptor preparation."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shlex
import stat
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np


STANDARD_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
ONE_TO_THREE = {value: key for key, value in THREE_TO_ONE.items()}
DEFAULT_ROSETTA_MAIN_CANDIDATES = (
    "/mnt/l/WSL/softwares/rosetta.source.release-408/main",
    "~/rosetta/source",
    "~/rosetta/main",
)


@dataclass(frozen=True)
class AtomRecord:
    line: str
    name: str
    resname: str
    chain: str
    resseq: int
    icode: str
    coord: tuple[float, float, float]
    element: str


@dataclass(frozen=True)
class MissingResidue:
    chain: str
    resseq: int
    resname: str
    region: str
    left_observed: int | None
    right_observed: int | None


@dataclass(frozen=True)
class SimpleGap:
    chain: str
    resseq: int
    resname: str
    left_resseq: int
    left_resname: str
    right_resseq: int
    right_resname: str


def add_repair_loops_parser(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--project-dir", type=Path, help="Project root; defaults to current directory")
    parser.add_argument("--state-dir", default="ple_project", help="Workflow output folder inside project-dir")
    parser.add_argument("--pdb", type=Path, help="Receptor PDB; required when project has multiple PDB files")
    parser.add_argument("--chain", help="Repair only this chain, for example A")
    parser.add_argument("--resseq", type=int, help="Repair only this missing residue number, for example 108")
    parser.add_argument("--outdir", type=Path, help="Output loop-repair workspace")
    parser.add_argument("--rosetta-main", type=Path, help="Rosetta main directory containing source/bin and database")
    parser.add_argument("--rosetta-np", type=int, default=1, help="MPI process count written into the run script")
    parser.add_argument("--rosetta-nstruct", type=int, default=1, help="Number of loopmodel structures to generate")
    parser.add_argument("--loop-flank", type=int, default=4, help="Residues on each side of the inserted residue included in the loop file")
    parser.add_argument("--run-rosetta", action="store_true", help="Run generated loopmodel scripts immediately")
    parser.add_argument("--force", action="store_true", help="Overwrite existing generated files")
    parser.add_argument("--dry-run", action="store_true", help="Print planned work without writing or running")


def repair_loops_command(args: argparse.Namespace) -> None:
    project_dir = (args.project_dir or Path.cwd()).expanduser().resolve()
    pdb = detect_pdb(project_dir, args.pdb)
    state_dir = project_dir / args.state_dir
    outdir = args.outdir.expanduser() if args.outdir else state_dir / "loop_repair" / pdb.stem
    if not outdir.is_absolute():
        outdir = project_dir / outdir
    outdir = outdir.resolve()

    seqres, atoms_by_chain = parse_pdb(pdb)
    missing = find_missing_residues(seqres, atoms_by_chain)
    targets = select_simple_gaps(missing, atoms_by_chain, args.chain, args.resseq)
    rosetta_main = resolve_rosetta_main(args.rosetta_main)

    print(f"pdb={pdb}")
    print(f"outdir={outdir}")
    print(f"missing_internal={sum(1 for item in missing if item.region == 'internal')}")
    if not targets:
        print("repair_targets=none")
        print("reason=no single-residue internal gap matched --chain/--resseq")
        return
    print("repair_targets=" + ",".join(f"{item.chain}:{item.resseq}{item.resname}" for item in targets))
    if args.dry_run:
        return

    outdir.mkdir(parents=True, exist_ok=True)
    write_missing_csv(outdir / "missing_residues.csv", missing)
    write_summary(outdir / "README.md", pdb, missing, targets, rosetta_main)

    generated_scripts: list[Path] = []
    for gap in targets:
        workspace = outdir / f"{clean_token(gap.chain)}_{gap.resseq}_{gap.resname}"
        if workspace.exists() and not args.force:
            print(f"workspace_exists={workspace}")
            generated_scripts.append(workspace / "run_loopmodel.sh")
            continue
        workspace.mkdir(parents=True, exist_ok=True)
        inserted = write_template_inserted_chain(
            pdb=pdb,
            workspace=workspace,
            gap=gap,
            atoms_by_chain=atoms_by_chain,
        )
        loop_file = workspace / "repair.loop"
        write_loop_file(loop_file, gap, args.loop_flank, atoms_by_chain)
        script = write_loopmodel_script(
            workspace=workspace,
            input_pdb=inserted,
            loop_file=loop_file,
            rosetta_main=rosetta_main,
            np_count=args.rosetta_np,
            nstruct=args.rosetta_nstruct,
        )
        write_audit_script(workspace / "audit_loopmodel.py", gap)
        generated_scripts.append(script)
        print(f"workspace={workspace}")
        print(f"inserted_pdb={inserted}")
        print(f"loop_file={loop_file}")
        print(f"run_script={script}")

    if args.run_rosetta:
        for script in generated_scripts:
            if not script.exists():
                raise FileNotFoundError(script)
            subprocess.run([str(script)], check=True)


def detect_pdb(project_dir: Path, explicit: Path | None) -> Path:
    if explicit is not None:
        path = explicit.expanduser()
        if not path.is_absolute():
            path = project_dir / path
        if not path.exists():
            raise FileNotFoundError(path)
        if path.suffix.lower() != ".pdb":
            raise ValueError("loop repair needs a raw PDB file")
        return path.resolve()
    candidates = sorted(path for path in project_dir.iterdir() if path.is_file() and path.suffix.lower() == ".pdb")
    if not candidates:
        raise FileNotFoundError(f"no receptor *.pdb found in {project_dir}; pass --pdb")
    if len(candidates) > 1:
        names = ", ".join(path.name for path in candidates)
        raise ValueError(f"multiple PDB files found: {names}; pass --pdb")
    return candidates[0].resolve()


def parse_pdb(path: Path) -> tuple[dict[str, list[str]], dict[str, list[AtomRecord]]]:
    seqres: dict[str, list[str]] = {}
    atoms_by_chain: dict[str, list[AtomRecord]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("SEQRES"):
                chain = line[11].strip() or " "
                seqres.setdefault(chain, []).extend(line[19:70].split())
                continue
            if not line.startswith("ATOM"):
                continue
            resname = line[17:20].strip()
            if resname not in STANDARD_RESNAMES:
                continue
            try:
                resseq = int(line[22:26])
                coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            chain = line[21].strip() or " "
            atom_name = line[12:16].strip()
            element = (line[76:78].strip() or atom_name[:1]).upper()
            atoms_by_chain.setdefault(chain, []).append(
                AtomRecord(
                    line=line.rstrip("\n"),
                    name=atom_name,
                    resname=resname,
                    chain=chain,
                    resseq=resseq,
                    icode=line[26].strip(),
                    coord=coord,
                    element=element,
                )
            )
    return seqres, atoms_by_chain


def observed_residues(atoms: list[AtomRecord]) -> dict[int, str]:
    residues: dict[int, str] = {}
    for atom in atoms:
        if atom.icode:
            continue
        residues.setdefault(atom.resseq, atom.resname)
    return residues


def find_missing_residues(seqres: dict[str, list[str]], atoms_by_chain: dict[str, list[AtomRecord]]) -> list[MissingResidue]:
    rows: list[MissingResidue] = []
    for chain, sequence in sorted(seqres.items()):
        observed = observed_residues(atoms_by_chain.get(chain, []))
        if not observed:
            continue
        first = min(observed)
        last = max(observed)
        observed_numbers = sorted(observed)
        for resseq, resname in enumerate(sequence, start=1):
            if resseq in observed:
                continue
            if first < resseq < last:
                region = "internal"
            elif resseq < first:
                region = "n_term"
            else:
                region = "c_term"
            left = max((item for item in observed_numbers if item < resseq), default=None)
            right = min((item for item in observed_numbers if item > resseq), default=None)
            rows.append(MissingResidue(chain, resseq, resname, region, left, right))
    return rows


def select_simple_gaps(
    missing: list[MissingResidue],
    atoms_by_chain: dict[str, list[AtomRecord]],
    chain_filter: str | None,
    resseq_filter: int | None,
) -> list[SimpleGap]:
    by_chain: dict[str, list[MissingResidue]] = {}
    for item in missing:
        if item.region != "internal":
            continue
        if chain_filter and item.chain != chain_filter:
            continue
        if resseq_filter is not None and item.resseq != resseq_filter:
            continue
        by_chain.setdefault(item.chain, []).append(item)

    targets: list[SimpleGap] = []
    for chain, rows in by_chain.items():
        observed = observed_residues(atoms_by_chain.get(chain, []))
        for item in sorted(rows, key=lambda row: row.resseq):
            left = item.resseq - 1
            right = item.resseq + 1
            if left not in observed or right not in observed:
                continue
            # Only auto-repair isolated one-residue gaps. Multi-residue gaps need
            # explicit modelling choices and more sampling.
            if any(other.resseq in {item.resseq - 1, item.resseq + 1} for other in rows):
                continue
            targets.append(SimpleGap(chain, item.resseq, item.resname, left, observed[left], right, observed[right]))
    return targets


def write_missing_csv(path: Path, missing: list[MissingResidue]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["chain", "resseq", "seqres_resname", "region", "left_observed", "right_observed"],
        )
        writer.writeheader()
        for item in missing:
            writer.writerow(
                {
                    "chain": item.chain,
                    "resseq": item.resseq,
                    "seqres_resname": item.resname,
                    "region": item.region,
                    "left_observed": item.left_observed or "",
                    "right_observed": item.right_observed or "",
                }
            )


def write_summary(path: Path, pdb: Path, missing: list[MissingResidue], targets: list[SimpleGap], rosetta_main: Path) -> None:
    internal = [item for item in missing if item.region == "internal"]
    lines = [
        "# PLE Loop Repair",
        "",
        f"Input PDB: `{pdb}`",
        f"Rosetta main: `{rosetta_main}`",
        "",
        f"Internal missing residues: {len(internal)}",
        f"Auto repair targets: {', '.join(f'{item.chain}:{item.resseq}{item.resname}' for item in targets) or 'none'}",
        "",
        "This command auto-repairs only isolated one-residue internal gaps. Multi-residue gaps are reported in `missing_residues.csv` and should be handled with larger Rosetta/Modeller sampling or manual curation.",
        "",
        "For each repair target, run:",
        "",
        "```bash",
        "cd <target-workspace>",
        "./run_loopmodel.sh",
        "python3 audit_loopmodel.py",
        "```",
        "",
        "The generated run script treats Rosetta MPI finalization warnings as non-fatal when the log contains `reported success`.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def resolve_rosetta_main(explicit: Path | None) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit.expanduser())
    env = os.environ.get("ROSETTA_MAIN")
    if env:
        candidates.append(Path(env).expanduser())
    candidates.extend(Path(item).expanduser() for item in DEFAULT_ROSETTA_MAIN_CANDIDATES)
    for candidate in candidates:
        if (candidate / "source" / "bin").exists() and (candidate / "database").exists():
            return candidate.resolve()
    raise FileNotFoundError("Rosetta main directory not found; pass --rosetta-main or set ROSETTA_MAIN")


def write_template_inserted_chain(
    *,
    pdb: Path,
    workspace: Path,
    gap: SimpleGap,
    atoms_by_chain: dict[str, list[AtomRecord]],
) -> Path:
    template_atoms = build_tripeptide_template(gap)
    chain_atoms = [atom for atom in atoms_by_chain.get(gap.chain, []) if not atom.icode and atom.element != "H"]
    target_lookup = {(atom.resseq, atom.name): atom for atom in chain_atoms}
    fit_template: list[tuple[float, float, float]] = []
    fit_target: list[tuple[float, float, float]] = []
    for template_resseq, target_resseq in ((1, gap.left_resseq), (3, gap.right_resseq)):
        for name in ("N", "CA", "C", "O", "CB"):
            template_atom = template_atoms.get((template_resseq, name))
            target_atom = target_lookup.get((target_resseq, name))
            if template_atom is not None and target_atom is not None:
                fit_template.append(template_atom.coord)
                fit_target.append(target_atom.coord)
    if len(fit_template) < 6:
        raise ValueError(f"not enough flank atoms to place {gap.chain}:{gap.resseq}")
    transform = kabsch_transform(fit_template, fit_target)

    serial = 1
    lines: list[str] = []
    inserted = False
    for atom in chain_atoms:
        if atom.resseq == gap.right_resseq and not inserted:
            for template_atom in template_atoms.values():
                if template_atom.resseq != 2:
                    continue
                x, y, z = transform(template_atom.coord)
                lines.append(
                    format_atom_line(
                        serial=serial,
                        atom_name=template_atom.name,
                        resname=gap.resname,
                        chain=gap.chain,
                        resseq=gap.resseq,
                        coord=(x, y, z),
                        element=template_atom.element,
                        bfactor=50.0,
                    )
                )
                serial += 1
            inserted = True
        lines.append(f"{atom.line[:6]}{serial:5d}{atom.line[11:]}")
        serial += 1
    if not inserted:
        raise ValueError(f"right flank not found for {gap.chain}:{gap.resseq}")
    lines.extend(["TER", "END"])
    output = workspace / f"{pdb.stem}_chain{clean_token(gap.chain)}_{gap.resname}{gap.resseq}_template_insert.pdb"
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output


def build_tripeptide_template(gap: SimpleGap) -> dict[tuple[int, str], AtomRecord]:
    try:
        from Bio.PDB import PDBIO
        from Bio.PDB.PICIO import read_PIC_seq
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
    except Exception as exc:
        raise RuntimeError("Biopython is required for template loop insertion; install biopython in this environment") from exc

    sequence = "".join(THREE_TO_ONE[item] for item in (gap.left_resname, gap.resname, gap.right_resname))
    structure = read_PIC_seq(SeqRecord(Seq(sequence), id=sequence), pdbid="TPL", chain="A")
    chain = next(next(structure.get_models()).get_chains())
    chain.internal_coord.internal_to_atom_coordinates()
    temp = Path.cwd() / f".ple_{os.getpid()}_{clean_token(gap.chain)}_{gap.resseq}_template.pdb"
    try:
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(temp))
        _, atoms_by_chain = parse_pdb(temp)
    finally:
        temp.unlink(missing_ok=True)
    atoms: dict[tuple[int, str], AtomRecord] = {}
    for atom in atoms_by_chain.get("A", []):
        atoms[(atom.resseq, atom.name)] = atom
    return atoms


def kabsch_transform(
    source: list[tuple[float, float, float]],
    target: list[tuple[float, float, float]],
):
    p = np.asarray(source, dtype=float)
    q = np.asarray(target, dtype=float)
    pc = p.mean(axis=0)
    qc = q.mean(axis=0)
    h = (p - pc).T @ (q - qc)
    u, _, vt = np.linalg.svd(h)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T
    t = qc - pc @ r.T

    def transform(coord: tuple[float, float, float]) -> tuple[float, float, float]:
        value = np.asarray(coord, dtype=float) @ r.T + t
        return float(value[0]), float(value[1]), float(value[2])

    return transform


def format_atom_line(
    *,
    serial: int,
    atom_name: str,
    resname: str,
    chain: str,
    resseq: int,
    coord: tuple[float, float, float],
    element: str,
    bfactor: float,
) -> str:
    x, y, z = coord
    return (
        f"ATOM  {serial:5d} {atom_name:<4s} {resname:>3s} {chain:1s}{resseq:4d}"
        f"    {x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{bfactor:6.2f}          {element:>2s}"
    )


def write_loop_file(path: Path, gap: SimpleGap, flank: int, atoms_by_chain: dict[str, list[AtomRecord]]) -> None:
    observed = observed_residues(atoms_by_chain.get(gap.chain, []))
    start = max(min(observed), gap.resseq - flank)
    stop = min(max(observed), gap.resseq + flank)
    path.write_text(f"LOOP {start} {stop} {gap.resseq} 0.0 0\n", encoding="utf-8")


def write_loopmodel_script(
    *,
    workspace: Path,
    input_pdb: Path,
    loop_file: Path,
    rosetta_main: Path,
    np_count: int,
    nstruct: int,
) -> Path:
    executable = rosetta_main / "source" / "bin" / "loopmodel.mpi.linuxgccrelease"
    database = rosetta_main / "database"
    script = workspace / "run_loopmodel.sh"
    command = [
        "mpirun",
        "-np",
        '"${ROSETTA_NP:-' + str(np_count) + '}"',
        shlex.quote(str(executable)),
        "-database",
        shlex.quote(str(database)),
        "-in:file:s",
        shlex.quote(input_pdb.name),
        "-in:file:fullatom",
        "-out:file:fullatom",
        "-loops:loop_file",
        shlex.quote(loop_file.name),
        "-loops:remodel",
        "perturb_kic",
        "-loops:refine",
        "refine_kic",
        "-loops:max_kic_build_attempts",
        "200",
        "-nstruct",
        '"${ROSETTA_NSTRUCT:-' + str(nstruct) + '}"',
        "-score:weights",
        "ref2015",
        "-out:path:pdb",
        "loopmodel_out",
        "-out:prefix",
        "loopmodel_",
        "-overwrite",
    ]
    script.write_text(
        "#!/usr/bin/env bash\n"
        "set -uo pipefail\n"
        "cd \"$(dirname \"$0\")\"\n"
        "mkdir -p loopmodel_out\n"
        + " ".join(command)
        + " 2>&1 | tee loopmodel.log\n"
        "status=${PIPESTATUS[0]}\n"
        "if grep -q 'reported success' loopmodel.log; then\n"
        "  exit 0\n"
        "fi\n"
        "exit ${status}\n",
        encoding="utf-8",
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return script


def write_audit_script(path: Path, gap: SimpleGap) -> None:
    script = f'''#!/usr/bin/env python3
from pathlib import Path
import math
import re

THREE={{'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}}
expected_gap=({gap.resseq!r}, {gap.resname!r})
models=sorted(Path('loopmodel_out').glob('loopmodel_*.pdb'))
if not models:
    raise SystemExit('FAIL no loopmodel_*.pdb in loopmodel_out')
p=models[0]
rows=[]; seen=set(); atoms={{}}
for line in p.read_text(errors='replace').splitlines():
    if not line.startswith('ATOM'):
        continue
    key=(line[21].strip() or '_', int(line[22:26]), line[26].strip())
    resn=line[17:20].strip(); atom=line[12:16].strip()
    if key not in seen:
        seen.add(key); rows.append((key,resn))
    atoms[(key,atom)]=tuple(float(line[i:i+8]) for i in (30,38,46))
seq=''.join(THREE.get(resn,'X') for _,resn in rows)
gap_rows=[(idx,key,resn) for idx,(key,resn) in enumerate(rows,start=1) if key[1] == expected_gap[0]]
print('model', p)
print('residue_count', len(rows))
print('gap_rows', gap_rows)
ok=True
if not gap_rows or gap_rows[0][2] != expected_gap[1]:
    print('FAIL missing expected inserted residue', expected_gap)
    ok=False
for left,right in (({gap.left_resseq}, {gap.resseq}), ({gap.resseq}, {gap.right_resseq})):
    left_row=next(((key,resn) for key,resn in rows if key[1]==left), None)
    right_row=next(((key,resn) for key,resn in rows if key[1]==right), None)
    if not left_row or not right_row:
        print('FAIL missing flank', left, right)
        ok=False
        continue
    c=atoms.get((left_row[0],'C')); n=atoms.get((right_row[0],'N'))
    if not c or not n:
        print('FAIL missing peptide atoms', left, right)
        ok=False
        continue
    dist=math.dist(c,n)
    print(f'C-N {{left}}->{{right}} {{dist:.3f}}')
    if not 1.0 <= dist <= 1.5:
        ok=False
log=Path('loopmodel.log')
if log.exists():
    text=log.read_text(errors='replace')
    m=re.findall(r'chainbreak\\s+([-+0-9.eE]+)', text)
    if m:
        print('chainbreak', m[-1])
        try:
            if float(m[-1]) > 0.1:
                ok=False
        except ValueError:
            ok=False
    print('reported_success', 'reported success' in text)
    if 'reported success' not in text:
        ok=False
print('OK' if ok else 'FAIL')
raise SystemExit(0 if ok else 1)
'''
    path.write_text(script, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def clean_token(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip()) or "chain"
