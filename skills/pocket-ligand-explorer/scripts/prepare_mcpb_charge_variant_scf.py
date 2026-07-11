#!/usr/bin/env python3
"""Create an MCPB charge/spin variant and prepare a fixed-geometry SCF test."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
from pathlib import Path


ATOMIC_NUMBERS = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "S": 16,
    "FE": 26,
}


def parse_mol2_atoms(path: Path) -> tuple[list[str], list[int], float]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    in_atom = False
    atom_line_indices: list[int] = []
    charge = 0.0
    for idx, line in enumerate(lines):
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom = True
            continue
        if line.startswith("@<TRIPOS>"):
            in_atom = False
            continue
        if not in_atom or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 9:
            atom_line_indices.append(idx)
            charge += float(parts[8])
    return lines, atom_line_indices, charge


def write_charge_adjusted_mol2(
    src: Path,
    dst: Path,
    *,
    target_charge: float,
    charge_atoms: set[str],
) -> tuple[float, float, list[str]]:
    lines, atom_indices, old_charge = parse_mol2_atoms(src)
    selected: list[int] = []
    for idx in atom_indices:
        parts = lines[idx].split()
        if not charge_atoms or parts[1] in charge_atoms:
            selected.append(idx)
    if not selected:
        raise SystemExit(f"No selected charge atoms {sorted(charge_atoms)} found in {src}")
    delta_each = (target_charge - old_charge) / len(selected)
    adjusted_names: list[str] = []
    for idx in selected:
        parts = lines[idx].split()
        old = float(parts[8])
        new = old + delta_each
        adjusted_names.append(parts[1])
        lines[idx] = (
            f"{int(parts[0]):7d} {parts[1]:<8s}"
            f"{float(parts[2]):10.4f}{float(parts[3]):10.4f}{float(parts[4]):10.4f} "
            f"{parts[5]:<8s}{int(parts[6]):4d} {parts[7]:<8s}"
            f"{new:11.6f}"
        )
    dst.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    _, _, new_charge = parse_mol2_atoms(dst)
    return old_charge, new_charge, adjusted_names


def replace_line(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"^{re.escape(key)}\s+.*$", re.M)
    if pattern.search(text):
        return pattern.sub(f"{key} {value}", text)
    return text.rstrip() + f"\n{key} {value}\n"


def create_mcpb_in(src: Path, dst: Path, *, group: str, charge: int, mult: int) -> None:
    text = src.read_text(encoding="utf-8", errors="replace")
    text = replace_line(text, "group_name", group)
    text = replace_line(text, "smmodel_chg", str(charge))
    text = replace_line(text, "smmodel_spin", str(mult))
    text = replace_line(text, "lgmodel_chg", str(charge))
    text = replace_line(text, "lgmodel_spin", str(mult))
    dst.write_text(text.rstrip() + "\n", encoding="utf-8", newline="\n")


def remove_card(lines: list[str], card: str) -> list[str]:
    out: list[str] = []
    i = 0
    target = card.upper()
    while i < len(lines):
        if lines[i].strip().upper().startswith(target):
            if "$END" in lines[i].upper():
                i += 1
                continue
            i += 1
            while i < len(lines) and "$END" not in lines[i].upper():
                i += 1
            i += 1
            continue
        out.append(lines[i])
        i += 1
    return out


def replace_card(lines: list[str], card: str, replacement: str) -> list[str]:
    out: list[str] = []
    i = 0
    replaced = False
    target = card.upper()
    while i < len(lines):
        if lines[i].strip().upper().startswith(target):
            out.extend(replacement.splitlines())
            replaced = True
            if "$END" not in lines[i].upper():
                i += 1
                while i < len(lines) and "$END" not in lines[i].upper():
                    i += 1
            i += 1
            continue
        out.append(lines[i])
        i += 1
    if not replaced:
        out = replacement.splitlines() + out
    return out


def write_energy_input(mcpb_dir: Path, *, group: str, charge: int, mult: int) -> Path:
    source = mcpb_dir / f"{group}_small_opt.inp"
    if not source.exists():
        raise FileNotFoundError(source)
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    lines = remove_card(lines, "$STATPT")
    lines = replace_card(lines, "$SYSTEM", " $SYSTEM MEMDDI=160 MWORDS=80 $END")
    lines = replace_card(
        lines,
        "$CONTRL",
        (
            " $CONTRL\n"
            "  SCFTYP=UHF DFTTYP=B3LYP RUNTYP=ENERGY\n"
            f"  ICHARG={charge} MULT={mult} MAXIT=200\n"
            " $END"
        ),
    )
    lines = replace_card(
        lines,
        "$SCF",
        (
            " $SCF\n"
            "  DIRSCF=.T. DIIS=.F. SOSCF=.T. DAMP=.T. SHIFT=.T.\n"
            "  ETHRSH=10.0 MAXDII=30\n"
            " $END"
        ),
    )
    lines = replace_card(
        lines,
        "$DFT",
        (
            " $DFT\n"
            "  NRAD=96 NLEB=302 NRAD0=96 NLEB0=302 SWOFF=0.0\n"
            " $END"
        ),
    )
    lines = replace_card(lines, "$GUESS", " $GUESS GUESS=HUCKEL $END")
    out = mcpb_dir / f"{group}_scf_energy_soscf_finegrid.inp"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return out


def infer_element(atom_name: str, element_field: str = "") -> str:
    elem = "".join(ch for ch in element_field.strip() if ch.isalpha()).upper()
    if elem:
        return "FE" if elem == "FE" else elem[:1]
    letters = "".join(ch for ch in atom_name if ch.isalpha()).upper()
    if letters.startswith("FE"):
        return "FE"
    return letters[:1]


def electron_count(pdb: Path, charge: int) -> int:
    total = 0
    for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        elem = infer_element(line[12:16].strip(), line[76:78] if len(line) >= 78 else "")
        if elem not in ATOMIC_NUMBERS:
            raise SystemExit(f"Unknown element {elem!r} in {pdb}: {line}")
        total += ATOMIC_NUMBERS[elem]
    return total - charge


def write_scripts(
    mcpb_dir: Path,
    *,
    group: str,
    model_charge: int,
    mult: int,
    cores: int,
    gamess: str,
) -> None:
    step1 = mcpb_dir / "run_01_mcpb_step1.sh"
    step1.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
mkdir -p logs gamess_logs
set +u
source /mnt/l/WSL/conda_envs/AmberTools25/amber.sh
set -u
MCPB.py -i mcpb.in -s 1 | tee logs/mcpb_step1.log
/mnt/l/WSL/softwares/conda_envs/md/bin/python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/build_mcpb_visual_check_pdb.py \
  --small-pdb {group}_small.pdb \
  --mol2-dir "$PWD" \
  --out {group}_small_visual_check.pdb \
  --report logs/{group}_small_visual_check_report.tsv \
  --model-charge {model_charge} \
  --mult {mult}
/mnt/l/WSL/softwares/conda_envs/md/bin/python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/audit_mcpb_small_model_hydrogens.py \
  --small-pdb {group}_small.pdb \
  --mol2-dir "$PWD" \
  --out logs/mcpb_small_model_hydrogen_audit.tsv \
  --max-warnings "${{PLE_MCPB_H_WARNING_MAX:-2}}" \
  --model-charge {model_charge} \
  --mult {mult}
""",
        encoding="utf-8",
        newline="\n",
    )
    step1.chmod(0o755)

    prep = mcpb_dir / "run_02_prepare_energy_input.sh"
    prep.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
/mnt/l/WSL/softwares/conda_envs/md/bin/python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/prepare_mcpb_charge_variant_scf.py \\
  --source-dir /dev/null \\
  --outdir "$PWD/.." \\
  --group {group} \\
  --model-charge {model_charge} \\
  --mult {mult} \\
  --prepare-energy-only
""",
        encoding="utf-8",
        newline="\n",
    )
    prep.chmod(0o755)

    job = f"{group}_scf_energy_soscf_finegrid"
    foreground = mcpb_dir / f"run_{job}_foreground.sh"
    foreground.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={job}
NCORES=${{NCORES:-{cores}}}
GAMESS=${{GAMESS:-{gamess}}}
VERSION=${{VERSION:-00}}
mkdir -p gamess_logs logs
echo "[GAMESS] $JOB start $(date) NCORES=$NCORES" | tee "gamess_logs/${{JOB}}.status"
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

    tmux = mcpb_dir / f"run_{job}_tmux.sh"
    tmux.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
SESSION=${{SESSION:-{job.lower()}}}
NCORES=${{NCORES:-{cores}}}
if tmux has-session -t "$SESSION" 2>/dev/null; then
  echo "tmux_session_exists=$SESSION"
  echo "tail=tail -f $PWD/gamess_logs/{job}.log"
  exit 0
fi
tmux new-session -d -s "$SESSION" "cd '$PWD'; env NCORES='$NCORES' bash {foreground.name}; echo; echo '[tmux] finished'; exec bash"
echo "tmux_session=$SESSION"
echo "tail=tail -f $PWD/gamess_logs/{job}.log"
""",
        encoding="utf-8",
        newline="\n",
    )
    tmux.chmod(0o755)

    monitor = mcpb_dir / f"monitor_{job}.sh"
    monitor.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={job}
echo "===== status ====="
cat "gamess_logs/${{JOB}}.status" 2>/dev/null || true
echo
echo "===== running ====="
pgrep -af "${{JOB}}|ddikick|gamess.00.x" || true
echo
echo "===== markers ====="
grep -E 'RUNTYP=|ITER EX|DENSITY CONVERGED|SCF IS UNCONVERGED|SCF HAS NOT CONVERGED|FINAL U-B3LYP|S-SQUARED|TERMINATED' "gamess_logs/${{JOB}}.log" 2>/dev/null | tail -n 100 || true
echo
tail -n 60 "gamess_logs/${{JOB}}.log" 2>/dev/null || true
""",
        encoding="utf-8",
        newline="\n",
    )
    monitor.chmod(0o755)


def run_antechamber_ligand(
    input_dir: Path,
    mcpb_dir: Path,
    *,
    ligand: str,
    ligand_charge: float,
    ligand_mult: int,
    ambertools_bin: Path,
) -> tuple[Path, Path]:
    charge_int = int(round(ligand_charge))
    if abs(charge_int - ligand_charge) > 1e-6:
        raise SystemExit("antechamber -nc requires an integer formal ligand charge")
    env = os.environ.copy()
    env["PATH"] = f"{ambertools_bin}{os.pathsep}{env.get('PATH', '')}"
    env.setdefault("AMBERHOME", str(ambertools_bin.parent))
    mol2_in = input_dir / f"{ligand}.mol2"
    mol2_out = input_dir / f"{ligand}_antechamber_nc{charge_int}.mol2"
    frcmod_out = input_dir / f"{ligand}_antechamber_nc{charge_int}.frcmod"
    antechamber_log = input_dir.parent / f"antechamber_{ligand}_nc{charge_int}.log"
    parmchk_log = input_dir.parent / f"parmchk2_{ligand}_nc{charge_int}.log"
    with antechamber_log.open("w", encoding="utf-8", newline="\n") as log:
        command = [
            str(ambertools_bin / "antechamber"),
            "-i",
            str(mol2_in),
            "-fi",
            "mol2",
            "-o",
            str(mol2_out),
            "-fo",
            "mol2",
            "-rn",
            ligand,
            "-c",
            "bcc",
            "-s",
            "2",
            "-at",
            "gaff2",
            "-nc",
            str(charge_int),
            "-m",
            str(ligand_mult),
        ]
        if ligand_mult != 1:
            command.extend(["-ek", f"spin={ligand_mult}"])
        subprocess.run(
            command,
            cwd=input_dir,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    with parmchk_log.open("w", encoding="utf-8", newline="\n") as log:
        subprocess.run(
            [
                str(ambertools_bin / "parmchk2"),
                "-i",
                str(mol2_out),
                "-f",
                "mol2",
                "-o",
                str(frcmod_out),
                "-s",
                "gaff2",
            ],
            cwd=input_dir,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    shutil.copy2(mol2_out, input_dir / f"{ligand}.mol2")
    shutil.copy2(frcmod_out, input_dir / f"{ligand}.frcmod")
    shutil.copy2(input_dir / f"{ligand}.mol2", mcpb_dir / f"{ligand}.mol2")
    shutil.copy2(input_dir / f"{ligand}.frcmod", mcpb_dir / f"{ligand}.frcmod")
    return mol2_out, frcmod_out


def copy_variant_inputs(
    source_dir: Path,
    outdir: Path,
    *,
    group: str,
    ligand: str,
    ligand_charge: float,
    ligand_mult: int,
    charge_atoms: set[str],
    model_charge: int,
    mult: int,
    cores: int,
    gamess: str,
    rebuild_ligand_with_antechamber: bool,
    ambertools_bin: Path,
) -> None:
    src_mcpb = source_dir / "mcpb"
    src_input = source_dir / "input"
    if not src_mcpb.exists():
        raise FileNotFoundError(src_mcpb)
    mcpb_dir = outdir / "mcpb"
    input_dir = outdir / "input"
    mcpb_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)

    for name in ["mcpb_original.pdb", "FE.mol2"]:
        shutil.copy2(src_mcpb / name, mcpb_dir / name)
    for name in ["ACT.mol2", "ACT.frcmod", "UNK.frcmod"]:
        src = src_input / name if (src_input / name).exists() else src_mcpb / name
        shutil.copy2(src, input_dir / name)
        shutil.copy2(src, mcpb_dir / name)
    old, new, adjusted = write_charge_adjusted_mol2(
        src_input / f"{ligand}.mol2",
        input_dir / f"{ligand}.mol2",
        target_charge=ligand_charge,
        charge_atoms=charge_atoms,
    )
    shutil.copy2(input_dir / f"{ligand}.mol2", mcpb_dir / f"{ligand}.mol2")
    rebuilt: tuple[Path, Path] | None = None
    if rebuild_ligand_with_antechamber:
        rebuilt = run_antechamber_ligand(
            input_dir,
            mcpb_dir,
            ligand=ligand,
            ligand_charge=ligand_charge,
            ligand_mult=ligand_mult,
            ambertools_bin=ambertools_bin,
        )
    pdb_candidates = sorted(src_input.glob("*.pdb"))
    for pdb in pdb_candidates:
        shutil.copy2(pdb, input_dir / pdb.name)
    create_mcpb_in(src_mcpb / "mcpb.in", mcpb_dir / "mcpb.in", group=group, charge=model_charge, mult=mult)
    write_scripts(
        mcpb_dir,
        group=group,
        model_charge=model_charge,
        mult=mult,
        cores=cores,
        gamess=gamess,
    )
    (outdir / "README_charge_variant.md").write_text(
        f"""# MCPB Charge Variant

Source: `{source_dir.as_posix()}`

Variant:

- ligand `{ligand}` mol2 target charge: `{ligand_charge}`
- ligand antechamber multiplicity: `{ligand_mult}`
- original `{ligand}` charge sum: `{old:.6f}`
- new `{ligand}` charge sum: `{new:.6f}`
- adjusted atoms: `{', '.join(adjusted)}`
- antechamber rebuild: `{rebuilt[0].name if rebuilt else 'no'}`
- MCPB small/large charge: `{model_charge}`
- multiplicity: `{mult}`

This variant changes the ligand template charge and the MCPB/GAMESS total
charge. If `antechamber rebuild` is not `no`, the final ligand mol2/frcmod in
`input/` and `mcpb/` were regenerated with AmberTools using the requested
formal charge.
""",
        encoding="utf-8",
        newline="\n",
    )


def run_step1_and_prepare_energy(outdir: Path, *, group: str, model_charge: int, mult: int) -> tuple[int, str]:
    mcpb_dir = outdir / "mcpb"
    subprocess.run(["bash", str(mcpb_dir / "run_01_mcpb_step1.sh")], cwd=mcpb_dir, check=True)
    small = mcpb_dir / f"{group}_small.pdb"
    if not small.exists():
        raise FileNotFoundError(small)
    ne = electron_count(small, model_charge)
    compatible = (ne % 2 == 0 and mult % 2 == 1) or (ne % 2 == 1 and mult % 2 == 0)
    if not compatible:
        raise SystemExit(
            f"Electron/multiplicity parity mismatch: electrons={ne}, mult={mult}. "
            "Even electrons need odd multiplicity; odd electrons need even multiplicity."
        )
    energy = write_energy_input(mcpb_dir, group=group, charge=model_charge, mult=mult)
    return ne, str(energy)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--ligand", default="UNK")
    parser.add_argument("--ligand-charge", type=float, default=-1.0)
    parser.add_argument("--ligand-mult", type=int, default=1)
    parser.add_argument("--charge-atoms", default="")
    parser.add_argument("--model-charge", type=int, required=True)
    parser.add_argument("--mult", type=int, required=True)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--gamess", default="/home/qin/softwares/gamess/rungms")
    parser.add_argument("--rebuild-ligand-with-antechamber", action="store_true")
    parser.add_argument("--ambertools-bin", type=Path, default=Path("/mnt/l/WSL/conda_envs/AmberTools25/bin"))
    parser.add_argument("--run-step1", action="store_true")
    parser.add_argument("--prepare-energy-only", action="store_true")
    args = parser.parse_args()

    if args.prepare_energy_only:
        energy = write_energy_input(args.outdir / "mcpb", group=args.group, charge=args.model_charge, mult=args.mult)
        write_scripts(
            args.outdir / "mcpb",
            group=args.group,
            model_charge=args.model_charge,
            mult=args.mult,
            cores=args.cores,
            gamess=args.gamess,
        )
        print(f"energy_input={energy}")
        print(f"run={args.outdir / 'mcpb' / ('run_' + args.group + '_scf_energy_soscf_finegrid_tmux.sh')}")
        return 0

    charge_atoms = {atom.strip() for atom in args.charge_atoms.split(",") if atom.strip()}
    copy_variant_inputs(
        args.source_dir,
        args.outdir,
        group=args.group,
        ligand=args.ligand,
        ligand_charge=args.ligand_charge,
        ligand_mult=args.ligand_mult,
        charge_atoms=charge_atoms,
        model_charge=args.model_charge,
        mult=args.mult,
        cores=args.cores,
        gamess=args.gamess,
        rebuild_ligand_with_antechamber=args.rebuild_ligand_with_antechamber,
        ambertools_bin=args.ambertools_bin,
    )
    print(f"outdir={args.outdir}")
    print(f"mcpb_dir={args.outdir / 'mcpb'}")
    if args.run_step1:
        electrons, energy_input = run_step1_and_prepare_energy(
            args.outdir, group=args.group, model_charge=args.model_charge, mult=args.mult
        )
        print(f"small_electrons={electrons}")
        print(f"energy_input={energy_input}")
        print(f"run={args.outdir / 'mcpb' / ('run_' + args.group + '_scf_energy_soscf_finegrid_tmux.sh')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
