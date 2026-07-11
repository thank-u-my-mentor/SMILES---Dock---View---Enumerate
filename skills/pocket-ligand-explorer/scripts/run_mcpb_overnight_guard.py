#!/usr/bin/env python3
"""Watch a ZFY OPT job, then cautiously advance a TJ MCPB QM chain.

This is intentionally conservative.  It can run several fixed-geometry SCF
continuations, try a fallback multiplicity, and proceed through GRADIENT,
OPTIMIZE, and HESSIAN only when the previous stage is chemically usable.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
RESTART_DIR = Path("/home/qin/softwares/gamess/restart")
GAMESS = Path("/home/qin/softwares/gamess/rungms")


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def append_report(path: Path, line: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(f"- `{stamp}` {line}\n")


def has_any(text: str, needles: tuple[str, ...] | list[str]) -> bool:
    return any(needle in text for needle in needles)


def is_finished(log: Path, status: Path) -> bool:
    text = read_text(log) + "\n" + read_text(status)
    return has_any(
        text,
        (
            "EXECUTION OF GAMESS TERMINATED NORMALLY",
            "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-",
            "exit_code=",
            "status=normal",
            "status=abnormal",
            "status=not_normal",
        ),
    )


def is_normal(log: Path) -> bool:
    return "EXECUTION OF GAMESS TERMINATED NORMALLY" in read_text(log)


def scf_failed(log: Path) -> bool:
    text = read_text(log)
    return has_any(text, ("SCF IS UNCONVERGED", "SCF HAS NOT CONVERGED", "NO GRADIENT, SCF DID NOT CONVERGE"))


def scf_succeeded(log: Path) -> bool:
    text = read_text(log)
    return is_normal(log) and "DENSITY CONVERGED" in text and not scf_failed(log)


def gradient_succeeded(log: Path) -> bool:
    text = read_text(log)
    return is_normal(log) and not scf_failed(log) and "GRADIENT" in text


def opt_usable(log: Path) -> bool:
    text = read_text(log)
    if not is_normal(log) or scf_failed(log):
        return False
    bad = (
        "THE GEOMETRY SEARCH IS NOT CONVERGED",
        "FAILURE TO LOCATE STATIONARY POINT",
        "NOT CONVERGED",
    )
    if has_any(text, bad):
        return False
    return "END OF GEOMETRY SEARCH" in text


def summarize_log(log: Path) -> str:
    text = read_text(log)
    last_nserch = ""
    for line in text.splitlines():
        if "NSERCH:" in line and "GRAD. MAX=" in line:
            last_nserch = line.strip()
    energy = ""
    ssq = ""
    for line in reversed(text.splitlines()):
        if not ssq and "S-SQUARED" in line and "=" in line:
            ssq = line.strip()
        if not energy and "FINAL U-B3LYP ENERGY" in line:
            energy = line.strip()
        if ssq and energy:
            break
    markers: list[str] = []
    if is_normal(log):
        markers.append("normal")
    if scf_succeeded(log):
        markers.append("density_converged")
    if scf_failed(log):
        markers.append("scf_failed")
    if "THE GEOMETRY SEARCH IS NOT CONVERGED" in text:
        markers.append("geometry_not_converged")
    if "END OF GEOMETRY SEARCH" in text:
        markers.append("end_of_geometry_search")
    return "; ".join(part for part in [", ".join(markers), last_nserch, energy, ssq] if part)


def run_cmd(cmd: list[str], cwd: Path, report: Path, label: str) -> int:
    append_report(report, f"RUN {label}: `{' '.join(cmd)}` in `{cwd}`")
    proc = subprocess.run(cmd, cwd=str(cwd))
    append_report(report, f"DONE {label}: exit_code={proc.returncode}")
    return proc.returncode


def wait_for_job(log: Path, status: Path, report: Path, label: str, poll_seconds: int, timeout_hours: float) -> bool:
    append_report(report, f"WAIT {label}: log=`{log}`")
    deadline = time.time() + timeout_hours * 3600
    while time.time() < deadline:
        if is_finished(log, status):
            append_report(report, f"FINISHED {label}: {summarize_log(log)}")
            return True
        time.sleep(poll_seconds)
    append_report(report, f"TIMEOUT {label}: no termination marker after {timeout_hours:g} h")
    return False


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


def extract_vec(dat: Path) -> list[str]:
    lines = read_text(dat).splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip().upper() == "$VEC":
            start = i
    if start is None:
        raise RuntimeError(f"No $VEC in {dat}")
    end = None
    for i in range(start + 1, len(lines)):
        if lines[i].strip().upper() == "$END":
            end = i
            break
    if end is None:
        raise RuntimeError(f"No $END after $VEC in {dat}")
    return lines[start : end + 1]


def write_foreground_runner(root: Path, job: str, cores: int) -> None:
    script = root / f"run_{job}_foreground.sh"
    write_text(
        script,
        f"""#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "${{BASH_SOURCE[0]}}")"
JOB={job}
NCORES=${{NCORES:-{cores}}}
GAMESS=${{GAMESS:-{GAMESS}}}
VERSION=${{VERSION:-00}}
mkdir -p gamess_logs logs
echo "[GAMESS] $JOB start $(date) NCORES=$NCORES" | tee "gamess_logs/${{JOB}}.status"
"$GAMESS" "$JOB" "$VERSION" "$NCORES" > "gamess_logs/${{JOB}}.log" 2>&1
rc=$?
echo "[GAMESS] $JOB exit_code=$rc end $(date)" | tee -a "gamess_logs/${{JOB}}.status"
if grep -q 'DENSITY CONVERGED' "gamess_logs/${{JOB}}.log"; then echo "[GAMESS] $JOB scf=density_converged" | tee -a "gamess_logs/${{JOB}}.status"; fi
if grep -q 'SCF IS UNCONVERGED\\|SCF HAS NOT CONVERGED' "gamess_logs/${{JOB}}.log"; then echo "[GAMESS] $JOB scf=unconverged" | tee -a "gamess_logs/${{JOB}}.status"; fi
if grep -q 'EXECUTION OF GAMESS TERMINATED NORMALLY' "gamess_logs/${{JOB}}.log"; then echo "[GAMESS] $JOB status=normal" | tee -a "gamess_logs/${{JOB}}.status"; else echo "[GAMESS] $JOB status=not_normal" | tee -a "gamess_logs/${{JOB}}.status"; fi
exit "$rc"
""",
    )
    script.chmod(0o755)


def run_gamess_job(root: Path, job: str, cores: int, report: Path) -> Path:
    runner = root / f"run_{job}_foreground.sh"
    if not runner.exists():
        write_foreground_runner(root, job, cores)
        runner = root / f"run_{job}_foreground.sh"
    run_cmd(["env", f"NCORES={cores}", "bash", str(runner.name)], root, report, job)
    return root / "gamess_logs" / f"{job}.log"


def write_energy_huckel(root: Path, source_job: str, out_job: str, charge: int, mult: int, cores: int) -> None:
    source = root / f"{source_job}.inp"
    lines = remove_vec(read_text(source).splitlines())
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
    lines = replace_card(lines, "$GUESS", " $GUESS GUESS=HUCKEL $END")
    lines = replace_card(
        lines,
        "$DFT",
        (
            " $DFT\n"
            "  NRAD=96 NLEB=302 NRAD0=96 NLEB0=302 SWOFF=0.0\n"
            " $END"
        ),
    )
    write_text(root / f"{out_job}.inp", "\n".join(lines) + "\n")
    write_foreground_runner(root, out_job, cores)


def prepare_moread_energy(root: Path, source_job: str, out_job: str, dat: Path, charge: int, mult: int, norb: int, cores: int, report: Path) -> None:
    helper = SCRIPT_DIR / "prepare_gamess_moread_energy.py"
    cmd = [
        sys.executable,
        str(helper),
        "--mcpb-dir",
        str(root),
        "--source-job",
        source_job,
        "--out-job",
        out_job,
        "--dat",
        str(dat),
        "--charge",
        str(charge),
        "--mult",
        str(mult),
        "--norb",
        str(norb),
        "--cores",
        str(cores),
        "--session",
        out_job.lower(),
        "--dft-no-coarse-grid",
    ]
    rc = run_cmd(cmd, root, report, f"prepare {out_job}")
    if rc != 0:
        raise RuntimeError(f"prepare_moread_energy failed for {out_job}")


def prepare_gradient(root: Path, source_job: str, out_job: str, dat: Path, charge: int, mult: int, norb: int, cores: int, report: Path) -> None:
    helper = SCRIPT_DIR / "prepare_gamess_moread_gradient.py"
    cmd = [
        sys.executable,
        str(helper),
        "--mcpb-dir",
        str(root),
        "--source-job",
        source_job,
        "--out-job",
        out_job,
        "--dat",
        str(dat),
        "--charge",
        str(charge),
        "--mult",
        str(mult),
        "--norb",
        str(norb),
        "--cores",
        str(cores),
        "--session",
        out_job.lower(),
    ]
    rc = run_cmd(cmd, root, report, f"prepare {out_job}")
    if rc != 0:
        raise RuntimeError(f"prepare_gradient failed for {out_job}")


def prepare_opt(root: Path, source_job: str, out_job: str, dat: Path, charge: int, mult: int, norb: int, nstep: int, cores: int, report: Path) -> None:
    helper = SCRIPT_DIR / "prepare_gamess_moread_optimize.py"
    cmd = [
        sys.executable,
        str(helper),
        "--mcpb-dir",
        str(root),
        "--source-job",
        source_job,
        "--out-job",
        out_job,
        "--dat",
        str(dat),
        "--charge",
        str(charge),
        "--mult",
        str(mult),
        "--norb",
        str(norb),
        "--nstep",
        str(nstep),
        "--opttol",
        "0.0002",
        "--cores",
        str(cores),
        "--session",
        out_job.lower(),
    ]
    rc = run_cmd(cmd, root, report, f"prepare {out_job}")
    if rc != 0:
        raise RuntimeError(f"prepare_opt failed for {out_job}")


def extract_last_coords(log: Path) -> list[str]:
    lines = read_text(log).splitlines()
    indices = [i for i, line in enumerate(lines) if "COORDINATES OF ALL ATOMS ARE (ANGS)" in line]
    if not indices:
        raise RuntimeError(f"No coordinate block in {log}")
    coords: list[str] = []
    for line in lines[indices[-1] + 1 :]:
        stripped = line.strip()
        if not stripped:
            if coords:
                break
            continue
        parts = stripped.split()
        if len(parts) == 5 and parts[0][0].isalpha():
            coords.append(f"{parts[0]:<2s} {float(parts[1]):8.1f} {float(parts[2]):14.6f} {float(parts[3]):14.6f} {float(parts[4]):14.6f}")
        elif coords:
            break
    if not coords:
        raise RuntimeError(f"Could not parse coordinates from {log}")
    return coords


def fill_hessian_input(root: Path, template_job: str, out_job: str, opt_log: Path, charge: int, mult: int, cores: int) -> None:
    template = root / f"{template_job}.inp"
    coords = extract_last_coords(opt_log)
    lines = remove_vec(read_text(template).splitlines())
    lines = replace_card(lines, "$SYSTEM", " $SYSTEM MEMDDI=400 MWORDS=200 $END")
    lines = replace_card(
        lines,
        "$CONTRL",
        (
            " $CONTRL\n"
            "  SCFTYP=UHF DFTTYP=B3LYP RUNTYP=HESSIAN\n"
            f"  ICHARG={charge} MULT={mult} MAXIT=200 COORD=CART UNITS=ANGS\n"
            " $END"
        ),
    )
    lines = replace_card(
        lines,
        "$SCF",
        " $SCF DIRSCF=.T. DIIS=.T. DAMP=.T. SHIFT=.T. ETHRSH=2.0 MAXDII=20 $END",
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
    data_idx = next((i for i, line in enumerate(lines) if line.strip().upper() == "$DATA"), None)
    if data_idx is None:
        raise RuntimeError(f"No $DATA in {template}")
    end_idx = next((i for i in range(data_idx + 1, len(lines)) if lines[i].strip().upper() == "$END"), None)
    if end_idx is None:
        raise RuntimeError(f"No $END after $DATA in {template}")
    new_lines = lines[: data_idx + 3] + coords + [" $END"] + lines[end_idx + 1 :]
    write_text(root / f"{out_job}.inp", "\n".join(new_lines) + "\n")
    write_foreground_runner(root, out_job, cores)


def run_scf_branch(
    *,
    root: Path,
    branch: str,
    source_job: str,
    initial_job: str,
    charge: int,
    mult: int,
    norb: int,
    cores: int,
    max_moread: int,
    report: Path,
    initial_prepared: bool,
) -> str | None:
    if not initial_prepared:
        write_energy_huckel(root, source_job, initial_job, charge, mult, cores)
    job = initial_job
    for attempt in range(max_moread + 1):
        log = run_gamess_job(root, job, cores, report)
        append_report(report, f"{branch} SCF attempt {attempt}: {summarize_log(log)}")
        if scf_succeeded(log):
            append_report(report, f"{branch} SCF accepted: `{job}`")
            return job
        dat = RESTART_DIR / f"{job}.dat"
        if attempt >= max_moread:
            break
        if not dat.exists() or "$VEC" not in read_text(dat):
            append_report(report, f"{branch} SCF cannot continue: missing `$VEC` in `{dat}`")
            break
        next_job = f"{initial_job}_moread{attempt + 1}"
        prepare_moread_energy(root, job, next_job, dat, charge, mult, norb, cores, report)
        job = next_job
    append_report(report, f"{branch} SCF failed after {max_moread + 1} attempt(s)")
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--zfy-dir", type=Path, required=True)
    parser.add_argument("--zfy-job", required=True)
    parser.add_argument("--tj-dir", type=Path, required=True)
    parser.add_argument("--tj-source-job", default="TJ_HID_M7_FE_small_opt")
    parser.add_argument("--tj-prepared-m7-job", default="TJ_HID_M7_FE_scf_energy_moread1_finegrid")
    parser.add_argument("--tj-prepared-m7-is-ready", action="store_true")
    parser.add_argument("--charge", type=int, default=1)
    parser.add_argument("--primary-mult", type=int, default=7)
    parser.add_argument("--fallback-mult", type=int, default=5)
    parser.add_argument("--norb", type=int, default=658)
    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--poll-seconds", type=int, default=180)
    parser.add_argument("--zfy-timeout-hours", type=float, default=8.0)
    parser.add_argument("--tj-scf-moread-retries", type=int, default=2)
    parser.add_argument("--tj-opt-nstep", type=int, default=80)
    parser.add_argument("--skip-fallback-mult", action="store_true")
    args = parser.parse_args()

    report = args.tj_dir / "logs" / "overnight_zfy_then_tj_guard.md"
    write_text(
        report,
        "# Overnight ZFY -> TJ MCPB Guard\n\n"
        f"- Started: `{time.strftime('%Y-%m-%d %H:%M:%S')}`\n"
        f"- ZFY job: `{args.zfy_job}`\n"
        f"- TJ dir: `{args.tj_dir}`\n"
        f"- Primary MULT: `{args.primary_mult}`\n"
        f"- Fallback MULT: `{args.fallback_mult}`\n\n",
    )

    zfy_log = args.zfy_dir / "logs" / f"{args.zfy_job}.log"
    zfy_status = args.zfy_dir / "logs" / f"{args.zfy_job}.status"
    if not wait_for_job(zfy_log, zfy_status, report, "ZFY OPT", args.poll_seconds, args.zfy_timeout_hours):
        append_report(report, "STOP: ZFY OPT did not finish inside the guard window; TJ was not started.")
        return 1
    append_report(report, "ZFY Hessian is intentionally not launched by this guard.")

    append_report(report, "Starting TJ SCF/GRAD/OPT/Hessian chain.")
    m7_job = run_scf_branch(
        root=args.tj_dir,
        branch=f"M{args.primary_mult}",
        source_job=args.tj_source_job,
        initial_job=args.tj_prepared_m7_job,
        charge=args.charge,
        mult=args.primary_mult,
        norb=args.norb,
        cores=args.cores,
        max_moread=args.tj_scf_moread_retries,
        report=report,
        initial_prepared=args.tj_prepared_m7_is_ready,
    )
    chosen_job = m7_job
    chosen_mult = args.primary_mult

    if chosen_job is None and not args.skip_fallback_mult:
        m5_initial = f"TJ_HID_M{args.fallback_mult}_FE_scf_energy_soscf_finegrid"
        chosen_job = run_scf_branch(
            root=args.tj_dir,
            branch=f"M{args.fallback_mult}",
            source_job=args.tj_source_job,
            initial_job=m5_initial,
            charge=args.charge,
            mult=args.fallback_mult,
            norb=args.norb,
            cores=args.cores,
            max_moread=args.tj_scf_moread_retries,
            report=report,
            initial_prepared=False,
        )
        chosen_mult = args.fallback_mult

    if chosen_job is None:
        append_report(report, "STOP: no SCF branch produced `DENSITY CONVERGED`.")
        return 2

    scf_dat = RESTART_DIR / f"{chosen_job}.dat"
    grad_job = f"{chosen_job}_gradient"
    prepare_gradient(args.tj_dir, chosen_job, grad_job, scf_dat, args.charge, chosen_mult, args.norb, args.cores, report)
    grad_log = run_gamess_job(args.tj_dir, grad_job, args.cores, report)
    append_report(report, f"GRADIENT summary: {summarize_log(grad_log)}")
    if not gradient_succeeded(grad_log):
        append_report(report, "STOP: gradient did not produce a usable result.")
        return 3

    grad_dat = RESTART_DIR / f"{grad_job}.dat"
    opt_job = f"TJ_HID_M{chosen_mult}_FE_small_opt_guard_nstep{args.tj_opt_nstep}"
    prepare_opt(args.tj_dir, grad_job, opt_job, grad_dat, args.charge, chosen_mult, args.norb, args.tj_opt_nstep, args.cores, report)
    opt_log = run_gamess_job(args.tj_dir, opt_job, args.cores, report)
    append_report(report, f"OPT summary: {summarize_log(opt_log)}")
    if not opt_usable(opt_log):
        append_report(report, "STOP: OPT did not reach a Hessian-ready geometry.")
        return 4

    fc_job = f"TJ_HID_M{chosen_mult}_FE_small_fc_guard"
    fill_hessian_input(args.tj_dir, "TJ_HID_M7_FE_small_fc", fc_job, opt_log, args.charge, chosen_mult, args.cores)
    fc_log = run_gamess_job(args.tj_dir, fc_job, args.cores, report)
    append_report(report, f"HESSIAN summary: {summarize_log(fc_log)}")
    if not is_normal(fc_log) or scf_failed(fc_log):
        append_report(report, "STOP: Hessian failed or SCF was not converged.")
        return 5
    append_report(report, f"SUCCESS: TJ Hessian completed for MULT={chosen_mult}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
