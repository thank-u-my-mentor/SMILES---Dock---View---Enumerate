#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


LEAP_TEMPLATE = """source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams frcmod.ions234lm_126_tip3p
loadamberparams UNK.frcmod
loadamberparams ACT.frcmod
UNK = loadmol2 UNK.mol2
ACT = loadmol2 ACT.mol2
FE = loadmol2 FE.mol2
{hoh_line}
x = loadpdb mcpb_original.pdb
check x
charge x
quit
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Run a ZFY MCPB/tleap preflight.")
    parser.add_argument("--mcpb-dir", type=Path, required=True)
    parser.add_argument("--amber-sh", default="/mnt/l/WSL/conda_envs/AmberTools25/amber.sh")
    parser.add_argument(
        "--tip3p-hoh-template",
        action="store_true",
        help="Use Amber's standard HOH template instead of the MCPB helper HOH.mol2.",
    )
    args = parser.parse_args()

    root = args.mcpb_dir
    logs = root / "logs"
    logs.mkdir(exist_ok=True)
    leap_in = root / "tleap_preflight.in"
    leap_log = logs / "tleap_preflight.log"
    hoh_line = "# HOH uses Amber TIP3P template" if args.tip3p_hoh_template else "HOH = loadmol2 HOH.mol2"
    leap_in.write_text(LEAP_TEMPLATE.format(hoh_line=hoh_line), encoding="utf-8", newline="\n")

    cmd = (
        "set -euo pipefail; "
        f"cd {root}; "
        "set +u; "
        f"source {args.amber_sh}; "
        "set -u; "
        "tleap -f tleap_preflight.in > logs/tleap_preflight.log 2>&1"
    )
    rc = subprocess.run(["bash", "-lc", cmd], check=False).returncode
    print(f"tleap_exit={rc}")
    print(f"tleap_input={leap_in}")
    print(f"tleap_log={leap_log}")
    if leap_log.exists():
        text = leap_log.read_text(encoding="utf-8", errors="replace")
        interesting = [
            line
            for line in text.splitlines()
            if any(
                marker in line
                for marker in (
                    "Warning",
                    "ERROR",
                    "FATAL",
                    "ATOMS NOT BONDED",
                    "Maximum coordination",
                    "Unknown",
                    "missing",
                    "The unperturbed charge",
                    "Checking",
                )
            )
        ]
        print("important_lines:")
        for line in interesting[-120:]:
            print(line)
        print("tail:")
        for line in text.splitlines()[-50:]:
            print(line)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
