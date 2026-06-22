#!/usr/bin/env python3
"""Write a simple GAMESS runner for MCPB.py-generated .inp files."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mcpb-dir", type=Path, required=True)
    parser.add_argument("--gamess-dir", type=Path, default=Path("/home/qin/softwares/gamess"))
    parser.add_argument("--version", default="00")
    parser.add_argument("--cores", type=int, default=16)
    args = parser.parse_args()

    mcpb = args.mcpb_dir
    inputs = [
        "FE_A508_6707_small_opt.inp",
        "FE_A508_6707_small_fc.inp",
        "FE_A508_6707_large_mk.inp",
    ]
    missing = [name for name in inputs if not (mcpb / name).exists()]
    if missing:
        raise FileNotFoundError(", ".join(missing))
    script = mcpb / "run_gamess_mcpb_qm.sh"
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        'ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        'cd "$ROOT"',
        "mkdir -p gamess_logs",
        f"GAMESS_DIR={args.gamess_dir}",
        f"VER={args.version}",
        f"CORES={args.cores}",
        "",
        "run_job() {",
        "  local inp=\"$1\"",
        "  local stem=\"${inp%.inp}\"",
        "  echo \"[GAMESS] $inp\"",
        "  \"$GAMESS_DIR/rungms\" \"$stem\" \"$VER\" \"$CORES\" > \"gamess_logs/${stem}.log\"",
        "  grep -qi 'terminated normally' \"gamess_logs/${stem}.log\"",
        "}",
        "",
    ]
    for name in inputs:
        lines.append(f"run_job {name!r}")
    lines.extend(
        [
            "",
            "echo all_gamess_jobs_finished",
            "echo next='MCPB.py -i mcpb.in -s 2'",
        ]
    )
    script.write_text("\n".join(lines) + "\n", encoding="utf-8")
    script.chmod(0o755)
    print(f"runner={script}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
