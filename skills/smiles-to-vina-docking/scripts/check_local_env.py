#!/usr/bin/env python
"""Check whether the current local environment can run the docking workflow."""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path


def module_status(module: str) -> tuple[bool, str]:
    try:
        spec = importlib.util.find_spec(module)
    except Exception as exc:
        return False, str(exc)
    if spec is None:
        return False, "not found"
    return True, "ok"


def command_status(command: str) -> tuple[bool, str]:
    expanded = Path(command).expanduser()
    if expanded.exists() and not expanded.is_dir():
        return True, str(expanded)
    path = shutil.which(command)
    if path:
        return True, path
    return False, "not found on PATH"


def run_version(command: str, args: list[str]) -> str:
    expanded = Path(command).expanduser()
    path = str(expanded) if expanded.exists() and not expanded.is_dir() else shutil.which(command)
    if not path:
        return ""
    try:
        proc = subprocess.run([path, *args], check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=8)
    except Exception:
        return ""
    text = (proc.stdout or "").strip().splitlines()
    return text[0] if text else ""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--docking-engine", default="vina", help="Docking command name or executable path")
    parser.add_argument("--strict", action="store_true", help="Exit nonzero if a required module/command is missing")
    args = parser.parse_args()

    print(f"python={sys.executable}")
    print(f"python_version={sys.version.split()[0]}")
    required_modules = ["rdkit", "meeko"]
    optional_modules = ["numpy", "sklearn", "pymol"]
    missing_required = []
    for module in required_modules:
        ok, detail = module_status(module)
        print(f"module:{module}={detail}")
        if not ok:
            missing_required.append(f"python module {module}")
    for module in optional_modules:
        ok, detail = module_status(module)
        label = "ok" if ok else detail
        print(f"module:{module}={label}")

    required_commands = [args.docking_engine]
    optional_commands = ["obabel", "mk_prepare_ligand.py", "pymol"]
    missing_commands = []
    for command in required_commands:
        ok, detail = command_status(command)
        print(f"command:{command}={detail}")
        version = run_version(command, ["--version"])
        if version:
            print(f"version:{command}={version}")
        if not ok:
            missing_commands.append(f"command {command}")
    for command in optional_commands:
        ok, detail = command_status(command)
        print(f"command:{command}={detail if ok else 'optional_missing'}")

    if missing_required or missing_commands:
        print("status=missing_required " + "; ".join([*missing_required, *missing_commands]))
        if args.strict:
            raise SystemExit(1)
    else:
        print("status=ok")


if __name__ == "__main__":
    main()
