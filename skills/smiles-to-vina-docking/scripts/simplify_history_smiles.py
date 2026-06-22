#!/usr/bin/env python
"""Rewrite dock_history.csv with one compact `smiles` column.

Old histories with `input_smiles` and `canonical_smiles` are read safely, backed
up, and rewritten using the current compact ledger schema.
"""

from __future__ import annotations

import argparse
import shutil
from datetime import datetime
from pathlib import Path

import dock_utils as dl


def history_path(args: argparse.Namespace) -> Path:
    if args.history_csv:
        return args.history_csv.expanduser().resolve()
    history_dir = dl.normalize_ledger_dir(args.history_dir)
    return dl.ledger_path(history_dir)


def migrate(path: Path, backup: bool = True) -> Path | None:
    if not path.exists():
        raise FileNotFoundError(path)
    rows = dl.read_csv(path)
    backup_path: Path | None = None
    if backup:
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        backup_path = path.with_name(f"{path.name}.bak-{stamp}")
        shutil.copy2(path, backup_path)
    dl.write_csv(path, rows)
    return backup_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--history-dir", "--ledger-dir", dest="history_dir", type=Path, default=dl.DEFAULT_LEDGER_DIR)
    parser.add_argument("--history-csv", type=Path)
    parser.add_argument("--no-backup", action="store_true")
    args = parser.parse_args()
    path = history_path(args)
    backup_path = migrate(path, backup=not args.no_backup)
    backup_text = f" backup={backup_path}" if backup_path else ""
    print(f"rewrote={path} rows={len(dl.read_csv(path))}{backup_text}", flush=True)


if __name__ == "__main__":
    main()
