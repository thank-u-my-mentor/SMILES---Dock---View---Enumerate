from __future__ import annotations

import csv
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: inspect_pse_ligands.py input.pse outdir")
        return 2

    pse = Path(sys.argv[1])
    outdir = Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)

    import pymol

    pymol.finish_launching(["pymol", "-cq"])
    from pymol import cmd

    cmd.load(str(pse))
    objects = cmd.get_names("objects")
    print("objects=" + ",".join(objects))

    with (outdir / "objects.txt").open("w", encoding="utf-8") as handle:
        for name in objects:
            object_type = cmd.get_type(name)
            atom_count = ""
            if object_type == "object:molecule":
                try:
                    atom_count = str(cmd.count_atoms("%" + name))
                except Exception:
                    atom_count = "count_failed"
            handle.write(f"{name}\t{object_type}\t{atom_count}\n")

    selections = {
        "ACE": "resn ACE",
        "UNL": "resn UNL",
        "LIG": "resn LIG",
        "organic": "organic",
        "hetero": "(organic or inorganic) and not solvent",
    }
    with (outdir / "selection_counts.txt").open("w", encoding="utf-8") as handle:
        for label, selection in selections.items():
            count = cmd.count_atoms(selection)
            handle.write(f"{label}\t{selection}\t{count}\n")
            print(f"{label}={count}")

    for label, selection in selections.items():
        if cmd.count_atoms(selection) == 0:
            continue
        cmd.save(str(outdir / f"{label}.pdb"), selection)
        try:
            cmd.save(str(outdir / f"{label}.sdf"), selection)
        except Exception as exc:
            print(f"warn: could not save {label}.sdf: {exc}")

        rows = []
        model = cmd.get_model(selection, state=1)
        for atom in model.atom:
            rows.append(
                {
                    "index": atom.index,
                    "model": atom.model,
                    "chain": atom.chain,
                    "resn": atom.resn,
                    "resi": atom.resi,
                    "name": atom.name,
                    "elem": atom.symbol,
                    "x": atom.coord[0],
                    "y": atom.coord[1],
                    "z": atom.coord[2],
                    "b": atom.b,
                    "q": atom.q,
                }
            )
        with (outdir / f"{label}_atoms.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    cmd.quit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
