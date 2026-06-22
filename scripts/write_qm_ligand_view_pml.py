from __future__ import annotations

import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: write_qm_ligand_view_pml.py qm_dir")
        return 2
    qm_dir = Path(sys.argv[1])
    opt_dir = qm_dir / "optimized"
    pml_wsl = opt_dir / "view_qm_optimized_ligands_wsl.pml"
    pml_windows = opt_dir / "view_qm_optimized_ligands_windows.pml"
    act = opt_dir / "ACT_current_acetate_opt.pdb"
    unl = opt_dir / "UNL_current_radical_opt.pdb"

    def write_pml(path: Path, act_path: str, unl_path: str) -> None:
        path.write_text(
            "\n".join(
                [
                    "delete all",
                    f"load {act_path}, ACT_qm_opt",
                    f"load {unl_path}, UNL_qm_opt",
                    "hide everything",
                    "show sticks, ACT_qm_opt or UNL_qm_opt",
                    "color palegreen, elem C",
                    "color red, elem O",
                    "color blue, elem N",
                    "color white, elem H",
                    "set stick_radius, 0.18",
                    "set sphere_scale, 0.20",
                    "show spheres, elem N+O",
                    "orient",
                    "zoom all, 4",
                    "",
                ]
            ),
            encoding="utf-8",
        )

    write_pml(pml_wsl, str(act), str(unl))
    write_pml(
        pml_windows,
        str(act).replace("/mnt/e/", "E:/"),
        str(unl).replace("/mnt/e/", "E:/"),
    )
    # Backward-compatible default: WSL path, because the user's active PyMOL is WSL PyMOL.
    (opt_dir / "view_qm_optimized_ligands.pml").write_text(
        "\n".join(
            [
                "# WSL PyMOL version. For Windows PyMOL use view_qm_optimized_ligands_windows.pml",
                f"run {pml_wsl}",
            ]
        ),
        encoding="utf-8",
    )
    print(f"pml_wsl={pml_wsl}")
    print(f"pml_windows={pml_windows}")
    print(f"pml_default={opt_dir / 'view_qm_optimized_ligands.pml'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
