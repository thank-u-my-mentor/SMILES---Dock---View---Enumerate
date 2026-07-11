#!/usr/bin/env python3
from pathlib import Path

import pandas as pd


ROOT = Path("/mnt/e/TJ/260706_mcpb_m2_200ns")
OUT = ROOT / "amber_200ns_pmemd_compare/analysis_vs_gromacs_manual_20260710"

PATHS = {
    "gromacs": ROOT
    / "md_200ns/analysis_reactive_C4O1_20260708/methods_code/reproduced_ul1_face_pseudo_rs/ul1_face_and_pseudo_rs_timeseries.csv",
    "amber_pmemd": OUT / "amber_ul1_face_pseudo_rs/ul1_face_and_pseudo_rs_timeseries.csv",
}


def main() -> None:
    rows = []
    for dataset, path in PATHS.items():
        df = pd.read_csv(path)
        subsets = []
        for cutoff in [3.50, 3.00, 2.90, 2.80, 2.75]:
            subsets.append((f"C4_O1_le_{cutoff:.2f}A", df[df.C4_O1_A <= cutoff]))
        subsets.append(("top_100_shortest", df.sort_values("C4_O1_A").head(100)))
        subsets.append(("all_frames", df))

        for name, sub in subsets:
            total = len(sub)
            for face in ["same", "opposite"]:
                face_total = int((sub.O1_H02_face_relation == face).sum())
                for rs in ["R_like_by_product_reference", "S_like_by_product_reference"]:
                    rs_total = int((sub.pseudo_RS_by_product_reference == rs).sum())
                    n = int(
                        (
                            (sub.O1_H02_face_relation == face)
                            & (sub.pseudo_RS_by_product_reference == rs)
                        ).sum()
                    )
                    rows.append(
                        {
                            "dataset": dataset,
                            "subset": name,
                            "n_total": total,
                            "face": face,
                            "pseudo_RS": rs,
                            "n": n,
                            "fraction_of_subset": n / total if total else float("nan"),
                            "fraction_within_face": n / face_total if face_total else float("nan"),
                            "fraction_within_RS": n / rs_total if rs_total else float("nan"),
                            "face_total": face_total,
                            "rs_total": rs_total,
                        }
                    )

    res = pd.DataFrame(rows)
    out_csv = OUT / "face_vs_pseudo_RS_crosstab.csv"
    res.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")
    print(
        res[res.subset.isin(["C4_O1_le_3.00A", "C4_O1_le_2.80A", "top_100_shortest"])].to_string(
            index=False
        )
    )


if __name__ == "__main__":
    main()
