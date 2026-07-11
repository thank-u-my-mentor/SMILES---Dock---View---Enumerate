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
        subsets = [
            ("C4_O1_le_3.00A", df[df.C4_O1_A <= 3.0]),
            ("C4_O1_le_2.80A", df[df.C4_O1_A <= 2.8]),
            ("top_100_shortest", df.sort_values("C4_O1_A").head(100)),
        ]
        for subset_name, sub in subsets:
            sub = sub.copy()
            sub["O1_side"] = sub.O1_height_to_C3C4C5_plane_A.apply(
                lambda x: "O1_positive_side" if x > 0 else "O1_negative_side"
            )
            total = len(sub)
            for side in ["O1_positive_side", "O1_negative_side"]:
                side_total = int((sub.O1_side == side).sum())
                for rs in ["R_like_by_product_reference", "S_like_by_product_reference"]:
                    rs_total = int((sub.pseudo_RS_by_product_reference == rs).sum())
                    n = int(((sub.O1_side == side) & (sub.pseudo_RS_by_product_reference == rs)).sum())
                    rows.append(
                        {
                            "dataset": dataset,
                            "subset": subset_name,
                            "n_total": total,
                            "O1_side": side,
                            "pseudo_RS": rs,
                            "n": n,
                            "fraction_of_subset": n / total if total else float("nan"),
                            "fraction_within_O1_side": n / side_total if side_total else float("nan"),
                            "fraction_within_RS": n / rs_total if rs_total else float("nan"),
                            "O1_side_total": side_total,
                            "rs_total": rs_total,
                        }
                    )
    res = pd.DataFrame(rows)
    out_csv = OUT / "O1_side_vs_pseudo_RS_crosstab.csv"
    res.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
