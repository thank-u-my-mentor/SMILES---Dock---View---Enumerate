#!/usr/bin/env python3
"""Run a small-model M2 RESP fit and audit mapping back to MCPB standard atoms."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path

from pymsmt.lib import get_lib_dict
from pymsmt.mcpb.resp_fitting import gene_resp_input_file, read_resp_file
from pymsmt.mol.gmsio import get_esp_from_gms


def pdb_key(line: str) -> str:
    atom = line[12:16].strip()
    res = line[17:20].strip()
    resid = int(line[22:26])
    return f"{resid}-{res}-{atom}"


def read_pdb_keys(path: Path) -> list[str]:
    keys: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            keys.append(pdb_key(line))
    return keys


def read_standard_fp(path: Path) -> tuple[list[str], dict[str, str]]:
    keys: list[str] = []
    atomtypes: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("LINK"):
            continue
        parts = line.split()
        key = parts[0]
        keys.append(key)
        atomtypes[key] = parts[-1]
    return keys, atomtypes


def lib_charge(libdict: dict, key: str) -> float:
    _resid, resname, atom = key.split("-", 2)
    lib_key = f"{resname}-{atom}"
    if lib_key not in libdict:
        raise KeyError(f"No library charge for {lib_key} while filling {key}")
    return float(libdict[lib_key][1])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-pdb", required=True)
    parser.add_argument("--standard-fp", required=True)
    parser.add_argument("--small-mk-log", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--mol2", nargs="*", default=[])
    parser.add_argument("--ion-id", type=int, default=2797)
    parser.add_argument("--charge", type=int, default=1)
    parser.add_argument("--ff", default="ff19SB")
    args = parser.parse_args()

    small_pdb = Path(args.small_pdb).resolve()
    standard_fp = Path(args.standard_fp).resolve()
    small_mk_log = Path(args.small_mk_log).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    local_small_pdb = outdir / "TJ2unlNH_small_for_resp.pdb"
    local_standard_fp = outdir / "TJ2unlNH_standard.fingerprint"
    shutil.copy2(small_pdb, local_small_pdb)
    shutil.copy2(standard_fp, local_standard_fp)

    small_keys = read_pdb_keys(local_small_pdb)
    standard_keys, atomtypes = read_standard_fp(local_standard_fp)
    small_fp = outdir / "TJ2unlNH_small_resp.fingerprint"
    small_fp.write_text("\n".join(small_keys) + "\n")

    mol2s = [str(Path(x).resolve()) for x in args.mol2]
    old_cwd = Path.cwd()
    os.chdir(outdir)
    try:
        gene_resp_input_file(
            str(local_small_pdb),
            [args.ion_id],
            str(local_standard_fp),
            args.ff,
            mol2s,
            1,
            [],
            args.charge,
        )
        # Amber resp is brittle with long absolute paths, so keep the ESP file
        # name local to outdir after chdir.
        esp_path = Path(small_mk_log.name.removesuffix(".log") + ".esp")
        get_esp_from_gms(str(small_mk_log), str(esp_path))

        subprocess.run(
            ["resp", "-O", "-i", "resp1.in", "-o", "resp1.out", "-p", "resp1.pch",
             "-t", "resp1.chg", "-e", str(esp_path), "-s", "resp1_calc.esp"],
            check=True,
        )
        subprocess.run(
            ["resp", "-O", "-i", "resp2.in", "-o", "resp2.out", "-p", "resp2.pch",
             "-q", "resp1.chg", "-t", "resp2.chg", "-e", str(esp_path),
             "-s", "resp2_calc.esp"],
            check=True,
        )
    finally:
        os.chdir(old_cwd)

    chgs = read_resp_file(str(outdir / "resp2.chg"))
    if len(chgs) != len(small_keys):
        raise RuntimeError(f"RESP charges {len(chgs)} != small atoms {len(small_keys)}")

    small_charge = dict(zip(small_keys, chgs))
    libdict, _chargedict = get_lib_dict(args.ff)
    for mol2 in mol2s:
        mol2_libdict, _mol2_chargedict = get_lib_dict(mol2)
        libdict.update(mol2_libdict)

    mapped: list[tuple[str, str, float, str]] = []
    for key in standard_keys:
        if key in small_charge:
            mapped.append((key, "small_M2_RESP", small_charge[key], atomtypes[key]))
        else:
            mapped.append((key, "ff19SB_missing_backbone", lib_charge(libdict, key), atomtypes[key]))

    tsv = outdir / "small_M2_RESP_to_standard_charge_map.tsv"
    with tsv.open("w", newline="") as handle:
        handle.write("key\tsource\tcharge\tatomtype\n")
        for key, source, charge, atomtype in mapped:
            handle.write(f"{key}\t{source}\t{charge:.8f}\t{atomtype}\n")

    report = outdir / "small_M2_RESP_mapping_report.txt"
    overlap = [k for k in standard_keys if k in small_charge]
    missing = [k for k in standard_keys if k not in small_charge]
    extra = [k for k in small_keys if k not in set(standard_keys)]
    mapped_total = sum(x[2] for x in mapped)
    small_total = sum(chgs)
    cap_total = sum(small_charge[k] for k in extra)
    overlap_total = sum(small_charge[k] for k in overlap)
    missing_total = sum(x[2] for x in mapped if x[1] == "ff19SB_missing_backbone")
    correction_to_charge = args.charge - mapped_total

    def resid_of(key: str) -> str:
        resid, resname, _atom = key.split("-", 2)
        return f"{resid}-{resname}"

    corrected: list[tuple[str, str, float, str]] = []
    missing_by_res: dict[str, list[str]] = {}
    extra_by_res: dict[str, list[str]] = {}
    for key in missing:
        missing_by_res.setdefault(resid_of(key), []).append(key)
    for key in extra:
        extra_by_res.setdefault(resid_of(key), []).append(key)

    residue_lines: list[str] = []
    for resid in sorted(set(missing_by_res) | set(extra_by_res)):
        miss_keys = missing_by_res.get(resid, [])
        extra_keys = extra_by_res.get(resid, [])
        miss_ff = sum(lib_charge(libdict, k) for k in miss_keys)
        cap_resp = sum(small_charge[k] for k in extra_keys)
        delta = cap_resp - miss_ff
        per_atom = delta / len(miss_keys) if miss_keys else 0.0
        residue_lines.append(
            f"{resid}: missing_ff={miss_ff:.8f}, cap_RESP={cap_resp:.8f}, "
            f"delta={delta:.8f}, per_missing_atom_shift={per_atom:.8f}"
        )

    for key in standard_keys:
        if key in small_charge:
            corrected.append((key, "small_M2_RESP", small_charge[key], atomtypes[key]))
        else:
            resid = resid_of(key)
            miss_keys = missing_by_res[resid]
            miss_ff = sum(lib_charge(libdict, k) for k in miss_keys)
            cap_resp = sum(small_charge[k] for k in extra_by_res.get(resid, []))
            shift = (cap_resp - miss_ff) / len(miss_keys)
            corrected.append(
                (key, "ff19SB_shifted_to_small_cap_group", lib_charge(libdict, key) + shift, atomtypes[key])
            )

    corrected_tsv = outdir / "small_M2_RESP_to_standard_charge_map_capgroup_corrected.tsv"
    with corrected_tsv.open("w", newline="") as handle:
        handle.write("key\tsource\tcharge\tatomtype\n")
        for key, source, charge, atomtype in corrected:
            handle.write(f"{key}\t{source}\t{charge:.8f}\t{atomtype}\n")

    corrected_total = sum(x[2] for x in corrected)

    global_shift = correction_to_charge / len(overlap)
    global_corrected: list[tuple[str, str, float, str]] = []
    for key in standard_keys:
        if key in small_charge:
            global_corrected.append(
                (key, "small_M2_RESP_plus_global_overlap_shift", small_charge[key] + global_shift, atomtypes[key])
            )
        else:
            global_corrected.append((key, "ff19SB_missing_backbone", lib_charge(libdict, key), atomtypes[key]))
    global_tsv = outdir / "small_M2_RESP_to_standard_charge_map_global_overlap_corrected.tsv"
    with global_tsv.open("w", newline="") as handle:
        handle.write("key\tsource\tcharge\tatomtype\n")
        for key, source, charge, atomtype in global_corrected:
            handle.write(f"{key}\t{source}\t{charge:.8f}\t{atomtype}\n")
    global_total = sum(x[2] for x in global_corrected)

    report.write_text(
        "\n".join(
            [
                "small-model M2 RESP mapping audit",
                f"small atoms: {len(small_keys)}",
                f"standard atoms: {len(standard_keys)}",
                f"overlap standard<-small: {len(overlap)}",
                f"missing standard filled by ff19SB: {len(missing)}",
                f"extra small cap atoms discarded: {len(extra)}",
                f"RESP small total charge: {small_total:.8f}",
                f"overlap RESP charge sum: {overlap_total:.8f}",
                f"discarded cap charge sum: {cap_total:.8f}",
                f"ff19SB missing-backbone charge sum: {missing_total:.8f}",
                f"mapped standard total before correction: {mapped_total:.8f}",
                f"correction needed to reach charge {args.charge}: {correction_to_charge:.8f}",
                f"cap-group corrected standard total: {corrected_total:.8f}",
                f"global-overlap corrected standard total: {global_total:.8f}",
                f"global shift on each overlapped small-RESP atom: {global_shift:.8f}",
                "",
                "per-residue cap-group correction:",
                *residue_lines,
                "",
                "missing standard atoms:",
                *missing,
                "",
                "extra small cap atoms:",
                *extra,
            ]
        )
        + "\n"
    )
    print(report.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
