"""Utilities for importing H++ pK output into PLE/Amber workflows."""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path


PKOUT_RE = re.compile(r"^\s*([A-Za-z]+)-?(\d+)\s+([<>]?\s*[-0-9.]+)\s+([<>]?\s*[-0-9.]+)")


@dataclass(frozen=True)
class HppRecord:
    source_resname: str
    resseq: int
    pkint: str
    pkhalf: str


@dataclass(frozen=True)
class MetalContact:
    chain: str
    resseq: int
    resname: str
    atom: str
    element: str
    distance_a: float


def parse_pkout(path: Path) -> list[HppRecord]:
    records: list[HppRecord] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = PKOUT_RE.match(line)
        if not match:
            continue
        raw_resname, resseq, pkint, pkhalf = match.groups()
        resname = raw_resname.upper()
        # H++ can emit pseudo terminal entries such as NTASP-16, aspNT-1, CTLYS-377.
        if resname.startswith("NT") or resname.endswith("NT") or resname.startswith("CT") or resname.endswith("CT"):
            continue
        records.append(HppRecord(resname, int(resseq), pkint.replace(" ", ""), pkhalf.replace(" ", "")))
    return records


def parse_pdb_atom(line: str) -> tuple[str, int, str, str, str, tuple[float, float, float]] | None:
    if not line.startswith(("ATOM  ", "HETATM")):
        return None
    try:
        atom = line[12:16].strip()
        resname = line[17:20].strip().upper()
        chain = line[21:22].strip() or "-"
        resseq = int(line[22:26])
        element = line[76:78].strip().upper() or atom[:1].upper()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except Exception:
        return None
    return chain, resseq, resname, atom, element, xyz


def find_metal_contacts(pdb: Path, *, chain: str = "A", cutoff_a: float = 2.8) -> list[MetalContact]:
    atoms = []
    metals = []
    for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        parsed = parse_pdb_atom(line)
        if not parsed:
            continue
        atom_chain, resseq, resname, atom, element, xyz = parsed
        if atom_chain != chain:
            continue
        record = (atom_chain, resseq, resname, atom, element, xyz)
        atoms.append(record)
        if element == "FE":
            metals.append(record)

    contacts: dict[tuple[str, int, str], MetalContact] = {}
    for metal in metals:
        for atom_chain, resseq, resname, atom, element, xyz in atoms:
            if element not in {"N", "O", "S"}:
                continue
            distance = math.dist(metal[5], xyz)
            if distance > cutoff_a:
                continue
            key = (atom_chain, resseq, atom)
            old = contacts.get(key)
            if old is None or distance < old.distance_a:
                contacts[key] = MetalContact(atom_chain, resseq, resname, atom, element, distance)
    return sorted(contacts.values(), key=lambda item: item.distance_a)


def hpp_resname_for_amber(record: HppRecord, *, ph: float = 7.4) -> tuple[str | None, str]:
    resname = record.source_resname
    if resname in {"HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYM"}:
        return resname, "explicit_hplusplus_resname"
    if resname == "HIS":
        pkhalf = numeric_pk(record.pkhalf)
        if pkhalf is not None and pkhalf > ph + 0.3:
            return "HIP", "heuristic_histidine_pkhalf_above_ph"
        return "HIE", "heuristic_histidine_default_epsilon"
    if resname == "ASP":
        pkhalf = numeric_pk(record.pkhalf)
        if pkhalf is not None and pkhalf > ph:
            return "ASH", "heuristic_acid_pkhalf_above_ph"
        return None, "standard_asp"
    if resname == "GLU":
        pkhalf = numeric_pk(record.pkhalf)
        if pkhalf is not None and pkhalf > ph:
            return "GLH", "heuristic_acid_pkhalf_above_ph"
        return None, "standard_glu"
    if resname == "CYS":
        pkhalf = numeric_pk(record.pkhalf)
        if pkhalf is not None and pkhalf < ph:
            return "CYM", "heuristic_cys_pkhalf_below_ph"
        return None, "standard_cys"
    if resname == "LYS":
        pkhalf = numeric_pk(record.pkhalf)
        if pkhalf is not None and pkhalf < ph:
            return "LYN", "heuristic_lys_pkhalf_below_ph"
        return None, "standard_lys"
    return None, "not_titration_name_for_amber"


def numeric_pk(text: str) -> float | None:
    value = text.strip()
    if value.startswith(">"):
        return 99.0
    if value.startswith("<"):
        return -99.0
    try:
        return float(value)
    except ValueError:
        return None


def write_protonation_map(
    *,
    pkout: Path,
    pdb: Path,
    out_csv: Path,
    chain: str = "A",
    ph: float = 7.4,
    metal_cutoff_a: float = 2.8,
) -> tuple[int, int]:
    records = parse_pkout(pkout)
    contacts = find_metal_contacts(pdb, chain=chain, cutoff_a=metal_cutoff_a)
    contact_by_residue = {(item.chain, item.resseq) for item in contacts}
    rows: list[dict[str, object]] = []
    for record in records:
        amber_resname, reason = hpp_resname_for_amber(record, ph=ph)
        if amber_resname is None:
            continue
        if (chain, record.resseq) in contact_by_residue and record.source_resname.startswith("HI"):
            reason = f"{reason};metal_contact_within_{metal_cutoff_a:g}A"
        rows.append(
            {
                "chain": chain,
                "resseq": record.resseq,
                "resname": amber_resname,
                "source_resname": record.source_resname,
                "pkhalf": record.pkhalf,
                "pkint": record.pkint,
                "reason": reason,
            }
        )

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["chain", "resseq", "resname", "source_resname", "pkhalf", "pkint", "reason"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    contacts_csv = out_csv.with_name(out_csv.stem + "_metal_contacts.csv")
    with contacts_csv.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["chain", "resseq", "resname", "atom", "element", "distance_a"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in contacts:
            writer.writerow(
                {
                    "chain": item.chain,
                    "resseq": item.resseq,
                    "resname": item.resname,
                    "atom": item.atom,
                    "element": item.element,
                    "distance_a": f"{item.distance_a:.3f}",
                }
            )
    return len(rows), len(contacts)
