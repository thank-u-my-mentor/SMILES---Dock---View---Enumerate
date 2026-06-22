#!/usr/bin/env python3
"""Element-level Metal-F ligand rescan.

This workflow does not modify the older curated XLSX/FASTA files. It builds a
new sidecar result set by:

1. Downloading/caching the wwPDB Chemical Component Dictionary.
2. Finding every CCD component whose formula contains F.
3. Finding every CCD component whose formula contains Fe/Co/Ni/Mn/Cr/Cu.
4. Using RCSB Search API comp_id searches to get candidate PDB entries.
5. Downloading mmCIF coordinates and computing ligand-F to metal distances.
6. Writing a human-readable XLSX plus protein/CDS FASTA for non-complex hits.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import shlex
import sys
import textwrap
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


PROJECT = Path("/home/qin/Metal-F_project")
INTERACTIONS = PROJECT / "Metal-F_interactions"
OLD_XLSX = INTERACTIONS / "Metal-F_geometry_fast_hits.xlsx"
CACHE = PROJECT / ".cache" / "metal_f_element_ligand_rescan"
CCD_GZ = CACHE / "components.cif.gz"
CCD_URL = "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
TARGET_METALS = {"FE", "CO", "NI", "MN", "CR", "CU"}
SOLVENTS = {"HOH", "WAT", "DOD"}

sys.path.insert(0, str(PROJECT))
import enrich_metal_f_hits_with_sequences as seqmeta  # noqa: E402


def read_url(url: str, *, data: bytes | None = None, headers: dict | None = None, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, data=data, headers=headers or {"User-Agent": "metal-f-rescan/1.0"})
    last_error = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError, ConnectionResetError) as exc:
            last_error = exc
            time.sleep(1.5 * (attempt + 1))
    raise last_error


def download_ccd() -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    if CCD_GZ.exists() and CCD_GZ.stat().st_size > 10_000_000:
        return CCD_GZ
    print(f"Downloading CCD: {CCD_URL}", flush=True)
    payload = read_url(CCD_URL, timeout=240)
    CCD_GZ.write_bytes(payload)
    return CCD_GZ


def parse_value(line: str) -> str:
    try:
        parts = shlex.split(line, posix=True)
    except ValueError:
        parts = line.split(None, 1)
    if len(parts) < 2:
        return ""
    return parts[1].strip()


def formula_elements(formula: str) -> set[str]:
    # CCD formulas are commonly like "C 19 H 17 F3 O2 S" or "Cu 2".
    return {item.upper() for item in re.findall(r"([A-Z][a-z]?)(?=[0-9\s+\-.]|$)", formula or "")}


def parse_ccd(path: Path) -> dict[str, dict]:
    components: dict[str, dict] = {}
    current: dict[str, str] | None = None

    def finish() -> None:
        if not current or not current.get("comp_id"):
            return
        formula = current.get("formula", "")
        elements = formula_elements(formula)
        comp_id = current["comp_id"].upper()
        components[comp_id] = {
            "comp_id": comp_id,
            "name": current.get("name", ""),
            "type": current.get("type", ""),
            "formula": formula,
            "elements": ";".join(sorted(elements)),
            "has_f": "F" in elements,
            "target_metal_elements": ";".join(sorted(elements & TARGET_METALS)),
        }

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("data_"):
                finish()
                current = {"comp_id": line[5:].strip().upper()}
                continue
            if current is None:
                continue
            if line.startswith("_chem_comp.id"):
                current["comp_id"] = parse_value(line).upper()
            elif line.startswith("_chem_comp.name"):
                current["name"] = parse_value(line)
            elif line.startswith("_chem_comp.type"):
                current["type"] = parse_value(line)
            elif line.startswith("_chem_comp.formula "):
                current["formula"] = parse_value(line)
        finish()
    return components


def chunks(values: list[str], size: int) -> list[list[str]]:
    return [values[i : i + size] for i in range(0, len(values), size)]


def search_comp_ids(comp_ids: list[str], *, chunk_size: int = 100, quiet: bool = False) -> set[str]:
    out: set[str] = set()
    comp_ids = sorted(set(comp_ids))
    for idx, chunk in enumerate(chunks(comp_ids, chunk_size), 1):
        if not quiet:
            print(f"RCSB comp_id query chunk {idx}/{math.ceil(len(comp_ids) / chunk_size)} ({len(chunk)} IDs)", flush=True)
        start = 0
        while True:
            payload = {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
                        "operator": "in",
                        "value": chunk,
                    },
                },
                "return_type": "entry",
                "request_options": {
                    "paginate": {"start": start, "rows": 1000},
                    "results_content_type": ["experimental"],
                },
            }
            data = read_url(
                SEARCH_URL,
                data=json.dumps(payload).encode("utf-8"),
                headers={"Content-Type": "application/json", "User-Agent": "metal-f-rescan/1.0"},
                timeout=120,
            )
            result = json.loads(data.decode("utf-8"))
            rows = result.get("result_set") or []
            out.update(row["identifier"].upper() for row in rows if row.get("identifier"))
            if start + len(rows) >= int(result.get("total_count") or 0) or not rows:
                break
            start += len(rows)
    return out


def download_cif(pdb_id: str) -> tuple[Path, bool]:
    cif_dir = CACHE / "cif"
    cif_dir.mkdir(parents=True, exist_ok=True)
    path = cif_dir / f"{pdb_id.upper()}.cif"
    if path.exists() and path.stat().st_size > 0:
        return path, False
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    path.write_bytes(read_url(url, timeout=180))
    return path, True


def as_list(data: dict, key: str) -> list[str]:
    value = data.get(key, [])
    if isinstance(value, str):
        return [value]
    return list(value)


def atom_site(cif_path: Path) -> list[dict]:
    data = MMCIF2Dict(str(cif_path))
    xs = as_list(data, "_atom_site.Cartn_x")
    ys = as_list(data, "_atom_site.Cartn_y")
    zs = as_list(data, "_atom_site.Cartn_z")
    elems = as_list(data, "_atom_site.type_symbol")
    comp_label = as_list(data, "_atom_site.label_comp_id")
    comp_auth = as_list(data, "_atom_site.auth_comp_id") or comp_label
    atom_label = as_list(data, "_atom_site.label_atom_id")
    atom_auth = as_list(data, "_atom_site.auth_atom_id") or atom_label
    chain = as_list(data, "_atom_site.auth_asym_id") or as_list(data, "_atom_site.label_asym_id")
    seq = as_list(data, "_atom_site.auth_seq_id") or as_list(data, "_atom_site.label_seq_id")
    model = as_list(data, "_atom_site.pdbx_PDB_model_num") or ["1"] * len(xs)
    group = as_list(data, "_atom_site.group_PDB") or [""] * len(xs)
    atoms = []
    for i in range(len(xs)):
        if str(model[i]) not in {"1", "."}:
            continue
        try:
            xyz = (float(xs[i]), float(ys[i]), float(zs[i]))
        except Exception:
            continue
        atoms.append(
            {
                "elem": str(elems[i]).strip().upper(),
                "comp": str(comp_auth[i] if i < len(comp_auth) else comp_label[i]).strip().upper(),
                "atom": str(atom_auth[i] if i < len(atom_auth) else atom_label[i]).strip(),
                "chain": str(chain[i]).strip() if i < len(chain) else "",
                "seq": str(seq[i]).strip() if i < len(seq) else "",
                "group": str(group[i]).strip().upper() if i < len(group) else "",
                "xyz": xyz,
            }
        )
    return atoms


def atom_label(atom: dict) -> str:
    return f"{atom['comp']}:{atom['chain']}:{atom['seq']}:{atom['atom']}:{atom['elem']}"


def distance(a: dict, b: dict) -> float:
    return math.sqrt(sum((a["xyz"][i] - b["xyz"][i]) ** 2 for i in range(3)))


def geometry_for_pdb(pdb_id: str, components: dict[str, dict], cutoff: float) -> dict:
    cif_path, cif_downloaded = download_cif(pdb_id)
    atoms = atom_site(cif_path)
    f_atoms = [
        atom
        for atom in atoms
        if atom["elem"] == "F"
        and atom["group"] == "HETATM"
        and atom["comp"] not in SOLVENTS
        and components.get(atom["comp"], {}).get("has_f")
    ]
    metal_atoms = [
        atom
        for atom in atoms
        if atom["elem"] in TARGET_METALS and atom["group"] == "HETATM" and atom["comp"] not in SOLVENTS
    ]
    nearest = None
    pairs = []
    for f_atom in f_atoms:
        for metal in metal_atoms:
            dist = distance(f_atom, metal)
            if nearest is None or dist < nearest[0]:
                nearest = (dist, f_atom, metal)
            if dist <= cutoff:
                pairs.append((dist, f_atom, metal))
    pairs.sort(key=lambda row: row[0])
    if not nearest:
        return {
            "status": "no_f_or_metal_atoms",
            "pdb_id": pdb_id,
            "cutoff_A": cutoff,
            "nearest_F_metal_A": "",
            "nearest_F_atom": "",
            "nearest_metal_atom": "",
            "passing_pairs_le_cutoff_A": "",
            "fluorine_comp_ids": ";".join(sorted({atom["comp"] for atom in f_atoms})),
            "metal_comp_ids": ";".join(sorted({atom["comp"] for atom in metal_atoms})),
            "metal_elements": ";".join(sorted({atom["elem"] for atom in metal_atoms})),
            "f_atom_count": len(f_atoms),
            "metal_atom_count": len(metal_atoms),
            "pair_count_le_cutoff": 0,
            "cif_path": str(cif_path),
            "cif_downloaded_this_run": "YES" if cif_downloaded else "NO",
        }
    f_pair_comps = sorted({pair[1]["comp"] for pair in pairs})
    m_pair_comps = sorted({pair[2]["comp"] for pair in pairs})
    f_all_comps = sorted({atom["comp"] for atom in f_atoms})
    m_all_comps = sorted({atom["comp"] for atom in metal_atoms})
    return {
        "status": "geometry_hit" if pairs else "no_geometry_hit",
        "pdb_id": pdb_id,
        "cutoff_A": cutoff,
        "nearest_F_metal_A": f"{nearest[0]:.3f}",
        "nearest_F_atom": atom_label(nearest[1]),
        "nearest_metal_atom": atom_label(nearest[2]),
        "passing_pairs_le_cutoff_A": ";".join(
            f"{dist:.3f}:{atom_label(f_atom)}--{atom_label(metal)}" for dist, f_atom, metal in pairs[:50]
        ),
        "fluorine_comp_ids": ";".join(f_pair_comps or f_all_comps),
        "fluorine_ligand_names": ";".join(components.get(comp, {}).get("name", "") for comp in (f_pair_comps or f_all_comps)),
        "fluorine_ligand_formulas": ";".join(
            components.get(comp, {}).get("formula", "") for comp in (f_pair_comps or f_all_comps)
        ),
        "metal_comp_ids": ";".join(m_pair_comps or m_all_comps),
        "metal_ligand_names": ";".join(components.get(comp, {}).get("name", "") for comp in (m_pair_comps or m_all_comps)),
        "metal_ligand_formulas": ";".join(components.get(comp, {}).get("formula", "") for comp in (m_pair_comps or m_all_comps)),
        "metal_elements": ";".join(sorted({pair[2]["elem"] for pair in pairs} or {atom["elem"] for atom in metal_atoms})),
        "f_atom_count": len(f_atoms),
        "metal_atom_count": len(metal_atoms),
        "pair_count_le_cutoff": len(pairs),
        "cif_path": str(cif_path),
        "cif_downloaded_this_run": "YES" if cif_downloaded else "NO",
    }


def load_old_status(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    frame = pd.read_excel(path)
    out = {}
    for _, row in frame.iterrows():
        pdb_id = str(row.get("pdb_id", "")).upper()
        if not pdb_id or pdb_id == "NAN":
            continue
        nearest = row.get("nearest_F_metal_Angstrom", row.get("nearest_F_metal_A", ""))
        out[pdb_id] = {
            "original_row_present": True,
            "original_geometry_filled": bool(pd.notna(nearest) and str(nearest).strip() and str(nearest).lower() != "nan"),
        }
    return out


def read_pdb_ids_from_xlsx(path: Path) -> list[str]:
    if not path.exists():
        return []
    frame = pd.read_excel(path)
    if "pdb_id" not in frame.columns:
        return []
    ids = []
    for value in frame["pdb_id"]:
        pdb_id = str(value or "").strip().upper()
        if re.fullmatch(r"[0-9][A-Z0-9]{3}", pdb_id):
            ids.append(pdb_id)
    return sorted(dict.fromkeys(ids))


def parse_pdb_ids(text: str) -> list[str]:
    ids = re.findall(r"\b[0-9][A-Za-z0-9]{3}\b", text or "")
    return sorted(dict.fromkeys(item.upper() for item in ids))


def write_candidate_cache(path: Path, candidates: list[str], f_entries: set[str], metal_entries: set[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pdb_id", "in_f_entries", "in_metal_entries"])
        writer.writeheader()
        for pdb_id in candidates:
            writer.writerow(
                {
                    "pdb_id": pdb_id,
                    "in_f_entries": "YES" if pdb_id in f_entries else "NO",
                    "in_metal_entries": "YES" if pdb_id in metal_entries else "NO",
                }
            )


def read_candidate_cache(path: Path) -> list[str]:
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8", errors="replace") as handle:
        rows = list(csv.DictReader(handle))
    ids = []
    for row in rows:
        pdb_id = str(row.get("pdb_id", "")).strip().upper()
        if re.fullmatch(r"[0-9][A-Z0-9]{3}", pdb_id):
            ids.append(pdb_id)
    return sorted(dict.fromkeys(ids))


def read_done_ledger(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    done = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle):
            pdb_id = str(row.get("pdb_id", "")).strip().upper()
            if pdb_id:
                done[pdb_id] = row
    return done


def append_done(path: Path, row: dict) -> None:
    fields = [
        "pdb_id",
        "status",
        "nearest_F_metal_A",
        "protein_entity_count",
        "complex_warning",
        "fasta_seed_included",
        "cif_path",
        "cif_downloaded_this_run",
        "cif_deleted",
        "cif_delete_error",
        "error",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists() and path.stat().st_size > 0
    with open(path, "a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerow({field: row.get(field, "") for field in fields})


def print_status(candidates: list[str], done: dict[str, dict]) -> None:
    candidate_set = set(candidates)
    done_set = candidate_set & set(done)
    hit_count = sum(1 for pdb_id in done_set if done[pdb_id].get("status") == "geometry_hit")
    error_count = sum(1 for pdb_id in done_set if done[pdb_id].get("status") == "error")
    remaining = len(candidate_set - done_set)
    print(f"candidate_total: {len(candidate_set)}")
    print(f"processed: {len(done_set)}")
    print(f"geometry_hits: {hit_count}")
    print(f"errors: {error_count}")
    print(f"remaining: {remaining}")
    if remaining:
        preview = sorted(candidate_set - done_set)[:20]
        print("next_unprocessed_preview:", ",".join(preview))


def get_entry_doi(pdb_id: str) -> str:
    try:
        data = seqmeta.get_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
        citation = data.get("rcsb_primary_citation") or {}
        return citation.get("pdbx_database_id_DOI") or ""
    except Exception:
        return ""


def enrich_hit(row: dict, old_status: dict[str, dict]) -> tuple[dict, list[dict], list[dict]]:
    pdb_id = row["pdb_id"].upper()
    entities = seqmeta.get_polymer_entities(pdb_id)
    protein_entity_count = len(entities)
    complex_warning = protein_entity_count > 1
    expressions = sorted({e["expression_system"] for e in entities if e.get("expression_system") and e["expression_system"] != "unknown"})
    organisms = sorted({e["source_organism"] for e in entities if e.get("source_organism") and e["source_organism"] != "unknown"})
    kingdoms = sorted({e["kingdom"] for e in entities if e.get("kingdom")})
    uniprots = sorted({u for e in entities for u in e.get("uniprot_ids", "").split(";") if u})
    ecoli = any("escherichia coli" in e.get("expression_system", "").lower() for e in entities)
    fasta_seed_included = protein_entity_count == 1
    protein_records = []
    dna_records = []

    if fasta_seed_included:
        for entity in entities:
            if not entity.get("protein_sequence"):
                continue
            header = (
                f"{pdb_id}_{entity['entity_id']}|chains={entity['chains']}"
                f"|uniprot={entity['uniprot_ids']}|organism={entity['source_organism']}"
                f"|kingdom={entity['kingdom']}|expression_system={entity['expression_system']}"
            )
            protein_records.append({"pdb_id": pdb_id, "header": header, "sequence": entity["protein_sequence"]})
            for uniprot_id in [item for item in entity.get("uniprot_ids", "").split(";") if item]:
                try:
                    refs = seqmeta.uniprot_embl_refs(uniprot_id)
                except Exception:
                    refs = []
                fetched = False
                for ref in refs:
                    if fetched:
                        break
                    coding_rows = seqmeta.ena_coding_rows_by_protein_id(ref.get("protein_id", ""))
                    coding = coding_rows[0].get("accession", "") if coding_rows else ""
                    records = seqmeta.fetch_ena_fasta(coding or ref.get("embl_accession", ""))
                    if not records:
                        continue
                    ena_header, dna_seq = records[0]
                    dna_header = (
                        f"{pdb_id}_{entity['entity_id']}|chains={entity['chains']}"
                        f"|uniprot={uniprot_id}|embl={ref.get('embl_accession','')}"
                        f"|coding={coding}|protein_id={ref.get('protein_id','')}"
                        f"|organism={entity['source_organism']}|kingdom={entity['kingdom']}"
                    )
                    dna_records.append(
                        {
                            "pdb_id": pdb_id,
                            "header": dna_header,
                            "sequence": dna_seq,
                            "ena_header": ena_header,
                            "uniprot_id": uniprot_id,
                            "embl_accession": ref.get("embl_accession", ""),
                            "coding_accession": coding,
                        }
                    )
                    fetched = True

    original = old_status.get(pdb_id, {"original_row_present": False, "original_geometry_filled": False})
    row.update(
        {
            **original,
            "supplement_status": (
                "already_filled"
                if original["original_geometry_filled"]
                else "fill_existing_blank"
                if original["original_row_present"]
                else "new_candidate"
            ),
            "doi": get_entry_doi(pdb_id),
            "protein_entity_count": protein_entity_count,
            "complex_warning": "YES" if complex_warning else "NO",
            "fasta_seed_included": "YES" if fasta_seed_included else "NO",
            "ecoli_expression_flag": "YES" if ecoli else "NO",
            "source_organism": ";".join(organisms) if organisms else "unknown",
            "kingdom": ";".join(kingdoms) if kingdoms else "Unknown",
            "expression_system": ";".join(expressions) if expressions else "unknown",
            "uniprot_ids": ";".join(uniprots),
            "protein_fasta_headers": ";".join(item["header"] for item in protein_records),
            "dna_fasta_headers": ";".join(item["header"] for item in dna_records),
        }
    )
    return row, protein_records, dna_records


def write_fasta(path: Path, records: list[dict]) -> None:
    seen = set()
    with open(path, "w", encoding="utf-8") as handle:
        for record in records:
            key = (record["header"], record["sequence"])
            if key in seen:
                continue
            seen.add(key)
            handle.write(f">{record['header']}\n")
            handle.write("\n".join(textwrap.wrap(record["sequence"], 80)) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", default=str(INTERACTIONS / "Metal-F_element_ligand_rescan"))
    parser.add_argument("--cutoff", type=float, default=3.5)
    parser.add_argument("--chunk-size", type=int, default=100)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--max-candidates", type=int, default=0)
    parser.add_argument("--pdb-ids", default="", help="Comma/space separated PDB IDs to process directly.")
    parser.add_argument("--pdb-ids-from-xlsx", default="", help="Read PDB IDs from an XLSX file with a pdb_id column.")
    parser.add_argument("--candidate-cache", default=str(CACHE / "Metal-F_element_ligand_candidates.csv"))
    parser.add_argument("--reuse-candidates", action="store_true", help="Read candidate-cache instead of querying RCSB.")
    parser.add_argument("--done-ledger", default=str(CACHE / "Metal-F_element_ligand_done.csv"))
    parser.add_argument("--resume", action="store_true", help="Skip PDB IDs already present in done-ledger.")
    parser.add_argument("--status-only", action="store_true", help="Print candidate/done/remaining counts and exit.")
    parser.add_argument(
        "--delete-cif-nonhits",
        action="store_true",
        help="Delete CIF files for candidates that do not satisfy the distance cutoff.",
    )
    parser.add_argument(
        "--delete-cif-nonhits-even-if-cached",
        action="store_true",
        help="Also delete non-hit CIF files that existed before this run.",
    )
    parser.add_argument("--quiet-candidates", action="store_true", help="Suppress comp_id candidate-query chunk messages.")
    parser.add_argument("--skip-dna", action="store_true")
    args = parser.parse_args()

    output_prefix = Path(args.prefix)
    ccd = parse_ccd(download_ccd())
    f_components = {cid: row for cid, row in ccd.items() if row["has_f"]}
    metal_components = {cid: row for cid, row in ccd.items() if row["target_metal_elements"]}
    print(f"CCD F components: {len(f_components)}", flush=True)
    print(f"CCD target-metal components: {len(metal_components)}", flush=True)

    direct_ids = parse_pdb_ids(args.pdb_ids)
    xlsx_ids = read_pdb_ids_from_xlsx(Path(args.pdb_ids_from_xlsx)) if args.pdb_ids_from_xlsx else []
    if direct_ids or xlsx_ids:
        f_entries = set()
        metal_entries = set()
        candidates = sorted(dict.fromkeys(direct_ids + xlsx_ids))
    elif args.reuse_candidates:
        f_entries = set()
        metal_entries = set()
        candidates = read_candidate_cache(Path(args.candidate_cache))
        if not candidates:
            raise SystemExit(f"No candidates found in {args.candidate_cache}. Run once without --reuse-candidates.")
    else:
        f_entries = search_comp_ids(list(f_components), chunk_size=args.chunk_size, quiet=args.quiet_candidates)
        metal_entries = search_comp_ids(list(metal_components), chunk_size=args.chunk_size, quiet=args.quiet_candidates)
        candidates = sorted(f_entries & metal_entries)
        write_candidate_cache(Path(args.candidate_cache), candidates, f_entries, metal_entries)
    if args.limit:
        candidates = candidates[: args.limit]
    if args.start or args.max_candidates:
        end = None if not args.max_candidates else args.start + args.max_candidates
        candidates = candidates[args.start:end]
    print(f"F-entry count: {len(f_entries)}", flush=True)
    print(f"metal-entry count: {len(metal_entries)}", flush=True)
    print(f"intersection candidate count: {len(candidates)}", flush=True)

    done_ledger = Path(args.done_ledger)
    done = read_done_ledger(done_ledger)
    if args.status_only:
        print_status(candidates, done)
        return
    if args.resume:
        before = len(candidates)
        candidates = [pdb_id for pdb_id in candidates if pdb_id not in done]
        print(f"resume_skip_already_processed: {before - len(candidates)}", flush=True)
        print(f"resume_remaining_this_selection: {len(candidates)}", flush=True)

    old_status = load_old_status(OLD_XLSX)
    hit_rows = []
    processed_rows = []
    protein_records_all: list[dict] = []
    dna_records_all: list[dict] = []
    errors = []

    for index, pdb_id in enumerate(candidates, 1):
        print(f"[{index}/{len(candidates)}] geometry {pdb_id}", flush=True)
        try:
            row = geometry_for_pdb(pdb_id, ccd, args.cutoff)
            ledger_row = dict(row)
            processed_rows.append(ledger_row)
            if row.get("status") != "geometry_hit":
                should_delete_cif = args.delete_cif_nonhits and row.get("cif_path") and (
                    row.get("cif_downloaded_this_run") == "YES" or args.delete_cif_nonhits_even_if_cached
                )
                if should_delete_cif:
                    try:
                        Path(row["cif_path"]).unlink(missing_ok=True)
                        ledger_row["cif_deleted"] = "YES"
                        processed_rows[-1]["cif_deleted"] = "YES"
                    except Exception as exc:
                        ledger_row["cif_delete_error"] = str(exc)
                        processed_rows[-1]["cif_delete_error"] = str(exc)
                append_done(done_ledger, ledger_row)
                continue
            enriched, protein_records, dna_records = enrich_hit(row, old_status)
            if args.skip_dna:
                dna_records = []
                enriched["dna_fasta_headers"] = ""
            hit_rows.append(enriched)
            protein_records_all.extend(protein_records)
            dna_records_all.extend(dna_records)
            append_done(done_ledger, enriched)
        except Exception as exc:
            error_row = {"pdb_id": pdb_id, "status": "error", "error": str(exc)}
            errors.append(error_row)
            processed_rows.append(error_row)
            append_done(done_ledger, error_row)

    hit_rows.sort(key=lambda row: (float(row.get("nearest_F_metal_A") or 999), row["pdb_id"]))
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    protein_fasta = output_prefix.with_name(output_prefix.name + "_protein.fasta")
    dna_fasta = output_prefix.with_name(output_prefix.name + "_cds_dna.fasta")
    xlsx = output_prefix.with_suffix(".xlsx")
    write_fasta(protein_fasta, protein_records_all)
    write_fasta(dna_fasta, dna_records_all)

    summary = [
        {"metric": "ccd_f_component_count", "value": len(f_components)},
        {"metric": "ccd_target_metal_component_count", "value": len(metal_components)},
        {"metric": "f_entry_count", "value": len(f_entries)},
        {"metric": "metal_entry_count", "value": len(metal_entries)},
        {"metric": "intersection_candidate_count", "value": len(candidates)},
        {"metric": "geometry_hit_count", "value": len(hit_rows)},
        {"metric": "protein_fasta_records", "value": len(protein_records_all)},
        {"metric": "dna_fasta_records", "value": len(dna_records_all)},
        {"metric": "cutoff_A", "value": args.cutoff},
        {"metric": "target_metals", "value": ";".join(sorted(TARGET_METALS))},
        {"metric": "protein_fasta", "value": str(protein_fasta)},
        {"metric": "dna_fasta", "value": str(dna_fasta)},
    ]

    with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
        pd.DataFrame(hit_rows).to_excel(writer, sheet_name="geometry_hits", index=False)
        pd.DataFrame(processed_rows).to_excel(writer, sheet_name="processed_candidates", index=False)
        pd.DataFrame(protein_records_all).to_excel(writer, sheet_name="protein_fasta", index=False)
        pd.DataFrame(dna_records_all).to_excel(writer, sheet_name="dna_fasta", index=False)
        pd.DataFrame(f_components.values()).sort_values("comp_id").to_excel(writer, sheet_name="ccd_f_ligands", index=False)
        pd.DataFrame(metal_components.values()).sort_values("comp_id").to_excel(
            writer, sheet_name="ccd_metal_ligands", index=False
        )
        pd.DataFrame(summary).to_excel(writer, sheet_name="run_summary", index=False)
        if errors:
            pd.DataFrame(errors).to_excel(writer, sheet_name="errors", index=False)

    print(f"XLSX: {xlsx}")
    print(f"Protein FASTA: {protein_fasta}")
    print(f"DNA FASTA: {dna_fasta}")


if __name__ == "__main__":
    main()
