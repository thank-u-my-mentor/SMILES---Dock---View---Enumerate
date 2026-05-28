#!/usr/bin/env python3
"""Filter fluorinated organic ligand candidates by ligand-F to metal distance.

Inputs are the Metal-Name candidate directories made by metal_fluoro_ligand_search.py.
The hard geometric rule is min organic-ligand fluorine atom to selected metal atom
distance <= 3.5 Angstrom.
"""

import argparse
import csv
import json
import math
import os
import re
import shlex
import time
import textwrap
import urllib.error
import urllib.request


PROJECT = "/home/qin/Metal-F_project"
METALS = {"FE", "CO", "NI", "CU", "MN"}
SKIP_COMP_IDS = {"F", "HF", "BF4", "PF6", "ALF", "FLC"}


def read_url(url, data=None, timeout=120):
    headers = {"User-Agent": "metal-fluoro-geometry-pipeline/1.0"}
    if data is not None:
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers)
    last_error = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return response.read()
        except (urllib.error.URLError, ConnectionResetError, TimeoutError) as exc:
            last_error = exc
            time.sleep(1.5 * (attempt + 1))
    raise last_error


def download_cif(pdb_id, cif_dir):
    os.makedirs(cif_dir, exist_ok=True)
    pdb_id = pdb_id.upper()
    path = os.path.join(cif_dir, f"{pdb_id}.cif")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    content = read_url(f"https://files.rcsb.org/download/{pdb_id}.cif")
    with open(path, "wb") as handle:
        handle.write(content)
    return path


def parse_loop(lines, start):
    tags = []
    i = start + 1
    while i < len(lines) and lines[i].strip().startswith("_"):
        tags.append(lines[i].strip())
        i += 1
    rows = []
    while i < len(lines):
        line = lines[i].rstrip("\n")
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            i += 1
            if rows:
                break
            continue
        if stripped == "loop_" or stripped.startswith("_"):
            break
        lexer = shlex.shlex(line, posix=True)
        lexer.whitespace_split = True
        lexer.commenters = ""
        values = list(lexer)
        if len(values) == len(tags):
            rows.append(dict(zip(tags, values)))
        i += 1
    return tags, rows, i


def mmcif_atom_site_rows(path):
    with open(path, encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            j = i + 1
            tags = []
            while j < len(lines) and lines[j].strip().startswith("_"):
                tags.append(lines[j].strip())
                j += 1
            if tags and tags[0].startswith("_atom_site."):
                rows = []
                while j < len(lines):
                    stripped = lines[j].strip()
                    if not stripped or stripped.startswith("#"):
                        j += 1
                        if rows:
                            break
                        continue
                    if stripped == "loop_" or stripped.startswith("_"):
                        break
                    try:
                        lexer = shlex.shlex(lines[j], posix=True)
                        lexer.whitespace_split = True
                        lexer.commenters = ""
                        values = list(lexer)
                    except ValueError:
                        values = []
                    if len(values) == len(tags):
                        rows.append(dict(zip(tags, values)))
                    j += 1
                return rows
            i = j
        else:
            i += 1
    return []


def atom_site_records(path):
    atoms = []
    for row in mmcif_atom_site_rows(path):
        try:
            x = float(row.get("_atom_site.Cartn_x", "nan"))
            y = float(row.get("_atom_site.Cartn_y", "nan"))
            z = float(row.get("_atom_site.Cartn_z", "nan"))
        except ValueError:
            continue
        if not all(math.isfinite(v) for v in (x, y, z)):
            continue
        model = row.get("_atom_site.pdbx_PDB_model_num", "1")
        if model not in ("1", ".", "?"):
            continue
        atoms.append(
            {
                "group": row.get("_atom_site.group_PDB", ""),
                "elem": row.get("_atom_site.type_symbol", "").upper(),
                "atom": row.get("_atom_site.label_atom_id", ""),
                "comp": row.get("_atom_site.label_comp_id", "").upper(),
                "chain": row.get("_atom_site.auth_asym_id")
                or row.get("_atom_site.label_asym_id", ""),
                "resi": row.get("_atom_site.auth_seq_id")
                or row.get("_atom_site.label_seq_id", ""),
                "x": x,
                "y": y,
                "z": z,
            }
        )
    return atoms


def distance(a, b):
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def is_organic_fluorine(atom):
    if atom["group"] != "HETATM":
        return False
    if atom["elem"] != "F":
        return False
    if atom["comp"] in SKIP_COMP_IDS:
        return False
    return True


def geometry_hits(path, cutoff):
    atoms = atom_site_records(path)
    ligand_f = [a for a in atoms if is_organic_fluorine(a)]
    metals = [a for a in atoms if a["group"] == "HETATM" and a["comp"] in METALS and a["elem"] in METALS]
    pairs = []
    for f_atom in ligand_f:
        for metal in metals:
            d = distance(f_atom, metal)
            if d <= cutoff:
                pairs.append((d, f_atom, metal))
    pairs.sort(key=lambda row: row[0])
    nearest = None
    for f_atom in ligand_f:
        for metal in metals:
            d = distance(f_atom, metal)
            if nearest is None or d < nearest[0]:
                nearest = (d, f_atom, metal)
    return pairs, nearest, ligand_f, metals


def atom_label(a):
    if not a:
        return ""
    return f"{a['comp']}:{a['chain']}:{a['resi']}:{a['atom']}"


def rcsb_entry_metadata(pdb_id):
    data = json.loads(read_url(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}").decode("utf-8"))
    citation = data.get("rcsb_primary_citation") or {}
    if not citation and data.get("citation"):
        primary = [row for row in data["citation"] if row.get("rcsb_is_primary") == "Y"]
        citation = primary[0] if primary else data["citation"][0]
    return {
        "pdb_id": pdb_id,
        "pdb_url": f"https://www.rcsb.org/structure/{pdb_id}",
        "doi": citation.get("pdbx_database_id_DOI", ""),
        "pubmed_id": str(citation.get("pdbx_database_id_PubMed", "") or ""),
        "citation_title": citation.get("title", ""),
        "citation_year": str(citation.get("year", "") or ""),
        "journal": citation.get("journal_abbrev") or citation.get("rcsb_journal_abbrev", ""),
        "authors": "; ".join(citation.get("rcsb_authors", [])[:8]),
        "structure_title": data.get("struct", {}).get("title", ""),
        "polymer_entity_ids": data.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", []),
    }


def rcsb_polymer_entities(pdb_id):
    meta = rcsb_entry_metadata(pdb_id)
    rows = []
    for entity_id in meta.get("polymer_entity_ids", []):
        try:
            data = json.loads(
                read_url(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}").decode("utf-8")
            )
        except Exception:
            continue
        entity_poly = data.get("entity_poly", {})
        if entity_poly.get("rcsb_entity_polymer_type") != "Protein":
            continue
        ids = data.get("rcsb_polymer_entity_container_identifiers", {})
        source = {}
        sources = data.get("rcsb_entity_source_organism") or data.get("entity_src_nat") or []
        if sources:
            source = sources[0]
        sequence = entity_poly.get("pdbx_seq_one_letter_code_can") or entity_poly.get("pdbx_seq_one_letter_code") or ""
        sequence = re.sub(r"\s+", "", sequence)
        uniprot_ids = ids.get("uniprot_ids", [])
        rows.append(
            {
                "pdb_id": pdb_id,
                "entity_id": str(entity_id),
                "chains": ",".join(ids.get("auth_asym_ids") or ids.get("asym_ids") or []),
                "sequence": sequence,
                "sequence_length": len(sequence),
                "uniprot_ids": ";".join(uniprot_ids),
                "organism": source.get("scientific_name")
                or source.get("ncbi_scientific_name")
                or source.get("pdbx_organism_scientific", ""),
            }
        )
    return rows


def write_csv(path, rows, fields):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    seen = set()
    with open(path, "w") as handle:
        for row in rows:
            key = (row["uniprot_ids"], row["sequence"])
            if key in seen or not row["sequence"]:
                continue
            seen.add(key)
            header = (
                f">{row['pdb_id']}_{row['entity_id']}|chains={row['chains']}"
                f"|uniprot={row['uniprot_ids']}|organism={row['organism']}"
                f"|length={row['sequence_length']}"
            )
            handle.write(header + "\n")
            handle.write("\n".join(textwrap.wrap(row["sequence"], 80)) + "\n")


def update_readme(path, cutoff, min_length):
    block = (
        "## Geometry Filter\n\n"
        f"Hard rule: keep an entry only if an organic ligand fluorine atom is within <= {cutoff} A "
        "of a FE/CO/NI/CU/MN metal atom in the mmCIF coordinates. The RCSB API result is only a "
        "metadata/substructure candidate set; this geometry CSV is the evidence layer.\n\n"
        f"FASTA rule: write only protein polymer entities with sequence length >= {min_length} aa to "
        "`*_hmmer_full_length_seed.fasta`; shorter fragments/domains are recorded in "
        "`*_protein_entities_excluded.csv` and should not be used as HMMER seeds by default.\n"
    )
    if os.path.exists(path):
        text = open(path, encoding="utf-8", errors="replace").read()
        text = re.split(r"\n## Geometry Filter\n", text, maxsplit=1)[0].rstrip() + "\n\n"
    else:
        text = ""
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text + block)


def candidate_ids(name_dir):
    path = os.path.join(name_dir, f"{os.path.basename(name_dir)}_pdb_ids.csv")
    ids = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            ids.append(row["pdb_id"].upper())
    return ids


def process_name(name, project, cutoff, min_length):
    outdir = os.path.join(project, f"Metal-{name}")
    cif_dir = os.path.join(outdir, "cif")
    ids = candidate_ids(outdir)
    summary = []
    literature = []
    fasta_rows = []
    excluded_fasta = []
    passed_ids = []
    metadata_cache = {}

    for pdb_id in ids:
        cif = download_cif(pdb_id, cif_dir)
        pairs, nearest, ligand_f, metals = geometry_hits(cif, cutoff)
        if not pairs:
            continue
        passed_ids.append(pdb_id)
        meta = metadata_cache.get(pdb_id)
        if meta is None:
            meta = rcsb_entry_metadata(pdb_id)
            metadata_cache[pdb_id] = meta
        entities = rcsb_polymer_entities(pdb_id)
        good_entities = [row for row in entities if row["sequence_length"] >= min_length]
        excluded_fasta.extend(
            [{**row, "exclude_reason": f"sequence_length<{min_length}"} for row in entities if row["sequence_length"] < min_length]
        )
        fasta_rows.extend(good_entities)
        uniprot_ids = sorted({u for row in good_entities for u in row["uniprot_ids"].split(";") if u})
        organisms = sorted({row["organism"] for row in good_entities if row["organism"]})
        pair_text = ";".join(f"{d:.3f}:{atom_label(f_atom)}--{atom_label(metal)}" for d, f_atom, metal in pairs[:20])
        nearest_text = nearest[0] if nearest else ""
        summary.append(
            {
                **meta,
                "motif_name": name,
                "cutoff_A": cutoff,
                "nearest_ligandF_metal_A": f"{nearest_text:.3f}" if nearest else "",
                "passing_pairs_le_3p5A": pair_text,
                "organic_f_atom_count": len(ligand_f),
                "metal_atom_count": len(metals),
                "uniprot_ids": ";".join(uniprot_ids),
                "organisms": ";".join(organisms),
                "kept_fasta_entity_count": len(good_entities),
                "excluded_short_entity_count": len(entities) - len(good_entities),
            }
        )
        literature.append(
            {
                "Title": meta.get("citation_title", ""),
                "Author": meta.get("authors", ""),
                "Publication Year": meta.get("citation_year", ""),
                "DOI": meta.get("doi", ""),
                "Abstract Note": f"PDB {pdb_id}; {name}; ligand-F to transition metal <= {cutoff} A; {meta.get('structure_title', '')}",
                "Publication Title": meta.get("journal", ""),
                "Url": meta.get("pdb_url", ""),
                "Tags": f"PDB;fluorinated-ligand;transition-metal;{name};distance<={cutoff}A",
                "PDB ID": pdb_id,
                "UniProt IDs": ";".join(uniprot_ids),
            }
        )

    write_csv(
        os.path.join(outdir, f"Metal-{name}_geometry_summary.csv"),
        summary,
        [
            "motif_name",
            "pdb_id",
            "cutoff_A",
            "nearest_ligandF_metal_A",
            "passing_pairs_le_3p5A",
            "organic_f_atom_count",
            "metal_atom_count",
            "doi",
            "pubmed_id",
            "citation_title",
            "citation_year",
            "journal",
            "authors",
            "structure_title",
            "pdb_url",
            "uniprot_ids",
            "organisms",
            "kept_fasta_entity_count",
            "excluded_short_entity_count",
        ],
    )
    write_csv(
        os.path.join(outdir, f"Metal-{name}_literature_seed.csv"),
        literature,
        ["Title", "Author", "Publication Year", "DOI", "Abstract Note", "Publication Title", "Url", "Tags", "PDB ID", "UniProt IDs"],
    )
    write_csv(
        os.path.join(outdir, f"Metal-{name}_protein_entities_kept.csv"),
        fasta_rows,
        ["pdb_id", "entity_id", "chains", "sequence_length", "uniprot_ids", "organism", "sequence"],
    )
    write_csv(
        os.path.join(outdir, f"Metal-{name}_protein_entities_excluded.csv"),
        excluded_fasta,
        ["pdb_id", "entity_id", "chains", "sequence_length", "uniprot_ids", "organism", "exclude_reason", "sequence"],
    )
    write_fasta(os.path.join(outdir, f"Metal-{name}_hmmer_full_length_seed.fasta"), fasta_rows)
    with open(os.path.join(outdir, f"Metal-{name}_geometry_pass_pdb_ids.txt"), "w") as handle:
        handle.write("\n".join(passed_ids) + ("\n" if passed_ids else ""))
    update_readme(os.path.join(outdir, "README.md"), cutoff, min_length)
    print(f"Metal-{name}: {len(passed_ids)} / {len(ids)} pass <= {cutoff} A")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", default="CF3")
    parser.add_argument("--project", default=PROJECT)
    parser.add_argument("--cutoff", type=float, default=3.5)
    parser.add_argument("--min-length", type=int, default=150)
    args = parser.parse_args()
    process_name(args.name, args.project, args.cutoff, args.min_length)


if __name__ == "__main__":
    main()
