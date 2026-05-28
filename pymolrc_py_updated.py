from pymol import cmd, stored

import csv
import json
import os
import re
import textwrap
import urllib.parse
import urllib.request


METAL_SELECTION = "resn FE+FE2+FE3+ZN+CU+MN+MG+CA+CO+NI"
FLUORIDE_SELECTION = "resn F and elem F"
DONOR_SELECTION = (
    "name NE2+ND1+OE1+OE2+OD1+OD2+SG+SD+OG+OG1+O"
    " and elem N+O+S"
    " and not resn HOH+WAT"
)
SINGLE_ATOM_ION_SELECTION = (
    "resn F+CL+BR+IOD+NA+K+LI+CS+MG+CA+ZN+FE+FE2+FE3+CU+MN+CO+NI"
)
PROJECT_DIR = "/home/qin/Metal-F_project"


def _has_atoms(selection):
    try:
        return cmd.count_atoms(selection) > 0
    except Exception:
        return False


def _as_bool(value):
    if isinstance(value, str):
        return value.strip().lower() not in ("0", "false", "no", "off", "")
    return bool(value)


def _object_selection(name):
    return f"%{name}"


def _atom_selection(model, index):
    return f"({_object_selection(model)} and index {int(index)})"


def _safe_name(text):
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", str(text)).strip("_")
    if not safe:
        safe = "site"
    if safe[0].isdigit():
        safe = "pdb_" + safe
    return safe


def _resolve_metal_name(name):
    raw = str(name or "CF3").strip().strip("\"'")
    if raw.startswith("Metal-"):
        raw = raw[6:]
    aliases = {
        "cf3": "CF3",
        "trifluoromethyl": "CF3",
        "c-f": "C-F",
        "cf": "C-F",
        "c_f": "C-F",
        "fluoroaryl": "fluoroaryl",
        "arylfluoride": "fluoroaryl",
        "fluoroaryl_ligand": "fluoroaryl",
    }
    return aliases.get(raw.lower(), raw)


def _molecule_objects(objects=""):
    if objects:
        names = [x.strip() for x in str(objects).replace("+", ",").split(",")]
        return [x for x in names if x]

    names = []
    for name in cmd.get_names("objects"):
        try:
            if cmd.get_type(name) == "object:molecule":
                names.append(name)
        except Exception:
            if _has_atoms(_object_selection(name)):
                names.append(name)
    return names


def _clear_scenes():
    for scene_name in list(cmd.get_scene_list()):
        try:
            cmd.scene(scene_name, "clear")
        except Exception:
            pass


def _resolve_scope(scope):
    scope = (scope or "auto").strip()
    lowered = scope.lower()

    if lowered in ("auto", "current", "selection"):
        if _has_atoms("sele"):
            return "sele"

        enabled = [
            name
            for name in _molecule_objects()
            if name in cmd.get_names("objects", enabled_only=1)
        ]
        if len(enabled) == 1:
            return _object_selection(enabled[0])

        if _has_atoms("visible"):
            return "visible"
        return "all"

    if lowered in ("camera", "view", "visible"):
        return "visible"

    object_names = set(_molecule_objects())
    if scope in object_names:
        return _object_selection(scope)

    return scope


def _style_distance_object(name):
    # Keep PyMOL's familiar distance style: yellow dashed line plus label.
    cmd.set("dash_color", "yellow", name)
    cmd.set("dash_radius", 0.03, name)
    cmd.set("dash_gap", 0.35, name)
    cmd.set("label_color", "yellow", name)
    cmd.show("dashes", name)
    cmd.show("labels", name)


def _atom_records(selection):
    stored._pymolrc_atom_records = []
    cmd.iterate_state(
        1,
        selection,
        "stored._pymolrc_atom_records.append((model, index, resn, resi, name))",
    )
    return list(stored._pymolrc_atom_records)


def _residue_records(selection):
    stored._pymolrc_residue_records = []
    cmd.iterate(
        selection,
        "stored._pymolrc_residue_records.append((model, chain, resn, resi, segi))",
    )
    seen = set()
    records = []
    for record in stored._pymolrc_residue_records:
        if record in seen:
            continue
        seen.add(record)
        records.append(record)
    return records


def _residue_label(record):
    model, chain, resn, resi, segi = record
    chain_label = chain or segi or "-"
    return f"{model}/{chain_label}/{resn}{resi}"


def _complete_sidechain_selection(residue_scope):
    # Keep CA as an anchor so long side chains do not look chopped off at the backbone.
    return f"({residue_scope}) and not name N+C+O+OXT"


def _nearby_residue_selection(seed_selection, scope, around):
    return f"byres ((polymer.protein and ({scope})) within {float(around)} of ({seed_selection}))"


def _close_pairs(selection_a, selection_b, cutoff):
    cutoff = float(cutoff)
    pairs = []
    atoms_a = _atom_records(selection_a)
    atoms_b = _atom_records(selection_b)

    for atom_a in atoms_a:
        sel_a = _atom_selection(atom_a[0], atom_a[1])
        for atom_b in atoms_b:
            if atom_a[0] == atom_b[0] and atom_a[1] == atom_b[1]:
                continue
            sel_b = _atom_selection(atom_b[0], atom_b[1])
            try:
                distance = cmd.get_distance(sel_a, sel_b)
            except Exception:
                continue
            if distance <= cutoff:
                pairs.append((distance, atom_a, atom_b))

    pairs.sort(key=lambda row: row[0])
    return pairs


def _nearest_pair(selection_a, selection_b):
    best = None
    atoms_a = _atom_records(selection_a)
    atoms_b = _atom_records(selection_b)

    for atom_a in atoms_a:
        sel_a = _atom_selection(atom_a[0], atom_a[1])
        for atom_b in atoms_b:
            if atom_a[0] == atom_b[0] and atom_a[1] == atom_b[1]:
                continue
            sel_b = _atom_selection(atom_b[0], atom_b[1])
            try:
                distance = cmd.get_distance(sel_a, sel_b)
            except Exception:
                continue
            if best is None or distance < best[0]:
                best = (distance, atom_a, atom_b)

    return best


def _atom_label(atom):
    if not atom:
        return ""
    model, index, resn, resi, name = atom
    return f"{model}/{resn}{resi}/{name}:{index}"


def _read_url(url, data=None, timeout=120):
    headers = {"User-Agent": "pymol-metal-f-scene-builder/1.0"}
    if data is not None:
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(url, data=data, headers=headers)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def _payload_from_rcsb_search_url(search_url):
    parsed = urllib.parse.urlparse(search_url)
    params = urllib.parse.parse_qs(parsed.query)
    if "request" not in params:
        raise ValueError("RCSB search URL does not contain a request= parameter.")

    payload = json.loads(params["request"][0])
    payload.pop("request_info", None)
    payload.setdefault("request_options", {})
    payload["request_options"]["results_content_type"] = ["experimental"]
    payload["request_options"]["paginate"] = {"start": 0, "rows": 1000}
    payload["return_type"] = "entry"
    return payload


def _rcsb_search_ids(search_url):
    payload = _payload_from_rcsb_search_url(search_url)
    endpoint = "https://search.rcsb.org/rcsbsearch/v2/query"
    rows = 1000
    start = 0
    pdb_ids = []

    while True:
        payload["request_options"]["paginate"] = {"start": start, "rows": rows}
        data = json.dumps(payload).encode("utf-8")
        result = json.loads(_read_url(endpoint, data=data).decode("utf-8"))
        result_set = result.get("result_set", [])
        pdb_ids.extend(row["identifier"] for row in result_set)

        total = int(result.get("total_count", len(pdb_ids)))
        start += rows
        if start >= total or not result_set:
            break

    return pdb_ids


def _download_rcsb_cif(pdb_id, cif_dir, overwrite=False):
    cif_dir = os.path.expanduser(cif_dir)
    os.makedirs(cif_dir, exist_ok=True)
    pdb_id = pdb_id.upper()
    path = os.path.join(cif_dir, f"{pdb_id}.cif")

    if os.path.exists(path) and os.path.getsize(path) > 0 and not overwrite:
        return path

    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    content = _read_url(url)
    with open(path, "wb") as handle:
        handle.write(content)
    return path


def _safe_get(mapping, *keys, default=""):
    current = mapping
    for key in keys:
        if not isinstance(current, dict):
            return default
        current = current.get(key)
        if current is None:
            return default
    return current


def _rcsb_entry_metadata(pdb_id):
    pdb_id = pdb_id.upper()
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    data = json.loads(_read_url(url).decode("utf-8"))
    citation = data.get("rcsb_primary_citation") or {}
    if not citation and data.get("citation"):
        primary = [row for row in data["citation"] if row.get("rcsb_is_primary") == "Y"]
        citation = primary[0] if primary else data["citation"][0]

    return {
        "pdb_id": pdb_id,
        "pdb_url": f"https://www.rcsb.org/structure/{pdb_id}",
        "structure_title": _safe_get(data, "struct", "title"),
        "doi": citation.get("pdbx_database_id_DOI", ""),
        "pubmed_id": str(citation.get("pdbx_database_id_PubMed", "") or ""),
        "citation_title": citation.get("title", ""),
        "citation_year": str(citation.get("year", "") or ""),
        "journal": citation.get("journal_abbrev") or citation.get("rcsb_journal_abbrev", ""),
        "authors": "; ".join(citation.get("rcsb_authors", [])[:8]),
    }


def _rcsb_polymer_entities(pdb_id):
    pdb_id = pdb_id.upper()
    entry = json.loads(
        _read_url(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}").decode("utf-8")
    )
    entity_ids = _safe_get(
        entry,
        "rcsb_entry_container_identifiers",
        "polymer_entity_ids",
        default=[],
    )
    rows = []
    for entity_id in entity_ids:
        try:
            data = json.loads(
                _read_url(
                    f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
                ).decode("utf-8")
            )
        except Exception:
            continue

        entity_poly = data.get("entity_poly", {})
        if entity_poly.get("rcsb_entity_polymer_type") != "Protein":
            continue

        identifiers = data.get("rcsb_polymer_entity_container_identifiers", {})
        source = {}
        sources = data.get("rcsb_entity_source_organism") or data.get("entity_src_nat") or []
        if sources:
            source = sources[0]

        sequence = entity_poly.get("pdbx_seq_one_letter_code_can") or entity_poly.get(
            "pdbx_seq_one_letter_code", ""
        )
        sequence = re.sub(r"\s+", "", sequence)
        uniprot_ids = identifiers.get("uniprot_ids", [])
        if not uniprot_ids:
            refs = identifiers.get("reference_sequence_identifiers", [])
            uniprot_ids = [
                row.get("database_accession", "")
                for row in refs
                if row.get("database_name") == "UniProt" and row.get("database_accession")
            ]

        rows.append(
            {
                "pdb_id": pdb_id,
                "entity_id": str(entity_id),
                "chains": ",".join(identifiers.get("auth_asym_ids") or identifiers.get("asym_ids") or []),
                "sequence": sequence,
                "sequence_length": str(len(sequence)) if sequence else "",
                "uniprot_ids": ";".join(uniprot_ids),
                "organism": source.get("scientific_name")
                or source.get("ncbi_scientific_name")
                or source.get("pdbx_organism_scientific", ""),
                "taxonomy_id": str(
                    source.get("ncbi_taxonomy_id") or source.get("pdbx_ncbi_taxonomy_id") or ""
                ),
                "embl_ena_query": (
                    "https://www.ebi.ac.uk/ena/browser/search?query="
                    + urllib.parse.quote(" OR ".join(uniprot_ids))
                    if uniprot_ids
                    else ""
                ),
            }
        )
    return rows


def _rcsb_metadata_bundle(pdb_id):
    meta = _rcsb_entry_metadata(pdb_id)
    entities = _rcsb_polymer_entities(pdb_id)
    meta["protein_entities"] = entities
    meta["protein_entity_ids"] = ";".join(row["entity_id"] for row in entities)
    meta["uniprot_ids"] = ";".join(
        sorted(
            {
                token
                for row in entities
                for token in row.get("uniprot_ids", "").split(";")
                if token
            }
        )
    )
    meta["organisms"] = ";".join(
        sorted({row.get("organism", "") for row in entities if row.get("organism")})
    )
    meta["embl_ena_queries"] = ";".join(
        sorted({row.get("embl_ena_query", "") for row in entities if row.get("embl_ena_query")})
    )
    return meta


def _write_fasta(rows, path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as handle:
        for row in rows:
            sequence = row.get("sequence", "")
            if not sequence:
                continue
            header = (
                f">{row['pdb_id']}_{row['entity_id']}|chains={row.get('chains', '')}"
                f"|uniprot={row.get('uniprot_ids', '')}|organism={row.get('organism', '')}"
            )
            handle.write(header + "\n")
            handle.write("\n".join(textwrap.wrap(sequence, 80)) + "\n")


def _write_csv(rows, path, fields):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


@cmd.extend
def simple_hetero(selection="organic or inorganic"):
    """Show multi-atom hetero groups as preset-simple-like sticks."""
    target = f"({selection}) and not solvent"
    sticks = f"({target}) and not ({SINGLE_ATOM_ION_SELECTION})"
    ions = f"({target}) and ({SINGLE_ATOM_ION_SELECTION})"

    if _has_atoms(sticks):
        cmd.hide("spheres", sticks)
        cmd.show("sticks", sticks)
    if _has_atoms(ions):
        cmd.show("spheres", ions)
        cmd.set("sphere_scale", 0.25, ions)

    print(f"[simple_hetero] Styled hetero atoms in: {selection}")


@cmd.extend
def metal_coord(metal_sele=METAL_SELECTION, cutoff=2.8, prefix="coord", scope="auto", clear=1):
    """Draw yellow labelled dashed metal-donor distances in the selected scope."""
    cutoff = float(cutoff)
    clear = _as_bool(clear)
    scope = _resolve_scope(scope)
    scoped_metals = f"({metal_sele}) and ({scope})"

    if clear:
        cmd.delete(f"{prefix}_*")

    metals = _atom_records(scoped_metals)
    if not metals:
        print(f"[metal_coord] No metal ions found in scope: {scope}")
        return 0

    count = 0
    for metal in metals:
        metal_atom = _atom_selection(metal[0], metal[1])
        donor_query = (
            f"(({DONOR_SELECTION}) and ({scope}) and not ({metal_sele}))"
            f" within {cutoff} of {metal_atom}"
        )
        donors = _atom_records(donor_query)

        for donor in donors:
            donor_atom = _atom_selection(donor[0], donor[1])
            distance_name = f"{prefix}_{count}"
            distance = cmd.distance(distance_name, metal_atom, donor_atom, cutoff)
            if distance > 0:
                _style_distance_object(distance_name)
                count += 1

    print(
        f"[metal_coord] {count} labelled distance(s) drawn for "
        f"{len(metals)} metal ion(s) in scope: {scope}"
    )
    return count


@cmd.extend
def metal_f_dist(
    f_sele=FLUORIDE_SELECTION,
    metal_sele=METAL_SELECTION,
    cutoff=4.0,
    prefix="f_metal",
    scope="auto",
    clear=1,
):
    """Draw F-metal distances below cutoff in the selected scope."""
    cutoff = float(cutoff)
    clear = _as_bool(clear)
    scope = _resolve_scope(scope)

    if clear:
        cmd.delete(f"{prefix}_*")

    f_scoped = f"({f_sele}) and ({scope})"
    metal_scoped = f"({metal_sele}) and ({scope})"
    pairs = _close_pairs(f_scoped, metal_scoped, cutoff)

    for i, (distance, fluorine, metal) in enumerate(pairs):
        name = f"{prefix}_{i}"
        cmd.distance(name, _atom_selection(fluorine[0], fluorine[1]), _atom_selection(metal[0], metal[1]))
        _style_distance_object(name)

    print(f"[metal_f_dist] {len(pairs)} F-metal distance(s) <= {cutoff:.2f} A in scope: {scope}")
    return len(pairs)


@cmd.extend
def show_metal(scope="auto", around=4.0, cutoff=2.8, prefix="coord"):
    """Show scoped metal sites with nearby side chains and labelled distances."""
    around = float(around)
    scope = _resolve_scope(scope)
    scoped_metals = f"({METAL_SELECTION}) and ({scope})"

    if not _has_atoms(scoped_metals):
        print(f"[show_metal] No metal ions found in scope: {scope}")
        return 0

    cmd.show("spheres", scoped_metals)
    cmd.set("sphere_scale", 0.25, scoped_metals)

    nearby_residues = _nearby_residue_selection(scoped_metals, scope, around)
    sidechains = _complete_sidechain_selection(nearby_residues)
    if _has_atoms(sidechains):
        cmd.show("sticks", sidechains)

    nearby_hetero = f"byres (((organic or inorganic) and not solvent) within {around} of ({scoped_metals}))"
    simple_hetero(nearby_hetero)

    count = metal_coord(scoped_metals, cutoff, prefix, scope, 1)
    print(f"[show_metal] Shown metals plus side chains/hetero groups within {around:.1f} A.")
    return count


@cmd.extend
def show_metal_f(scope="auto", around=4.0, f_cutoff=4.0, coord_cutoff=2.8, prefix="mf"):
    """Show F-metal site, nearby side chains, and F-metal/metal-donor distances."""
    scope = _resolve_scope(scope)
    show_metal(scope, around, coord_cutoff, f"{prefix}_coord")

    f_scoped = f"({FLUORIDE_SELECTION}) and ({scope})"
    if _has_atoms(f_scoped):
        cmd.show("spheres", f_scoped)
        cmd.set("sphere_scale", 0.25, f_scoped)

    return metal_f_dist(
        FLUORIDE_SELECTION,
        METAL_SELECTION,
        f_cutoff,
        f"{prefix}_fmetal",
        scope,
        1,
    )


@cmd.extend
def build_metal_f_scenes(
    session="/home/qin/Metal-F_project/Metal-F.pse",
    output="/home/qin/Metal-F_project/Metal-F_auto_scenes.pse",
    f_cutoff=4.0,
    coord_cutoff=2.8,
    around=4.0,
    objects="",
    only_close=0,
    png_dir="",
    ray=0,
    csv_path="",
    metadata=None,
    fasta_path="",
    entity_csv_path="",
    literature_csv_path="",
):
    """Build one scene per molecule object containing close F-metal pairs."""
    f_cutoff = float(f_cutoff)
    coord_cutoff = float(coord_cutoff)
    around = float(around)
    only_close = _as_bool(only_close)
    ray = int(ray)

    if session and str(session).lower() not in ("0", "none", "loaded"):
        cmd.delete("all")
        _clear_scenes()
        cmd.load(session)
    else:
        _clear_scenes()

    cmd.delete("mf_*")
    names = _molecule_objects(objects)
    made = []
    summary = []
    metadata = metadata or {}
    entity_rows_all = []
    literature_rows = []

    if png_dir:
        os.makedirs(os.path.expanduser(png_dir), exist_ok=True)

    for obj in names:
        obj_sel = _object_selection(obj)
        if not _has_atoms(obj_sel):
            continue
        pdb_id = obj.replace("pdb_", "", 1).upper()
        entry_meta = metadata.get(pdb_id, {})

        f_sel = f"({FLUORIDE_SELECTION}) and ({obj_sel})"
        metal_sel = f"({METAL_SELECTION}) and ({obj_sel})"
        pairs = _close_pairs(f_sel, metal_sel, f_cutoff)
        nearest = pairs[0] if pairs else _nearest_pair(f_sel, metal_sel)
        site_core = f"(({f_sel}) or ({metal_sel}))"
        site_residues = _nearby_residue_selection(site_core, obj_sel, around)
        metal_residues = _nearby_residue_selection(metal_sel, obj_sel, around)
        f_residues = _nearby_residue_selection(f_sel, obj_sel, around)
        coord_residues = _nearby_residue_selection(metal_sel, obj_sel, coord_cutoff)
        site_residue_labels = [_residue_label(row) for row in _residue_records(site_residues)]
        metal_residue_labels = [_residue_label(row) for row in _residue_records(metal_residues)]
        f_residue_labels = [_residue_label(row) for row in _residue_records(f_residues)]
        coord_residue_labels = [_residue_label(row) for row in _residue_records(coord_residues)]
        pair_labels = [
            f"{distance:.3f}:{_atom_label(fluorine)}--{_atom_label(metal)}"
            for distance, fluorine, metal in pairs
        ]
        if only_close and not pairs:
            print(f"[build_metal_f_scenes] skip {obj}: no F-metal pair <= {f_cutoff:.2f} A")
            if nearest:
                summary.append(
                    {
                        "object": obj,
                        "scene": "",
                        "nearest_f_metal_A": f"{nearest[0]:.3f}",
                        "within_cutoff": "0",
                        "fluoride_atom": _atom_label(nearest[1]),
                        "metal_atom": _atom_label(nearest[2]),
                        "site_residues_4A": ";".join(site_residue_labels),
                        "metal_neighbor_residues_4A": ";".join(metal_residue_labels),
                        "fluoride_neighbor_residues_4A": ";".join(f_residue_labels),
                        "metal_coord_residues": ";".join(coord_residue_labels),
                        "f_metal_pairs_within_cutoff": ";".join(pair_labels),
                        **entry_meta,
                    }
                )
            continue

        safe = _safe_name(obj)
        cmd.disable("all")
        cmd.enable(obj)
        cmd.hide("everything", "all")
        cmd.show("cartoon", f"({obj_sel}) and polymer.protein")

        sidechains = _complete_sidechain_selection(site_residues)
        if _has_atoms(sidechains):
            cmd.show("sticks", sidechains)
            cmd.color("gray70", f"({sidechains}) and elem C")

        hetero_near_site = (
            f"byres (((organic or inorganic) and not solvent and {obj_sel}) "
            f"within {around} of ({site_core}))"
        )
        simple_hetero(hetero_near_site)

        if _has_atoms(f_sel):
            cmd.show("spheres", f_sel)
            cmd.set("sphere_scale", 0.28, f_sel)
        if _has_atoms(metal_sel):
            cmd.show("spheres", metal_sel)
            cmd.set("sphere_scale", 0.28, metal_sel)

        show_metal(obj_sel, around, coord_cutoff, f"mf_{safe}_coord")
        metal_f_dist(FLUORIDE_SELECTION, METAL_SELECTION, f_cutoff, f"mf_{safe}_fmetal", obj_sel, 0)

        focus = f"({site_core}) expand 6"
        if _has_atoms(focus):
            cmd.orient(focus)
            cmd.zoom(focus, 3.5)

        scene_name = f"MF_{safe}"
        cmd.scene(scene_name, "store")
        made.append(scene_name)
        if nearest:
            summary.append(
                {
                    "object": obj,
                    "scene": scene_name,
                    "nearest_f_metal_A": f"{nearest[0]:.3f}",
                    "within_cutoff": "1" if pairs else "0",
                    "fluoride_atom": _atom_label(nearest[1]),
                    "metal_atom": _atom_label(nearest[2]),
                    "site_residues_4A": ";".join(site_residue_labels),
                    "metal_neighbor_residues_4A": ";".join(metal_residue_labels),
                    "fluoride_neighbor_residues_4A": ";".join(f_residue_labels),
                    "metal_coord_residues": ";".join(coord_residue_labels),
                    "f_metal_pairs_within_cutoff": ";".join(pair_labels),
                    **entry_meta,
                }
            )

        if png_dir:
            png_path = os.path.join(os.path.expanduser(png_dir), f"{scene_name}.png")
            cmd.png(png_path, width=1800, height=1400, dpi=300, ray=ray)

        min_distance = nearest[0] if nearest else None
        if min_distance is None:
            print(f"[build_metal_f_scenes] scene {scene_name}: no close pair filter applied")
        else:
            print(f"[build_metal_f_scenes] scene {scene_name}: nearest F-metal {min_distance:.2f} A")

    if made:
        cmd.set("scene_buttons", 1)
        cmd.disable("all")
        for name in cmd.get_names("objects"):
            cmd.disable(name)
        cmd.scene(made[0], "recall")
        if output:
            cmd.save(output)

    if csv_path and summary:
        csv_path = os.path.expanduser(csv_path)
        os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)
        fields = [
            "pdb_id",
            "object",
            "scene",
            "nearest_f_metal_A",
            "within_cutoff",
            "fluoride_atom",
            "metal_atom",
            "f_metal_pairs_within_cutoff",
            "site_residues_4A",
            "metal_neighbor_residues_4A",
            "fluoride_neighbor_residues_4A",
            "metal_coord_residues",
            "doi",
            "pubmed_id",
            "citation_title",
            "citation_year",
            "journal",
            "authors",
            "structure_title",
            "pdb_url",
            "protein_entity_ids",
            "uniprot_ids",
            "organisms",
            "embl_ena_queries",
        ]
        _write_csv(summary, csv_path, fields)
        print(f"[build_metal_f_scenes] summary CSV: {csv_path}")

    for pdb_id, meta in sorted(metadata.items()):
        entity_rows_all.extend(meta.get("protein_entities", []))
        literature_rows.append(
            {
                "Title": meta.get("citation_title", ""),
                "Author": meta.get("authors", ""),
                "Publication Year": meta.get("citation_year", ""),
                "DOI": meta.get("doi", ""),
                "Abstract Note": (
                    f"PDB {pdb_id}; F-metal structural candidate; "
                    f"{meta.get('structure_title', '')}"
                ),
                "Publication Title": meta.get("journal", ""),
                "Url": meta.get("pdb_url", ""),
                "Tags": "PDB;metal-fluoride;structure-seed",
                "PDB ID": pdb_id,
                "UniProt IDs": meta.get("uniprot_ids", ""),
            }
        )

    if entity_rows_all and fasta_path:
        _write_fasta(entity_rows_all, os.path.expanduser(fasta_path))
        print(f"[build_metal_f_scenes] FASTA: {fasta_path}")

    if entity_rows_all and entity_csv_path:
        _write_csv(
            entity_rows_all,
            os.path.expanduser(entity_csv_path),
            [
                "pdb_id",
                "entity_id",
                "chains",
                "sequence_length",
                "uniprot_ids",
                "organism",
                "taxonomy_id",
                "embl_ena_query",
                "sequence",
            ],
        )
        print(f"[build_metal_f_scenes] protein entity CSV: {entity_csv_path}")

    if literature_rows and literature_csv_path:
        _write_csv(
            literature_rows,
            os.path.expanduser(literature_csv_path),
            [
                "Title",
                "Author",
                "Publication Year",
                "DOI",
                "Abstract Note",
                "Publication Title",
                "Url",
                "Tags",
                "PDB ID",
                "UniProt IDs",
            ],
        )
        print(f"[build_metal_f_scenes] literature seed CSV: {literature_csv_path}")

    print(f"[build_metal_f_scenes] built {len(made)} scene(s): {', '.join(made)}")
    if output:
        print(f"[build_metal_f_scenes] saved: {output}")
    return len(made)


@cmd.extend
def build_metal_f_from_rcsb(
    search_url,
    output=PROJECT_DIR + "/Metal-F_auto_scenes.pse",
    cif_dir=PROJECT_DIR + "/cif",
    f_cutoff=4.0,
    coord_cutoff=2.8,
    around=4.0,
    only_close=0,
    png_dir="",
    csv_path=PROJECT_DIR + "/Metal-F_scene_summary.csv",
    fasta_path=PROJECT_DIR + "/Metal-F_protein_entities.fasta",
    entity_csv_path=PROJECT_DIR + "/Metal-F_protein_entities.csv",
    literature_csv_path=PROJECT_DIR + "/Metal-F_literature_seed.csv",
    overwrite=0,
    ray=0,
):
    """Download RCSB hits, load them, and build F-metal scenes automatically."""
    overwrite = _as_bool(overwrite)
    pdb_ids = _rcsb_search_ids(search_url)
    print(f"[build_metal_f_from_rcsb] RCSB returned {len(pdb_ids)} entrie(s): {', '.join(pdb_ids)}")

    cmd.delete("all")
    _clear_scenes()
    loaded = []
    metadata = {}
    for pdb_id in pdb_ids:
        try:
            cif_path = _download_rcsb_cif(pdb_id, cif_dir, overwrite)
            object_name = f"pdb_{pdb_id.lower()}"
            cmd.load(cif_path, object_name)
            loaded.append(object_name)
            metadata[pdb_id.upper()] = _rcsb_metadata_bundle(pdb_id)
            print(f"[build_metal_f_from_rcsb] loaded {pdb_id}: {cif_path}")
        except Exception as exc:
            print(f"[build_metal_f_from_rcsb] failed {pdb_id}: {exc}")

    if not loaded:
        print("[build_metal_f_from_rcsb] No structures loaded.")
        return 0

    return build_metal_f_scenes(
        session="loaded",
        output=output,
        f_cutoff=f_cutoff,
        coord_cutoff=coord_cutoff,
        around=around,
        objects=",".join(loaded),
        only_close=only_close,
        png_dir=png_dir,
        ray=ray,
        csv_path=csv_path,
        metadata=metadata,
        fasta_path=fasta_path,
        entity_csv_path=entity_csv_path,
        literature_csv_path=literature_csv_path,
    )


@cmd.extend
def build_metal_name_scenes(
    name="CF3",
    project=PROJECT,
    cutoff=3.5,
    coord_cutoff=2.8,
    around=4.0,
    output="",
):
    """Build scenes for a Metal-Name geometry-pass set."""
    name = _resolve_metal_name(name)
    project = os.path.expanduser(project)
    outdir = os.path.join(project, f"Metal-{name}")
    ids_path = os.path.join(outdir, f"Metal-{name}_geometry_pass_pdb_ids.txt")
    cif_dir = os.path.join(outdir, "cif")
    if not output:
        output = os.path.join(outdir, f"Metal-{name}_scenes.pse")

    if not os.path.exists(ids_path):
        print(f"[build_metal_name_scenes] Missing geometry pass list: {ids_path}")
        return 0

    with open(ids_path) as handle:
        pdb_ids = [line.strip().upper() for line in handle if line.strip()]
    if not pdb_ids:
        print(f"[build_metal_name_scenes] No geometry-pass PDB IDs for Metal-{name}.")
        return 0

    cmd.delete("all")
    _clear_scenes()
    loaded = []
    for pdb_id in pdb_ids:
        cif_path = os.path.join(cif_dir, f"{pdb_id}.cif")
        if not os.path.exists(cif_path):
            print(f"[build_metal_name_scenes] Missing CIF: {cif_path}")
            continue
        obj = f"pdb_{pdb_id.lower()}"
        cmd.load(cif_path, obj)
        loaded.append(obj)

    if not loaded:
        print("[build_metal_name_scenes] No CIF files loaded.")
        return 0

    return build_metal_f_scenes(
        session="loaded",
        output=output,
        f_cutoff=cutoff,
        coord_cutoff=coord_cutoff,
        around=around,
        objects=",".join(loaded),
        only_close=0,
        png_dir="",
        ray=0,
        csv_path=os.path.join(outdir, f"Metal-{name}_pymol_scene_summary.csv"),
        metadata={},
        fasta_path="",
        entity_csv_path="",
        literature_csv_path="",
    )


@cmd.extend
def metal_cf3():
    """Shortcut: build Metal-CF3 scenes."""
    return build_metal_name_scenes("CF3")


@cmd.extend
def metal_cf():
    """Shortcut: build Metal-C-F scenes."""
    return build_metal_name_scenes("C-F")


@cmd.extend
def metal_fluoroaryl():
    """Shortcut: build Metal-fluoroaryl scenes."""
    return build_metal_name_scenes("fluoroaryl")


cmd.set("dash_color", "yellow")
cmd.set("dash_radius", 0.03)
cmd.set("dash_gap", 0.35)
cmd.set("label_color", "yellow")

print("=" * 50)
print("pymolrc.py loaded | cmds: metal_coord, metal_f_dist, show_metal, show_metal_f, build_metal_f_scenes, build_metal_f_from_rcsb, build_metal_name_scenes, metal_cf3, metal_cf, metal_fluoroaryl")
print("=" * 50)
