from pymol import cmd, stored


METAL_SELECTION = "resn FE+FE2+FE3+ZN+CU+MN+MG+CA+CO+NI"
DONOR_SELECTION = (
    "name NE2+ND1+OE1+OE2+OD1+OD2+SG+SD+OG+OG1+O"
    " and not resn HOH+WAT"
)


def _has_atoms(selection):
    try:
        return cmd.count_atoms(selection) > 0
    except Exception:
        return False


def _resolve_scope(scope):
    scope = (scope or "auto").strip()
    if scope.lower() in ("auto", "current"):
        if _has_atoms("sele"):
            return "sele"
        return "visible"
    return scope


def _atom_selection(model, index):
    return f"({model} and index {index})"


def _style_distance_object(name):
    cmd.set("dash_color", "yellow", name)
    cmd.set("dash_radius", 0.05, name)
    cmd.set("dash_gap", 0.35, name)
    cmd.set("label_color", "yellow", name)
    cmd.show("dashes", name)
    cmd.show("labels", name)


def _simple_nearby_hetero(selection):
    hetero = f"(({selection}) and (organic or inorganic) and not solvent and not ({METAL_SELECTION}))"
    if not _has_atoms(hetero):
        return
    cmd.hide("spheres", hetero)
    cmd.show("sticks", hetero)


@cmd.extend
def simple_hetero(selection="organic or inorganic"):
    """Show phosphate/sulfate/ligands as preset-simple-like sticks, not spheres."""
    target = f"({selection}) and not solvent and not ({METAL_SELECTION})"
    if not _has_atoms(target):
        print("[simple_hetero] No matching hetero atoms found.")
        return
    cmd.hide("spheres", target)
    cmd.show("sticks", target)
    print(f"[simple_hetero] Displayed {cmd.count_atoms(target)} atom(s) as sticks.")


@cmd.extend
def metal_coord(metal_sele=METAL_SELECTION, cutoff=2.8, prefix="coord", scope="auto", clear=1):
    """Draw yellow labelled dashed distances from scoped metals to donor atoms."""
    cutoff = float(cutoff)
    clear = int(clear)
    scope = _resolve_scope(scope)
    scoped_metals = f"({metal_sele}) and ({scope})"

    if clear:
        cmd.delete(f"{prefix}_*")

    stored.metals = []
    cmd.iterate_state(
        1,
        scoped_metals,
        "stored.metals.append((model, index, resi, resn, name))",
    )

    if not stored.metals:
        print(f"[metal_coord] No metal ions found in scope: {scope}")
        return

    count = 0
    for model, m_idx, resi, resn, mname in stored.metals:
        metal_atom = _atom_selection(model, m_idx)
        donor_query = (
            f"(({DONOR_SELECTION}) and ({scope}) and not ({metal_sele}))"
            f" within {cutoff} of {metal_atom}"
        )
        stored.donors = []
        cmd.iterate_state(1, donor_query, "stored.donors.append((model, index))")

        for donor_model, d_idx in stored.donors:
            donor_atom = _atom_selection(donor_model, d_idx)
            distance_name = f"{prefix}_{count}"
            distance = cmd.distance(distance_name, metal_atom, donor_atom, cutoff)
            if distance > 0:
                _style_distance_object(distance_name)
                count += 1

    print(
        f"[metal_coord] {count} labelled distance(s) drawn for "
        f"{len(stored.metals)} metal ion(s) in scope: {scope}"
    )


@cmd.extend
def show_metal(scope="auto", around=4.0, cutoff=2.8, prefix="coord"):
    """Show scoped metal sites with nearby side chains and labelled distances."""
    around = float(around)
    scope = _resolve_scope(scope)
    scoped_metals = f"({METAL_SELECTION}) and ({scope})"

    if not _has_atoms(scoped_metals):
        print(f"[show_metal] No metal ions found in scope: {scope}")
        return

    cmd.show("spheres", scoped_metals)
    cmd.set("sphere_scale", 0.25, scoped_metals)

    nearby_residues = f"byres ((polymer.protein) within {around} of ({scoped_metals}))"
    sidechains = f"({nearby_residues}) and not name N+C+O+CA+OXT"
    if _has_atoms(sidechains):
        cmd.show("sticks", sidechains)

    nearby_hetero = f"byres (((organic or inorganic) and not solvent) within {around} of ({scoped_metals}))"
    _simple_nearby_hetero(nearby_hetero)

    metal_coord(scoped_metals, cutoff, prefix, scope, 1)
    print(f"[show_metal] Shown metals plus side chains/hetero groups within {around} A.")


cmd.set("dash_color", "yellow")
cmd.set("dash_radius", 0.05)
cmd.set("label_color", "yellow")

print("=" * 50)
print("pymolrc.py loaded | cmds: metal_coord, show_metal, simple_hetero")
print("=" * 50)
