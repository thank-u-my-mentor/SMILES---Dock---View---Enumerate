"""Clean PyMOL startup helpers for metal-site visualization.

This file intentionally does not include the old Metal-F/RCSB scene builders.
It keeps only small, robust commands for catalytic metal inspection.
"""

from pymol import cmd


METAL_SELECTION = (
    "(resn FE+FE2+FE3+ZN+CU+MN+MG+CA+CO+NI) "
    "or elem Fe or elem Zn or elem Cu or elem Mn or elem Mg "
    "or elem Ca or elem Co or elem Ni"
)

SINGLE_ATOM_ION_SELECTION = (
    "resn F+CL+BR+IOD+NA+K+LI+CS+MG+CA+ZN+FE+FE2+FE3+CU+MN+CO+NI"
)

SOLVENT_SELECTION = "solvent or resn HOH+WAT+SOL+DOD"

DONOR_SELECTION = (
    "((polymer.protein and "
    "name ND1+NE2+OE1+OE2+OD1+OD2+OG+OG1+SG+SD+OH+OXT "
    "and elem N+O+S) "
    "or (((organic or inorganic) and not solvent) and elem N+O+S)) "
    "and not resn HOH+WAT+SOL+DOD"
)


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


def _molecule_objects():
    names = []
    for name in cmd.get_names("objects"):
        try:
            if cmd.get_type(name) == "object:molecule":
                names.append(name)
        except Exception:
            if _has_atoms(_object_selection(name)):
                names.append(name)
    return names


def _resolve_scope(scope):
    scope = str(scope or "auto").strip()
    lowered = scope.lower()

    if lowered in ("auto", "current", "selection"):
        if _has_atoms(f"({METAL_SELECTION}) and sele"):
            return "sele"

        enabled = [
            name
            for name in _molecule_objects()
            if name in cmd.get_names("objects", enabled_only=1)
        ]
        if len(enabled) == 1:
            return _object_selection(enabled[0])

        if _has_atoms(f"({METAL_SELECTION}) and visible"):
            return "visible"
        return "all"

    if lowered in ("camera", "view", "visible"):
        return "visible"

    if scope in set(_molecule_objects()):
        return _object_selection(scope)

    return scope


def _atom_records(selection):
    try:
        model = cmd.get_model(selection, state=1)
    except Exception:
        return []
    return list(model.atom)


def _atom_selection(atom):
    return f"({_object_selection(atom.model)} and index {int(atom.index)})"


def _distance(sel_a, sel_b):
    try:
        return float(cmd.get_distance(sel_a, sel_b, state=1))
    except Exception:
        return 0.0


def _style_distance_object(name):
    cmd.set("dash_color", "yellow", name)
    cmd.set("dash_radius", 0.03, name)
    cmd.set("dash_gap", 0.35, name)
    cmd.set("label_color", "yellow", name)
    cmd.show("dashes", name)
    cmd.show("labels", name)


def _color_carbons(selection, carbon_color):
    try:
        cmd.color(carbon_color, f"({selection}) and elem C")
    except Exception:
        cmd.color("green", f"({selection}) and elem C")


def simple_hetero(selection="organic or inorganic", carbon_color="palegreen"):
    """Show hetero groups as sticks, while leaving N/O/S/metal element colors alone."""
    target = f"({selection}) and not ({SOLVENT_SELECTION})"
    sticks = f"({target}) and not ({SINGLE_ATOM_ION_SELECTION})"
    ions = f"({target}) and ({SINGLE_ATOM_ION_SELECTION})"

    if _has_atoms(sticks):
        cmd.hide("spheres", sticks)
        cmd.show("sticks", sticks)
        _color_carbons(sticks, carbon_color)

    if _has_atoms(ions):
        cmd.show("spheres", ions)
        cmd.set("sphere_scale", 0.35, ions)

    print(f"[simple_hetero] shown={cmd.count_atoms(target)} selection=({selection})")
    return cmd.count_atoms(target)


def metal_coord(metal_selection=METAL_SELECTION, cutoff=2.8, prefix="coord", scope="auto", clear=1):
    """Draw yellow Fe/Zn/etc.-donor distances without restricting donors to sele."""
    cutoff = float(cutoff)
    scope_sel = _resolve_scope(scope)
    metals = _atom_records(f"({metal_selection}) and ({scope_sel})")

    if clear:
        cmd.delete(f"{prefix}_*")

    if not metals and scope_sel == "sele":
        scope_sel = "visible" if _has_atoms(f"({METAL_SELECTION}) and visible") else "all"
        metals = _atom_records(f"({metal_selection}) and ({scope_sel})")

    if not metals:
        print(f"[metal_coord] no metal atoms found in scope: {scope_sel}")
        return 0

    count = 0
    for metal in metals:
        metal_atom = _atom_selection(metal)
        object_scope = _object_selection(metal.model)
        donor_query = (
            f"(({DONOR_SELECTION}) and ({object_scope}) and not ({METAL_SELECTION})) "
            f"within {cutoff} of ({metal_atom})"
        )
        donors = _atom_records(donor_query)

        for donor in donors:
            donor_atom = _atom_selection(donor)
            distance = _distance(metal_atom, donor_atom)
            if distance <= 0 or distance > cutoff:
                continue
            name = f"{prefix}_{count}_{donor.resn}{donor.resi}_{donor.name}"
            cmd.distance(name, metal_atom, donor_atom)
            _style_distance_object(name)
            count += 1

    print(f"[metal_coord] distances={count} metals={len(metals)} scope={scope_sel}")
    return count


def show_metal(scope="auto", around=4.0, cutoff=2.8, prefix="coord", carbon_color="palegreen"):
    """Show a metal center, nearby side chains/hetero groups, and donor distances.

    Usage examples:
        show_metal
        show_metal sele
        show_metal all, 4.0, 2.8
    """
    around = float(around)
    cutoff = float(cutoff)
    scope_sel = _resolve_scope(scope)
    metal_sel = f"({METAL_SELECTION}) and ({scope_sel})"
    metals = _atom_records(metal_sel)

    if not metals and scope_sel == "sele":
        scope_sel = "visible" if _has_atoms(f"({METAL_SELECTION}) and visible") else "all"
        metal_sel = f"({METAL_SELECTION}) and ({scope_sel})"
        metals = _atom_records(metal_sel)

    if not metals:
        print(f"[show_metal] no metal atoms found in scope: {scope_sel}")
        return 0

    object_scope = " or ".join(sorted({_object_selection(atom.model) for atom in metals}))
    nearby_residues = (
        f"byres ((polymer.protein and ({object_scope})) within {around} of ({metal_sel}))"
    )
    sidechains = f"({nearby_residues}) and not name N+C+O+OXT"
    nearby_hetero = (
        f"byres (((organic or inorganic) and not ({SOLVENT_SELECTION}) "
        f"and ({object_scope}) and not ({METAL_SELECTION})) within {around} of ({metal_sel}))"
    )

    cmd.show("spheres", metal_sel)
    cmd.set("sphere_scale", 0.45, metal_sel)

    if _has_atoms(sidechains):
        cmd.show("sticks", sidechains)
        _color_carbons(sidechains, carbon_color)

    simple_hetero(nearby_hetero, carbon_color)
    distance_count = metal_coord(metal_sel, cutoff, prefix, object_scope, 1)

    cmd.zoom(f"({metal_sel}) or ({sidechains}) or ({nearby_hetero})", 8)
    cmd.deselect()

    print(
        f"[show_metal] metals={len(metals)} sidechain_atoms={cmd.count_atoms(sidechains)} "
        f"hetero_atoms={cmd.count_atoms(nearby_hetero)} distances={distance_count}"
    )
    return distance_count


cmd.extend("simple_hetero", simple_hetero)
cmd.extend("metal_coord", metal_coord)
cmd.extend("show_metal", show_metal)
cmd.extend("ple_show_metal", show_metal)

cmd.set("dash_color", "yellow")
cmd.set("dash_radius", 0.03)
cmd.set("dash_gap", 0.35)
cmd.set("label_color", "yellow")

print("[pymolrc] clean metal helpers loaded: show_metal, ple_show_metal, metal_coord, simple_hetero")
