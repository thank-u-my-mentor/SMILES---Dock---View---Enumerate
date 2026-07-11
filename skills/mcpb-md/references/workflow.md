# MCPB.py To MD Workflow

## Overall Flow

```text
manual/literature active-site structure
-> clean PDB and ligand names
-> define Fe donors, His protonation, ligand charge/radical/spin
-> prepare mol2/frcmod
-> MCPB.py -s 1 to generate small/large/standard models
-> QM small_opt and small_fc/Hessian
-> QM large_mk or reduced charge model
-> MCPB.py -s 2/3/4
-> tleap/Amber
-> GROMACS conversion
-> EM/NVT/NPT/production
-> PBC repair and Fe-donor/RDC/distance audits
```

## MCPB Models

- `small model`: capped minimal QM model for OPT/Hessian/Seminario bonded terms.
- `large model`: larger capped QM model for MK/RESP electrostatic potential.
- `standard model`: MCPB mapping model for final metal-site mol2 residues; it is not usually run as QM.

The models can differ in atom count because they serve different purposes. Do not directly substitute one model's log into another model's fingerprint without auditing atom order and charge closure.

## Required Checks

After `MCPB.py -s 1`:

```text
small atom count and donor atoms
large atom count and donor atoms
standard fingerprint
charge and multiplicity
His protonation
ligand/radical hydrogens
Fe coordination geometry
```

After OPT:

```text
true OPT convergence, not just normal program termination
GRAD.MAX and RMS below requested tolerance
stable S-SQUARED
no large geometry jumps near Fe donors
```

After Hessian:

```text
all displacement points converged
no SCF IS UNCONVERGED bad point
imaginary modes checked; small peripheral modes may be acceptable for MCPB bonded terms
```

After large_mk or reduced charge model:

```text
SCF converged
ESP/MK points generated
RESP total charge closes
per-residue and per-atom charges are chemically sane
no extreme neighboring +/- charge cancellation
```

After `MCPB.py -s 4` and tleap:

```text
tleap Errors = 0
dry and solv prmtop/inpcrd generated
Fe-donor bonds present
residue names map correctly
internal residue indices may differ from original PDB numbering
```

## Charge Strategy

Prefer final spin-consistent parameters:

```text
M2 Hessian bonded terms + M2 RESP charges
```

Avoid final mixed-spin parameters:

```text
M2 Hessian bonded terms + M6 RESP charges
```

If M2 large_mk SCF does not converge:

1. Try conservative SCF strategies and good `$VEC` history.
2. Try chemically recapped reduced charge models rather than naked atom deletion.
3. Use ChgModD-like fixed backbone/CB logic when mapping small/reduced RESP charges.
4. Document nonstandard charge normalization explicitly.

ChgModB/C/D affects RESP constraints after ESP generation; it does not solve GAMESS SCF.
