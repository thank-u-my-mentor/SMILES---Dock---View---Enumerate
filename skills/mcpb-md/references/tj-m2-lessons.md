# TJ Fe/M2 Project Lessons

## Successful Evidence Chain

The successful M2 small-model OPT/Hessian chain is:

```text
TJ2unlNH_M2_opt_from_smallstep60_punchhess_nstep40
NSERCH=12
E=-3001.8042149339
GRAD.MAX=0.0001972
RMS=0.0000509
S-SQUARED about 0.778
```

Successful Hessian:

```text
TJ2unlNH_M2_hessian_from_punchhess_opt12_diis_shift
No bad displacement points
Two very small peripheral imaginary modes
Accepted for MCPB/Seminario bonded terms
```

## Important Lessons

- Reusing good `$VEC`, `$GRAD`, and `$HESS` history can rescue difficult OPT; do not discard good optimizer history casually.
- A single low GRAD.MAX point can be a false local step if subsequent steps rebound. Prefer smooth trajectories such as old NSERCH 57-60 when they show stable energy/gradient behavior.
- Normal GAMESS termination does not mean OPT convergence.
- Hessian bad points are displacement points where SCF fails or gives zero/garbage energy; they can contaminate force constants.
- For Hessian displacement points, stronger DIIS/damping/shift is often safer than aggressive SOSCF.
- Use full/fine DFT grid consistently when possible for final-quality metal-center calculations.
- `mcpb.in` `smmodel_spin` and `lgmodel_spin` must be checked every time. A stale M6 setting can silently regenerate mixed-spin large_mk inputs.

## Large MK Findings

The M6 large_mk completed but is mixed-spin relative to M2 Hessian:

```text
TJ2unlNH_large_mk_uhf_direct_fullgrid_8core.log
Normal termination
S-SQUARED about 11.533
Use only as provisional or constraint-behavior reference, not final M2 charge.
```

Multiple M2 109-atom large_mk and manual-trim attempts failed or became unstable. Naked deletion of ACE caps did not solve the M2 SCF problem.

Small-model M2 MK/RESP succeeded cleanly:

```text
SCF in 4 iterations
E=-3001.8042149377
S-SQUARED=0.778
MK points=2776
RESP total charge about +1
```

The current best candidate route is:

```text
small-M2 RESP + ChgModD-like fixed backbone/CB normalization
```

It generated sane mol2 charges and passed MCPB.py `-s 4` and tleap using the internal-index tleap input.
