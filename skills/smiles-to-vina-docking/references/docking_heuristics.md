# Docking Heuristics For Lead Analog Search

## Why a polar substituent can lower Vina affinity

- The added donor or acceptor may not satisfy hydrogen-bond geometry. Distance alone is not enough; angle and orientation matter.
- A polar group buried without a partner pays a desolvation penalty.
- Extra atoms can introduce steric strain even in an apparently empty pocket because Vina optimizes a simplified flexible ligand against a mostly rigid receptor.
- Extra rotatable bonds add entropic and scoring penalties.
- A substituent can pull the ligand into a different pose that scores worse or loses the original anchor interaction.
- PDBQT protonation, tautomer, charge assignment, or aromaticity perception can differ from the intended chemistry.
- A single Vina run can overfit stochastic search noise; require repeatability for small improvements.

## Priority order for small edits

1. Halogen scan on solvent-facing aromatic C-H sites: F before Cl unless a hydrophobic pocket is obvious.
2. Methyl scan for shape filling with minimal polarity change.
3. H-bond acceptor scan: OMe, CN, carbonyl/formyl when geometry points toward donor residues.
4. H-bond donor scan: OH or NH2 only when the group can remain exposed or point at a strong acceptor.
5. Tiny polar handles: hydroxymethyl or aminomethyl only when there is enough room and one added rotatable bond is acceptable.
6. Paired edit: remove a substituent that causes clash or torsion, then add one new group elsewhere.

## Stop conditions

- Stop a local edit family if all members score worse than baseline by more than 0.5 kcal/mol.
- Promote a candidate only if it beats the redocked baseline by at least 0.3 kcal/mol in first pass.
- Confirm promoted candidates with at least three seeds or higher exhaustiveness.
- Treat a hit as suspicious if the score improves but the pose loses the known ligand-like orientation or creates buried unsatisfied donors/acceptors.

## Drug-likeness limits for first-pass filtering

- Heavy atoms: 12 to 45.
- Molecular weight: 180 to 650 Da.
- cLogP: -1 to 6.
- HBD: at most 5.
- HBA: at most 10.
- TPSA: at most 150.
- Rotatable bonds: at most 10, ideally at most 7 for kinase-like pockets.
- Formal charge: prefer -1, 0, or +1; avoid zwitterions unless the binding model supports them.

## Record every negative result

Use a CSV or SQLite table with: parent SMILES, candidate SMILES, edit label, changed atom indices if known, RDKit validity, prep status, Vina status, affinity, pose path, log path, properties, and rejection reason. A useful optimizer learns where not to continue.
