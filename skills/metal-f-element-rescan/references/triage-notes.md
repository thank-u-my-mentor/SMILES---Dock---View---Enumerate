# Metal-F Rescan Triage Notes

## Why CCD Element Enumeration Over-Recalls

CCD formulas answer "does this chemical component contain an element", not "is this a biological catalytic metal center." A formula-level metal component set will include:

- free metal ions (`FE`, `CU`, `NI`, `MN`)
- small salts/counterions
- hemes and Fe-S clusters
- synthetic coordination complexes
- metal-containing inhibitors, probes, crystallization additives, and model compounds

Therefore, do not interpret the raw count of metal-containing components as enzyme-center diversity.

## Evidence Hierarchy

1. CCD formula element membership: useful for broad recall only.
2. RCSB comp_id entry search: candidate entry recall only.
3. mmCIF atom coordinates: required for F-to-metal distance.
4. RCSB polymer entity metadata: required for complex warnings, expression host, organism, UniProt IDs.
5. Manual review: required for whether a hit is a practical recombinant enzyme target.

## Practical Recombinant-Enzyme Filter

Prefer:

- one protein entity
- E. coli or E. coli BL21/BL21(DE3) expression host
- clear UniProt ID
- protein length compatible with downstream HMMER/tree work
- fluorinated ligand and target metal within 3.5 Å

Deprioritize:

- `protein_entity_count > 1`
- native purified multi-subunit complexes
- unknown expression host
- metal component classified as synthetic/counterion unless it directly supports the desired chemistry.
