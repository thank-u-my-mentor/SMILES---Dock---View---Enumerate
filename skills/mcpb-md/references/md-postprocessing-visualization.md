# MD Postprocessing And Visualization

This reference captures the TJ Fe/M2 postprocessing logic used after a successful MCPB.py -> Amber/GROMACS setup. It is meant to make the plotting workflow reproducible when Codex is not available.

## Core Reaction Descriptors

- `C4-O1 distance`: the main near-attack distance. It is descriptive for classical MD; no covalent bond is formed.
- `signed triple product around C4`: `V = det[(O1-C4), (C5-C4), (C3-C4)]`. Its sign is calibrated against a product reference before assigning `R-like` or `S-like`.
- `R-like/S-like`: a product-reference-calibrated frame label, not proof that an unreacted MD frame has formed a product.
- `Group 1/2/3`: unsupervised active-site geometry clusters. These are not R/S labels; they must be checked against R/S by composition or same-axis scatter plots.
- `Stride 10 frame`: a postprocessing subsample of every 10th saved trajectory frame. In the TJ 200 ns production trajectory saved every 5 ps/frame, stride 10 means one analyzed frame per 50 ps. It is not the MD integration timestep and not a separate simulation.

Use Angstrom as `\u00c5` in plain labels and `$\AA^3$` for scalar-triple-product volume axes.

## Recommended Plot Hierarchy

1. Stability plots: C-alpha RMSD and C-alpha RMSF.
2. Reactive distance plots: C4-O1 histograms for each engine and engine comparison.
3. R/S descriptor plots: triple-product-vs-C4-O1 scatter for GROMACS and Amber on matched axes.
4. Cluster sanity plots: the same GROMACS strided frames colored once by R/S and once by Group 1/2/3.
5. Pocket interpretation plots: R/S pocket-distance heatmaps and selected feature-difference plots.

For presentations, prefer readable square plots with large axis labels, tick labels, and visible legend markers. Do not rely on tiny scatter markers in legends.

If Group 1/2/3 are mostly separated by C4-O1 distance or pocket basin, report them as auxiliary pocket states rather than pro-R/pro-S predictors. The product-calibrated R/S label should remain the primary stereochemical analysis axis.

## Pocket Contact Features

Current selected contact features are deliberately interpretable:

- F336/PHE336 minimum distance to C4/O1.
- F347/PHE347 minimum distance to C4/O1.
- H270/HD2 minimum distance to C4/O1.
- Q269/GLN255 minimum distance to C4/O1.
- F336 ring-centroid to UL1 C4 distance.

`min distance` means a per-frame contact feature: within one frame, take the closest relevant atom pair between the residue and C4/O1. It does not mean the minimum over the entire trajectory. Across frames, report median and mean separately:

- Median is robust to occasional excursions and is good for heatmaps.
- Mean is a conventional average and is more sensitive to tail behavior.

## Random-Forest Feature Importance

The current Random Forest is a lightweight, interpretable, GNN-like contact-feature classifier. Residues/features are treated as contact nodes and distances are treated as edge-like descriptors.

Feature importance means mean decrease in Gini impurity across trees. It is not an interaction energy and does not prove causality. `Combined` means the GROMACS and Amber near-attack feature rows were pooled to test whether a signal persists across engines; it is not a third physical trajectory.

Use the random forest as an exploratory audit, not as a final mechanism model. Each frame is one row, the R/S label comes from the already-defined scalar triple product, and the current inputs are a small hand-selected feature set. MD frames are temporally correlated, so cross-validation over individual frames can overstate generalization. A convincing predictive model needs replicate/mutant-level validation and a mutation-compatible UL1-centered pocket feature set.

For mutant comparisons, do not depend only on F336 ring-centroid features, because F336R removes the aromatic ring. Use mutation-compatible UL1-centered features:

- residue or side-chain centroid to C4, O1, and substrate phenyl centroid;
- residue min distance to C4/O1;
- contact occupancy under thresholds such as 3.5, 4.0, 4.5, and 5.0 \u00c5;
- residue class descriptors: aromatic, cationic, anionic, polar, hydrophobic, H-bond donor/acceptor;
- optional atom-level pocket graph/contact matrix for all atoms or residues within a UL1-centered cutoff.

## Future Pocket-Graph Direction

The next robust cross-mutant model should extract all residues or atoms within a UL1-centered cutoff, analogous to a PyMOL `around 4` selection but computed directly from trajectory coordinates/topology. Represent the pocket as a frame-wise graph:

- nodes: pocket residues or atoms with residue class and atom element/type metadata;
- edges: distances/contact occupancies to UL1 atoms and between nearby pocket nodes;
- global labels: C4-O1 distance, R-like/S-like triple-product label, active-site group, engine/replicate/mutant.

Do not present this as a trained GNN until the full feature extraction and validation are actually implemented.

## Full-Frame UL1-Centered Pocket Scan

Use `scripts/ul1_pocket_residue_scan.py` when the goal is to avoid hand-picking only a few familiar residues such as F336/F347. The script streams every saved frame and identifies all non-solvent residues within a UL1-centered cutoff, then writes:

- a long feature table: one row per contacted residue per frame;
- a residue summary table: `seen_fraction` is the real 6 Å pocket occupancy;
- an R/S-split summary table: R/S labels are posterior scalar-triple-product labels, not clustering inputs;
- a fixed-width frame matrix with PCA, UMAP, and HDBSCAN pocket-state labels.

Recommended features include side-chain centroid distances to UL1 C4, O1, and the UL1 phenyl centroid; per-frame closest heavy-atom distances; contact flags; closest atom identities; residue class flags; and aromatic ring-plane angles. Do not treat one per-frame minimum distance as a golden standard. Use it as one descriptor among several.

Use `scripts/summarize_ul1_pocket_scan.py` to combine GROMACS and Amber/pmemd scans. It reports all-frame and near-attack (`C4-O1 <= 3 Å`) residue presence, R-like minus S-like median feature shifts, and posterior R/S composition of unsupervised pocket states.

For model building, first inspect PCA/UMAP/HDBSCAN pocket states with R/S as a posterior label. Only then consider Random Forest, XGBoost, or a graph model. Any supervised model should use time-block or replicate/mutant-level validation, because individual MD frames are autocorrelated.

## Deep Pocket-State Search When Shallow Features Fail

If residue min-distance/contact flags do not separate R-like and S-like frames, run `scripts/ul1_pocket_deep_state_search.py` before moving to Random Forest. It evaluates several richer, coordinate-based feature loops:

- `local_coords`: residue side-chain/backbone centroids, closest atoms, ring centroids, and ring normals in a UL1-fixed local 3D coordinate system.
- `ul1_distance_matrix`: residue-to-every-UL1-heavy-atom distance fingerprints. This is more informative than one min distance because it preserves the substrate-facing geometry.
- `pair_network`: pairwise residue-centroid distances plus pocket atom-cloud moments.
- `shape_grid`: coarse 3D pocket atom-density grids around UL1 in the UL1-fixed local frame.

Cluster each feature loop by PCA/UMAP plus HDBSCAN parameter sweeps, then overlay R/S labels only afterward. Compare loops by coverage, purity, NMI/ARI, weighted R/S enrichment, and time-block stability. A good state definition should not merely isolate a tiny number of frames or reproduce only the C4-O1 distance axis.

For the TJ M2 200 ns analysis, the first deep-state run found that GROMACS separated best with `local_coords`, while Amber/pmemd separated best with `ul1_distance_matrix`; `shape_grid` was also strong, and `pair_network` was weaker. This suggests that R/S-like sampling is encoded in local 3D pocket/substrate geometry, not in a single residue contact.
