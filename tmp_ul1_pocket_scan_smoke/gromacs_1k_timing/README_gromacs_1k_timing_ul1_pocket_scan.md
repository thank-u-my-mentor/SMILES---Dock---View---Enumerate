# gromacs_1k_timing UL1-centered pocket residue scan

- Frames analyzed: 1000
- Frame stride: 1 (`1` means every saved frame; no postprocessing stride).
- UL1 around cutoff: 6.0 A.
- R/S labels were merged after feature extraction. They were not used to define PCA/UMAP/HDBSCAN states.
- `min_to_UL1_heavy_A` and `min_to_C4O1_A` are per-frame closest-pair distances, not trajectory-wide minima.
- Side-chain centroid distances, ring-plane angles, residue class flags, closest atom identities, and contact occupancies are included to avoid relying on a single hand-picked feature.
- Pocket clusters are unsupervised density states from the residue-feature matrix. Use `*_cluster_RS_summary.csv` to check posterior R/S enrichment.

## Key outputs
- `gromacs_1k_timing_ul1_pocket_frame_residue_features.csv.gz`: long table; one row per contacted residue per frame.
- `gromacs_1k_timing_ul1_pocket_residue_summary.csv`: residue occupancy and median/mean distance summary.
- `gromacs_1k_timing_ul1_pocket_residue_summary_by_RS.csv`: same summary split by R-like/S-like.
- `gromacs_1k_timing_ul1_pocket_frame_matrix_embeddings.csv.gz`: fixed-width frame matrix with PCA/UMAP coordinates and pocket cluster labels.
