# gromacs_full deep UL1-pocket state search

This run evaluates four feature loops:

1. `local_coords`: residue centroids and closest atoms in a UL1-fixed local 3D frame.
2. `ul1_distance_matrix`: residue-to-every-UL1-heavy-atom distance fingerprints.
3. `pair_network`: pairwise residue-centroid distances and pocket atom-cloud moments.
4. `shape_grid`: coarse 3D atom-density grid around UL1.

R/S labels are merged only after feature extraction and are used for posterior enrichment scoring.
The best-state score is not a mechanism proof; it ranks which unsupervised state definition most enriches R-like vs S-like frames.
Blocked logistic AUC on PCA20 is included only as a diagnostic for whether the feature space contains R/S information under time-block validation.

Main files:
- `gromacs_full_all_deep_feature_hdbscan_scores.csv`
- `gromacs_full_all_deep_feature_best_state_details.csv`
- `gromacs_full_deep_feature_loop_score_comparison.png`
