---
name: mcpb-md
description: Run, monitor, debug, and audit metal-center MCPB.py workflows that use GAMESS/Gaussian-style QM, AmberTools/tleap, ParmEd, or GROMACS for metalloprotein MD. Use for Fe/heme/nonheme active-site parameterization, MCPB.py small/large/standard model questions, OPT/Hessian/large_mk SCF failures, RESP/ChgMod charge fitting, tleap conversion, GROMACS EM/NVT/NPT/production setup, and cleanup of MCPB/GAMESS project files.
---

# MCPB MD

Use this skill for metal-center parameterization and MD setup workflows, especially when a project involves MCPB.py, GAMESS logs, AmberTools, Fe coordination, open-shell spin states, RESP charges, or GROMACS continuation.

## Required Reading

Read only what is relevant:

- For the end-to-end workflow and decision points, read `references/workflow.md`.
- For this user's local paths and command environment, read `references/local-paths.md`.
- For the TJ Fe/M2 project-specific lessons, read `references/tj-m2-lessons.md`.
- For MD postprocessing, R/S descriptors, Amber-vs-GROMACS comparison, and publication-style visualization, read `references/md-postprocessing-visualization.md`.

## Operating Rules

- Preserve successful QM evidence chains. Do not delete OPT/Hessian/MK inputs or logs that prove convergence unless the user explicitly approves deletion.
- Prefer archive/quarantine folders over deletion for failed GAMESS attempts; failed paths often encode useful SCF strategy lessons.
- Treat `MCPB.py -s 3` charge fitting as separate from GAMESS SCF. ChgModB/C/D constrains RESP after ESP generation; it does not make `large_mk` SCF converge.
- Verify `smmodel_spin` and `lgmodel_spin` in `mcpb.in` before generating inputs. Do not rely on memory.
- For final/publishable workflows, avoid mixed-spin parameters unless clearly labeled provisional.
- For open-shell Fe systems, keep `S-SQUARED`, charge, multiplicity, donor geometry, OPT convergence, Hessian bad points, and RESP charge sanity checks in the audit trail.
- Before starting MD, confirm tleap produced Amber files, then convert/check GROMACS files and run EM before NVT/NPT/production.

## Bundled Scripts

- `scripts/clean_mcpb_pdb_ter_by_geometry.py`: remove only false internal `TER` records whose neighboring C/N atoms form normal peptide geometry.
- `scripts/check_mcpb_topology_bonds.py`: audit tleap/Amber topology connectivity around MCPB residues.
- `scripts/check_gro_fe_distances.py`: report Fe-donor distances from a GROMACS `.gro`.
- `scripts/check_pdb_fe_distances.py`: report Fe-donor distances from a PDB.
- `scripts/diagnose_gro_overlap.py`: quick overlap/geometry diagnosis for GROMACS structures.
- `scripts/postprocess_reactive_md.py`: compute RMSD, reactive atom-pair distances, closest frames, RMSD-metric clusters, and residue short-distance enrichment.
- `scripts/plot_reactive_md_outputs.py`: plot the CSV outputs from `postprocess_reactive_md.py`.
- `scripts/scan_rmsd_cluster_cutoffs.py`: scan RMSD-metric clustering thresholds to choose a useful population cutoff.
- `scripts/cluster_reactive_geometry_features.py`: cluster frames by mechanistic active-site features instead of whole-protein RMSD.
- `scripts/reactive_fel_cv_analysis.py`: compute reaction-coordinate CVs and 2D free-energy landscapes, including PHE ring centroid/plane diagnostics.
- `scripts/plot_reactive_cv_density_heatmaps.py`: plot raw occupancy/log-count heatmaps for reactive CVs without free-energy conversion.
- `scripts/make_reactive_pymol_bundle.py`: merge closest reactive-distance frames into PyMOL states, add Fe-donor bonds/styling, and optionally write a small XTC subset.
- `scripts/make_md_publication_style_summary.py`: build publication-style GROMACS MD plots for RMSD, C4-O1 distance, pro-R/pro-S scalar triple product, active-site groups, and pocket-state enrichment.
- `scripts/make_amber_gromacs_comparison_summary.py`: compare GROMACS and Amber/pmemd trajectories with matched RMSD/RMSF/C4-O1/R-S plots and selected pocket contact features.
- `scripts/make_integrated_proRS_visuals.py`: integrate R/S composition, Group 1/2/3 enrichment, pocket contact heatmaps, and feature-importance summaries.
- `scripts/make_engine_distance_rs_numeric_audit.py`: write exact all-frame C4-O1 distance/R-like/S-like counts, threshold counts, and Amber-minus-GROMACS bin differences.
- `scripts/refresh_pro_rs_scatter.py`: refresh only the product-calibrated R/S triple-product scatter plot when the main summary already exists.
- `scripts/ul1_pocket_residue_scan.py`: stream full MD trajectories and build UL1-centered residue features for every frame, including residue presence, side-chain centroid distances, contact flags, closest atom identities, aromatic ring-plane angles, PCA/UMAP coordinates, and HDBSCAN pocket states. R/S labels are merged only after feature extraction.
- `scripts/summarize_ul1_pocket_scan.py`: combine GROMACS and Amber/pmemd UL1-pocket scans into engine comparison tables and plots for all frames, near-attack frames, R/S feature deltas, and posterior cluster enrichment.
- `scripts/ul1_pocket_deep_state_search.py`: build deeper coordinate-based UL1-pocket feature loops, including local 3D residue coordinates, residue-to-UL1 atom distance matrices, residue-pair networks, and local 3D pocket shape grids; evaluate PCA/UMAP/HDBSCAN state definitions against posterior R/S enrichment.
- `scripts/summarize_ul1_deep_state_search.py`: combine GROMACS and Amber/pmemd deep state-search outputs, rank feature loops, list R/S-enriched states, and audit time-block stability.
- `scripts/small_m2_resp_map.py`: audit small-model M2 RESP charges against MCPB standard-model atoms.
- `scripts/make_smallM2_chgmod3like_mol2.py`: build ChgModD-like mol2 files using small-M2 RESP charges while fixing backbone/CB atoms.

When scripts need project-specific paths, inspect the current project directory and adjust arguments instead of hard-coding old paths.
