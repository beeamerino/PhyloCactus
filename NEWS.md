# PhyloCactus 0.4.5

* Corrected the default rooting pattern of `resolve_rooting_outgroup()` to the Talinaceae terminals that form a clade of the reference topology; the previous default spanned Portulacaceae and Talinaceae, which is not a clade of that tree and could not root it.
* Resolved marker alias collisions in `integrate_and_clean_markers()` with a declared rule and an exported record, and made `run_concatenation_pipeline()` abort on duplicated terminal names instead of silently keeping the first row.
* `integrate_publication_tree()` now reads chronograms written in NEXUS, which is the format `TreeAnnotator` produces.
* `run_joint_realignment()` accepts `mafft_exec`, matching `run_alignment_pipeline()`.
* The saturation screen reports a degenerate regression as such instead of flagging it as saturation.
* `report_marker_group_coverage()` now runs after `exclude_markers`, so its table describes the partitions the supermatrix actually contains.
* Corrected published counts in the README, the vignettes and the function documentation: manual exclusions, `phyC` outgroup coverage, the scope of the topological constraint scaffold, the scope of the `ASTRAL-III` step, and the name of the distributed checklist file.
* Declared the root of the chronogram as uncalibrated and documented the deactivated `ACPT_root` bound.
* Regenerated the reference run and republished every distributed artifact with 1,024 terminal tips and 12,809 sites across 11 partitions.
* Divergence times estimated under penalized likelihood in `treePL` using an optimal rate smoothing of `smooth = 1e-12` identified via single-threaded deterministic cross-validation across 18 orders of magnitude (1e3 to 1e-14).
* Added `send_run_notification()` and a `notify` argument to `automate_treePL()`, which emails the status, timings and output paths when a local dating run ends. Credentials live outside the package and a failed notification never fails the run.
* The `checklist_path` argument of `infer_gene_trees()` is retained and remains reserved for the molecular diagnostic section, which is still under construction.

# PhyloCactus 0.4.4

* Updated outgroup sampling in data preparation and tree rooting to include Talinaceae (*Talinum*, *Talinella*) alongside Portulacaceae and Anacampserotaceae.
* Standardized primary clade support metric to Felsenstein Bootstrap Proportions (FBP) across documentation, functions, and collapsing workflows (`integrate_publication_tree()`), retaining Transfer Bootstrap Expectation (TBE) as a complementary diagnostic metric.
* Added support for identifying and distinguishing topologically constrained nodes in tree visualization (`classify_constrained_nodes()`).
* Replaced shell wrapper execution of treePL with direct system calls (`run_treePL_direct()`) and native R cross-validation optimization (`run_treePL_cv()`).
* Retired legacy `run_treePL()` wrapper function.
* Improved documentation, metadata manifests, and error handling across all analytical modules.
