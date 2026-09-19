# Changelog

## PhyloCactus 0.4.5

- Corrected the default rooting pattern of
  [`resolve_rooting_outgroup()`](https://beeamerino.github.io/PhyloCactus/reference/resolve_rooting_outgroup.md)
  to the Talinaceae terminals that form a clade of the reference
  topology; the previous default spanned Portulacaceae and Talinaceae,
  which is not a clade of that tree and could not root it.
- Resolved marker alias collisions in
  [`integrate_and_clean_markers()`](https://beeamerino.github.io/PhyloCactus/reference/integrate_and_clean_markers.md)
  with a declared rule and an exported record, and made
  [`run_concatenation_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_concatenation_pipeline.md)
  abort on duplicated terminal names instead of silently keeping the
  first row.
- [`integrate_publication_tree()`](https://beeamerino.github.io/PhyloCactus/reference/integrate_publication_tree.md)
  now reads chronograms written in NEXUS, which is the format
  `TreeAnnotator` produces.
- [`run_joint_realignment()`](https://beeamerino.github.io/PhyloCactus/reference/run_joint_realignment.md)
  accepts `mafft_exec`, matching
  [`run_alignment_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_alignment_pipeline.md).
- The saturation screen reports a degenerate regression as such instead
  of flagging it as saturation.
- [`report_marker_group_coverage()`](https://beeamerino.github.io/PhyloCactus/reference/report_marker_group_coverage.md)
  now runs after `exclude_markers`, so its table describes the
  partitions the supermatrix actually contains.
- Corrected published counts in the README, the vignettes and the
  function documentation: manual exclusions, `phyC` outgroup coverage,
  the scope of the topological constraint scaffold, the scope of the
  `ASTRAL-III` step, and the name of the distributed checklist file.
- Declared the root of the chronogram as uncalibrated and documented the
  deactivated `ACPT_root` bound.
- Regenerated the reference run and republished every distributed
  artifact with 1,024 terminal tips and 12,809 sites across 11
  partitions.
- The reference chronogram is dated at `smooth = 1e-12`, the lowest
  cross-validation error on a grid from 1e3 to 1e-14. The error curve
  levels off below about 1e-06, so this value is one of several that fit
  the data about equally well.
  [`run_treePL_cv()`](https://beeamerino.github.io/PhyloCactus/reference/run_treePL_cv.md)
  now classifies a cross-validation curve as an interior minimum, an
  edge of the grid or a plateau, and reports the result; the smoothing
  sensitivity table of Tutorial 2 spans the plateau.
- [`run_treePL_cv()`](https://beeamerino.github.io/PhyloCactus/reference/run_treePL_cv.md)
  and
  [`automate_treePL()`](https://beeamerino.github.io/PhyloCactus/reference/automate_treePL.md)
  gain `cv_nthreads`, the number of threads used for cross-validation
  (default `1L`), and their default `cvstop` is lowered from 1e-8 to
  1e-14.
- The tree-space figure of
  [`validate_phylogenies()`](https://beeamerino.github.io/PhyloCactus/reference/validate_phylogenies.md)
  identifies each tree by color and shape in a single legend labelled by
  source, with *et al.* in italics, and no longer draws point labels
  inside the panels. Labels are set with the new `tree_labels` argument.
  `ggrepel` is no longer a dependency.
- Fixed the edge test of the cross-validation curve. It compared the
  selected smoothing value with the ends of the grid through
  [`all.equal()`](https://rdrr.io/r/base/all.equal.html), whose absolute
  tolerance treats any two values below about 1e-8 as equal, so a
  minimum at 1e-12 on a grid ending at 1e-14 was reported as an edge.
  The position of the minimum is now read from the sorted grid.
- When
  [`automate_treePL()`](https://beeamerino.github.io/PhyloCactus/reference/automate_treePL.md)
  resumes from a cached maximum-likelihood run, it classifies the
  cross-validation curve on disk and reports its shape; a file it cannot
  classify is reported and does not stop the run.
- `phylocactus_table_tree_distances.csv` as first released in 0.4.5 was
  the 0.4.4 table: Tutorial 4 read the focal tree from the installed
  package, which still held the 0.4.4 chronogram. The table is replaced
  with distances computed on the 0.4.5 chronogram, and Tutorial 4 now
  takes the focal tree from the run (`8_Dating/dated_summary_hpd.tree`),
  falling back to the distributed chronogram when that file is absent.
- The full-size tree and chronogram in `extended_trees/` are replaced
  with those of the 0.4.5 run (1,024 tips; the previous files, from
  0.4.4, had 1,023).
- `phylocactus_table_species.csv` is regenerated from a new IUCN Red
  List query (version 2026-1). One record changes: the assessment
  details of *Stenocactus coptonogonus*, previously empty, are now
  filled; categories and counts are unchanged.
- Documentation corrections: the IUCN token instructions pointed to the
  retired API v3 and now use
  [`rredlist::rl_use_iucn()`](https://docs.ropensci.org/rredlist/reference/rl_use_iucn.html);
  the description of `ASTRAL-III` and of the effect of rooting in
  [`validate_phylogenies()`](https://beeamerino.github.io/PhyloCactus/reference/validate_phylogenies.md)
  (distances and tree space are computed on unrooted trees); the
  reference for the Different Phylogenetic Information distance (Smith
  2020, *Bioinformatics*) and the DOI of Wright (2024); the page size of
  the full-size tree figures; the links to `extended_trees/`, which
  pointed to a branch that does not exist. The tutorials were also
  edited for concision.
- `DESCRIPTION` declares `biocViews:`, so that `remotes` resolves the
  Bioconductor dependencies.
- Added
  [`send_run_notification()`](https://beeamerino.github.io/PhyloCactus/reference/send_run_notification.md)
  and a `notify` argument to
  [`automate_treePL()`](https://beeamerino.github.io/PhyloCactus/reference/automate_treePL.md),
  which emails the status, timings and output paths when a local dating
  run ends. Credentials live outside the package and a failed
  notification never fails the run.
- [`run_alignment_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_alignment_pipeline.md)
  puts every sequence of a marker on the same strand before alignment
  (`fix_strand`), because some GenBank accessions are deposited as
  reverse complements.
- [`run_marker_screening()`](https://beeamerino.github.io/PhyloCactus/reference/run_marker_screening.md)
  flags markers whose ingroup and outgroup sequences share almost no
  k-mers (`no_detectable_homology`): two sequences can carry the same
  gene name and come from different regions.
- Clusters whose GenBank description matches no rule in `genes_map.csv`
  are reported with a warning instead of being dropped without notice.
- The bootstrap replicates are written without a support metric; FBP or
  TBE is chosen when the replicates are summarized with
  [`map_branch_supports()`](https://beeamerino.github.io/PhyloCactus/reference/map_branch_supports.md).
- `treePL` runs report when the gradient optimization did not improve on
  the simulated annealing, in which case the chronogram is not an
  optimized solution.
- The documentation of
  [`rescale_tree()`](https://beeamerino.github.io/PhyloCactus/reference/rescale_tree.md)
  states what the scaling factor changes: the scale of the rate
  parameters and the numerical conditioning of the optimization, not
  which branches `treePL` shortens.
- [`evaluate_dna_barcoding()`](https://beeamerino.github.io/PhyloCactus/reference/evaluate_dna_barcoding.md)
  is a placeholder that signals an error. The DNA barcoding section is
  under construction; its scope and open questions are set out in
  Tutorial 5.
- The `checklist_path` argument of
  [`infer_gene_trees()`](https://beeamerino.github.io/PhyloCactus/reference/infer_gene_trees.md)
  is retained and remains reserved for the molecular diagnostic section,
  which is still under construction.

## PhyloCactus 0.4.4

- Updated outgroup sampling in data preparation and tree rooting to
  include Talinaceae (*Talinum*, *Talinella*) alongside Portulacaceae
  and Anacampserotaceae.
- Standardized primary clade support metric to Felsenstein Bootstrap
  Proportions (FBP) across documentation, functions, and collapsing
  workflows
  ([`integrate_publication_tree()`](https://beeamerino.github.io/PhyloCactus/reference/integrate_publication_tree.md)),
  retaining Transfer Bootstrap Expectation (TBE) as a complementary
  diagnostic metric.
- Added support for identifying and distinguishing topologically
  constrained nodes in tree visualization
  ([`classify_constrained_nodes()`](https://beeamerino.github.io/PhyloCactus/reference/classify_constrained_nodes.md)).
- Replaced shell wrapper execution of treePL with direct system calls
  ([`run_treePL_direct()`](https://beeamerino.github.io/PhyloCactus/reference/run_treePL_direct.md))
  and native R cross-validation optimization
  ([`run_treePL_cv()`](https://beeamerino.github.io/PhyloCactus/reference/run_treePL_cv.md)).
- Retired legacy `run_treePL()` wrapper function.
- Improved documentation, metadata manifests, and error handling across
  all analytical modules.
