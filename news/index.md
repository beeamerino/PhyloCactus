# Changelog

## PhyloCactus 0.5.2

- Step 7 in parallel (`classify_barcoding_folds(workers = )`) stops with
  an error that names the worker when a worker process ends without
  writing its status (memory exhausted, a killed process, a crash of R),
  instead of waiting without end. A worker that has written no process
  id after `getOption("PhyloCactus.worker_grace_seconds")` seconds (by
  default three polls, never under 60) counts as ended. The analyses and
  the distributed results are those of 0.5.0.

## PhyloCactus 0.5.1

Documentation, metadata and the website. The analyses and the
distributed results are those of 0.5.0.

- `DESCRIPTION` describes both branches of the package and lists the
  external programs it calls (`SystemRequirements`);
  `citation("PhyloCactus")` takes the new title.
- The help pages of
  [`assemble_ingroup_phylotar()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_ingroup_phylotar.md)
  and
  [`assemble_outgroup_phylotar()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_outgroup_phylotar.md),
  and Tutorial 1, document the cut of records over 2000 bases introduced
  in 0.5.0, with the tables that list them.
- Tutorial 2 reads the weak support of the backbone: soft polytomies
  from missing data, not evidence of a simultaneous divergence.
- Tutorial 6 explains the two version numbers of its reports.
- The homology check of
  [`integrate_and_clean_markers()`](https://beeamerino.github.io/PhyloCactus/reference/integrate_and_clean_markers.md)
  states, in its help page and its warning, the interval of the
  reference run of 0.5.0 (0.26 to 0.58 at k = 20).
- The website no longer serves the pages of the functions removed in
  0.5.0: each deployment replaces the previous one.
- The checklist of Korotkova et al. (2021) is cited wherever it is
  named; the comments of the source code describe the code.
- The 0.5.0 entry below is completed with changes it left out.

## PhyloCactus 0.5.0

### Molecular diagnostic section (DNA barcoding)

- A new branch builds a reference library of Cactaceae from the Sanger
  marker records of the phylotaR mining of the phylogeny, and optionally
  the user’s own sequences
  ([`assemble_barcoding_dataset()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_barcoding_dataset.md),
  [`curate_barcoding_markers()`](https://beeamerino.github.io/PhyloCactus/reference/curate_barcoding_markers.md),
  [`screen_barcoding_markers()`](https://beeamerino.github.io/PhyloCactus/reference/screen_barcoding_markers.md),
  [`finalize_barcoding_library()`](https://beeamerino.github.io/PhyloCactus/reference/finalize_barcoding_library.md)).
  Names are reconciled against the checklist of Caryophyllales.org
  (Korotkova et al. 2021); records over 2000 bases are cut to the extent
  of the shorter records of their locus; identical sequences of a
  species are collapsed. The library of this version holds 15 loci.
  [`extract_barcoding_plastome_loci()`](https://beeamerino.github.io/PhyloCactus/reference/extract_barcoding_plastome_loci.md)
  crops the loci of the library from GenBank plastomes as extra records;
  [`barcoding_outgroup_taxids()`](https://beeamerino.github.io/PhyloCactus/reference/barcoding_outgroup_taxids.md)
  gives the six NCBI Taxonomy IDs of the outgroup the controls query the
  library with.
- The library is validated before it is used: folds that ask for the
  species (scheme E) or, with the species left out, for the genus
  (scheme G)
  ([`build_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/build_barcoding_folds.md));
  the barcode gap per locus
  ([`analyze_barcode_gap()`](https://beeamerino.github.io/PhyloCactus/reference/analyze_barcode_gap.md));
  classification of the folds by nearest neighbour and by IdTaxa, on one
  machine in parallel processes (`classify_barcoding_folds(workers = )`)
  or on a computing cluster
  ([`classify_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/classify_barcoding_folds.md),
  [`generate_barcoding_job_scripts()`](https://beeamerino.github.io/PhyloCactus/reference/generate_barcoding_job_scripts.md),
  [`merge_barcoding_chunks()`](https://beeamerino.github.io/PhyloCactus/reference/merge_barcoding_chunks.md));
  negative controls for leakage, outgroup queries and resubstitution
  ([`run_barcoding_controls()`](https://beeamerino.github.io/PhyloCactus/reference/run_barcoding_controls.md));
  thresholds
  ([`sweep_barcoding_threshold()`](https://beeamerino.github.io/PhyloCactus/reference/sweep_barcoding_threshold.md));
  metrics per class
  ([`summarise_barcoding_metrics()`](https://beeamerino.github.io/PhyloCactus/reference/summarise_barcoding_metrics.md),
  [`summarise_barcoding_genus_discrimination()`](https://beeamerino.github.io/PhyloCactus/reference/summarise_barcoding_genus_discrimination.md)).
- A sample is identified locus by locus, in three states (species,
  genus, no answer), never combined into one verdict, from sequences or
  a GenBank flat file
  ([`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md)),
  a sample sheet
  ([`identify_barcoding_samples()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_samples.md)),
  an assembly
  ([`identify_barcoding_assembly()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_assembly.md),
  with `minimap2` and `samtools`) or paired Illumina reads
  ([`identify_barcoding_reads()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_reads.md),
  with `fastp` and GetOrganelle). A species is named only above a
  species threshold per locus, set from the folds of scheme G; `rbcL`
  answers at genus rank only (argument `genus_only_loci`), because its
  species answers were right in 11 of 34 on the library of this version.
  Each run writes a self-contained HTML report
  ([`report_barcoding_identification()`](https://beeamerino.github.io/PhyloCactus/reference/report_barcoding_identification.md)).
- Tutorial 5 is now the tutorial of the branch: the script runs steps 1
  to 11 and draws the figures of the library and its validation, with
  its limitations. The plastome MN517611.1 is distributed in `extdata`
  as the example query.
- Tutorial 6 identifies public samples of each kind of input (a Sanger
  record, specimens of several loci, GenBank plastomes, a sample sheet,
  a genome assembly and Illumina reads) with the library of this
  version, and gives the rates measured with it; Tutorial 5 keeps the
  library and its validation. Two new functions draw an identification
  from the tables of its report:
  [`plot_barcoding_identification()`](https://beeamerino.github.io/PhyloCactus/reference/plot_barcoding_identification.md)
  (the answer per locus and the sampling of the named genus) and
  [`plot_barcoding_identification_summary()`](https://beeamerino.github.io/PhyloCactus/reference/plot_barcoding_identification_summary.md)
  (several identifications side by side). The function reference is now
  Tutorial 7; its former address is redirected.

### Changes that also reach the phylogeny

- The reference run is regenerated and every distributed artifact
  republished: a supermatrix of 1218 terminals, 13 700 sites and 13
  partitions, analysed as 1208 terminals once identical sequences are
  collapsed; the chronogram is dated at a smoothing of 1e-05, an
  interior minimum of the cross-validation, and its root is 45.53 Ma.
  `MANIFEST.csv` gains `expected_tips`, the number of terminals each
  distributed tree must have.
- [`assemble_ingroup_phylotar()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_ingroup_phylotar.md)
  and
  [`assemble_outgroup_phylotar()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_outgroup_phylotar.md)
  cut every record over 2000 bases to the extent the shorter records of
  its marker cover, as the identification branch does, and list the long
  records in `TABLE_LONG_RECORDS_CUT_INGROUP.csv` and
  `TABLE_LONG_RECORDS_CUT_OUTGROUP.csv`.
- The family filter of
  [`validate_phylogenies()`](https://beeamerino.github.io/PhyloCactus/reference/validate_phylogenies.md)
  keeps the species of Cactaceae in the checklist (Korotkova et
  al. 2021) instead of the species of the constraint table; the tree
  distances of the reference run do not change.
- 31 species of the new mining are added to `cactus_constraints.csv`, so
  that the constraint places every terminal of the supermatrix.
- One style for the figures of every tutorial:
  [`theme_phylocactus()`](https://beeamerino.github.io/PhyloCactus/reference/theme_phylocactus.md),
  [`phylocactus_colours()`](https://beeamerino.github.io/PhyloCactus/reference/phylocactus_colours.md)
  and
  [`phylocactus_panel_tags()`](https://beeamerino.github.io/PhyloCactus/reference/phylocactus_panel_tags.md).
- The phylotaR mining of both branches keeps sequences of 100 to 5000
  bases, and the outgroup of the phylogeny is mined through the same
  function as the ingroup.
- The manual exclusions are keyed by locus and accession instead of
  cluster number, which a new mining renumbers; the logs report the
  exclusions actually applied.
- 21 `ITS` records annotated in GenBank as pseudogenes are added to
  `manual_exclusions_ingroup.csv`, read by both branches.
- `genes_map.csv` gains variant names of loci found in the new mining.
- The subject of every email notification of the package carries a
  cactus emoji.

### Removals

- Removed `infer_gene_trees()`. No stage of the pipeline, no tutorial
  script and no test called it, and on its documented input (one
  sequence per species) its species monophyly table was always empty.
  Its `checklist_path` argument had no effect. Per-locus monophyly is
  planned for version 0.6.
  [`extract_species_binomial()`](https://beeamerino.github.io/PhyloCactus/reference/extract_species_binomial.md),
  which it shared with the marker integration step, is kept.
- Removed the placeholder `evaluate_dna_barcoding()` and the `barcoding`
  documentation topic. The placeholder only signalled an error, and its
  arguments belonged to the implementation withdrawn in 0.4.2. The
  molecular diagnostic section was rebuilt as a separate branch with its
  own functions (above).

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
- `evaluate_dna_barcoding()` is a placeholder that signals an error. The
  DNA barcoding section is under construction; its scope and open
  questions are set out in Tutorial 5.
- The `checklist_path` argument of `infer_gene_trees()` is retained and
  remains reserved for the molecular diagnostic section, which is still
  under construction.

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
