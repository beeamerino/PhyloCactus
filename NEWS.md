# PhyloCactus 0.5.5

This version closes the package part of the final closure review of 0.5.4 (25 November 2026). It corrects the help page of `build_constraint_scaffold()`, which described relationships the constraint imposes as estimated from the data, and fixes the other findings whose destination was the package. No analysis was run again: the supermatrix, the constraint tree, the maximum-likelihood tree, the chronograms and the identification library are those of 0.5.3. The finding of the review that each change closes is given in brackets.

## What the constraint imposes

* The help of `build_constraint_scaffold()` said that the constrained search estimates every relationship except the monophyly of higher ranks; the function writes a fixed topology. One description, in the help page, README (Stage 2) and Tutorial 2 (Module 8), lists what the 34 constrained bipartitions of the reference run fix: the monophyly of 19 groups (14 in Cactaceae, 5 outgroup genera), 6 bipartitions of the backbone of Cactaceae (*Leuenbergeria*, *Pereskia*, Opuntioideae, *Maihuenia*, Cactoideae, and *Blossfeldia* sister to the rest of Cactoideae), 6 within Opuntioideae and Cactoideae, and 3 of the outgroup families (Anacampserotaceae and its internal order, Talinaceae). It names the relationships left to the data (the trichotomy of Anacampserotaceae, Cactaceae and Portulacaceae, the base of core Cactoideae, everything within the groups), says that the topology is written in step 7 of the function in `R/inference.R`, and names the imposed relationships that published analyses leave unresolved or contest (the tribes of Opuntioideae, core Cactoideae II). Tutorial 2 no longer calls the relationships among the major lineages "established"; Tutorial 7 describes the function the same way [AP-01, AM-04, AB-07].
* `build_constraint_scaffold()` loses an unreachable genus-constraint block and requires only `Specie_name` and `Clade`; its help says that the 21 groups and the topology are fixed in the code. The scaffold of the reference run is identical to `phylocactus_constraint_tree.tree` [AP-18].

## Dating

* `run_treePL_cv()`: each priming repeat runs with its own seed, s to s + `n_prime` - 1 from the seed of the configuration, recorded in `prime_<label>.csv`. Up to 0.5.4 the repeats used one seed, gave the same parameters, and `prime_rule` could not change the selection. The distributed chronogram was dated with the earlier behaviour [AP-02].
* `automate_treePL()` reuses results already in `results_dir` only when they were produced with the same seed, `rescale_factor`, scaled `numsites`, calibration file, maximum-likelihood tree, priming and cross-validation settings and, for the replicates, replicate file, number of replicates, smoothing and optimization lines; each stage records them in `run_parameters_ML_tree.csv` or `run_parameters_BS_tree.csv`. A call that differs stops; a folder with results and no such record (written before 0.5.5) stops too, unless `reuse_unverified = TRUE`, which reuses it with a warning. Only the folders `BS_1` to `BS_<num_bs>` of the call are dated and collected [AP-03].
* `automate_treePL()` restores the random state of the session on exit, reads indented calibration lines, and says when it dates an already rooted tree without applying `outgroup` [AP-20].
* `generate_temporal_bootstrap_script()` takes `blopt`, with the default of `calculate_temporal_bootstraps()`, `nr_safe`; the temporal replicates of the reference run were estimated with `NR-FAST`, the default of `RAxML-NG` 1.2.2-master, as their log records, and Tutorial 2 passes `blopt = "nr_fast"` to reproduce them [AP-04]. The script reads its threads from `SLURM_CPUS_PER_TASK` and is made executable [AP-22].
* Row 8 of `calibrations_bounds.csv` and Tutorial 2 name the reference node of Wu *et al*. (2026) as the most recent common ancestor of Cactaceae and Portulacaceae in a sampling without Anacampserotaceae, not the crown of the ACP clade; row 9 (ACPT) is unchanged [AP-12, C-01].
* Tutorial 2 names the node dated as "Cactoideae crown" in its two age tables, the MRCA of the terminals of Cactoideae without *Blossfeldia*, as Cacteae plus core Cactoideae (Cactoideae without *Blossfeldia*), as do the smoothing-sensitivity chunk and `inst/scripts` [AB-02].
* Tutorial 2 names 41.82 Ma without "median", says that core Cactoideae II was checked, not tested, by the search that released it, and no longer calls Rhipsalideae monophyletic in that search, where its monophyly stayed imposed (BMM's answers 3, 4 and 5 to the review) [AB-13, AB-05, AB-04].

## Inference and data assembly

* `preprocess_partitions()` removes the reduced files of an earlier check before `RAxML-NG --check`, so the reduction it reports is that of the matrix checked in the call, and takes `threads` [AP-05].
* `run_modeltest_ng()` takes `output_dir` and `seed` (`--rngseed`); both default to the earlier behaviour [AP-23].
* `run_marker_screening()` takes `gamma_shape`, the shape of the gamma distribution of rates in the corrected K80 distance of the saturation test, default 1, the value `ape` used for `gamma = TRUE` in earlier versions; the help, the log, Tutorial 1 and Tutorial 7 say that the thresholds 0.3 and 0.5 were set with it [AP-06].
* A batch of GenBank metadata whose retries all fail is logged and counted, its accessions are left out of the cache so that a new call retries them, and the call stops when they exceed `max_failed_fraction` (default 0.05) of those requested; an error of `phylotaR::run()` is logged and given as a warning [AP-07]. `fetch_genbank_metadata()` deprecates `force_download`, which did nothing [AP-30].
* `cp_normalize_marker()` turns typographic dashes into hyphens; its pattern replaced a hyphen with a hyphen [AP-25].
* `run_concatenation_pipeline()` takes `registry_file`; `NULL` keeps the lookup in the tutorial layout. Hyphenated epithets are written with an underscore in terminal labels and by `extract_species_binomial()`, as `standardize_taxon()`, the constraint table and the distributed trees write them [AP-26].
* `calculate_ml_tree()`, `map_branch_supports()`, `calculate_temporal_bootstraps()`, `calculate_rf_distances()` and `run_modeltest_ng()` check their input files and name each missing one; `calculate_rf_distances()` checks its binary; `map_branch_supports()` reports its own name when `RAxML-NG` fails [AP-19].
* The help of every SLURM generator says that the default partition and modules are those of the NLHPC cluster of the reference run [AP-21]; the help of `run_local_bootstraps()` says that performance cores are detected only when `threads = NULL` [AP-28], and that of `generate_ml_search_script()` names the `ml_search/` subfolder [AP-29].
* `validate_phylogenies()` takes `family` and roots through `root_on_clade()`; distances are computed on unrooted trees and do not change [AP-24].

## Distributed data

* The full-size Figures 7 and 9 (`extended_trees/` and `inst/extdata/phylocactus_figure_*.pdf`) and the annotated trees of `extended_trees/` carry the 28 groups of 0.5.4 instead of those of 0.5.3; their trees, supports and ages are unchanged [AP-09].
* `PROVENANCE.csv` was republished with `publish_run_outputs.R` on the reference run: its row of `cactus_constraints.csv` described the file of 0.5.3, and the rows of the two figures and of `calibrations_bounds.csv` now match the files. A test compares every md5 and size of `PROVENANCE.csv` with the distributed file [AP-08].

## Documentation

* README, the overview and Tutorial 1 attributed errors of orthology assessment and alignment to incomplete lineage sorting and reticulation. They now say that ILS and reticulation make the histories of loci discordant, which concatenation does not model, and that orthology errors come from paralogous copies and from homeologous copies left by hybridization or polyploidy [AP-32].
* The overview and Tutorials 2 to 6 open with the problem they address [AP-31].
* README and the workflow figure of the package give five published trees, as Tutorial 4 compares since 0.5.4 [AP-10, AM-03]; the line under the README title names molecular identification [AP-38]; the lifecycle badge is stable (BMM's answer 7) [AP-43].
* Tutorial 3 restricts the comparison of threatened species inside and outside the tree to assessed species, and counts 265 assessed species without sequence data (38.5% threatened), leaving out the one Not Evaluated as the manuscript does [AB-08, AB-09].
* Tutorial 4 records the licence of the source of each reference tree, and a chunk computes its table of RF with and without the constrained bipartitions from the outputs of `validate_phylogenies()` [AP-11, AP-13].
* Tutorials 5 to 7 follow the section structure of the other tutorials and have headings in sentence case (BMM's answer 8) [AP-45].
* Software names are set in backticks, without bold, in the vignettes [AP-34]; Cactaceae is bold once in Tutorial 2 [AP-35]; Tutorials 3 and 7 cite Korotkova *et al*. (2021) once [AP-33].
* Every exported function has examples; 16 of the 26 that had none run on distributed data [AP-27]. `add_outgroup_bar()` is in the reference index of the site and the internal `extract_species_binomial()` is not [AP-14, AP-30].

## Package metadata and tests

* `LICENSE.md` holds the full MIT text, excluded from the build, so that GitHub detects the licence [AP-16]; `parallel` is in Imports [AP-17]; the empty `biocViews` field is gone [AP-36].
* The two tests of the cactus emoji in the classifier banner and the report run only in a UTF-8 locale [AP-39].

# PhyloCactus 0.5.4

This version adds the dated tree of Carruthers *et al*. (2026) to the comparison with published phylogenies (Module 13, Tutorial 4), corrects four points of the documentation, names groups of the chronogram that the figures of Tutorial 3 left without a label, and adds `add_outgroup_bar()`, which marks the outgroup in two of those figures. No analysis was run again: the supermatrix, the constraint tree, the maximum-likelihood tree, the chronograms and the identification library are those of 0.5.3.

## Comparison with published phylogenies

* `validate_phylogenies()` and `plot_tree_space()` give the tree of Carruthers *et al*. (2026), named `Carruthers` in `trees_mapping_list`, the sixth colour of the categorical palette, its own shape and the default legend label `"Carruthers et al. (2026)"`.
* Tutorial 4 compares the focal tree with five published trees. The new one is `best.tre_dated` of the OSF deposit of Carruthers *et al*. (2026) (<https://doi.org/10.17605/OSF.IO/9TBHA>); the tutorial gives its md5, prunes it to its clade of Cactaceae, shortens its labels from `Genus_species__<NCBI taxid>` to `Genus_species`, and gives the md5 of the pruned file the reference run read. Like the other reference trees, it is not distributed with the package.
* Under `pruning_strategy = "common_set"`, the species shared by every tree of set A3 are those of A1, so panel A3 of the tree-space figure has the points of panel A1. Its label in `SUPP_tree_space_coordinates.csv` now says so; the pairwise distances of A3 keep all the species of the shared genera and their label.
* A sixth tree changes the species shared by all trees, so `phylocactus_table_tree_distances.csv` and `phylocactus_table_tree_space.csv` were regenerated for every tree, not only for the new one; the distances of set A2 (the focal tree and the trees of Thompson *et al*. 2024 and Amaral *et al*. 2022) do not involve the new tree and are unchanged.

## Text and data

* Tutorial 2 no longer gives a running time of the maximum-likelihood search on an Apple M2 Pro, which was not measured on the current supermatrix; it says that the module `raxml-ng/1.1.0-mpi-zen4-n` of the cluster runs RAxML-NG 1.2.2-master, the version in the logs of the run; and it no longer attributes the maximum clade credibility target tree of `TreeAnnotator` to Maurin (2020), who used 0% burn-in and mean node heights.
* `cactus_constraints.csv` places Rhipsalideae at level 3 (`Annot_level_3`) in Core Cactoideae II, where `build_constraint_scaffold()` places the tribe. The constraint tree built from the table is unchanged.
* `cactus_constraints.csv` gives level 4 (`Annot_level_4`) names to lineages of Cacteae that had none, with the names of the lineages in Figure 2 of Guerrero *et al*. (2019): Mammilloid clade (*Mammillaria*, *Coryphantha*, *Cochemiea*, *Pelecyphora*, *Cumarinia*), *Turbinicarpus* (with *Ariocarpus* and *Kadenicarpus*), *Rapicactus* (with *Lophophora*, *Acharagma* and *Obregonia*), *Ferocactus* (with *Thelocactus*, *Stenocactus*, *Kroenleinia* and *Leuchtenbergia*) and *Astrophytum*, *Echinocactus* (with two species of *Homalocephala*); *Eulychnia* takes *Austrocactus*. Each group is a clade of the chronogram, so the chronogram and the maximum-likelihood tree of Tutorial 3 now label it; *Coryphantha*, *Cochemiea*, *Ariocarpus*, *Astrophytum* and *Stenocactus* are no longer labelled on their own. Only `Clade` builds the constraint tree, and the tree built from the table is unchanged.
* In Core Cactoideae I, *Selenicereus grandiflorus* no longer carries the level 4 name Hylocereeae, so that the tribe, a clade of the chronogram and of the maximum-likelihood tree without that species, is labelled in Tutorial 3. *Leptocereus* and three species of *Armatocereus* form the level 4 group *Leptocereus*, *Armatocereus*, a clade of both trees; it is part of the lineage that Figure 2 of Guerrero *et al*. (2019) names after five genera, which is not monophyletic in these trees, and its name is not taken from that figure.
* `add_clade_labels_by_level()` sets a label that ends in "clade" in roman type, as it does for suprageneric names.
* `add_outgroup_bar()` (new) draws a bar and a label over each run of consecutive terminals outside a given ingroup, at the position of the aligned bars of `add_clade_labels_by_level()`. Figures 6 and 8 of Tutorial 3 use it to mark the outgroup, which is not a clade of the rooted trees. In the same figures a level 4 label of two names is set on two lines, and the scale bar of Figure 6 is placed below the tree.

# PhyloCactus 0.5.3

This version closes the audit of 26 October 2026. The trees and tables of the reference run of 13 to 15 October 2026 are those distributed since 0.5.0, except that the registry of species loses columns and three tables are added. The chronogram was dated again (Module 10, 30 October 2026) under three calibration schemes, each one fixed age on the crown of the ACP clade and no calibration inside Cactaceae; the maximum-likelihood search, the bootstrap and the temporal replicates were not run again.

## Dating

* `calibrations_bounds.csv` holds three calibration schemes for `ACP_root`, each a fixed age from a published chronogram: scheme 1, the young tree of Zuntini *et al*. (2024), 53.22 Ma, the stem of Cactaceae, which on this tree is the ACP crown (the active row); scheme 2, Ramírez-Barahona *et al*. (2020), 41.82 Ma, the median height of the stem of Portulacaceae in the BEAST analysis `RC_complete`; scheme 3, Arakaki *et al*. (2011), 37.0 Ma. The ACP and ACPT crowns of Wu *et al*. (2026) are added as references.
* Module 10 was run for each scheme from the configuration of the distributed run with only the calibration lines changed: priming, cross-validation, the maximum-likelihood tree and the same 100 temporal replicates, at the smoothing each cross-validation selects (1e-05 in all three), with seed 79992967, and TreeAnnotator v2.7.7 of BEAST 2. The distributed chronogram and its summary are those of scheme 1: the crown of Cactaceae is 27.33 Ma (23.51 to 30.37 Ma over the replicates). `PROVENANCE.csv` records the scheme and the smoothing. Tutorial 2 gives the three schemes, why scheme 1 is the main one, and their ages.

## What left `inst/extdata`, and how to obtain it

* The content of the IUCN Red List assessment of each species. `phylocactus_table_species.csv` keeps 3 of its 36 `iucn_*` columns (`iucn_found`, `iucn_category` and `iucn_category_name`); the other 33, with the texts, ranges, threats and habitats of each assessment, are no longer distributed, as sections 3 and 4 of the Terms of Use of the IUCN Red List require. A user obtains them with their own token through `get_iucn_data()`, as in Tutorial 3, which writes them to `9_Visualization/tables/TABLE_species_all_data.csv`. The new `summarise_iucn_registry()` reduces that table to the category of each species and to two count tables, distributed as `phylocactus_table_iucn_endemic_by_country.csv` (`country`, `iucn_category`, `n`) and `phylocactus_table_iucn_classes.csv` (`field`, `class`, `n`: threats, habitats and endemic species in the tree), both with the columns `redlist_version` and `access_date` (version 2026-1, read on 15 October 2026). Tutorial 3 draws Figures 14 to 17 and the endemism row of Table 1 from them, with the counts of 0.5.2.
* The checklist of Caryophyllales.org (Korotkova *et al*. 2021), which may be used but not redistributed. A user requests it from the team of Caryophyllales.org and gives its path in `~/.Renviron` as `PHYLOCACTUS_CHECKLIST`. The new `phylocactus_checklist()` returns that path, reports the md5 of the file and warns when it is not the file of the reference run, because species counts then differ from those of the tutorials. The tutorials and the default `checklist_path` of `assemble_barcoding_dataset()` and `extract_barcoding_plastome_loci()` call it; the functions that read the checklist take a path as before, so no analysis changes.
* The four published trees that Tutorial 4 compares with the focal tree, `reference_amaral_2022.tree`, `reference_thompson_2024.tree`, `reference_zuntini_2024.tree` and `reference_devos_2025_astral.tree`. Tutorial 4 lists, for each, where the file is, its name there and the md5 of the file the reference run read, and gives the pruning of the Kew Tree of Life 4.0 to Cactaceae; a user downloads them into a folder named by `PHYLOCACTUS_REFERENCE_TREES`. Figure 1 of Tutorial 4 is drawn by the new `plot_tree_space()` from the tree-space coordinates of the reference run, distributed as `phylocactus_table_tree_space.csv`, and its table of distances from `phylocactus_table_tree_distances.csv`, so the page is built without the trees and without the checklist.
* The authors, ranks, taxa and identifiers of the checklist in `phylocactus_table_species.csv` (`checklist_fullName`, `checklist_author`, `checklist_rank`, `checklist_taxon` and `checklist_uuid`), which the registry carried for 1965 species; `summarise_iucn_registry()` now leaves them out, and the registry keeps the names and the result of the match to the checklist.
* These files remain in the earlier tags of the repository and in their archives; the history is not rewritten.

## Provenance of the reference run

* `PROVENANCE.csv` records the seeds of the reference run, read from its logs: 1111 for the maximum-likelihood search; 2150, 4300, 6450 and 8600 for the four bootstrap tasks; 87786 for the temporal bootstrap; 79992967 for treePL. It records the versions of MAFFT, ModelTest-NG, phylotaR, BLAST+, DECIPHER and TreeAnnotator beside those of RAxML-NG, treePL and R; the version of the package installed when the run was first published, 0.5.0 (`run_package_version`), beside the one that republished it (`package_version`); and the name and md5 of the checklist the run read.

## Code

* `run_treePL_cv()` counts the random subsamples of the cross-validation that failed, read from `cv_<label>.log` (an increment of the running chi-square of 1e+30 or more), and warns when one or more failed at a smoothing value of 0.1 or below. The warning names each such value with the number of subsamples that failed there, says whether the selected value is among them, and names the value with the lowest chi-square among those at which no subsample failed; the result returns `cv_failures` and `cv_smoothing_without_failures`. `automate_treePL()` repeats the check on a run read from disk.
* `run_treePL_cv()` and `automate_treePL()` gain `cv_outlier_factor` (default 5). A random subsample of the cross-validation that fails, or adds more than that many times the median increment of the subsamples, at any smoothing value of 0.1 or below, is left out of the selection, and the values are compared on one set of subsamples. A warning names the subsamples left out; both curves are written to `cv_<label>_selection.csv` and the account to `cv_<label>_selection.txt`. `cv_outlier_factor = NULL` selects on every subsample. In the re-dating, one subsample of each scheme had moved the selection to 1e-09 to 1e-11. The default of 5 was set after the curves of the re-dating had been read. When the failed subsample is left out, the warning on failed subsamples says so, and the factor must be `NULL` or one number above 1.
* `report_smoothing_sensitivity()` gains `nodes`, a named list of terminal sets whose ages are reported beside those of the calibrated nodes. Under one fixed-age calibration the table otherwise holds one constant age; Module 10 passes the root and the crowns of Cactaceae, Cactoideae, Cacteae and Opuntioideae.
* New exported functions: `phylocactus_checklist()`, `summarise_iucn_registry()` and `plot_tree_space()` (above).

## Text

* The README, the tutorials, their scripts and the help pages state the numbers of the distributed run, recomputed from `inst/extdata`: among them the bipartitions the constraint fixes and those the data estimate, the support of the nodes of the maximum-likelihood tree, and the counts of the curation in Tutorial 1, which its script now states as the vignette does. Help pages no longer carry numbers of earlier runs.
* The intervals of the dated replicates are called HPD only where TreeAnnotator names its column. Tutorial 2 states how many bootstrap replicates were dated, which branches treePL lengthens, and that the fixed age of scheme 2 is the median height of its source.
* The help page of `assemble_ingroup_phylotar()` says which of the distributed exclusions state a ground, and, for those that read only "manual conflict resolution", what the metadata of the mining show: 43 of the 44 found there were placed in the clusters of more than one locus. Tutorial 1 no longer states that other groups can be studied without modifying the package.
* The README lists the external programs with the versions of the reference run; `SystemRequirements` names TreeAnnotator.

# PhyloCactus 0.5.2

* Step 7 in parallel (`classify_barcoding_folds(workers = )`) stops with an error that names the worker when a worker process ends without writing its status (memory exhausted, a killed process, a crash of R), instead of waiting without end. A worker that has written no process id after `getOption("PhyloCactus.worker_grace_seconds")` seconds (by default three polls, never under 60) counts as ended. The analyses and the distributed results are those of 0.5.0.

# PhyloCactus 0.5.1

Documentation, metadata and the website. The analyses and the distributed results are those of 0.5.0.

* `DESCRIPTION` describes both branches of the package and lists the external programs it calls (`SystemRequirements`); `citation("PhyloCactus")` takes the new title.
* The help pages of `assemble_ingroup_phylotar()` and `assemble_outgroup_phylotar()`, and Tutorial 1, document the cut of records over 2000 bases introduced in 0.5.0, with the tables that list them.
* Tutorial 2 reads the weak support of the backbone: soft polytomies from missing data, not evidence of a simultaneous divergence.
* Tutorial 6 explains the two version numbers of its reports.
* The homology check of `integrate_and_clean_markers()` states, in its help page and its warning, the interval of the reference run of 0.5.0 (0.26 to 0.58 at k = 20).
* The website no longer serves the pages of the functions removed in 0.5.0: each deployment replaces the previous one.
* The checklist of Korotkova et al. (2021) is cited wherever it is named; the comments of the source code describe the code.
* The 0.5.0 entry below is completed with changes it left out.

# PhyloCactus 0.5.0

## Molecular diagnostic section (DNA barcoding)

* A new branch builds a reference library of Cactaceae from the Sanger marker records of the phylotaR mining of the phylogeny, and optionally the user's own sequences (`assemble_barcoding_dataset()`, `curate_barcoding_markers()`, `screen_barcoding_markers()`, `finalize_barcoding_library()`). Names are reconciled against the checklist of Caryophyllales.org (Korotkova et al. 2021); records over 2000 bases are cut to the extent of the shorter records of their locus; identical sequences of a species are collapsed. The library of this version holds 15 loci. `extract_barcoding_plastome_loci()` crops the loci of the library from GenBank plastomes as extra records; `barcoding_outgroup_taxids()` gives the six NCBI Taxonomy IDs of the outgroup the controls query the library with.
* The library is validated before it is used: folds that ask for the species (scheme E) or, with the species left out, for the genus (scheme G) (`build_barcoding_folds()`); the barcode gap per locus (`analyze_barcode_gap()`); classification of the folds by nearest neighbour and by IdTaxa, on one machine in parallel processes (`classify_barcoding_folds(workers = )`) or on a computing cluster (`classify_barcoding_folds()`, `generate_barcoding_job_scripts()`, `merge_barcoding_chunks()`); negative controls for leakage, outgroup queries and resubstitution (`run_barcoding_controls()`); thresholds (`sweep_barcoding_threshold()`); metrics per class (`summarise_barcoding_metrics()`, `summarise_barcoding_genus_discrimination()`).
* A sample is identified locus by locus, in three states (species, genus, no answer), never combined into one verdict, from sequences or a GenBank flat file (`identify_barcoding_query()`), a sample sheet (`identify_barcoding_samples()`), an assembly (`identify_barcoding_assembly()`, with `minimap2` and `samtools`) or paired Illumina reads (`identify_barcoding_reads()`, with `fastp` and GetOrganelle). A species is named only above a species threshold per locus, set from the folds of scheme G; `rbcL` answers at genus rank only (argument `genus_only_loci`), because its species answers were right in 11 of 34 on the library of this version. Each run writes a self-contained HTML report (`report_barcoding_identification()`).
* Tutorial 5 is now the tutorial of the branch: the script runs steps 1 to 11 and draws the figures of the library and its validation, with its limitations. The plastome MN517611.1 is distributed in `extdata` as the example query.
* Tutorial 6 identifies public samples of each kind of input (a Sanger record, specimens of several loci, GenBank plastomes, a sample sheet, a genome assembly and Illumina reads) with the library of this version, and gives the rates measured with it; Tutorial 5 keeps the library and its validation. Two new functions draw an identification from the tables of its report: `plot_barcoding_identification()` (the answer per locus and the sampling of the named genus) and `plot_barcoding_identification_summary()` (several identifications side by side). The function reference is now Tutorial 7; its former address is redirected.

## Changes that also reach the phylogeny

* The reference run is regenerated and every distributed artifact republished: a supermatrix of 1218 terminals, 13 700 sites and 13 partitions, analysed as 1208 terminals once identical sequences are collapsed; the chronogram is dated at a smoothing of 1e-05, an interior minimum of the cross-validation, and its root is 45.53 Ma. `MANIFEST.csv` gains `expected_tips`, the number of terminals each distributed tree must have.
* `assemble_ingroup_phylotar()` and `assemble_outgroup_phylotar()` cut every record over 2000 bases to the extent the shorter records of its marker cover, as the identification branch does, and list the long records in `TABLE_LONG_RECORDS_CUT_INGROUP.csv` and `TABLE_LONG_RECORDS_CUT_OUTGROUP.csv`.
* The family filter of `validate_phylogenies()` keeps the species of Cactaceae in the checklist (Korotkova et al. 2021) instead of the species of the constraint table; the tree distances of the reference run do not change.
* 31 species of the new mining are added to `cactus_constraints.csv`, so that the constraint places every terminal of the supermatrix.
* One style for the figures of every tutorial: `theme_phylocactus()`, `phylocactus_colours()` and `phylocactus_panel_tags()`.
* The phylotaR mining of both branches keeps sequences of 100 to 5000 bases, and the outgroup of the phylogeny is mined through the same function as the ingroup.
* The manual exclusions are keyed by locus and accession instead of cluster number, which a new mining renumbers; the logs report the exclusions actually applied.
* 21 `ITS` records annotated in GenBank as pseudogenes are added to `manual_exclusions_ingroup.csv`, read by both branches.
* `genes_map.csv` gains variant names of loci found in the new mining.
* The subject of every email notification of the package carries a cactus emoji.

## Removals

* Removed `infer_gene_trees()`. No stage of the pipeline, no tutorial script and no test called it, and on its documented input (one sequence per species) its species monophyly table was always empty. Its `checklist_path` argument had no effect. Per-locus monophyly is planned work, with no version set. `extract_species_binomial()`, which it shared with the marker integration step, is kept.
* Removed the placeholder `evaluate_dna_barcoding()` and the `barcoding` documentation topic. The placeholder only signalled an error, and its arguments belonged to the implementation withdrawn in 0.4.2. The molecular diagnostic section was rebuilt as a separate branch with its own functions (above).

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
* The reference chronogram is dated at `smooth = 1e-12`, the lowest cross-validation error on a grid from 1e3 to 1e-14. The error curve levels off below about 1e-06, so this value is one of several that fit the data about equally well. `run_treePL_cv()` now classifies a cross-validation curve as an interior minimum, an edge of the grid or a plateau, and reports the result; the smoothing sensitivity table of Tutorial 2 spans the plateau.
* `run_treePL_cv()` and `automate_treePL()` gain `cv_nthreads`, the number of threads used for cross-validation (default `1L`), and their default `cvstop` is lowered from 1e-8 to 1e-14.
* The tree-space figure of `validate_phylogenies()` identifies each tree by color and shape in a single legend labelled by source, with *et al.* in italics, and no longer draws point labels inside the panels. Labels are set with the new `tree_labels` argument. `ggrepel` is no longer a dependency.
* Fixed the edge test of the cross-validation curve. It compared the selected smoothing value with the ends of the grid through `all.equal()`, whose absolute tolerance treats any two values below about 1e-8 as equal, so a minimum at 1e-12 on a grid ending at 1e-14 was reported as an edge. The position of the minimum is now read from the sorted grid.
* When `automate_treePL()` resumes from a cached maximum-likelihood run, it classifies the cross-validation curve on disk and reports its shape; a file it cannot classify is reported and does not stop the run.
* `phylocactus_table_tree_distances.csv` as first released in 0.4.5 was the 0.4.4 table: Tutorial 4 read the focal tree from the installed package, which still held the 0.4.4 chronogram. The table is replaced with distances computed on the 0.4.5 chronogram, and Tutorial 4 now takes the focal tree from the run (`8_Dating/dated_summary_hpd.tree`), falling back to the distributed chronogram when that file is absent.
* The full-size tree and chronogram in `extended_trees/` are replaced with those of the 0.4.5 run (1,024 tips; the previous files, from 0.4.4, had 1,023).
* `phylocactus_table_species.csv` is regenerated from a new IUCN Red List query (version 2026-1). One record changes: the assessment details of *Stenocactus coptonogonus*, previously empty, are now filled; categories and counts are unchanged.
* Documentation corrections: the IUCN token instructions pointed to the retired API v3 and now use `rredlist::rl_use_iucn()`; the description of `ASTRAL-III` and of the effect of rooting in `validate_phylogenies()` (distances and tree space are computed on unrooted trees); the reference for the Different Phylogenetic Information distance (Smith 2020, *Bioinformatics*) and the DOI of Wright (2024); the page size of the full-size tree figures; the links to `extended_trees/`, which pointed to a branch that does not exist. The tutorials were also edited for concision.
* `DESCRIPTION` declares `biocViews:`, so that `remotes` resolves the Bioconductor dependencies.
* Added `send_run_notification()` and a `notify` argument to `automate_treePL()`, which emails the status, timings and output paths when a local dating run ends. Credentials live outside the package and a failed notification never fails the run.
* `run_alignment_pipeline()` puts every sequence of a marker on the same strand before alignment (`fix_strand`), because some GenBank accessions are deposited as reverse complements.
* `run_marker_screening()` flags markers whose ingroup and outgroup sequences share almost no k-mers (`no_detectable_homology`): two sequences can carry the same gene name and come from different regions.
* Clusters whose GenBank description matches no rule in `genes_map.csv` are reported with a warning instead of being dropped without notice.
* The bootstrap replicates are written without a support metric; FBP or TBE is chosen when the replicates are summarized with `map_branch_supports()`.
* `treePL` runs report when the gradient optimization did not improve on the simulated annealing, in which case the chronogram is not an optimized solution.
* The documentation of `rescale_tree()` states what the scaling factor changes: the scale of the rate parameters and the numerical conditioning of the optimization, not which branches `treePL` shortens.
* `evaluate_dna_barcoding()` is a placeholder that signals an error. The DNA barcoding section is under construction; its scope and open questions are set out in Tutorial 5.
* The `checklist_path` argument of `infer_gene_trees()` is retained and remains reserved for the molecular diagnostic section, which is still under construction.

# PhyloCactus 0.4.4

* Updated outgroup sampling in data preparation and tree rooting to include Talinaceae (*Talinum*, *Talinella*) alongside Portulacaceae and Anacampserotaceae.
* Standardized primary clade support metric to Felsenstein Bootstrap Proportions (FBP) across documentation, functions, and collapsing workflows (`integrate_publication_tree()`), retaining Transfer Bootstrap Expectation (TBE) as a complementary diagnostic metric.
* Added support for identifying and distinguishing topologically constrained nodes in tree visualization (`classify_constrained_nodes()`).
* Replaced shell wrapper execution of treePL with direct system calls (`run_treePL_direct()`) and native R cross-validation optimization (`run_treePL_cv()`).
* Retired legacy `run_treePL()` wrapper function.
* Improved documentation, metadata manifests, and error handling across all analytical modules.
