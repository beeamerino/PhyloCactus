#' Collapse a Rooting Outgroup Specification into a RAxML-NG Argument
#'
#' `RAxML-NG` accepts `--outgroup` as a comma-separated list of terminals. Accepting a vector here
#' and joining it at the call site keeps every function in the pipeline taking the same argument
#' shape as `ape::root()`, so a single `rooting_outgroup` object can be declared once and passed
#' through the maximum-likelihood search, the bootstraps and the dating step without reshaping.
#'
#' @param outgroup Character vector, or `NULL`.
#' @return A single comma-separated string, or `NULL` when nothing usable was supplied.
#' @noRd
.format_outgroup <- function(outgroup) {
  if (is.null(outgroup)) return(NULL)
  outgroup <- outgroup[!is.na(outgroup) & nzchar(outgroup)]
  if (length(outgroup) == 0L) return(NULL)
  paste(outgroup, collapse = ",")
}

#' @noRd
.check_cli_exit <- function(status, tool_name, args = NULL) {
  exit_code <- if (is.integer(status)) status else attr(status, "status")
  if (is.null(exit_code)) exit_code <- 0L
  if (!identical(exit_code, 0L)) {
    stderr_msg <- if (is.character(status)) paste(utils::tail(status, 20), collapse = "\n") else ""
    stop(tool_name, " failed (exit ", exit_code, ").",
         if (nzchar(stderr_msg)) paste0("\nOutput:\n", stderr_msg) else "",
         call. = FALSE)
  }
  invisible(status)
}

#' Preprocess Partitions and Validate Alignment Syntax
#'
#' Validates PHYLIP alignment syntax and partition file coordinates using `RAxML-NG` (Kozlov *et al.*, 2019).
#' Verifies site ranges, formatting compatibility, and data integrity prior to substitution model evaluation.
#'
#' @param phy_matrix Character. Path to input PHYLIP supermatrix file.
#' @param part_file Character. Path to input partition mapping text file in RAxML-style format
#'   (`PARTITION_raxml_style.txt`, exported by [run_concatenation_pipeline()]).
#' @param raxml_path Character. System command or full path to executable `RAxML-NG` binary.
#' @param output_dir Character. Output directory for validated partition outputs. Defaults to `dirname(phy_matrix)`.
#' @param prefix Character. Basename prefix for every file written by this step. Defaults to `"cactus"`.
#' @param force_check Logical. Re-run validation even if the validated partition file already exists? Defaults to `FALSE`.
#' @param model_handling Character. Model field written to the validated partition map. `"force_dna"`
#'   (default) writes the datatype token `DNA`, leaving the substitution model to be selected by
#'   [run_modeltest_ng()]. `"preserve"` keeps whatever model string `RAxML-NG` emitted, which is only
#'   appropriate when the partition map is passed straight to `RAxML-NG` without model selection.
#'   `ModelTest-NG` cannot parse the `RAxML-NG` model syntax and aborts on it.
#' @return Character path to the validated partition map file ready for model evaluation.
#' @seealso [run_concatenation_pipeline()] for the partition map this function consumes,
#'   and [run_modeltest_ng()] for the model selection step that consumes its output.
#' @references
#' Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and
#' user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics*, 35(21), 4453-4455.
#' \doi{10.1093/bioinformatics/btz305}
#' @export
preprocess_partitions <- function(phy_matrix, part_file, raxml_path, output_dir = dirname(phy_matrix), prefix = "cactus", force_check = FALSE, model_handling = c("force_dna", "preserve")) {
  model_handling <- match.arg(model_handling)

  if (Sys.which(raxml_path) == "") {
    stop("Executable '", raxml_path, "' not found in your system's PATH.\n",
         "Please ensure RAxML-NG is installed and available, or provide the full absolute path.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # The filename records what the file is (a validated partition map), not which internal RAxML-NG
  # prefix happened to produce it. Whether terminals were collapsed is reported by the return
  # attributes and by the presence of the reduced matrix, not asserted by the name.
  clean_part   <- file.path(output_dir, paste0(prefix, "_partitions_validated.txt"))
  prefix_check <- file.path(output_dir, prefix)

  if (!force_check && file.exists(clean_part)) {
    message("CACHE: Partition map already validated, loading: ", clean_part)
    cached_reduced_phy <- paste0(prefix_check, ".raxml.reduced.phy")
    attr(clean_part, "was_reduced")  <- file.exists(cached_reduced_phy)
    attr(clean_part, "analysed_phy") <- if (file.exists(cached_reduced_phy)) cached_reduced_phy else phy_matrix
    return(clean_part)
  }

  # Run system raxmlcheck validation
  aln_check <- system2(
    raxml_path,
    args = c(
      "--check",
      "--msa", shQuote(phy_matrix),
      "--model", shQuote(part_file),
      "--prefix", shQuote(prefix_check),
      "--threads", "4"
    ),
    stdout = "",
    stderr = ""
  )
  .check_cli_exit(aln_check, "RAxML-NG --check")

  # RAxML-NG writes *.raxml.reduced.* only when it actually drops duplicate terminals or gap-only
  # columns. Both branches are legitimate outcomes and the caller needs to know which one occurred,
  # because a reduced matrix has fewer tips than the supermatrix exported by Module 6.
  reduced_part_file <- paste0(prefix_check, ".raxml.reduced.partition")
  reduced_phy_file  <- paste0(prefix_check, ".raxml.reduced.phy")
  was_reduced <- file.exists(reduced_part_file)

  if (was_reduced) {
    lines <- readLines(reduced_part_file, warn = FALSE)
    message("RAxML-NG reduced the matrix. Analysed alignment: ", reduced_phy_file)
  } else {
    message("No reduction applied. Analysed alignment: ", phy_matrix)
    lines <- readLines(part_file, warn = FALSE)
  }

  # RAxML-NG expands a bare "DNA" datatype token into its own default model (GTR+FC+G4m+B) and
  # writes that expansion into the reduced partition file. Passing that expansion to ModelTest-NG
  # would both pre-empt the model selection ModelTest-NG exists to perform and crash its partition
  # parser, so the datatype token is restored here.
  if (model_handling == "force_dna") {
    out_lines <- sub("^[^,]+,", "DNA,", lines)
  } else {
    out_lines <- lines
  }

  writeLines(out_lines, clean_part)

  attr(clean_part, "was_reduced")   <- was_reduced
  attr(clean_part, "analysed_phy")  <- if (was_reduced) reduced_phy_file else phy_matrix
  return(clean_part)
}

#' Evaluate Nucleotide Substitution Models via ModelTest-NG
#'
#' Evaluates nucleotide substitution model fit per predefined supermatrix partition using `ModelTest-NG` (Darriba *et al.*, 2020).
#' Selecting optimal substitution models under the AICc criterion controls for mutational rate heterogeneity across molecular locus alignments,
#' mitigating systematic long-branch attraction (LBA) bias during maximum-likelihood inference.
#'
#' @param modeltest_exec_path Character. System command or full path to executable `ModelTest-NG` binary.
#' @param aln_file Character. Path to validated PHYLIP supermatrix file.
#' @param part_file Character. Path to the validated partition map written by [preprocess_partitions()].
#'   Its model field must carry the datatype token `DNA`; `ModelTest-NG` cannot parse `RAxML-NG`
#'   model syntax such as `GTR+FC+G4m+B` and aborts with a segmentation fault on it.
#' @param prefix Character. Output filename prefix. Defaults to `"MODELTEST_cactus_phylo"`.
#' @param threads Integer. Number of processing threads. Defaults to `4`.
#' @return Character path to the resulting partition file containing selected model parameters (`.part.aicc`).
#' @references
#' Darriba, D., Posada, D., Kozlov, A. M., Stamatakis, A., Morel, B., & Flouri, T. (2020). ModelTest-NG: a new and
#' scalable tool for the selection of DNA and protein evolutionary models. *Molecular Biology and Evolution*, 37(1), 291-294.
#' \doi{10.1093/molbev/msz189}
#' @examples
#' \dontrun{
#' run_modeltest_ng(
#'   modeltest_exec_path = "modeltest-ng",
#'   aln_file = "ALIGNMENT_supermatrix.phy",
#'   part_file = "cactus_partitions_validated.txt"
#' )
#' }
#' @export
run_modeltest_ng <- function(modeltest_exec_path, aln_file, part_file, prefix = "MODELTEST_cactus_phylo", threads = 4) {
  
  if (Sys.which(modeltest_exec_path) == "") {
    stop("Executable '", modeltest_exec_path, "' not found in your system's PATH.\n",
         "Please ensure ModelTest-NG is installed and available, or provide the full absolute path.")
  }

  if (!file.exists(part_file)) {
    stop("Partition map not found: ", part_file, call. = FALSE)
  }

  # ModelTest-NG 0.1.7 reports an unparseable model field and then crashes while unwinding, so the
  # exit status reads as a segmentation fault rather than a format error. Rejecting the input here
  # turns that into a message that names the actual problem.
  part_lines <- readLines(part_file, warn = FALSE)
  part_lines <- part_lines[nzchar(trimws(part_lines))]
  bad_models <- part_lines[!grepl("^\\s*(DNA|BIN|PROT|MORPH)\\s*,", part_lines, ignore.case = TRUE)]
  if (length(bad_models) > 0) {
    stop("ModelTest-NG only accepts a datatype token in the model field of a partition map.\n",
         "These lines in '", part_file, "' carry something else:\n  ",
         paste(utils::head(bad_models, 5), collapse = "\n  "),
         "\nRe-run preprocess_partitions() with model_handling = \"force_dna\".",
         call. = FALSE)
  }

  # Invocation parameters setup
  args <- c(
    "--datatype", "nt",
    "--input", shQuote(aln_file),
    "--partitions", shQuote(part_file),
    "--output", shQuote(prefix),
    "--processes", as.character(threads),
    "--template", "raxml"
  )
  
  message("Executing ModelTest-NG over partitions database...")
  exit_status <- system2(
    command = modeltest_exec_path,
    args = args,
    stdout = TRUE,
    stderr = TRUE
  )
  .check_cli_exit(exit_status, "ModelTest-NG")
  
  # Returning an unverified path would let a silently failed run propagate into the ML search as a
  # missing-file error several stages later.
  expected_out <- paste0(prefix, ".part.aicc")
  if (!file.exists(expected_out)) {
    stop("ModelTest-NG exited cleanly but did not write '", expected_out, "'.\n",
         "Check the ModelTest-NG log at '", paste0(prefix, ".log"), "'.",
         call. = FALSE)
  }
  return(expected_out)
}

#' Synthesize Multifurcating Monophyly Constraint Scaffold
#'
#' Constructs a Newick multifurcating constraint tree enforcing monophyly of established higher taxonomic ranks (e.g., subfamilies, tribes).
#' Constrained maximum-likelihood searches restrict branch topology exploration to scientifically verified monophyletic backbone clades,
#' preventing aberrant tree topologies when analyzing sparse supermatrices.
#'
#' @param alignment_path Character. Path to input PHYLIP supermatrix alignment file.
#' @param constraints_csv_path Character. Path to taxonomy CSV table mapping species binomials to taxonomic ranks.
#' @param output_dir Character. Directory path to save generated constraint scaffold file. Defaults to `dirname(alignment_path)`.
#' @param prefix Character. Run prefix for the generated files. Defaults to `"cactus"`, giving
#'   `cactus_constraints.tree`. Keep it equal to the prefix used by the rest of the run so that
#'   [resolve_run_paths()] locates the scaffold.
#' @return Character string path to the saved Newick constraint tree file (`<prefix>_constraints.tree`).
#' @examples
#' \dontrun{
#' build_constraint_scaffold(
#'   alignment_path = "ALIGNMENT_supermatrix.phy",
#'   constraints_csv_path = "Cactaceae_taxonomy.csv"
#' )
#' }
#' @export
build_constraint_scaffold <- function(alignment_path, constraints_csv_path, output_dir = dirname(alignment_path),
                                      prefix = "cactus") {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  out_missing     <- file.path(output_dir, paste0("TABLE_", prefix, "_alignment_taxa_missing_from_constraints.csv"))
  out_dup_clade   <- file.path(output_dir, paste0("TABLE_", prefix, "_constraint_duplicate_clade_membership.csv"))
  out_clade_sizes <- file.path(output_dir, paste0("TABLE_", prefix, "_constraint_clade_sizes.csv"))
  out_tree        <- file.path(output_dir, paste0(prefix, "_constraints.tree"))
  out_log         <- file.path(output_dir, paste0("LOG_", prefix, "_constraint_tree.txt"))
  
  log_message <- function(...) {
    msg <- paste0(...)
    timestamped <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
    message(timestamped)
    cat(timestamped, "\n", file = out_log, append = TRUE)
  }
  
  if (file.exists(out_log)) file.remove(out_log)
  log_message("Starting constraint tree construction.")
  
  # 1. VALIDATE AND READ ALIGNMENT
  if (!file.exists(alignment_path)) stop("Missing alignment: ", alignment_path)
  if (!file.exists(constraints_csv_path)) stop("Missing constraints table: ", constraints_csv_path)
  
  aln <- tryCatch(
    ape::read.dna(alignment_path, format = "sequential"),
    error = function(e) stop("Could not read alignment at ", alignment_path, " | ", conditionMessage(e))
  )
  
  species_alignment <- rownames(aln)
  if (is.null(species_alignment) || length(species_alignment) == 0) {
    stop("No taxa found in alignment: ", alignment_path)
  }
  species_alignment <- sort(unique(species_alignment))
  log_message("Alignment taxa loaded: ", length(species_alignment))
  
  # 2. READ CONSTRAINTS
  cf_all <- tryCatch(
    utils::read.csv(constraints_csv_path, stringsAsFactors = FALSE),
    error = function(e) stop("Could not read constraints table at ", constraints_csv_path, " | ", conditionMessage(e))
  )
  
  required_cols <- c("Specie_name", "Clade", "Subfam", "Family")
  missing_cols <- setdiff(required_cols, colnames(cf_all))
  if (length(missing_cols) > 0) stop("Missing required columns in constraints table: ", paste(missing_cols, collapse = ", "))
  
  cf_all$Specie_name <- trimws(cf_all$Specie_name)
  cf_all$Clade       <- trimws(cf_all$Clade)
  cf_all$Subfam      <- trimws(cf_all$Subfam)
  cf_all$Family      <- trimws(cf_all$Family)
  
  if (any(is.na(cf_all$Specie_name)) || any(cf_all$Specie_name == "")) stop("Constraints table contains empty or NA values in Specie_name.")
  if (any(is.na(cf_all$Clade)) || any(cf_all$Clade == "")) stop("Constraints table contains empty or NA values in Clade.")
  
  log_message("Constraint rows loaded: ", nrow(cf_all))
  
  # 3. COVERAGE AND CONSISTENCY AUDITS
  missing_taxa <- setdiff(species_alignment, cf_all$Specie_name)
  if (length(missing_taxa) > 0) {
    utils::write.csv(data.frame(Species_missing = missing_taxa), out_missing, row.names = FALSE)
    stop(length(missing_taxa), " taxa in the alignment are missing from constraints table. ",
         "Missing list written to: ", out_missing)
  }
  
  log_message("All alignment taxa are represented in constraints table.")
  
  # Check duplicates
  cf_align <- cf_all[cf_all$Specie_name %in% species_alignment, ]
  cf_unique <- unique(cf_align[, c("Specie_name", "Clade")])
  dup_counts <- table(cf_unique$Specie_name)
  true_dups <- names(dup_counts[dup_counts > 1])
  
  if (length(true_dups) > 0) {
    dup_detail <- cf_unique[cf_unique$Specie_name %in% true_dups, ]
    dup_detail <- dup_detail[order(dup_detail$Specie_name, dup_detail$Clade), ]
    utils::write.csv(dup_detail, out_dup_clade, row.names = FALSE)
    stop("Species assigned to multiple clades detected. Details written to: ", out_dup_clade)
  }
  
  log_message("No duplicate clade membership detected among alignment taxa.")
  
  # 4. FILTER TO ALIGNMENT TAXA
  cf <- unique(cf_align)
  cf$Genus <- sub("_.*$", "", cf$Specie_name)
  
  genus_to_constraint <- character(0)
  target_genera <- sort(intersect(unique(cf$Genus), genus_to_constraint))
  log_message("Genera with additional genus-level constraints: ",
              ifelse(length(target_genera) == 0, "none", paste(target_genera, collapse = ", ")))
  
  genus_newick_lut <- NULL
  if (length(target_genera) > 0) {
    # Emulate group_by collapse
    spp_by_genus <- split(cf$Specie_name, cf$Genus)
    spp_by_genus <- spp_by_genus[names(spp_by_genus) %in% target_genera]
    genus_newick_lut <- sapply(spp_by_genus, function(spp) {
      paste0("(", paste(sort(unique(spp)), collapse = ","), ")")
    })
  }
  
  # Helper
  extract_required_clade <- function(df, clade_name) {
    res <- sort(unique(df$Specie_name[df$Clade == clade_name]))
    if (length(res) == 0) stop("Required clade missing or empty in constraints table: ", clade_name)
    return(res)
  }

  # Outgroup families are optional: a run may or may not include Talinaceae. Returning an empty
  # vector rather than stopping keeps the scaffold usable with either outgroup sampling.
  extract_optional_clade <- function(df, clade_name) {
    sort(unique(df$Specie_name[df$Clade == clade_name]))
  }
  
  collapse_clade <- function(x, genus_newick_lut = NULL) {
    if (length(x) == 0) return("")
    
    x <- sort(unique(x))
    genus_vec <- sub("_.*$", "", x)
    
    if (!is.null(genus_newick_lut) && length(genus_newick_lut) > 0) {
      grouped <- split(x, genus_vec)
      new_elements <- character(0)
      
      for (g in sort(names(grouped))) {
        if (g %in% names(genus_newick_lut)) {
          new_elements <- c(new_elements, genus_newick_lut[[g]])
        } else {
          new_elements <- c(new_elements, sort(grouped[[g]]))
        }
      }
      x <- sort(unique(new_elements))
    }
    
    if (length(x) == 1) return(x)
    return(paste0("(", paste(x, collapse = ","), ")"))
  }
  
  # 5. EXTRACT REQUIRED CLADES
  clades <- list(
    talinopsis    = extract_required_clade(cf, "Talinopsis"),
    grahamia      = extract_required_clade(cf, "Grahamia"),
    anacampseros  = extract_required_clade(cf, "Anacampseros"),
    portulaca     = extract_required_clade(cf, "Portulaca"),
    talinum       = extract_optional_clade(cf, "Talinum"),
    talinella     = extract_optional_clade(cf, "Talinella"),
    leuenbergeria = extract_required_clade(cf, "Leuenbergeria"),
    pereskia      = extract_required_clade(cf, "Pereskia"),
    tephrocacteae = extract_required_clade(cf, "Tephrocacteae"),
    cylindropuntieae = extract_required_clade(cf, "Cylindropuntieae"),
    opuntieae     = extract_required_clade(cf, "Opuntieae"),
    maihuenia     = extract_required_clade(cf, "Maihuenia"),
    blossfeldia   = extract_required_clade(cf, "Blossfeldia"),
    core_I        = extract_required_clade(cf, "Core I"),
    rhipsalideae  = extract_required_clade(cf, "Rhipsalideae"),
    notocacteae   = extract_required_clade(cf, "Notocacteae"),
    bct_core      = extract_required_clade(cf, "BCT"),
    calymmanthium = extract_required_clade(cf, "Calymmanthium"),
    copiapoa      = extract_required_clade(cf, "Copiapoa"),
    frailea       = extract_required_clade(cf, "Frailea"),
    cacteae       = extract_required_clade(cf, "Cacteae")
  )
  
  clade_sizes <- data.frame(
    clade = names(clades),
    n_taxa = vapply(clades, length, integer(1)),
    stringsAsFactors = FALSE
  )
  clade_sizes <- clade_sizes[order(clade_sizes$clade), ]
  utils::write.csv(clade_sizes, out_clade_sizes, row.names = FALSE)
  log_message("Clade size summary written to: ", out_clade_sizes)
  
  # 6. COLLAPSE CLADES TO NEWICK STRINGS
  s_tal   <- collapse_clade(clades$talinopsis, genus_newick_lut)
  s_gra   <- collapse_clade(clades$grahamia, genus_newick_lut)
  s_ana   <- collapse_clade(clades$anacampseros, genus_newick_lut)
  s_por   <- collapse_clade(clades$portulaca, genus_newick_lut)
  
  s_leu   <- collapse_clade(clades$leuenbergeria, genus_newick_lut)
  s_per   <- collapse_clade(clades$pereskia, genus_newick_lut)
  s_teph  <- collapse_clade(clades$tephrocacteae, genus_newick_lut)
  s_cyl   <- collapse_clade(clades$cylindropuntieae, genus_newick_lut)
  s_opu   <- collapse_clade(clades$opuntieae, genus_newick_lut)
  s_mai   <- collapse_clade(clades$maihuenia, genus_newick_lut)
  s_blo   <- collapse_clade(clades$blossfeldia, genus_newick_lut)
  s_coreI  <- collapse_clade(clades$core_I, genus_newick_lut)
  s_rhip <- collapse_clade(clades$rhipsalideae, genus_newick_lut)
  s_noto <- collapse_clade(clades$notocacteae, genus_newick_lut)
  s_bct <- collapse_clade(clades$bct_core, genus_newick_lut)
  s_caly  <- collapse_clade(clades$calymmanthium, genus_newick_lut)
  s_copi  <- collapse_clade(clades$copiapoa, genus_newick_lut)
  s_frai  <- collapse_clade(clades$frailea, genus_newick_lut)
  s_cact  <- collapse_clade(clades$cacteae, genus_newick_lut)
  
  s_talinum   <- collapse_clade(clades$talinum, genus_newick_lut)
  s_talinella <- collapse_clade(clades$talinella, genus_newick_lut)

  # 7. BUILD FIXED TOPOLOGY
  #
  # The outgroup families enter as a basal polytomy, one clade each, and their inter-family
  # relationships are left to the likelihood. With Portulacaceae and Anacampserotaceae alone this is
  # equivalent to the previous nested form: an unrooted constraint carries the same bipartitions
  # either way, so the scaffold is unchanged for runs without Talinaceae. It stops being equivalent
  # the moment a third outgroup family is present, because a nested form would then impose one of
  # the three possible resolutions of the Cactaceae-Anacampserotaceae-Portulacaceae-Talinaceae
  # quartet, which is the relationship the analysis is meant to test. Kew Tree of Life release 4.0
  # gives that node a quartet support of 0.47 with 0.33 on an alternative, so it is not settled and
  # must not be assumed here.
  s_anacampserotaceae <- paste0("(", s_tal, ",(", s_gra, ",", s_ana, "))")
  s_talinaceae <- if (nzchar(s_talinum) && nzchar(s_talinella)) {
    paste0("(", s_talinum, ",", s_talinella, ")")
  } else if (nzchar(s_talinum)) {
    s_talinum
  } else {
    s_talinella
  }
  outgroup_families <- c(s_por, s_anacampserotaceae, s_talinaceae)
  outgroup_families <- outgroup_families[nzchar(outgroup_families)]
  log_message("Outgroup families entering the scaffold as a basal polytomy: ",
              length(outgroup_families))
  s_opuntioideae <- paste0("(", s_cyl, ",(", s_teph, ",", s_opu, "))")
  s_coreII <- paste0("(", s_rhip, ",(", s_noto, ",", s_bct, "))")
  s_core_cactoideae <- paste0("(", s_frai, ",", s_caly, ",", s_copi, ",(", s_coreII, ",", s_coreI, "))")
  s_cactoideae <- paste0("(", s_cact, ",", s_core_cactoideae, ")")
  s_cactaceae <- paste0("(", s_leu, ",(", s_per, ",(", s_opuntioideae, ",(", s_mai, ",(", s_blo, ",", s_cactoideae, ")))))")
  
  constraint_tree <- paste0("(", paste(c(outgroup_families, s_cactaceae), collapse = ","), ");")
  constraint_tree <- gsub("\\s+", "", constraint_tree)
  log_message("Constraint tree string constructed.")
  
  # 8. VALIDATE TREE
  tree_obj <- tryCatch(
    ape::read.tree(text = constraint_tree),
    error = function(e) stop("Generated Newick is invalid: ", conditionMessage(e))
  )
  
  if (is.null(tree_obj$tip.label) || length(tree_obj$tip.label) == 0) {
    stop("Generated tree has no tip labels.")
  }
  
  tree_tips <- sort(unique(tree_obj$tip.label))
  if (!setequal(tree_tips, species_alignment)) {
    missing_in_tree <- setdiff(species_alignment, tree_tips)
    extra_in_tree <- setdiff(tree_tips, species_alignment)
    
    detail_msg <- paste0(
      "Final constraint tree tips do not match alignment taxa.",
      if (length(missing_in_tree) > 0) paste0(" Missing in tree: ", paste(missing_in_tree, collapse = ", "), ".") else "",
      if (length(extra_in_tree) > 0) paste0(" Extra in tree: ", paste(extra_in_tree, collapse = ", "), ".") else ""
    )
    stop(detail_msg)
  }
  
  if (anyDuplicated(tree_obj$tip.label) > 0) stop("Generated tree contains duplicated tip labels.")
  
  log_message("Tree parsed successfully.")
  log_message("Tree tips: ", length(tree_obj$tip.label))
  log_message("Tree is binary: ", ape::is.binary(tree_obj))
  
  # 9. EXPORT TREE
  writeLines(constraint_tree, con = out_tree)
  log_message("Constraint tree written to: ", out_tree)
  log_message("Constraint tree newick length: ", nchar(constraint_tree))
  log_message("Run finished successfully.")
  
  return(out_tree)
}

#' Infer Maximum-Likelihood Phylogeny under Constrained Search
#'
#' Infers the maximum-likelihood evolutionary hypothesis explaining the concatenated supermatrix under specified partition models
#' and topological constraint scaffolds using `RAxML-NG` (Kozlov *et al.*, 2019).
#' Executes multiple independent tree searches starting from randomized and parsimony starting trees to avoid local likelihood Optima.
#'
#' @param raxml_bin_path Character. System command or full path to executable `RAxML-NG` binary.
#' @param aln_file Character. Path to input PHYLIP supermatrix alignment file.
#' @param part_file Character. Path to partition file specifying substitution models per partition.
#' @param constraint_file Character. Path to Newick topological constraint scaffold file.
#' @param outgroup Character vector of terminals passed to `RAxML-NG --outgroup`, or `NULL` (default); multiple terminals are joined with commas. `RAxML-NG` writes an unrooted topology with these terminals placed first, so this argument orders the output rather than rooting the tree: the root is imposed downstream by `automate_treePL()` via `ape::root(..., resolve.root = TRUE)`. Declaring the same set at every stage keeps the output ordering consistent across the maximum-likelihood search, the bootstrap replicates and the temporal bootstraps. Derive it with `resolve_rooting_outgroup()` rather than naming a terminal by hand.
#' @param n_init_trees Character. Initial starting tree specifications. Defaults to \verb{"rand{25},pars{25}"} (25 random + 25 parsimony trees).
#' @param seed Integer. Random seed for reproducible tree search initialization. Defaults to `NULL` (random).
#' @param n_workers Integer. Parallel worker process count for RAxML-NG. If `NULL` (default), derived
#'   from `threads`, `n_init_trees` and `min_threads_per_worker` so that workers divide the starting
#'   tree count evenly and no worker idles in the final round. Set to `1` to reproduce a strictly
#'   sequential search.
#' @param min_threads_per_worker Integer. Lower bound on threads per worker when `n_workers` is derived. Defaults to `4L`.
#' @param threads Integer. Total number of CPU threads available to the search. Defaults to `4`.
#' @param output_dir Character. Directory path to save resulting maximum-likelihood tree files. Defaults to `dirname(aln_file)`.
#' @param prefix Character. Output file prefix. Defaults to `"cactus_search"`.
#' @return A named list containing paths to the best ML tree (`.raxml.bestTree`) and all evaluated trees (`.raxml.mlTrees`).
#' @references
#' Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and
#' user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics*, 35(21), 4453-4455.
#' \doi{10.1093/bioinformatics/btz305}
#' @examples
#' \dontrun{
#' calculate_ml_tree(
#'   raxml_bin_path = "raxml-ng",
#'   aln_file = "ALIGNMENT_supermatrix.phy",
#'   part_file = "MODELTEST_cactus_phylo.part.aicc",
#'   constraint_file = "cactus_constraints.tree"
#' )
#' }
#' @export
calculate_ml_tree <- function(raxml_bin_path, aln_file, part_file, constraint_file,
                              outgroup = NULL, n_init_trees = "rand{25},pars{25}",
                              seed = NULL, n_workers = NULL, min_threads_per_worker = 4L,
                              threads = 4,
                              output_dir = dirname(aln_file), prefix = "cactus_search") {

  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
    message("Using generated seed: ", seed)
  }

  # Wall-clock is set by how many sequential rounds each worker runs, not by the thread count
  # alone. Leaving n_workers at 1 puts every starting tree in its own round, which is why a
  # 50-tree search can take a full day on a workstation that could have run it in a quarter of
  # the rounds at the same total core-work.
  n_trees <- .parse_n_start_trees(n_init_trees)
  if (is.null(n_workers)) {
    plan <- .plan_ml_workers(n_trees, threads, min_threads_per_worker)
  } else {
    n_workers <- max(1L, as.integer(n_workers))
    plan <- list(workers = n_workers,
                 threads_per_worker = max(1L, as.integer(threads) %/% n_workers),
                 threads_used = n_workers * max(1L, as.integer(threads) %/% n_workers),
                 rounds = as.integer(ceiling(n_trees / n_workers)))
  }
  message("ML search plan: ", n_trees, " starting trees, ", plan$workers, " worker(s) x ~",
          plan$threads_per_worker, " threads, ", plan$rounds, " sequential round(s).")

  if (Sys.which(raxml_bin_path) == "") {
    stop("Executable '", raxml_bin_path, "' not found in your system's PATH.\n",
         "Please ensure RAxML-NG is installed and available, or provide the full absolute path.")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  out_prefix <- file.path(output_dir, prefix)
  
  args <- c(
    "--msa", shQuote(aln_file),
    "--model", shQuote(part_file),
    "--tree-constraint", shQuote(constraint_file),
    "--tree", shQuote(n_init_trees),
    "--seed", as.character(seed),
    "--workers", as.character(plan$workers),
    "--threads", as.character(threads),
    "--prefix", shQuote(out_prefix)
  )
  
  og <- .format_outgroup(outgroup)
  if (!is.null(og)) {
    args <- c(args, "--outgroup", shQuote(og))
  }
  
  ml_status <- system2(command = raxml_bin_path, args = args)
  .check_cli_exit(ml_status, "RAxML-NG ML")
  
  message("\nMaximum-likelihood tree calculation completed. \U0001f335")
  return(list(
    bestTree = paste0(out_prefix, ".raxml.bestTree"),
    mlTrees  = paste0(out_prefix, ".raxml.mlTrees")
  ))
}

#' Count Starting Trees Declared in a RAxML-NG `--tree` Specification
#'
#' Parses specifications such as \verb{"rand{25},pars{25}"} into the total number of independent
#' tree searches RAxML-NG will perform. A term without an explicit count contributes one tree.
#' @param n_init_trees Character or numeric. Starting tree specification.
#' @return Integer count of independent starting trees.
#' @noRd
.parse_n_start_trees <- function(n_init_trees) {
  if (is.numeric(n_init_trees)) return(max(1L, as.integer(n_init_trees)))
  terms <- trimws(strsplit(as.character(n_init_trees), ",")[[1]])
  terms <- terms[nzchar(terms)]
  if (length(terms) == 0) return(1L)
  counts <- vapply(terms, function(tm) {
    m <- regmatches(tm, regexpr("[0-9]+", tm))
    if (length(m) == 1L && nzchar(m)) as.integer(m) else 1L
  }, integer(1))
  as.integer(sum(counts))
}

#' Plan the Worker and Thread Split for a Parallel ML Tree Search
#'
#' RAxML-NG distributes independent starting-tree searches across `--workers`, each worker using
#' `threads / workers` threads on one tree at a time. Total core-work is fixed, so wall-clock is
#' governed by the number of sequential rounds, `ceiling(n_trees / workers)`. A worker count that
#' does not divide `n_trees` leaves workers idle in the final round: 16 workers over 50 trees runs
#' four rounds with only two workers busy in the last one, the same wall-clock as 25 workers over
#' two full rounds but using more cores to get there. This planner therefore restricts the worker
#' count to divisors of `n_trees`.
#'
#' The lower bound on threads per worker exists because per-worker scaling saturates once each
#' thread holds too few site patterns. As a working floor, keep at least a few hundred patterns per
#' thread and re-validate on new hardware rather than trusting the default.
#'
#' @param n_trees Integer. Total independent starting trees.
#' @param threads Integer. Total threads available to the job.
#' @param min_threads_per_worker Integer. Lower bound on threads assigned to each worker.
#' @return A list with `workers`, `threads_per_worker`, `threads_used` and `rounds`.
#' @noRd
.plan_ml_workers <- function(n_trees, threads, min_threads_per_worker = 4L) {
  n_trees <- max(1L, as.integer(n_trees))
  threads <- max(1L, as.integer(threads))
  min_tpw <- max(1L, as.integer(min_threads_per_worker))

  max_workers <- max(1L, threads %/% min_tpw)
  divisors <- seq_len(n_trees)[n_trees %% seq_len(n_trees) == 0L]
  feasible <- divisors[divisors <= max_workers]
  workers <- if (length(feasible) > 0L) max(feasible) else 1L

  tpw <- max(1L, threads %/% workers)
  list(
    workers = as.integer(workers),
    threads_per_worker = as.integer(tpw),
    threads_used = as.integer(workers * tpw),
    rounds = as.integer(ceiling(n_trees / workers))
  )
}

#' Generate a SLURM Batch Script for the Constrained Maximum-Likelihood Search
#'
#' Writes a single SLURM batch script running the same constrained `RAxML-NG` search as
#' [calculate_ml_tree()], sized for a compute node instead of a workstation. The search is one job
#' rather than a job array: independent starting-tree searches are distributed across `--workers`
#' inside the job, and RAxML-NG writes one `.raxml.bestTree` directly, so no collection step is
#' needed to compare log-likelihoods across tasks.
#'
#' @details
#' Wall-clock is governed by the number of sequential rounds each worker must run,
#' `ceiling(n_trees / workers)`, not by the thread count alone. Total core-work is fixed, so adding
#' threads to a single worker shortens each tree while adding workers shortens the number of rounds,
#' and the second lever is the one that matters once a search is already using an efficient thread
#' count per worker. Running 50 starting trees on one worker executes 50 rounds; the same 50 trees
#' across 25 workers execute 2.
#'
#' Reference timing illustrating worker parallelization: on an Apple M2 Pro (8 threads, `--workers 1`,
#' `RAxML-NG` 1.2.2, SSE3 kernels), searching 50 starting trees sequentially required 47115 s (~938 s per tree).
#' In contrast, the production run on an HPC cluster node (AMD EPYC 9754, 75 threads, `--workers 25`) over the
#' full supermatrix (1023 terminals, 12806 sites, 5954 patterns, 11 partitions) completed 50 starting trees in
#' 2393 s (~40 minutes), executing two parallel rounds.
#'
#' `workers` is derived automatically and constrained to a divisor of the starting tree
#' count, so that no worker sits idle in the final round.
#'
#' @param alignment_file Character. Path to input PHYLIP supermatrix alignment file.
#' @param partition_file Character. Path to partition file specifying substitution models per partition.
#' @param constraint_file Character. Path to Newick topological constraint scaffold file.
#' @param outgroup Character vector of terminals passed to `RAxML-NG --outgroup`, or `NULL` (default); multiple terminals are joined with commas. `RAxML-NG` writes an unrooted topology with these terminals placed first, so this argument orders the output rather than rooting the tree: the root is imposed downstream by `automate_treePL()` via `ape::root(..., resolve.root = TRUE)`. Declaring the same set at every stage keeps the output ordering consistent across the maximum-likelihood search, the bootstrap replicates and the temporal bootstraps. Derive it with `resolve_rooting_outgroup()` rather than naming a terminal by hand.
#' @param n_init_trees Character. Starting tree specification passed to `--tree`. Defaults to \verb{"rand{25},pars{25}"}.
#' @param seed Integer. Random seed for reproducible tree search initialization. Defaults to `NULL` (random).
#' @param threads Integer. CPU cores requested per SLURM task (`--cpus-per-task`). Defaults to `75`.
#' @param workers Integer. Worker process count. If `NULL` (default), derived from `threads`,
#'   `n_init_trees` and `min_threads_per_worker` so that workers divide the starting tree count evenly.
#' @param min_threads_per_worker Integer. Lower bound on threads per worker when `workers` is derived. Defaults to `3L`.
#' @param preparse Logical. Pre-parse the alignment into compressed binary `.rba` format before the search? Defaults to `TRUE`.
#' @param output_dir Character. Destination directory for the script and its outputs. Defaults to `getwd()`.
#' @param script_name Character. Name of the generated Bash script. Defaults to `"run_ml_search.sh"`.
#' @param prefix Character. RAxML-NG output prefix inside the job. Defaults to `"cactus_search"`.
#' @param cluster_job_name Character. SLURM job name. Defaults to `"cactus_ml"`.
#' @param cluster_partition Character. SLURM partition. Defaults to `"main"`.
#' @param cluster_nodes Integer. Number of compute nodes requested (\verb{--nodes}). Defaults to \code{1L}.
#' @param cluster_mem Character. Memory allocation (`--mem`). Defaults to `"16G"`.
#' @param cluster_time Character. Time limit (`--time`). Defaults to `"02:00:00"`.
#' @param cluster_queue Character. Optional SLURM queue / QoS. Defaults to `NULL`.
#' @param cluster_mail_user Character. Notification recipient email address for SLURM (\verb{--mail-user}). Defaults to \code{Sys.getenv("MY_EMAIL", "")}.
#' @param load_module Character vector. Environment modules to load. Defaults to the NLHPC Leftraru
#'   `RAxML-NG` toolchain.
#' @param raxml_exec Character. `RAxML-NG` executable invoked inside the job. Defaults to `"raxml-ng-mpi"`.
#' @return Character path to the generated SLURM batch script, carrying the resolved parallel plan
#'   as the `ml_plan` attribute.
#' @seealso [calculate_ml_tree()] for the equivalent local run.
#' @references
#' Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and
#' user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics*, 35(21), 4453-4455.
#' \doi{10.1093/bioinformatics/btz305}
#' @examples
#' \dontrun{
#' generate_ml_search_script(
#'   alignment_file = "cactus.raxml.reduced.phy",
#'   partition_file = "cactus_modeltest.part.aicc",
#'   constraint_file = "cactus_constraints.tree",
#'   outgroup = resolve_rooting_outgroup(ape::read.tree("bestTree.tree")$tip.label),
#'   output_dir = "7_Phylogenetics"
#' )
#' }
#' @export
generate_ml_search_script <- function(alignment_file, partition_file, constraint_file,
                                      outgroup = NULL, n_init_trees = "rand{25},pars{25}",
                                      seed = NULL,
                                      threads = 75, workers = NULL, min_threads_per_worker = 3L,
                                      preparse = TRUE,
                                      output_dir = getwd(),
                                      script_name = "run_ml_search.sh",
                                      prefix = "cactus_search",
                                      cluster_job_name = "cactus_ml",
                                      cluster_partition = "main",
                                      cluster_nodes = 1L,
                                      cluster_mem = "16G", cluster_time = "02:00:00",
                                      cluster_queue = NULL,
                                      cluster_mail_user = Sys.getenv("MY_EMAIL", ""),
                                      load_module = c("gcc/14.2.0-nlhpc", "openmpi/5.0.3-o", "raxml-ng/1.1.0-mpi-zen4-n"),
                                      raxml_exec = "raxml-ng-mpi") {

  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
    message("Using generated seed: ", seed)
  }

  n_trees <- .parse_n_start_trees(n_init_trees)
  threads <- as.integer(threads)

  if (is.null(workers)) {
    plan <- .plan_ml_workers(n_trees, threads, min_threads_per_worker)
  } else {
    workers_val <- max(1L, as.integer(workers))
    plan <- list(
      workers = workers_val,
      threads_per_worker = max(1L, threads %/% workers_val),
      threads_used = workers_val * max(1L, threads %/% workers_val),
      rounds = as.integer(ceiling(n_trees / workers_val))
    )
  }

  ml_dir <- file.path(output_dir, "ml_search")
  if (!dir.exists(ml_dir)) dir.create(ml_dir, recursive = TRUE, showWarnings = FALSE)

  bash_script <- file.path(ml_dir, script_name)
  out_prefix_base <- prefix
  rba_prefix_base <- paste0(prefix, "_parsed")

  aln_basename  <- basename(alignment_file)
  part_basename <- basename(partition_file)
  cons_basename <- basename(constraint_file)

  L <- function(...) cat(paste0(..., "\n"), file = bash_script, append = TRUE)

  cat("#!/bin/bash\n", file = bash_script)
  L("#SBATCH -J ", cluster_job_name)
  L("#SBATCH -p ", cluster_partition)
  L("#SBATCH --nodes=", cluster_nodes)
  L("#SBATCH --ntasks-per-node=1")
  L("#SBATCH --cpus-per-task=", threads)
  L("#SBATCH --mem=", cluster_mem)
  L("#SBATCH -t ", cluster_time)
  L("#SBATCH -o slurm_%j.out")
  L("#SBATCH -e slurm_%j.err")
  if (!is.null(cluster_queue) && nzchar(cluster_queue)) L("#SBATCH -q ", cluster_queue)
  L("#SBATCH --mail-type=ALL")
  if (!is.null(cluster_mail_user) && nzchar(trimws(cluster_mail_user))) {
    L("#SBATCH --mail-user=", trimws(cluster_mail_user))
  }
  L("#")
  L("")

  if (!is.null(load_module) && length(load_module) > 0) {
    modules_str <- paste(load_module, collapse = " ")
    if (nzchar(trimws(modules_str))) {
      L("module load ", modules_str)
      L("")
    }
  }

  L("echo \"=== RAxML-NG version used by this job ===\"")
  L(raxml_exec, " --version | head -n 2")
  L("echo \"=========================================\"")
  L("")
  L("# Move to submit directory and dynamically resolve paths")
  L("cd \"${SLURM_SUBMIT_DIR:-$PWD}\"")
  L("SCRIPT_DIR=\"$PWD\"")
  L("RUN_DIR=\"$(cd \"${SCRIPT_DIR}/..\" && pwd)\"")
  L("")
  L("ALIGNMENT_FILE=\"${RUN_DIR}/", aln_basename, "\"")
  L("PARTITION_FILE=\"${RUN_DIR}/", part_basename, "\"")
  L("CONSTRAINT_FILE=\"${RUN_DIR}/", cons_basename, "\"")
  L("OUT_PREFIX=\"${SCRIPT_DIR}/", out_prefix_base, "\"")
  L("RBA_PREFIX=\"${SCRIPT_DIR}/", rba_prefix_base, "\"")
  L("RBA_FILE=\"${RBA_PREFIX}.raxml.rba\"")
  L("")

  L("TOTAL_THREADS=${SLURM_CPUS_PER_TASK:-", threads, "}")
  L("NUM_WORKERS=", plan$workers)
  L("SEED=", seed)
  L("")
  L("export OMP_PROC_BIND=false")
  L("export OMPI_MCA_hwloc_base_binding_policy=none")
  L("")
  L("echo \"ML search: ", n_trees, " starting trees over ${NUM_WORKERS} workers",
    " (~", plan$threads_per_worker, " threads each), ", plan$rounds, " sequential round(s)\"")
  L("")

  if (preparse) {
    L("if [ ! -f \"$RBA_FILE\" ]; then")
    L("    echo \"Pre-parsing alignment to compressed binary RBA format...\"")
    L("    ", raxml_exec, " --parse --msa \"$ALIGNMENT_FILE\"",
      " --model \"$PARTITION_FILE\" --prefix \"$RBA_PREFIX\"")
    L("fi")
    L("")
  }

  cmd <- paste0(raxml_exec, " ")
  if (preparse) {
    cmd <- paste0(cmd, "--msa \"$RBA_FILE\" ")
  } else {
    cmd <- paste0(cmd, "--msa \"$ALIGNMENT_FILE\" --model \"$PARTITION_FILE\" ")
  }
  cmd <- paste0(cmd, "--tree-constraint \"$CONSTRAINT_FILE\" ")
  cmd <- paste0(cmd, "--tree ", shQuote(n_init_trees, type = "sh"), " ")
  og <- .format_outgroup(outgroup)
  if (!is.null(og)) {
    cmd <- paste0(cmd, "--outgroup ", shQuote(og, type = "sh"), " ")
  }
  cmd <- paste0(
    cmd,
    "--seed ${SEED} ",
    "--threads ${TOTAL_THREADS} ",
    "--workers ${NUM_WORKERS} ",
    "--force perf_threads ",
    "--extra thread-nopin ",
    "--prefix \"$OUT_PREFIX\""
  )
  L(cmd)

  if (.Platform$OS.type == "unix") {
    chmod_status <- system(paste("chmod +x", shQuote(bash_script)))
    .check_cli_exit(chmod_status, "chmod")
  } else {
    Sys.chmod(bash_script, mode = "0755")
  }

  message("SLURM ML search script generated: ", bash_script)
  message("  Starting trees : ", n_trees)
  message("  Workers        : ", plan$workers, " x ~", plan$threads_per_worker, " threads",
          " (", plan$threads_used, " of ", threads, " cores in use)")
  message("  Sequential rounds: ", plan$rounds)

  attr(bash_script, "ml_plan") <- plan
  return(bash_script)
}

#' Map Bootstrap Support Values onto Reference Phylogeny
#'
#' Maps clade support derived from non-parametric bootstrap replicates onto the best maximum-likelihood
#' tree topology. The replicates themselves carry no metric: `RAxML-NG --bootstrap` writes plain
#' topologies with branch lengths, and the metric is chosen here, at the summarising step. The same
#' replicate file can therefore be summarised under both metrics without recomputation.
#'
#' Defaults to Felsenstein's Bootstrap Percentage (FBP; Felsenstein, 1985), which is the metric the
#' Cactaceae and Caryophyllales dating literature reports and the only one against which this tree can
#' be compared. Transfer Bootstrap Expectation (TBE; Lemoine *et al.*, 2018) is available through
#' `metric = "tbe"` and belongs in a clearly labelled secondary column.
#'
#' The two are not on a common scale and TBE must never be reported as though it were a bootstrap
#' percentage. TBE is bounded below by FBP and its inflation grows with clade size. Measured on the
#' 986 unconstrained nodes of the 1023 terminal supermatrix tree (2026-09-06): median TBE 0.783
#' against median FBP 0.468, with the gap reaching 0.704 for clades of 51 to 200 terminals. Reporting
#' TBE would place 61 percent of nodes above 0.70; FBP places 29 percent.
#'
#' Support values are meaningless for any bipartition imposed through `--tree-constraint`, because
#' every replicate reproduces it by construction. Those nodes must be reported as constrained rather
#' than supported. Use `classify_constrained_nodes()` to separate them before tabulating support.
#'
#' @param raxml_bin Character. System command or full path to executable `RAxML-NG` binary.
#' @param best_tree Character. Path to reference maximum-likelihood tree file.
#' @param bootstraps_file Character. Path to concatenated non-parametric bootstrap trees file.
#' @param metric Character. Bootstrap support metric: `"fbp"` (Felsenstein's Bootstrap Percentage, the default and the comparable metric) or `"tbe"` (Transfer Bootstrap Expectation, secondary).
#' @param threads Integer. Number of CPU threads. Defaults to `4`.
#' @param output_dir Character. Output directory for annotated support tree. Defaults to `dirname(best_tree)`.
#' @param prefix Character. Output file prefix. Defaults to `"cactus_support"`.
#' @return Character path to the annotated support tree file (`.raxml.support`).
#' @references
#' Felsenstein, J. (1985). Confidence limits on phylogenies: an approach using the bootstrap.
#' *Evolution*, 39(4), 783-791. \doi{10.1111/j.1558-5646.1985.tb00420.x}
#'
#' Lemoine, F., Domelevo Entfellner, J. B., Wilkinson, E., Correia, D., Davila Felipe, M.,
#' De Oliveira, T., & Gascuel, O. (2018). Renewing Felsenstein's phylogenetic bootstrap in the era of
#' big data. *Nature*, 556(7702), 452-456. \doi{10.1038/s41586-018-0043-0}
#' @examples
#' \dontrun{
#' # Primary metric, comparable with the published literature.
#' map_branch_supports(
#'   raxml_bin = "raxml-ng",
#'   best_tree = "cactus_search.raxml.bestTree",
#'   bootstraps_file = "cactus_ALL_bootstraps.tree",
#'   metric = "fbp",
#'   prefix = "cactus_support_fbp"
#' )
#'
#' # Secondary metric, same replicates, no recomputation.
#' map_branch_supports(
#'   raxml_bin = "raxml-ng",
#'   best_tree = "cactus_search.raxml.bestTree",
#'   bootstraps_file = "cactus_ALL_bootstraps.tree",
#'   metric = "tbe",
#'   prefix = "cactus_support_tbe"
#' )
#' }
#' @export
map_branch_supports <- function(raxml_bin, best_tree, bootstraps_file, metric = c("fbp", "tbe"), threads = 4, output_dir = dirname(best_tree), prefix = "cactus_support") {
  metric <- match.arg(metric)
  
  if (Sys.which(raxml_bin) == "") {
    stop("Executable '", raxml_bin, "' not found in your system's PATH.\n",
         "Please ensure RAxML-NG is installed and available, or provide the full absolute path.")
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_prefix <- file.path(output_dir, prefix)
  
  args <- c(
    "--support",
    "--tree", shQuote(best_tree),
    "--bs-trees", shQuote(bootstraps_file),
    "--bs-metric", shQuote(metric),
    "--threads", as.character(threads),
    "--prefix", shQuote(out_prefix)
  )
  
  ml_status <- system2(command = raxml_bin, args = args)
  .check_cli_exit(ml_status, "RAxML-NG ML")
  support_file <- paste0(out_prefix, ".raxml.support")
  message("Branch support values (", toupper(metric), ") mapped onto best tree: ", support_file)
  return(support_file)
}

#' Compute Robinson-Foulds Distances Across Maximum-Likelihood Trees
#'
#' Computes pairwise Robinson-Foulds (RF) topological distances across tree topologies generated during maximum-likelihood search in `RAxML-NG`.
#' Quantifying topological variance evaluates whether independent search runs converged on identical tree topologies.
#'
#' @param raxml_bin_path Character. System command or full path to executable `RAxML-NG` binary.
#' @param ml_trees_file Character. Path to input `.raxml.mlTrees` file containing multiple ML tree search replicates.
#' @param output_dir Character. Output directory for RF distance calculations. Defaults to `dirname(ml_trees_file)`.
#' @param prefix Character. Output file prefix. Defaults to `"cactus_RF"`.
#' @return Character path to the resulting RF distance output file (`.raxml.rfdist`).
#' @export
calculate_rf_distances <- function(raxml_bin_path, ml_trees_file, output_dir = dirname(ml_trees_file), prefix = "cactus_RF") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_prefix <- file.path(output_dir, prefix)
  
  args <- c(
    "--rfdist",
    "--tree", shQuote(ml_trees_file),
    "--prefix", shQuote(out_prefix)
  )
  
  status <- system2(command = raxml_bin_path, args = args)
  .check_cli_exit(status, "RAxML-NG RF")
  rf_file <- paste0(out_prefix, ".raxml.rfdist")
  message("Pairwise Robinson-Foulds distance matrix saved to: ", rf_file)
  return(rf_file)
}

#' Detect Performance-Core Count for Local Thread Allocation
#'
#' Attempts to detect the number of Performance cores (P-cores) on Apple Silicon (macOS) via
#' `sysctl -n hw.perflevel0.logicalcpu`, so that local `RAxML-NG` runs are not scheduled onto
#' Efficiency cores (E-cores), which would otherwise create thread contention and reduce
#' throughput. This `sysctl` key only exists on Apple Silicon Macs with heterogeneous P/E core
#' topology (M1/M2/M3/M4 series); on Intel Macs, Linux, or Windows it fails silently and this
#' function falls back to `parallel::detectCores(logical = FALSE)`, and ultimately to a hardcoded
#' `8L` if physical core count cannot be determined either.
#' @noRd
.detect_pcores <- function() {
  if (identical(tolower(Sys.info()[["sysname"]]), "darwin")) {
    pcores_cmd <- suppressWarnings(system("sysctl -n hw.perflevel0.logicalcpu", intern = TRUE))
    if (length(pcores_cmd) > 0 && !is.na(as.integer(pcores_cmd[1])) && as.integer(pcores_cmd[1]) > 0) {
      return(as.integer(pcores_cmd[1]))
    }
  }
  phys <- parallel::detectCores(logical = FALSE)
  if (!is.na(phys) && phys > 0) {
    return(as.integer(phys))
  }
  return(8L)
}

#' Generate Non-Parametric Bootstrap Trees Locally
#'
#' Performs non-parametric bootstrap resampling over supermatrix site columns to infer a distribution of bootstrap tree topologies (`RAxML-NG`).
#' Evaluates topological variation under non-parametric resampling to quantify node support via Felsenstein Bootstrap Proportions (FBP) or Transfer Bootstrap Expectation (TBE).
#' Optimizes multi-threading for local workstations (e.g., Apple Silicon) by allocating threads to Performance cores (P-cores)
#' and configuring parallel workers according to a configurable thread-to-worker ratio, avoiding thread contention with
#' Efficiency cores (E-cores).
#'
#' Unlike `generate_bootstrap_script()`'s HPC-oriented default (empirically calibrated on a 128-core AMD EPYC node),
#' the `threads_per_worker = 4L` default here has not been empirically validated for bootstrap workloads on Apple
#' Silicon specifically (only single-threaded, non-coarse-grained ML tree search timings are currently available for
#' that hardware). Treat this default as a reasonable starting heuristic and re-validate empirically (e.g., time a
#' short run with `bs_trees` set low, at a couple of `threads_per_worker` values) before committing to a long local run.
#'
#' @param raxml_bin_path Character. System command or full path to executable `RAxML-NG` binary.
#' @param aln_file Character. Path to input PHYLIP supermatrix alignment file or pre-parsed `.rba` binary file.
#' @param part_file Character. Path to partition file specifying substitution models. Optional if `aln_file` is an `.rba` binary file. Defaults to `NULL`.
#' @param constraint_file Character. Path to Newick topological constraint scaffold file.
#' @param bs_trees Integer. Total number of non-parametric bootstrap trees to generate. Defaults to `500`.
#' @param outgroup Character vector of terminals passed to `RAxML-NG --outgroup`, or `NULL` (default); multiple terminals are joined with commas. `RAxML-NG` writes an unrooted topology with these terminals placed first, so this argument orders the output rather than rooting the tree: the root is imposed downstream by `automate_treePL()` via `ape::root(..., resolve.root = TRUE)`. Declaring the same set at every stage keeps the output ordering consistent across the maximum-likelihood search, the bootstrap replicates and the temporal bootstraps. Derive it with `resolve_rooting_outgroup()` rather than naming a terminal by hand.
#' @param seed Integer. Random seed for reproducible bootstrap initialization. Defaults to `NULL` (random).
#' @param threads Integer. Number of CPU threads. Defaults to `8` (or detected P-cores).
#' @param workers Integer. Parallel worker process count. If `NULL` (default), automatically calculated as `max(1L, as.integer(threads / threads_per_worker))`.
#' @param threads_per_worker Integer. Thread-to-worker ratio used to derive `workers` when `workers = NULL`. Defaults to `4L`
#'   (unvalidated heuristic for local bootstrap workloads; see Details). Ignored if `workers` is set explicitly.
#' @param output_dir Character. Output directory for generated bootstrap trees. Defaults to `dirname(aln_file)`.
#' @param prefix Character. Output file prefix. Defaults to `"cactus_bs"`.
#' @return Character path to the output bootstrap trees file (`.raxml.bootstraps`).
#' @references
#' Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and
#' user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics*, 35(21), 4453-4455.
#' \doi{10.1093/bioinformatics/btz305}
#' @export
run_local_bootstraps <- function(raxml_bin_path, aln_file, part_file = NULL, constraint_file,
                                 bs_trees = 500, outgroup = NULL, seed = NULL,
                                 threads = 8, workers = NULL, threads_per_worker = 4L,
                                 output_dir = dirname(aln_file),
                                 prefix = "cactus_bs") {
                                 
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
    message("Using generated seed: ", seed)
  }

  if (Sys.which(raxml_bin_path) == "") {
    stop("Executable '", raxml_bin_path, "' not found in your system's PATH.\n",
         "Please ensure RAxML-NG is installed and available, or provide the full absolute path.")
  }

  if (!file.exists(aln_file)) stop("Alignment file not found: ", aln_file)
  if (!file.exists(constraint_file)) stop("Constraint tree file not found: ", constraint_file)
  
  is_rba <- grepl("\\.rba$", aln_file, ignore.case = TRUE)
  if (!is_rba && (is.null(part_file) || !file.exists(part_file))) {
    stop("A valid part_file must be provided when aln_file is not a pre-parsed .rba binary file.")
  }

  if (is.null(threads) || is.na(threads)) {
    threads <- .detect_pcores()
  }
  threads <- as.integer(threads)

  if (is.null(workers) || is.na(workers)) {
    workers <- max(1L, as.integer(threads / as.integer(threads_per_worker)))
  }
  workers <- as.integer(workers)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_prefix <- file.path(output_dir, prefix)
  
  args <- c(
    "--bootstrap",
    "--msa", shQuote(aln_file)
  )
  
  if (!is_rba && !is.null(part_file)) {
    args <- c(args, "--model", shQuote(part_file))
  }
  
  # No `--bs-metric` here. In `--bootstrap` mode RAxML-NG resamples alignment columns and infers one
  # maximum-likelihood tree per replicate; the only output is `.raxml.bootstraps`, plain topologies
  # with branch lengths and no node labels. The metric is consumed later, by `--support`, which
  # summarises those replicates. Passing it here was accepted and silently dropped: the run of
  # 2026-09-04 carried `--bs-metric tbe` and its log never lists the option among the analysis
  # settings, reports `run mode: Bootstrapping` rather than naming a metric, and wrote no support
  # file. Removed because it reads as though the replicates were committed to one metric, which
  # cost a round of confusion on 2026-09-06. The replicates support any metric; choose it at
  # `map_branch_supports()`.
  args <- c(
    args,
    "--tree-constraint", shQuote(constraint_file),
    "--bs-trees", as.character(bs_trees),
    "--seed", as.character(seed),
    "--threads", as.character(threads),
    "--workers", as.character(workers),
    "--prefix", shQuote(out_prefix)
  )
  
  og <- .format_outgroup(outgroup)
  if (!is.null(og)) {
    args <- c(args, "--outgroup", shQuote(og))
  }
  
  status <- system2(command = raxml_bin_path, args = args)
  .check_cli_exit(status, "RAxML-NG Bootstraps")
  return(paste0(out_prefix, ".raxml.bootstraps"))
}

#' Generate HPC SLURM Batch Script for Parallel Bootstrapping
#'
#' Generates an executable Bash script with SLURM scheduler directives to parallelize non-parametric bootstrapping
#' across HPC compute nodes (optimized for NLHPC Leftraru Epu and general Slurm clusters).
#' Uses coarse-grained parallelization via Slurm Job Arrays (`#SBATCH --array=1-N`), pre-parsing to compressed binary `.rba` format
#' to optimize disk I/O, dynamic per-task seed multiplication for statistical independence, and a configurable thread-to-worker ratio.
#'
#' The default `threads_per_worker = 4L` is configured to optimize throughput on multi-core compute nodes. In a production run of this workload
#' (1023 taxa, 11 partitions), using `--threads 40 --workers 10` (4 threads per worker) on an AMD EPYC node provided efficient per-worker
#' memory and CPU allocation. This ratio is dataset- and hardware-dependent (it trades per-worker single-tree search speed against the number of trees searched in
#' parallel); if you migrate to different node hardware or a markedly different supermatrix size, re-validate it
#' empirically with a short trial run before committing a full job array to it.
#'
#' @param alignment_file Character. Path to input PHYLIP alignment file.
#' @param partition_file Character. Path to partition file.
#' @param constraint_file Character. Path to constraint scaffold tree file.
#' @param outgroup Character vector of terminals passed to `RAxML-NG --outgroup`, or `NULL` (default); multiple terminals are joined with commas. `RAxML-NG` writes an unrooted topology with these terminals placed first, so this argument orders the output rather than rooting the tree: the root is imposed downstream by `automate_treePL()` via `ape::root(..., resolve.root = TRUE)`. Declaring the same set at every stage keeps the output ordering consistent across the maximum-likelihood search, the bootstrap replicates and the temporal bootstraps. Derive it with `resolve_rooting_outgroup()` rather than naming a terminal by hand.
#' @param bs_per_rep Integer. Number of bootstrap trees generated per chunk replicate. Defaults to `500`.
#' @param max_reps Integer. Total number of parallel chunk replicates (tasks) to spawn in the SLURM array (`1-max_reps`). Defaults to `2`,
#'   giving a total array target of `bs_per_rep * max_reps = 1000` bootstrap trees. This provides a safety margin over the autoMRE
#'   convergence point observed in a production run of this exact dataset (1023 taxa, 11 partitions), which converged (`bs-cutoff = 0.03`,
#'   FBP) after 600 of 1000 collected trees, i.e. convergence is not guaranteed at a fixed replicate count for every dataset or taxon
#'   sampling scheme. Always confirm convergence with `check_bs_convergence()` on the collected trees (via `collect_bootstraps()`) rather
#'   than assuming `bs_per_rep * max_reps` is sufficient; increase `max_reps` and re-run `collect_bootstraps()` if it is not.
#' @param base_seed Integer. Base random seed for dynamic seed calculation (`SEED=$(( SLURM_ARRAY_TASK_ID * base_seed ))`). Defaults to `NULL` (random).
#' @param threads Integer. Number of CPU cores requested per SLURM task (`--cpus-per-task`). Defaults to `40`.
#' @param workers Integer. Number of RAxML-NG worker processes. If `NULL` (default), calculated dynamically as `max(1L, as.integer(threads / threads_per_worker))`.
#' @param threads_per_worker Integer. Thread-to-worker ratio used to derive `workers` when `workers = NULL`. Defaults to `4L`
#'   (empirically validated for this workload on an AMD EPYC node; see Details). Ignored if `workers` is set explicitly.
#' @param preparse Logical. Generate compressed binary `.rba` format with `raxml-ng --parse` prior to array execution to minimize disk I/O contention? Defaults to `TRUE`.
#' @param output_dir Character. Destination directory for script and chunk logs. Defaults to `getwd()`.
#' @param script_name Character. Name of output Bash script file. Defaults to `"run_bs_chunks.sh"`.
#' @param cluster_job_name Character. SLURM job name identifier. Defaults to `"cactus_bs"`.
#' @param cluster_partition Character. SLURM partition name. Defaults to `"main"` (Leftraru Epu AMD EPYC 9754 partition).
#' @param cluster_nodes Integer. Number of compute nodes requested (\verb{--nodes}). Defaults to \code{1L}.
#' @param cluster_mem Character. Memory allocation string for SLURM (`--mem`). Defaults to `"10G"`.
#' @param cluster_time Character. Time limit allocation string for SLURM (`--time`). Defaults to `"24:00:00"`.
#' @param cluster_queue Character. Optional SLURM queue / QoS name. Defaults to `NULL`.
#' @param cluster_mail_user Character. Notification recipient email address for SLURM (\verb{--mail-user}). Defaults to \code{Sys.getenv("MY_EMAIL", "")}.
#' @param load_module Character vector or string. Environment module(s) to load prior to execution. Defaults to `c("gcc/14.2.0-nlhpc", "openmpi/5.0.3-o", "raxml-ng/1.1.0-mpi-zen4-n")`.
#' @param raxml_exec Character. Executable `RAxML-NG` binary command. Defaults to `"raxml-ng-mpi"`.
#' @return Character path to the generated SLURM batch script file.
#' @export
generate_bootstrap_script <- function(alignment_file, partition_file, constraint_file,
                                      outgroup = NULL, bs_per_rep = 500, max_reps = 2,
                                      base_seed = NULL,
                                      threads = 40, workers = NULL, threads_per_worker = 4L,
                                      preparse = TRUE,
                                      output_dir = getwd(),
                                      script_name = "run_bs_chunks.sh",
                                      cluster_job_name = "cactus_bs", 
                                      cluster_partition = "main",
                                      cluster_nodes = 1L,
                                      cluster_mem = "10G", cluster_time = "24:00:00",
                                      cluster_queue = NULL,
                                      cluster_mail_user = Sys.getenv("MY_EMAIL", ""),
                                      load_module = c("gcc/14.2.0-nlhpc", "openmpi/5.0.3-o", "raxml-ng/1.1.0-mpi-zen4-n"),
                                      raxml_exec = "raxml-ng-mpi") {
  
  if (is.null(base_seed)) {
    base_seed <- sample.int(100000L, 1L)
    message("Using generated base seed: ", base_seed)
  }

  threads <- as.integer(threads)
  threads_per_worker <- as.integer(threads_per_worker)
  if (is.null(workers)) {
    workers_val <- max(1L, as.integer(threads / threads_per_worker))
  } else {
    workers_val <- as.integer(workers)
  }

  bs_dir <- file.path(output_dir, "bs_chunks")
  if (!dir.exists(bs_dir)) dir.create(bs_dir, recursive = TRUE)
  
  bash_script <- file.path(bs_dir, script_name)
  aln_basename  <- basename(alignment_file)
  part_basename <- basename(partition_file)
  cons_basename <- basename(constraint_file)

  cat("#!/bin/bash\n", file = bash_script)
  cat(paste0("#SBATCH -J ", cluster_job_name, "\n"), file = bash_script, append = TRUE)
  cat(paste0("#SBATCH -p ", cluster_partition, "\n"), file = bash_script, append = TRUE)
  cat(paste0("#SBATCH --nodes=", cluster_nodes, "\n"), file = bash_script, append = TRUE)
  cat("#SBATCH --ntasks-per-node=1\n", file = bash_script, append = TRUE)
  cat(paste0("#SBATCH --cpus-per-task=", threads, "\n"), file = bash_script, append = TRUE)
  cat(paste0("#SBATCH --mem=", cluster_mem, "\n"), file = bash_script, append = TRUE)
  cat(paste0("#SBATCH -t ", cluster_time, "\n"), file = bash_script, append = TRUE)
  cat(paste0("#SBATCH --array=1-", max_reps, "\n"), file = bash_script, append = TRUE)
  cat("#SBATCH -o slurm_%A_%a.out\n", file = bash_script, append = TRUE)
  cat("#SBATCH -e slurm_%A_%a.err\n", file = bash_script, append = TRUE)
  
  if (!is.null(cluster_queue) && cluster_queue != "") {
    cat(paste0("#SBATCH -q ", cluster_queue, "\n"), file = bash_script, append = TRUE)
  }
  cat("#SBATCH --mail-type=ALL\n", file = bash_script, append = TRUE)
  if (!is.null(cluster_mail_user) && nzchar(trimws(cluster_mail_user))) {
    cat(paste0("#SBATCH --mail-user=", trimws(cluster_mail_user), "\n"), file = bash_script, append = TRUE)
  }
  cat("#\n\n", file = bash_script, append = TRUE)
  
  # Module loading
  if (!is.null(load_module) && length(load_module) > 0) {
    modules_str <- paste(load_module, collapse = " ")
    if (nzchar(trimws(modules_str))) {
      cat(paste0("module load ", modules_str, "\n\n"), file = bash_script, append = TRUE)
    }
  }

  cat("echo \"=== RAxML-NG version used by this job ===\"\n", file = bash_script, append = TRUE)
  cat(paste0(raxml_exec, " --version | head -n 2\n"), file = bash_script, append = TRUE)
  cat("echo \"=========================================\"\n\n", file = bash_script, append = TRUE)

  # Dynamic directory resolution
  cat("# Move to submit directory and dynamically resolve paths\n", file = bash_script, append = TRUE)
  cat("cd \"${SLURM_SUBMIT_DIR:-$PWD}\"\n", file = bash_script, append = TRUE)
  cat("SCRIPT_DIR=\"$PWD\"\n", file = bash_script, append = TRUE)
  cat("RUN_DIR=\"$(cd \"${SCRIPT_DIR}/..\" && pwd)\"\n", file = bash_script, append = TRUE)
  cat("BS_DIR=\"${SCRIPT_DIR}\"\n\n", file = bash_script, append = TRUE)

  cat(paste0("ALIGNMENT_FILE=\"${RUN_DIR}/", aln_basename, "\"\n"), file = bash_script, append = TRUE)
  cat(paste0("PARTITION_FILE=\"${RUN_DIR}/", part_basename, "\"\n"), file = bash_script, append = TRUE)
  cat(paste0("CONSTRAINT_FILE=\"${RUN_DIR}/", cons_basename, "\"\n\n"), file = bash_script, append = TRUE)

  # Worker and thread calculation
  cat(paste0("TOTAL_THREADS=${SLURM_CPUS_PER_TASK:-", threads, "}\n"), file = bash_script, append = TRUE)
  if (is.null(workers)) {
    cat(paste0("NUM_WORKERS=$(( TOTAL_THREADS / ", threads_per_worker, " ))\n"), file = bash_script, append = TRUE)
    cat("if [ \"$NUM_WORKERS\" -lt 1 ]; then\n    NUM_WORKERS=1\nfi\n\n", file = bash_script, append = TRUE)
  } else {
    cat(paste0("NUM_WORKERS=", workers_val, "\n\n"), file = bash_script, append = TRUE)
  }

  # Dynamic seed and chunk prefix
  cat(paste0("BASE_SEED=", base_seed, "\n"), file = bash_script, append = TRUE)
  cat("SEED=$(( SLURM_ARRAY_TASK_ID * BASE_SEED ))\n\n", file = bash_script, append = TRUE)

  cat("export OMP_PROC_BIND=false\n", file = bash_script, append = TRUE)
  cat("export OMPI_MCA_hwloc_base_binding_policy=none\n\n", file = bash_script, append = TRUE)

  CHUNK_DIR="${BS_DIR}/bs_rep_${SLURM_ARRAY_TASK_ID}"
  cat(paste0("CHUNK_DIR=\"", CHUNK_DIR, "\"\n"), file = bash_script, append = TRUE)
  cat("mkdir -p \"$CHUNK_DIR\"\n", file = bash_script, append = TRUE)
  cat("PREFIX=\"${CHUNK_DIR}/cactus_bs_rep_${SLURM_ARRAY_TASK_ID}\"\n", file = bash_script, append = TRUE)
  cat("RBA_PREFIX=\"${CHUNK_DIR}/cactus_alignment\"\n", file = bash_script, append = TRUE)
  cat("RBA_FILE=\"${RBA_PREFIX}.raxml.rba\"\n\n", file = bash_script, append = TRUE)

  # Pre-parsing step
  if (preparse) {
    cat("if [ ! -f \"$RBA_FILE\" ]; then\n", file = bash_script, append = TRUE)
    cat("    echo \"Pre-parsing alignment to compressed binary RBA format...\"\n", file = bash_script, append = TRUE)
    cat(paste0("    ", raxml_exec, " --parse --msa \"$ALIGNMENT_FILE\" --model \"$PARTITION_FILE\" --prefix \"$RBA_PREFIX\"\n"), file = bash_script, append = TRUE)
    cat("fi\n\n", file = bash_script, append = TRUE)
  }

  cat("echo \"Running bootstrap replicate chunk ${SLURM_ARRAY_TASK_ID} on ${TOTAL_THREADS} threads (${NUM_WORKERS} workers), seed: ${SEED}\"\n\n", file = bash_script, append = TRUE)
  
  # RAxML-NG bootstrap execution
  cmd <- paste0(raxml_exec, " --bootstrap ")
  if (preparse) {
    cmd <- paste0(cmd, "--msa \"$RBA_FILE\" ")
  } else {
    cmd <- paste0(cmd, "--msa \"$ALIGNMENT_FILE\" --model \"$PARTITION_FILE\" ")
  }
  
  cmd <- paste0(
    cmd,
    "--tree-constraint \"$CONSTRAINT_FILE\" "
  )
  
  og <- .format_outgroup(outgroup)
  if (!is.null(og)) {
    cmd <- paste0(cmd, "--outgroup ", shQuote(og, type = "sh"), " ")
  }
  
  # `--bs-metric` is not written into the cluster script either: in `--bootstrap` mode RAxML-NG
  # ignores it and writes only the replicate topologies. Declaring it here suggested that the
  # replicates produced on the cluster were locked to one metric, which they are not.
  cmd <- paste0(
    cmd,
    "--bs-trees ", bs_per_rep, " ",
    "--threads ${TOTAL_THREADS} ",
    "--workers ${NUM_WORKERS} ",
    "--seed ${SEED} ",
    "--force perf_threads ",
    "--extra thread-nopin ",
    "--prefix \"${PREFIX}\"\n"
  )
  cat(cmd, file = bash_script, append = TRUE)
  
  if (.Platform$OS.type == "unix") {
    chmod_status <- system(paste("chmod +x", shQuote(bash_script)))
    .check_cli_exit(chmod_status, "chmod")
  } else {
    Sys.chmod(bash_script, mode = "0755")
  }
  message("Bash SLURM array script generated: ", bash_script)
  return(bash_script)
}

#' Collect and Concatenate Parallel Bootstrap Tree Outputs
#'
#' Scans a target directory for chunked bootstrap output files (`.raxml.bootstraps`) generated by SLURM Job Arrays or local runs,
#' verifies that all expected replicates are present, and concatenates them into a single Newick tree file.
#'
#' @param bs_dir Character. Directory path containing chunked bootstrap output files (e.g., `bs_chunks/`).
#' @param output_dir Character. Output directory path to save concatenated bootstrap file. Defaults to `bs_dir`.
#' @param prefix Character. File output prefix. Defaults to `"cactus_ALL_bootstraps"`.
#' @param expected_trees Integer. Optional expected total count of bootstrap tree replicates (e.g., `1000`). If provided, verifies that the concatenated tree count exactly matches this value. Defaults to `NULL`.
#' @return Character path to the concatenated bootstrap tree file (`.tree`).
#' @export
collect_bootstraps <- function(bs_dir, output_dir = bs_dir, prefix = "cactus_ALL_bootstraps", expected_trees = NULL) {
  
  if (!dir.exists(bs_dir)) {
    stop("Bootstrap directory does not exist: ", bs_dir)
  }

  bs_files <- list.files(bs_dir, pattern = "\\.raxml\\.bootstraps$", recursive = TRUE, full.names = TRUE)
  if (length(bs_files) == 0) {
    stop("No bootstrap tree files (*.raxml.bootstraps) found in directory: ", bs_dir)
  }
  
  # Natural numeric sorting of files
  bs_files <- bs_files[order(nchar(bs_files), bs_files)]
  message("Found ", length(bs_files), " bootstrap chunk file(s) in: ", bs_dir)
  
  tree_list <- lapply(bs_files, function(f) {
    tr <- tryCatch(
      ape::read.tree(f),
      error = function(e) {
        stop("Failed to read bootstrap trees from: ", f, " | Error: ", conditionMessage(e))
      }
    )
    if (inherits(tr, "phylo")) {
      tr <- list(tr)
      class(tr) <- "multiPhylo"
    }
    return(tr)
  })

  all_bs_trees <- do.call(c, tree_list)
  class(all_bs_trees) <- "multiPhylo"
  total_trees <- length(all_bs_trees)
  
  if (!is.null(expected_trees)) {
    expected_trees <- as.integer(expected_trees)
    if (total_trees != expected_trees) {
      stop("Bootstrap replicate integrity check failed: expected ", expected_trees,
           " bootstrap trees, but collected ", total_trees, " trees across ", length(bs_files),
           " chunk file(s). Please check SLURM job array logs for failed or incomplete tasks.")
    }
    message("Bootstrap integrity verified: ", total_trees, " of ", expected_trees, " expected trees collected.")
  } else {
    message("Collected ", total_trees, " bootstrap trees across ", length(bs_files), " chunk file(s).")
  }
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  bs_concat_file <- file.path(output_dir, paste0(prefix, ".tree"))
  ape::write.tree(all_bs_trees, file = bs_concat_file)
  
  message("Concatenated bootstrap trees saved to: ", bs_concat_file)
  return(bs_concat_file)
}

#' Parse RAxML-NG autoMRE Bootstrap Convergence Log Output
#'
#' Internal helper used by `check_bs_convergence()`. Parses the textual `.raxml.log` output of
#' `RAxML-NG --bsconverge` to extract whether the autoMRE bootstopping criterion converged, at
#' how many trees, and how many bootstrap trees were analyzed in total. Extracted as a standalone,
#' independently testable package-internal function (previously inline code inside
#' `check_bs_convergence()`), so the parsing logic itself can be validated against fixture log
#' text without requiring a real `RAxML-NG` installation.
#'
#' @param log_lines Character vector of lines from a `RAxML-NG --bsconverge` `.raxml.log` file.
#' @return A list with `converged` (logical), `trees_at_convergence` (integer or `NA`), and
#'   `trees_analyzed` (integer or `NA`).
#' @noRd
.parse_bs_convergence_log <- function(log_lines) {
  converged <- FALSE
  trees_at_convergence <- NA_integer_
  trees_analyzed <- NA_integer_

  # Check for "Bootstopping test converged after X trees"
  conv_match <- grep("Bootstopping test converged after\\s+(\\d+)\\s+trees", log_lines, value = TRUE, ignore.case = TRUE)
  if (length(conv_match) > 0) {
    converged <- TRUE
    trees_at_convergence <- as.integer(sub(".*converged after\\s+(\\d+)\\s+trees.*", "\\1", conv_match[1], ignore.case = TRUE))
  }

  # Check for "Bootstopping test did not converge after X trees"
  not_conv_match <- grep("Bootstopping test did not converge after\\s+(\\d+)\\s+trees", log_lines, value = TRUE, ignore.case = TRUE)
  if (length(not_conv_match) > 0) {
    converged <- FALSE
    trees_analyzed <- as.integer(sub(".*did not converge after\\s+(\\d+)\\s+trees.*", "\\1", not_conv_match[1], ignore.case = TRUE))
  }

  # Fallback: check table rows with YES / NO
  table_yes <- grep("^\\s*(\\d+)\\s+.*YES\\s*$", log_lines, value = TRUE)
  if (!converged && length(table_yes) > 0) {
    converged <- TRUE
    trees_at_convergence <- as.integer(sub("^\\s*(\\d+)\\s+.*", "\\1", table_yes[1]))
  }

  # Parse total loaded trees: "Loaded X trees with Y taxa"
  loaded_match <- grep("Loaded\\s+(\\d+)\\s+trees", log_lines, value = TRUE)
  if (length(loaded_match) > 0) {
    trees_analyzed <- as.integer(sub(".*Loaded\\s+(\\d+)\\s+trees.*", "\\1", loaded_match[1]))
  }

  if (is.na(trees_analyzed) && !is.na(trees_at_convergence)) {
    trees_analyzed <- trees_at_convergence
  }

  list(converged = converged, trees_at_convergence = trees_at_convergence, trees_analyzed = trees_analyzed)
}

#' Check Bootstrap Convergence Criterion in RAxML-NG
#'
#' Evaluates whether the generated pool of non-parametric bootstrap trees has achieved statistical
#' convergence (`RAxML-NG`), using the autoMRE bootstopping criterion at a permutation cutoff
#' (typically 0.03).
#'
#' The criterion is metric-specific and the two metrics do not converge at the same replicate count.
#' TBE is the more stable summary and reaches the cutoff earlier, so a pool declared converged under
#' TBE is not thereby converged under FBP. The default is FBP, matching `map_branch_supports()`, so
#' that the convergence statement and the reported support values refer to the same quantity. Run the
#' test under whichever metric is being reported, and under both when both are tabulated.
#'
#' @param raxml_bin_path Character. System command or full path to executable `RAxML-NG` binary.
#' @param bs_trees_file Character. Path to concatenated bootstrap tree file (e.g., `cactus_ALL_bootstraps.tree`).
#' @param bs_cutoff Numeric. Permutation cutoff threshold for convergence. Defaults to `0.03`.
#' @param bs_metric Character. Branch support metric for the bootstopping test: `"fbp"` (Felsenstein's Bootstrap Percentage, the default) or `"tbe"` (Transfer Bootstrap Expectation). Must match the metric being reported.
#' @param seed Integer. Random seed for reproducible convergence testing. Defaults to `NULL` (random).
#' @param threads Integer. Number of CPU threads. Defaults to `4`.
#' @param output_dir Character. Directory path to save convergence report logs. Defaults to `dirname(bs_trees_file)`.
#' @param prefix Character. Output file prefix. Defaults to `"cactus_bs_convergence"`.
#' @return A structured S3 object of class `cactus_bs_convergence` containing:
#' \itemize{
#'   \item `converged`: Logical indicating whether the autoMRE convergence criterion was satisfied.
#'   \item `trees_at_convergence`: Integer number of trees required to reach convergence (or `NA` if not converged).
#'   \item `trees_analyzed`: Total number of bootstrap trees evaluated.
#'   \item `cutoff`: Numeric convergence cutoff threshold.
#'   \item `metric`: Support metric used (`"fbp"` or `"tbe"`).
#'   \item `log_file`: Character path to the resulting RAxML-NG convergence log file (`.raxml.log`).
#' }
#' @export
check_bs_convergence <- function(raxml_bin_path, bs_trees_file, bs_cutoff = 0.03,
                                 bs_metric = c("fbp", "tbe"), seed = NULL, threads = 4,
                                 output_dir = dirname(bs_trees_file), prefix = "cactus_bs_convergence") {

  bs_metric <- match.arg(bs_metric)

  if (Sys.which(raxml_bin_path) == "") {
    stop("Executable '", raxml_bin_path, "' not found in your system's PATH.\n",
         "Please ensure RAxML-NG is installed and available, or provide the full absolute path.")
  }

  if (!file.exists(bs_trees_file)) {
    stop("Bootstrap tree file not found: ", bs_trees_file)
  }

  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
    message("Using generated seed: ", seed)
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_prefix <- file.path(output_dir, prefix)
  
  args <- c(
    "--bsconverge",
    "--bs-trees", shQuote(bs_trees_file),
    "--bs-metric", shQuote(bs_metric),
    "--bs-cutoff", as.character(bs_cutoff),
    "--seed", as.character(seed),
    "--threads", as.character(threads),
    "--prefix", shQuote(out_prefix)
  )
  
  status <- system2(command = raxml_bin_path, args = args)
  .check_cli_exit(status, "RAxML-NG Convergence")
  
  log_file <- paste0(out_prefix, ".raxml.log")
  if (!file.exists(log_file)) {
    stop("Convergence log file was not generated at: ", log_file)
  }
  
  log_lines <- readLines(log_file, warn = FALSE)
  parsed <- .parse_bs_convergence_log(log_lines)
  converged <- parsed$converged
  trees_at_convergence <- parsed$trees_at_convergence
  trees_analyzed <- parsed$trees_analyzed

  result <- structure(
    list(
      converged = converged,
      trees_at_convergence = trees_at_convergence,
      trees_analyzed = trees_analyzed,
      cutoff = bs_cutoff,
      metric = bs_metric,
      log_file = log_file
    ),
    class = c("cactus_bs_convergence", "list")
  )
  
  cat("\n====================================================\n")
  cat("  RAxML-NG Bootstrap Convergence Assessment (autoMRE)\n")
  cat("====================================================\n")
  cat("  Converged:             ", if (converged) "YES \U0001f335" else "NO \U0000274c", "\n")
  if (converged) {
    cat("  Trees at convergence:  ", trees_at_convergence, "\n")
  }
  cat("  Trees analyzed:        ", if (!is.na(trees_analyzed)) trees_analyzed else "N/A", "\n")
  cat("  Cutoff threshold:      ", bs_cutoff, "\n")
  cat("  Support metric:        ", toupper(bs_metric), "\n")
  cat("  Log file:              ", log_file, "\n")
  cat("====================================================\n\n")

  if (!converged) {
    warning("Bootstrap convergence NOT achieved with ", if (!is.na(trees_analyzed)) trees_analyzed else "current",
            " trees (cutoff <= ", bs_cutoff, "). Consider generating additional bootstrap replicates.", call. = FALSE)
  }
  
  return(invisible(result))
}

#' @export
print.cactus_bs_convergence <- function(x, ...) {
  cat("====================================================\n")
  cat("  RAxML-NG Bootstrap Convergence Assessment (autoMRE)\n")
  cat("====================================================\n")
  cat("  Converged:             ", if (x$converged) "YES \U0001f335" else "NO \U0000274c", "\n")
  if (x$converged) {
    cat("  Trees at convergence: ", x$trees_at_convergence, "\n")
  }
  cat("  Trees analyzed:        ", x$trees_analyzed, "\n")
  cat("  Cutoff threshold:      ", x$cutoff, "\n")
  cat("  Support metric:        ", toupper(x$metric), "\n")
  cat("  Log file:              ", x$log_file, "\n")
  cat("====================================================\n")
  invisible(x)
}

#' @export
as.character.cactus_bs_convergence <- function(x, ...) {
  return(x$log_file)
}

#' Estimate Temporal Bootstrap Replicates Constrained to Best ML Topology
#'
#' Generates non-parametric bootstrap tree replicates where branch lengths are re-estimated while holding the focal maximum-likelihood topology constrained.
#' Temporal bootstraps propagate branch length uncertainty into downstream penalized likelihood divergence time estimation (`treePL`),
#' providing empirical confidence intervals for node ages without introducing topological variance.
#'
#' @param raxml_bin_path Character. System command or full path to executable `RAxML-NG` binary.
#' @param aln_file Character. Path to input PHYLIP supermatrix alignment file.
#' @param part_file Character. Path to partition file specifying substitution models.
#' @param best_tree_file Character. Path to reference maximum-likelihood tree topology file (used as constraint).
#' @param bs_trees Integer. Total number of temporal bootstrap trees to generate. Defaults to `500`.
#' @param outgroup Character vector of terminals passed to `RAxML-NG --outgroup`, or `NULL` (default); multiple terminals are joined with commas. `RAxML-NG` writes an unrooted topology with these terminals placed first, so this argument orders the output rather than rooting the tree: the root is imposed downstream by `automate_treePL()` via `ape::root(..., resolve.root = TRUE)`. Declaring the same set at every stage keeps the output ordering consistent across the maximum-likelihood search, the bootstrap replicates and the temporal bootstraps. Derive it with `resolve_rooting_outgroup()` rather than naming a terminal by hand.
#' @param seed Integer. Random seed for reproducible temporal bootstrap initialization. Defaults to `NULL` (random).
#' @param threads Integer. Number of CPU threads. Defaults to `4`.
#' @param workers Integer or NULL. Number of parallel tree search workers. If `NULL`, auto-configured. Defaults to `NULL`.
#' @param blopt Character. Branch length optimization method (`"nr_safe"` or `"nr_fast"`). Defaults to `"nr_safe"` for robust numerical convergence across partitioned alignments.
#' @param output_dir Character. Directory path to save output temporal bootstrap trees. Defaults to `dirname(aln_file)`.
#' @param prefix Character. Output file prefix. Defaults to `"cactus_temporal"`.
#' @return Character path to the resulting temporal bootstrap trees file (`.raxml.bootstraps`).
#' @export
calculate_temporal_bootstraps <- function(raxml_bin_path, aln_file, part_file, best_tree_file,
                                          bs_trees = 500, outgroup = NULL, seed = NULL,
                                          threads = 4, workers = NULL, blopt = "nr_safe",
                                          output_dir = dirname(aln_file),
                                          prefix = "cactus_temporal") {
                                          
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
    message("Using generated seed: ", seed)
  }

  if (Sys.which(raxml_bin_path) == "") {
    stop("Executable '", raxml_bin_path, "' not found in your system's PATH.\n",
         "Please ensure RAxML-NG is installed and available, or provide the full absolute path.")
  }

  threads <- as.integer(threads)
  if (is.null(workers)) {
    workers_val <- max(1L, as.integer(threads / 4L))
  } else {
    workers_val <- as.integer(workers)
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_prefix <- file.path(output_dir, prefix)
  
  # No `--bs-metric`: these replicates feed `treePL`, which reads branch lengths, not node support.
  # In `--bootstrap` mode RAxML-NG drops the option in any case.
  args <- c(
    "--bootstrap",
    "--msa", shQuote(aln_file),
    "--model", shQuote(part_file),
    "--tree-constraint", shQuote(best_tree_file),
    "--bs-trees", as.character(bs_trees),
    "--seed", as.character(seed),
    "--threads", as.character(threads),
    "--workers", as.character(workers_val),
    "--blopt", blopt,
    "--force", "perf_threads",
    "--prefix", shQuote(out_prefix)
  )
  
  og <- .format_outgroup(outgroup)
  if (!is.null(og)) {
    args <- c(args, "--outgroup", shQuote(og))
  }
  
  status <- system2(command = raxml_bin_path, args = args)
  .check_cli_exit(status, "RAxML-NG Temporal")
  tb_file <- paste0(out_prefix, ".raxml.bootstraps")
  message("Temporal bootstrap trees saved to: ", tb_file)
  return(tb_file)
}


#' Generate High-Performance Computing (HPC/SLURM) Batch Script for Temporal Bootstrapping
#'
#' Writes a self-contained SLURM batch script for running high-throughput temporal bootstrap replicate inference
#' using MPI-enabled `RAxML-NG` (`raxml-ng-mpi`) on a compute cluster. Temporal bootstraps constrain the maximum-likelihood
#' tree topology and re-estimate branch lengths across bootstrap matrices, generating empirical distributions of branch
#' lengths for downstream divergence time estimation (`treePL`) without introducing topological discordance.
#'
#' @param alignment_file Character. Path to input PHYLIP alignment file.
#' @param partition_file Character. Path to partition file specifying substitution models.
#' @param best_tree_file Character. Path to best scoring maximum-likelihood tree topology file (used as constraint).
#' @param outgroup Character vector of terminals passed to `RAxML-NG --outgroup`, or `NULL` (default); multiple terminals are joined with commas. `RAxML-NG` writes an unrooted topology with these terminals placed first, so this argument orders the output rather than rooting the tree: the root is imposed downstream by `automate_treePL()` via `ape::root(..., resolve.root = TRUE)`. Declaring the same set at every stage keeps the output ordering consistent across the maximum-likelihood search, the bootstrap replicates and the temporal bootstraps. Derive it with `resolve_rooting_outgroup()` rather than naming a terminal by hand.
#' @param bs_trees Integer. Total number of temporal bootstrap trees to generate. Defaults to `500`.
#' @param base_seed Integer. Random seed for reproducible initialization. Defaults to `NULL` (random).
#' @param threads Integer. Number of CPU cores requested per SLURM task (`--cpus-per-task`). Defaults to `40`.
#' @param workers Integer or NULL. Number of RAxML-NG worker processes. If `NULL` (default), calculated dynamically as `max(1L, as.integer(threads / threads_per_worker))`.
#' @param threads_per_worker Integer. Thread-to-worker ratio used to derive `workers` when `workers = NULL`. Defaults to `4L`.
#' @param preparse Logical. Generate compressed binary `.rba` format with `raxml-ng --parse` prior to execution to minimize disk I/O contention? Defaults to `TRUE`.
#' @param output_dir Character. Destination directory for script and outputs. Defaults to `getwd()`.
#' @param script_name Character. Name of output Bash script file. Defaults to `"run_temporal_bs.sh"`.
#' @param cluster_job_name Character. SLURM job name identifier. Defaults to `"cactus_temp_bs"`.
#' @param cluster_partition Character. SLURM partition name. Defaults to `"main"`.
#' @param cluster_nodes Integer. Number of compute nodes requested (`--nodes`). Defaults to `1L`.
#' @param cluster_mem Character. Memory allocation string for SLURM (`--mem`). Defaults to `"16G"`.
#' @param cluster_time Character. Time limit allocation string for SLURM (`--time`). Defaults to `"04:00:00"`.
#' @param cluster_queue Character. Optional SLURM queue / QoS name. Defaults to `NULL`.
#' @param cluster_mail_user Character. Notification recipient email address for SLURM (`--mail-user`). Defaults to `Sys.getenv("MY_EMAIL", "")`.
#' @param load_module Character vector or string. Environment module(s) to load prior to execution. Defaults to `c("gcc/14.2.0-nlhpc", "openmpi/5.0.3-o", "raxml-ng/1.1.0-mpi-zen4-n")`.
#' @param raxml_exec Character. Executable `RAxML-NG` binary command. Defaults to `"raxml-ng-mpi"`.
#' @return Character path to the generated SLURM batch script file.
#' @export
generate_temporal_bootstrap_script <- function(alignment_file, partition_file, best_tree_file,
                                               outgroup = NULL, bs_trees = 500,
                                               base_seed = NULL,
                                               threads = 40, workers = NULL, threads_per_worker = 4L,
                                               preparse = TRUE,
                                               output_dir = getwd(),
                                               script_name = "run_temporal_bs.sh",
                                               cluster_job_name = "cactus_temp_bs",
                                               cluster_partition = "main",
                                               cluster_nodes = 1L,
                                               cluster_mem = "16G", cluster_time = "04:00:00",
                                               cluster_queue = NULL,
                                               cluster_mail_user = Sys.getenv("MY_EMAIL", ""),
                                               load_module = c("gcc/14.2.0-nlhpc", "openmpi/5.0.3-o", "raxml-ng/1.1.0-mpi-zen4-n"),
                                               raxml_exec = "raxml-ng-mpi") {
  
  if (is.null(base_seed)) {
    base_seed <- sample.int(100000L, 1L)
    message("Using generated base seed: ", base_seed)
  }

  threads <- as.integer(threads)
  threads_per_worker <- as.integer(threads_per_worker)
  if (is.null(workers)) {
    workers_val <- max(1L, as.integer(threads / threads_per_worker))
  } else {
    workers_val <- as.integer(workers)
  }

  temp_dir <- file.path(output_dir, "cactus_temporal_bs")
  if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
  
  bash_script <- file.path(temp_dir, script_name)
  aln_basename  <- basename(alignment_file)
  part_basename <- basename(partition_file)
  tree_basename <- basename(best_tree_file)

  cat("#!/bin/bash\n", file = bash_script)
  cat(paste0("#SBATCH -J ", cluster_job_name, "\n"), file = bash_script, append = TRUE)
  cat(paste0("#SBATCH -p ", cluster_partition, "\n"), file = bash_script, append = TRUE)
  cat(paste0("#SBATCH --nodes=", cluster_nodes, "\n"), file = bash_script, append = TRUE)
  cat("#SBATCH --ntasks-per-node=1\n", file = bash_script, append = TRUE)
  cat(paste0("#SBATCH --cpus-per-task=", threads, "\n"), file = bash_script, append = TRUE)
  cat(paste0("#SBATCH --mem=", cluster_mem, "\n"), file = bash_script, append = TRUE)
  cat(paste0("#SBATCH -t ", cluster_time, "\n"), file = bash_script, append = TRUE)
  cat("#SBATCH -o slurm_%j.out\n", file = bash_script, append = TRUE)
  cat("#SBATCH -e slurm_%j.err\n", file = bash_script, append = TRUE)
  
  if (!is.null(cluster_queue) && cluster_queue != "") {
    cat(paste0("#SBATCH -q ", cluster_queue, "\n"), file = bash_script, append = TRUE)
  }
  cat("#SBATCH --mail-type=ALL\n", file = bash_script, append = TRUE)
  if (!is.null(cluster_mail_user) && nzchar(trimws(cluster_mail_user))) {
    cat(paste0("#SBATCH --mail-user=", trimws(cluster_mail_user), "\n"), file = bash_script, append = TRUE)
  }
  cat("#\n\n", file = bash_script, append = TRUE)
  
  # Module loading
  if (!is.null(load_module) && length(load_module) > 0) {
    modules_str <- paste(load_module, collapse = " ")
    if (nzchar(trimws(modules_str))) {
      cat(paste0("module load ", modules_str, "\n\n"), file = bash_script, append = TRUE)
    }
  }

  cat("echo \"=== RAxML-NG version used by this job ===\"\n", file = bash_script, append = TRUE)
  cat(paste0(raxml_exec, " --version | head -n 2\n"), file = bash_script, append = TRUE)
  cat("echo \"=========================================\"\n\n", file = bash_script, append = TRUE)

  # Dynamic directory resolution
  cat("# Move to submit directory and dynamically resolve paths\n", file = bash_script, append = TRUE)
  cat("cd \"${SLURM_SUBMIT_DIR:-$PWD}\"\n", file = bash_script, append = TRUE)
  cat("SCRIPT_DIR=\"$PWD\"\n", file = bash_script, append = TRUE)
  cat("RUN_DIR=\"$(cd \"${SCRIPT_DIR}/..\" && pwd)\"\n", file = bash_script, append = TRUE)
  cat("TEMP_DIR=\"${SCRIPT_DIR}\"\n\n", file = bash_script, append = TRUE)

  cat(paste0("ALIGNMENT_FILE=\"${RUN_DIR}/", aln_basename, "\"\n"), file = bash_script, append = TRUE)
  cat(paste0("PARTITION_FILE=\"${RUN_DIR}/", part_basename, "\"\n"), file = bash_script, append = TRUE)
  cat(paste0("BEST_TREE_FILE=\"${RUN_DIR}/ml_search/", tree_basename, "\"\n\n"), file = bash_script, append = TRUE)

  cat("# If ml_search subdirectory does not contain the tree, check RUN_DIR directly\n", file = bash_script, append = TRUE)
  cat("if [ ! -f \"${BEST_TREE_FILE}\" ]; then\n", file = bash_script, append = TRUE)
  cat(paste0("  BEST_TREE_FILE=\"${RUN_DIR}/", tree_basename, "\"\n"), file = bash_script, append = TRUE)
  cat("fi\n\n", file = bash_script, append = TRUE)

  cat(paste0("TOTAL_THREADS=", threads, "\n"), file = bash_script, append = TRUE)
  cat(paste0("NUM_WORKERS=", workers_val, "\n"), file = bash_script, append = TRUE)
  cat(paste0("SEED=", base_seed, "\n\n"), file = bash_script, append = TRUE)

  if (preparse) {
    cat("# Pre-parse alignment to compressed binary RBA format\n", file = bash_script, append = TRUE)
    cat("PARSED_RBA=\"${TEMP_DIR}/cactus_temporal.raxml.rba\"\n", file = bash_script, append = TRUE)
    cat("if [ ! -f \"${PARSED_RBA}\" ]; then\n", file = bash_script, append = TRUE)
    cat("  echo \"Pre-parsing alignment to compressed binary RBA format...\"\n", file = bash_script, append = TRUE)
    cat(paste0("  ", raxml_exec, " --parse \\\n"), file = bash_script, append = TRUE)
    cat("    --msa \"${ALIGNMENT_FILE}\" \\\n", file = bash_script, append = TRUE)
    cat("    --model \"${PARTITION_FILE}\" \\\n", file = bash_script, append = TRUE)
    cat("    --prefix \"${TEMP_DIR}/cactus_temporal\"\n", file = bash_script, append = TRUE)
    cat("fi\n\n", file = bash_script, append = TRUE)
    cat("INPUT_MSA=\"${PARSED_RBA}\"\n\n", file = bash_script, append = TRUE)
  } else {
    cat("INPUT_MSA=\"${ALIGNMENT_FILE}\"\n\n", file = bash_script, append = TRUE)
  }

  # No `--bs-metric`, for the same reason as the local temporal bootstrap: these replicates are read
  # by `treePL` for their branch lengths, and `--bootstrap` mode ignores the option regardless.
  cmd <- paste0(
    raxml_exec, " ",
    "--bootstrap ",
    "--msa \"${INPUT_MSA}\" ",
    "--tree-constraint \"${BEST_TREE_FILE}\" ",
    "--bs-trees ", bs_trees, " ",
    "--threads ${TOTAL_THREADS} ",
    "--workers ${NUM_WORKERS} ",
    "--seed ${SEED} ",
    "--force perf_threads ",
    "--extra thread-nopin "
  )

  og <- .format_outgroup(outgroup)
  if (!is.null(og)) {
    cmd <- paste0(cmd, "--outgroup ", shQuote(og, type = "sh"), " ")
  }

  if (!preparse) {
    cmd <- paste0(cmd, "--model \"${PARTITION_FILE}\" ")
  }

  cmd <- paste0(cmd, "--prefix \"${TEMP_DIR}/cactus_temporal\"\n")

  cat("echo \"Starting temporal bootstrap tree inference...\"\n", file = bash_script, append = TRUE)
  cat(cmd, file = bash_script, append = TRUE)
  cat("\necho \"Temporal bootstrap inference complete. Output saved to: ${TEMP_DIR}/cactus_temporal.raxml.bootstraps\"\n", file = bash_script, append = TRUE)

  message("Generated temporal bootstrap HPC SLURM script: ", bash_script)
  return(bash_script)
}

                                          

# ------------------------------------------------------------------------------
# Single-locus gene tree inference. Relocated from R/barcoding.R in v0.4.2: the
# individual gene trees produced here feed the gene tree / species tree discordance
# analyses of Module 13 (ASTRAL-III, Stage 4), which are independent of the DNA
# barcoding strategy that remains under design.
# ------------------------------------------------------------------------------

#' Extract Standardized Species Binomial from Sequence Header
#'
#' Internal helper to parse species binomials from FASTA headers formatted as
#' `Genus_species|accession|sid` or `Genus species`.
#'
#' @param header Character string of the sequence header.
#' @return Standardized species binomial formatted with underscores (e.g., `Copiapoa_cinerea`).
#' @keywords internal
extract_species_binomial <- function(header) {
  vapply(header, function(h) {
    if (is.na(h) || !nzchar(h)) return(NA_character_)
    first_field <- stringr::str_split(h, "\\|", simplify = TRUE)[1]
    first_field <- stringr::str_trim(first_field)
    if (!nzchar(first_field)) return(NA_character_)
    clean_field <- stringr::str_replace_all(first_field, "[ -]+", "_")
    parts <- strsplit(clean_field, "_", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      paste0(parts[1], "_", parts[2])
    } else {
      clean_field
    }
  }, character(1), USE.NAMES = FALSE)
}

#' Infer Single-Locus Maximum-Likelihood Gene Trees and Assess Species Monophyly
#'
#' Reconstructs individual gene trees for all curated marker alignments using `RAxML-NG`
#' (with automated fallback to `phangorn` or `ape`), assesses reciprocal monophyly for
#' all multi-accession species, and prepares gene tree collections for coalescent analyses (e.g., `ASTRAL-III`).
#'
#' @param fasta_dir Character. Directory containing curated aligned FASTA sequence files (e.g., `4_Cleaned/cleaned_markers_ingroup` or `5_MAFFT_Cleaned/aligned_markers`).
#' @param output_dir Character. Root destination directory path to store Newick tree files and monophyly audit tables.
#' @param include_outgroup Logical. Include outgroup taxa in single-locus gene tree reconstruction? Defaults to `FALSE` to avoid Long-Branch Attraction (LBA) artifacts when evaluating species monophyly within the ingroup radiation. Set to `TRUE` when preparing unrooted gene trees for `ASTRAL-III`.
#' @param method Character. Phylogenetic inference method (`"auto"`, `"raxml"`, `"phangorn"`, or `"nj"`). Defaults to `"auto"`.
#' @param model Character. Nucleotide substitution model for maximum-likelihood search. Defaults to `"GTR+G"`.
#' @param threads Integer. Number of computational threads. Defaults to `2L`.
#' @param checklist_path Character. Path to accepted botanical checklist CSV (e.g., `CactaceaeFullList_accepted.csv`).
#' @return A list containing the species monophyly summary table and the per-marker tree paths.
#' @references
#' Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019).
#' RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference.
#' *Bioinformatics*, 35(21), 4453–4455. \doi{10.1093/bioinformatics/btz305}
#'
#' Schliep, K. P. (2011). phangorn: phylogenetic analysis in R.
#' *Bioinformatics*, 27(4), 592–593. \doi{10.1093/bioinformatics/btq706}
#' @examples
#' \dontrun{
#' infer_gene_trees(
#'   fasta_dir = "4_Cleaned/cleaned_markers_ingroup",
#'   output_dir = "4_Gene_Trees",
#'   include_outgroup = FALSE
#' )
#' }
#' @export
infer_gene_trees <- function(
  fasta_dir,
  output_dir,
  include_outgroup = FALSE,
  method = c("auto", "raxml", "phangorn", "nj"),
  model = "GTR+G",
  threads = 2L,
  checklist_path = NULL
) {
  
  method <- match.arg(method)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  trees_dir <- file.path(output_dir, "gene_trees")
  dir.create(trees_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Starting single-locus gene tree inference... \U0001f335")
  message("Input FASTA directory: ", fasta_dir)
  message("Include outgroup: ", include_outgroup)
  message("Output directory: ", output_dir)
  
  fasta_files <- sort(list.files(fasta_dir, pattern = "\\.fasta$", full.names = TRUE))
  if (length(fasta_files) == 0) {
    stop("No FASTA files found in: ", fasta_dir, call. = FALSE)
  }
  
  # Check if raxml-ng is available
  raxml_available <- FALSE
  if (method %in% c("auto", "raxml")) {
    chk <- suppressWarnings(system2("raxml-ng", args = "--version", stdout = FALSE, stderr = FALSE))
    raxml_available <- identical(chk, 0L)
    if (!raxml_available && method == "raxml") {
      warning("raxml-ng binary not found in $PATH. Falling back to phangorn ML.", call. = FALSE)
    }
  }
  
  base_engine <- if (method == "nj") {
    "nj"
  } else if (method == "phangorn" || (!raxml_available && method %in% c("auto", "raxml"))) {
    "phangorn"
  } else {
    "raxml"
  }
  message("Inference engine selected: ", base_engine)
  
  outgroup_genera <- c("Portulaca", "Anacampseros", "Talinopsis", "Grahamia", "Talinum", "Talinella")
  
  monophyly_rows <- list()
  tree_paths <- character(0)
  all_trees_list <- list()
  
  for (f in fasta_files) {
    marker_name <- tools::file_path_sans_ext(basename(f))
    message("Inferring gene tree for marker: ", marker_name)
    current_engine <- base_engine
    
    dna_in <- Biostrings::readDNAStringSet(f)
    if (length(dna_in) < 4) {
      message("Skipping ", marker_name, " (fewer than 4 sequences).")
      next
    }
    
    # Filter outgroups if include_outgroup is FALSE
    if (!include_outgroup) {
      sp_names <- vapply(names(dna_in), extract_species_binomial, character(1))
      genus_names <- sub("_.*", "", sp_names)
      is_out <- genus_names %in% outgroup_genera
      dna_in <- dna_in[!is_out]
    }
    
    if (length(dna_in) < 4) {
      message("Skipping ", marker_name, " (fewer than 4 sequences remaining after filtering).")
      next
    }
    
    # Standardize names to species binomials or unique tips
    sp_vec <- vapply(names(dna_in), extract_species_binomial, character(1))
    tip_names <- make.unique(sp_vec, sep = "_dup")
    names(dna_in) <- tip_names
    
    tree_out_file <- file.path(trees_dir, paste0(marker_name, ".tree"))
    tr <- NULL
    
    if (current_engine == "raxml") {
      temp_phy <- file.path(output_dir, paste0("TEMP_", marker_name, ".phy"))
      
      # Convert to PHYLIP
      seq_chars <- as.character(dna_in)
      taxa <- names(dna_in)
      header <- paste(length(taxa), nchar(seq_chars[1]))
      writeLines(c(header, paste(taxa, seq_chars)), temp_phy)
      
      cmd_args <- c(
        "--msa", temp_phy,
        "--model", model,
        "--threads", as.character(threads),
        "--prefix", file.path(output_dir, paste0("TEMP_RAXML_", marker_name)),
        "--seed", "111",
        "--redo"
      )
      
      status <- system2("raxml-ng", args = cmd_args, stdout = FALSE, stderr = FALSE)
      best_tree_file <- file.path(output_dir, paste0("TEMP_RAXML_", marker_name, ".raxml.bestTree"))
      
      if (identical(status, 0L) && file.exists(best_tree_file)) {
        tr <- ape::read.tree(best_tree_file)
        ape::write.tree(tr, file = tree_out_file)
      } else {
        message("RAxML-NG failed for ", marker_name, ". Falling back to phangorn.")
        current_engine <- "phangorn"
      }
      
      # Clean temp files
      temp_files <- list.files(output_dir, pattern = paste0("^TEMP_.*", marker_name), full.names = TRUE)
      unlink(temp_files)
    }
    
    if (current_engine == "phangorn" || (is.null(tr) && current_engine != "nj")) {
      dna_char <- as.character(dna_in)
      split_dna <- strsplit(toupper(dna_char), "", fixed = TRUE)
      dna_mat <- do.call(rbind, split_dna)
      dna_bin <- ape::as.DNAbin(dna_mat)
      rownames(dna_bin) <- names(dna_in)
      
      phy_dat <- phangorn::as.phyDat(dna_bin)
      dm <- phangorn::dist.ml(phy_dat)
      tr_nj <- ape::nj(dm)
      
      tr <- tryCatch({
        fit <- phangorn::pml(tr_nj, data = phy_dat)
        fit_opt <- phangorn::optim.pml(fit, model = "GTR", optGamma = TRUE, rearrangement = "NNI", control = phangorn::pml.control(trace = 0))
        fit_opt$tree
      }, error = function(e) {
        tr_nj
      })
      ape::write.tree(tr, file = tree_out_file)
    } else if (current_engine == "nj" && is.null(tr)) {
      dna_char <- as.character(dna_in)
      split_dna <- strsplit(toupper(dna_char), "", fixed = TRUE)
      dna_mat <- do.call(rbind, split_dna)
      dna_bin <- ape::as.DNAbin(dna_mat)
      rownames(dna_bin) <- names(dna_in)
      
      dm <- ape::dist.dna(dna_bin, model = "raw", pairwise.deletion = TRUE)
      tr <- ape::nj(dm)
      ape::write.tree(tr, file = tree_out_file)
    }
    
    if (!is.null(tr)) {
      tree_paths <- c(tree_paths, tree_out_file)
      all_trees_list[[marker_name]] <- tr
      
      # Assess species monophyly
      base_species <- sub("_dup.*", "", tr$tip.label)
      sp_counts <- table(base_species)
      multi_species <- names(sp_counts[sp_counts >= 2])
      
      n_mono <- 0L
      n_poly <- 0L
      mono_species_list <- character(0)
      non_mono_species_list <- character(0)
      
      if (length(multi_species) > 0) {
        for (sp in multi_species) {
          tips_sp <- tr$tip.label[base_species == sp]
          is_mono <- tryCatch(ape::is.monophyletic(tr, tips_sp), error = function(e) FALSE)
          if (isTRUE(is_mono)) {
            n_mono <- n_mono + 1L
            mono_species_list <- c(mono_species_list, sp)
          } else {
            n_poly <- n_poly + 1L
            non_mono_species_list <- c(non_mono_species_list, sp)
          }
        }
      }
      
      pct_mono <- if (length(multi_species) > 0) round((n_mono / length(multi_species)) * 100, 2) else NA_real_
      
      monophyly_rows[[marker_name]] <- dplyr::tibble(
        marker = marker_name,
        engine = current_engine,
        n_tips = length(tr$tip.label),
        n_species_total = length(unique(base_species)),
        n_multi_accession_species = length(multi_species),
        n_monophyletic_species = n_mono,
        n_non_monophyletic_species = n_poly,
        species_monophyly_rate_pct = pct_mono,
        monophyletic_species = paste(mono_species_list, collapse = ";"),
        non_monophyletic_species = paste(non_mono_species_list, collapse = ";")
      )
    }
  }
  
  # Concatenate all gene trees into one file for ASTRAL
  all_trees_file <- file.path(output_dir, "all_gene_trees.tree")
  if (length(all_trees_list) > 0) {
    all_trees_obj <- do.call(c, all_trees_list)
    class(all_trees_obj) <- "multiPhylo"
    ape::write.tree(all_trees_obj, file = all_trees_file)
  }
  
  compiled_monophyly <- dplyr::bind_rows(monophyly_rows)
  readr::write_csv(compiled_monophyly, file.path(output_dir, "TABLE_gene_tree_species_monophyly.csv"))
  
  message("\nSingle-locus gene tree inference finalized. \U0001f335")
  return(invisible(list(
    monophyly_summary = compiled_monophyly,
    tree_paths = tree_paths,
    all_trees_file = all_trees_file
  )))
}
