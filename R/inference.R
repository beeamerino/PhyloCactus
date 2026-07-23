#' Preprocess Partitions and Validate Alignment Syntax
#'
#' Validates PHYLIP alignment syntax and partition file coordinates using `RAxML-NG` (Kozlov *et al.*, 2019).
#' Verifies site ranges, formatting compatibility, and data integrity prior to substitution model evaluation.
#'
#' @param phy_matrix Character. Path to input PHYLIP supermatrix file.
#' @param part_file Character. Path to input partition mapping text file.
#' @param raxml_path Character. System command or full path to executable `RAxML-NG` binary.
#' @param output_dir Character. Output directory for validated partition outputs. Defaults to `dirname(phy_matrix)`.
#' @param force_check Logical. Bypass validation check if cleaned partition output file already exists? Defaults to `FALSE`.
#' @return Character path to the cleaned partition map file ready for model evaluation.
#' @references
#' Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and
#' user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics*, 35(21), 4453-4455.
#' \doi{10.1093/bioinformatics/btz305}
#' @export
preprocess_partitions <- function(phy_matrix, part_file, raxml_path, output_dir = dirname(phy_matrix), force_check = FALSE) {
  
  if (Sys.which(raxml_path) == "") {
    stop("Executable '", raxml_path, "' not found in your system's PATH.\n",
         "Please ensure RAxML-NG is installed and available, or provide the full absolute path.")
  }

  clean_part <- file.path(output_dir, "cactus_check.raxml.reduced.clean.partition")
  if (!force_check && file.exists(clean_part)) {
    message("CACHE: Partition validation checked already, loading: ", clean_part)
    return(clean_part)
  }
  
  prefix_check <- file.path(output_dir, "cactus_check")
  
  # Run system raxmlcheck validation
  aln_check <- system2(
    raxml_path,
    args = c(
      "--check",
      "--msa", shQuote(phy_matrix),
      "--model", shQuote(part_file),
      "--prefix", shQuote(prefix_check),
      "--threads 4"
    ),
    stdout = "",
    stderr = ""
  )
  
  # Read reduced output partition generated from check
  reduced_part_file <- file.path(output_dir, "cactus_check.raxml.reduced.partition")
  if (!file.exists(reduced_part_file)) {
    stop("Alignment validation check failed to create: ", reduced_part_file)
  }
  
  lines <- readLines(reduced_part_file, warn = FALSE)
  
  # Modeling to DNA replacements formatting
  dna_lines <- sub("^[^,]+,", "DNA,", lines)
  
  writeLines(dna_lines, clean_part)
  return(clean_part)
}

#' Evaluate Nucleotide Substitution Models via ModelTest-NG
#'
#' Evaluates nucleotide substitution model fit per predefined supermatrix partition using `ModelTest-NG` (Darriba *et al.*, 2020).
#' Selecting optimal substitution models under the AICc criterion controls for mutational rate heterogeneity across genomic regions,
#' mitigating systematic long-branch attraction (LBA) bias during maximum-likelihood inference.
#'
#' @param modeltest_exec_path Character. System command or full path to executable `ModelTest-NG` binary.
#' @param aln_file Character. Path to validated PHYLIP supermatrix file.
#' @param part_file Character. Path to cleaned partition mapping file.
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
#'   part_file = "PARTITION_raxml_ng.txt"
#' )
#' }
#' @export
run_modeltest_ng <- function(modeltest_exec_path, aln_file, part_file, prefix = "MODELTEST_cactus_phylo", threads = 4) {
  
  if (Sys.which(modeltest_exec_path) == "") {
    stop("Executable '", modeltest_exec_path, "' not found in your system's PATH.\n",
         "Please ensure ModelTest-NG is installed and available, or provide the full absolute path.")
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
  
  # Best model file returned path
  expected_out <- paste0(prefix, ".part.aicc")
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
#' @return Character string path to the saved Newick constraint tree file (`cactus_constraints.tree`).
#' @examples
#' \dontrun{
#' build_constraint_scaffold(
#'   alignment_path = "ALIGNMENT_supermatrix.phy",
#'   constraints_csv_path = "Cactaceae_taxonomy.csv"
#' )
#' }
#' @export
build_constraint_scaffold <- function(alignment_path, constraints_csv_path, output_dir = dirname(alignment_path)) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  out_missing     <- file.path(output_dir, "TABLE_cactus_alignment_taxa_missing_from_constraints.csv")
  out_dup_clade   <- file.path(output_dir, "TABLE_cactus_constraint_duplicate_clade_membership.csv")
  out_clade_sizes <- file.path(output_dir, "TABLE_cactus_constraint_clade_sizes.csv")
  out_tree        <- file.path(output_dir, "cactus_constraints.tree")
  out_log         <- file.path(output_dir, "LOG_cactus_constraint_tree.txt")
  
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
    read.csv(constraints_csv_path, stringsAsFactors = FALSE),
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
    write.csv(data.frame(Species_missing = missing_taxa), out_missing, row.names = FALSE)
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
    write.csv(dup_detail, out_dup_clade, row.names = FALSE)
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
  write.csv(clade_sizes, out_clade_sizes, row.names = FALSE)
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
  
  # 7. BUILD FIXED TOPOLOGY
  s_outgroup <- paste0("(", s_por, ",(", s_tal, ",(", s_gra, ",", s_ana, ")))")
  s_opuntioideae <- paste0("(", s_cyl, ",(", s_teph, ",", s_opu, "))")
  s_coreII <- paste0("(", s_rhip, ",(", s_noto, ",", s_bct, "))")
  s_core_cactoideae <- paste0("(", s_frai, ",", s_caly, ",", s_copi, ",(", s_coreII, ",", s_coreI, "))")
  s_cactoideae <- paste0("(", s_cact, ",", s_core_cactoideae, ")")
  s_cactaceae <- paste0("(", s_leu, ",(", s_per, ",(", s_opuntioideae, ",(", s_mai, ",(", s_blo, ",", s_cactoideae, ")))))")
  
  constraint_tree <- paste0("(", s_outgroup, ",", s_cactaceae, ");")
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
#' @param outgroup Character. Optional outgroup taxon binomial to root the resulting tree topology. Defaults to `NULL`.
#' @param n_init_trees Character. Initial starting tree specifications. Defaults to `"rand{25},pars{25}"` (25 random + 25 parsimony trees).
#' @param seed Integer. Random seed for reproducible tree search initialization. Defaults to `111`.
#' @param n_workers Integer. Parallel worker process count for RAxML-NG. Defaults to `1`.
#' @param threads Integer. Number of CPU threads per worker. Defaults to `4`.
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
                              seed = 111, n_workers = 1, threads = 4, 
                              output_dir = dirname(aln_file), prefix = "cactus_search") {
  
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
    "--workers", as.character(n_workers),
    "--threads", as.character(threads),
    "--prefix", shQuote(out_prefix)
  )
  
  if (!is.null(outgroup) && outgroup != "") {
    args <- c(args, "--outgroup", shQuote(outgroup))
  }
  
  system2(command = raxml_bin_path, args = args)
  
  message("\nMaximum-likelihood tree calculation completed. \U0001f335")
  return(list(
    bestTree = paste0(out_prefix, ".raxml.bestTree"),
    mlTrees  = paste0(out_prefix, ".raxml.mlTrees")
  ))
}

#' Map Transfer Bootstrap Expectation (TBE) Support Values onto Reference Phylogeny
#'
#' Maps statistical clade support metrics derived from non-parametric bootstrap replicates onto the best maximum-likelihood tree topology.
#' Implements Transfer Bootstrap Expectation (TBE; Lemoine *et al.*, 2018) as the primary support metric, which provides
#' robust support evaluation for large plant phylogenies without underestimating support for minor clade position shifts.
#'
#' @param raxml_bin Character. System command or full path to executable `RAxML-NG` binary.
#' @param best_tree Character. Path to reference maximum-likelihood tree file.
#' @param bootstraps_file Character. Path to concatenated non-parametric bootstrap trees file.
#' @param metric Character. Bootstrap support metric: `"tbe"` (Transfer Bootstrap Expectation) or `"fbp"` (Felsenstein's Bootstrap Percentage). Defaults to `"tbe"`.
#' @param threads Integer. Number of CPU threads. Defaults to `4`.
#' @param output_dir Character. Output directory for annotated support tree. Defaults to `dirname(best_tree)`.
#' @param prefix Character. Output file prefix. Defaults to `"cactus_support"`.
#' @return Character path to the annotated support tree file (`.raxml.support`).
#' @references
#' Lemoine, F., Entfellner, J. B., Gascuel, O., & Gascuel, O. (2018). Renewing Felsenstein’s phylogenetic
#' bootstrap in the era of big data. *Nature*, 556(7702), 452-456. \doi{10.1038/s41586-018-0043-0}
#' @examples
#' \dontrun{
#' map_branch_supports(
#'   raxml_bin = "raxml-ng",
#'   best_tree = "cactus_search.raxml.bestTree",
#'   bootstraps_file = "cactus_ALL_bootstraps.tree",
#'   metric = "tbe"
#' )
#' }
#' @export
map_branch_supports <- function(raxml_bin, best_tree, bootstraps_file, metric = "tbe", threads = 4, output_dir = dirname(best_tree), prefix = "cactus_support") {
  
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
  
  system2(command = raxml_bin, args = args)
  return(paste0(out_prefix, ".raxml.support"))
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
  
  system2(command = raxml_bin_path, args = args)
  return(paste0(out_prefix, ".raxml.rfdist"))
}

#' Generate Non-Parametric Bootstrap Trees Locally
#'
#' Performs non-parametric bootstrap resampling over supermatrix site columns to infer a distribution of bootstrap tree topologies (`RAxML-NG`).
#' Evaluates topological variation under non-parametric resampling to quantify node support via Transfer Bootstrap Expectation (TBE).
#'
#' @param raxml_bin_path Character. System command or full path to executable `RAxML-NG` binary.
#' @param aln_file Character. Path to input PHYLIP supermatrix alignment file.
#' @param part_file Character. Path to partition file specifying substitution models.
#' @param constraint_file Character. Path to Newick topological constraint scaffold file.
#' @param bs_trees Integer. Total number of non-parametric bootstrap trees to generate. Defaults to `500`.
#' @param outgroup Character. Optional outgroup taxon binomial to root bootstrap topologies. Defaults to `NULL`.
#' @param seed Integer. Random seed for reproducible bootstrap initialization. Defaults to `111`.
#' @param threads Integer. Number of CPU threads. Defaults to `8`.
#' @param workers Integer. Parallel worker process count. Defaults to `1`.
#' @param output_dir Character. Output directory for generated bootstrap trees. Defaults to `dirname(aln_file)`.
#' @param prefix Character. Output file prefix. Defaults to `"cactus_bs"`.
#' @return Character path to the output bootstrap trees file (`.raxml.bootstraps`).
#' @references
#' Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and
#' user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics*, 35(21), 4453-4455.
#' \doi{10.1093/bioinformatics/btz305}
#' @export
run_local_bootstraps <- function(raxml_bin_path, aln_file, part_file, constraint_file,
                                 bs_trees = 500, outgroup = NULL, seed = 111,
                                 threads = 8, workers = 1, output_dir = dirname(aln_file),
                                 prefix = "cactus_bs") {
                                 
  if (Sys.which(raxml_bin_path) == "") {
    stop("Executable '", raxml_bin_path, "' not found in your system's PATH.\n",
         "Please ensure RAxML-NG is installed and available, or provide the full absolute path.")
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_prefix <- file.path(output_dir, prefix)
  
  args <- c(
    "--bootstrap",
    "--msa", shQuote(aln_file),
    "--model", shQuote(part_file),
    "--tree-constraint", shQuote(constraint_file),
    "--bs-trees", as.character(bs_trees),
    "--bs-metric", "tbe",
    "--seed", as.character(seed),
    "--workers", as.character(workers),
    "--threads", as.character(threads),
    "--prefix", shQuote(out_prefix)
  )
  
  if (!is.null(outgroup) && outgroup != "") {
    args <- c(args, "--outgroup", shQuote(outgroup))
  }
  
  system2(command = raxml_bin_path, args = args)
  return(paste0(out_prefix, ".raxml.bootstraps"))
}

#' Generate HPC SLURM Batch Script for Parallel Bootstrapping
#'
#' Generates an executable Bash script with SLURM scheduler directives to parallelize non-parametric bootstrapping across HPC compute nodes.
#' Chunking bootstrap searches into parallel sub-jobs accelerates support estimation for large supermatrices.
#'
#' @param alignment_file Character. Path to input PHYLIP alignment file.
#' @param partition_file Character. Path to partition file.
#' @param constraint_file Character. Path to constraint scaffold tree file.
#' @param outgroup Character. Optional outgroup taxon binomial. Defaults to `NULL`.
#' @param bs_per_rep Integer. Number of bootstrap trees generated per chunk replicate. Defaults to `500`.
#' @param max_reps Integer. Total number of parallel chunk replicates to spawn. Defaults to `4`.
#' @param base_seed Integer. Base random seed. Defaults to `111`.
#' @param seed_step Integer. Seed increment value between consecutive chunks. Defaults to `1000`.
#' @param threads Integer. Number of CPU cores requested per SLURM task. Defaults to `120`.
#' @param workers Integer. Number of RAxML-NG worker processes. Defaults to `20`.
#' @param output_dir Character. Destination directory for script and chunk logs. Defaults to `getwd()`.
#' @param script_name Character. Name of output Bash script file. Defaults to `"run_bs_chunks.sh"`.
#' @param cluster_job_name Character. SLURM job name identifier. Defaults to `"C-raxml.bs"`.
#' @param cluster_mem Character. Memory allocation string for SLURM. Defaults to `"32G"`.
#' @param cluster_time Character. Time limit allocation string for SLURM. Defaults to `"7-00:00"`.
#' @param cluster_partition Character. SLURM partition name. Defaults to `"general"`.
#' @param cluster_queue Character. SLURM queue name. Defaults to `"public"`.
#' @param load_module Character. Environment module command to load prior to execution. Defaults to `"raxml-ng-1.1.0-gcc-11.2.0"`.
#' @param raxml_exec Character. Executable `RAxML-NG` binary name. Defaults to `"raxml-ng-mpi"`.
#' @return Character path to the generated SLURM batch script file.
#' @export
generate_bootstrap_script <- function(alignment_file, partition_file, constraint_file,
                                      outgroup = NULL, bs_per_rep = 500, max_reps = 4,
                                      base_seed = 111, seed_step = 1000, 
                                      threads = 120, workers = 20,
                                      output_dir = getwd(),
                                      script_name = "run_bs_chunks.sh",
                                      cluster_job_name = "C-raxml.bs", 
                                      cluster_mem = "32G", cluster_time = "7-00:00",
                                      cluster_partition = "general", cluster_queue = "public",
                                      load_module = "raxml-ng-1.1.0-gcc-11.2.0",
                                      raxml_exec = "raxml-ng-mpi") {
  
  bs_dir <- file.path(output_dir, "bs_chunks")
  if (!dir.exists(bs_dir)) dir.create(bs_dir, recursive = TRUE)
  
  bash_script <- file.path(bs_dir, script_name)
  cat("#!/bin/bash\n", file = bash_script)
  cat(paste0("#SBATCH -J ", cluster_job_name, "\n"), file = bash_script, append = TRUE)
  cat("#SBATCH --nodes=1\n", file = bash_script, append = TRUE)
  cat("#SBATCH --ntasks-per-node=1\n", file = bash_script, append = TRUE)
  cat(paste0("#SBATCH --cpus-per-task=", threads, "\n"), file = bash_script, append = TRUE)
  cat(paste0("#SBATCH --mem=", cluster_mem, "\n"), file = bash_script, append = TRUE)
  cat(paste0("#SBATCH -t ", cluster_time, "\n"), file = bash_script, append = TRUE)
  cat(paste0("#SBATCH -p ", cluster_partition, "\n"), file = bash_script, append = TRUE)
  cat(paste0("#SBATCH -q ", cluster_queue, "\n"), file = bash_script, append = TRUE)
  cat("#SBATCH --mail-type=ALL\n", file = bash_script, append = TRUE)
  cat("#\n\n", file = bash_script, append = TRUE)
  
  if (!is.null(load_module) && load_module != "") {
    cat(paste0("module load ", load_module, "\n\n"), file = bash_script, append = TRUE)
  }
  
  for (i in 1:max_reps) {
    seed_i <- base_seed + i * seed_step
    rep_dir <- file.path(bs_dir, sprintf("bs_rep_%02d", i))
    
    prefix_i <- file.path(rep_dir, sprintf("cactus_bs_rep_%02d", i))
    
    cmd <- paste0(
      "echo \"Running bootstrap replica ", i, " (", bs_per_rep, " BS, seed ", seed_i, ")\"\n",
      "mkdir -p ", shQuote(rep_dir), "\n",
      raxml_exec, " --bootstrap ",
      "--msa ", shQuote(alignment_file), " ",
      "--model ", shQuote(partition_file), " ",
      "--tree-constraint ", shQuote(constraint_file), " "
    )
    if (!is.null(outgroup) && outgroup != "") {
      cmd <- paste0(cmd, "--outgroup ", shQuote(outgroup), " ")
    }
    cmd <- paste0(cmd,
      "--bs-trees ", bs_per_rep, " ",
      "--bs-metric tbe ",
      "--threads ", threads, " ",
      "--workers ", workers, " ",
      "--seed ", seed_i, " ",
      "--prefix ", shQuote(prefix_i), "\n\n"
    )
    cat(cmd, file = bash_script, append = TRUE)
  }
  
  system(paste("chmod +x", shQuote(bash_script)))
  message("Bash script ready: ", bash_script)
  return(bash_script)
}

#' Collect and Concatenate Parallel Bootstrap Tree Outputs
#'
#' Scans a target directory for chunked bootstrap output files (`.raxml.bootstraps`) and concatenates them into a single Newick tree file.
#'
#' @param bs_dir Character. Directory path containing chunked bootstrap output files.
#' @param output_dir Character. Output directory path to save concatenated bootstrap file. Defaults to `bs_dir`.
#' @param prefix Character. File output prefix. Defaults to `"cactus_ALL_bootstraps"`.
#' @return Character path to the concatenated bootstrap tree file (`.tree`).
#' @export
collect_bootstraps <- function(bs_dir, output_dir = bs_dir, prefix = "cactus_ALL_bootstraps") {
  
  bs_files <- list.files(bs_dir, pattern = "\\.raxml\\.bootstraps$", recursive = TRUE, full.names = TRUE)
  if (length(bs_files) == 0) {
    stop("No BS files found in ", bs_dir)
  }
  
  message("Found ", length(bs_files), " BS files")
  
  all_bs_trees <- unlist(lapply(bs_files, ape::read.tree), recursive = FALSE)
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  bs_concat_file <- file.path(output_dir, paste0(prefix, ".tree"))
  ape::write.tree(all_bs_trees, file = bs_concat_file)
  
  message("BS trees (", length(all_bs_trees), " total) concatenated to: ", bs_concat_file)
  return(bs_concat_file)
}

#' Check Bootstrap Convergence Criterion in RAxML-NG
#'
#' Evaluates whether the generated pool of non-parametric bootstrap trees has achieved statistical convergence (`RAxML-NG`).
#' Convergence is assessed using the MRE-based cutoff criterion (typically <= 0.03), ensuring that a sufficient number of bootstrap replicates were sampled.
#'
#' @param raxml_bin_path Character. System command or full path to executable `RAxML-NG` binary.
#' @param bs_trees_file Character. Path to concatenated bootstrap tree file.
#' @param bs_cutoff Numeric. Permutation cutoff threshold for convergence. Defaults to `0.03`.
#' @param seed Integer. Random seed for reproducible convergence testing. Defaults to `111`.
#' @param threads Integer. Number of CPU threads. Defaults to `4`.
#' @param output_dir Character. Directory path to save convergence report logs. Defaults to `dirname(bs_trees_file)`.
#' @param prefix Character. Output file prefix. Defaults to `"cactus_bs_convergence"`.
#' @return Character path to the convergence log file (`.raxml.log`).
#' @export
check_bs_convergence <- function(raxml_bin_path, bs_trees_file, bs_cutoff = 0.03,
                                 seed = 111, threads = 4, 
                                 output_dir = dirname(bs_trees_file), prefix = "cactus_bs_convergence") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_prefix <- file.path(output_dir, prefix)
  
  args <- c(
    "--bsconverge",
    "--bs-trees", shQuote(bs_trees_file),
    "--bs-metric", "tbe",
    "--bs-cutoff", as.character(bs_cutoff),
    "--seed", as.character(seed),
    "--threads", as.character(threads),
    "--prefix", shQuote(out_prefix)
  )
  
  system2(command = raxml_bin_path, args = args)
  return(paste0(out_prefix, ".raxml.log"))
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
#' @param outgroup Character. Optional outgroup taxon binomial. Defaults to `NULL`.
#' @param seed Integer. Random seed for reproducible temporal bootstrap initialization. Defaults to `111`.
#' @param threads Integer. Number of CPU threads. Defaults to `4`.
#' @param output_dir Character. Directory path to save output temporal bootstrap trees. Defaults to `dirname(aln_file)`.
#' @param prefix Character. Output file prefix. Defaults to `"cactus_temporal"`.
#' @return Character path to the resulting temporal bootstrap trees file (`.raxml.bootstraps`).
#' @export
calculate_temporal_bootstraps <- function(raxml_bin_path, aln_file, part_file, best_tree_file,
                                          bs_trees = 500, outgroup = NULL, seed = 111,
                                          threads = 4, output_dir = dirname(aln_file),
                                          prefix = "cactus_temporal") {
                                          
  if (Sys.which(raxml_bin_path) == "") {
    stop("Executable '", raxml_bin_path, "' not found in your system's PATH.\n",
         "Please ensure RAxML-NG is installed and available, or provide the full absolute path.")
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_prefix <- file.path(output_dir, prefix)
  
  args <- c(
    "--bootstrap",
    "--msa", shQuote(aln_file),
    "--model", shQuote(part_file),
    "--tree-constraint", shQuote(best_tree_file),
    "--bs-trees", as.character(bs_trees),
    "--bs-metric", "tbe",
    "--seed", as.character(seed),
    "--threads", as.character(threads),
    "--prefix", shQuote(out_prefix)
  )
  
  if (!is.null(outgroup) && outgroup != "") {
    args <- c(args, "--outgroup", shQuote(outgroup))
  }
  
  system2(command = raxml_bin_path, args = args)
  return(paste0(out_prefix, ".raxml.bootstraps"))
}
                                          
