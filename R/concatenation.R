#' Concatenate Locus Alignments and Build Partition Coordinate Maps
#'
#' Concatenates individual orthologous locus alignments end-to-end into a unified multilocus supermatrix.
#' Combining independent molecular loci increases statistical power to resolve difficult ancestral nodes while
#' allowing partitioned substitution modeling to account for mutational rate heterogeneity across molecular locus alignments.
#' Exports one coordinate map per partition format rather than one per downstream program: a
#' RAxML-style map (`PARTITION_raxml_style.txt`) read by `RAxML-NG`, `ModelTest-NG` and `IQ-TREE`,
#' and a NEXUS SETS block (`PARTITION_nexus_charset.nex`). Configuration scripts for
#' `PartitionFinder2` and `MrBayes` are exported separately under their program names.
#'
#' @param input_dir Character. Path to directory containing curated, aligned locus FASTA files.
#' @param output_dir Character. Path to destination root directory for concatenated alignments and partition maps.
#' @param outgroup_pattern Character, named character vector, or `NULL`. Regular expression identifying outgroup terminals, or a named character vector of expressions (e.g. `c(Anacampserotaceae = "...", Portulacaceae = "...", Talinaceae = "...")`) to additionally compute per-group coverage columns. When supplied, [report_marker_group_coverage()] runs before concatenation and its table is written to `logs_and_qc/SUPP_TABLE_marker_group_coverage.csv`. Defaults to `NULL`.
#' @param min_coverage Numeric. Fraction of non-gap, non-missing characters at which a terminal counts as covered by a marker, passed to [report_marker_group_coverage()]. Defaults to `0.2`.
#' @param exclude_markers Character vector or `NULL`. Markers to leave out of the supermatrix, named
#'   as they appear in the alignment file names without the `Masked_` prefix and without the
#'   extension. The files are not modified or removed, so an excluded locus remains available for
#'   other analyses; only this supermatrix is built without it. The excluded markers are named in
#'   the run log and listed in `logs_and_qc/TABLE_markers_excluded.csv`. A name that matches no
#'   alignment raises a warning rather than failing, so a typo is visible instead of silent.
#'   Defaults to `NULL`.
#' @return A data frame containing supermatrix dimensions, taxon coverage, and locus partition bounds.
#' @examples
#' \dontrun{
#' run_concatenation_pipeline(
#'   input_dir = "5_MAFFT_Cleaned/aligned_markers",
#'   output_dir = "6_Concatenated"
#' )
#'
#' # A dating matrix without the nuclear locus that carries no outgroup coverage.
#' run_concatenation_pipeline(
#'   input_dir = "5_MAFFT_Cleaned/aligned_markers",
#'   output_dir = "6_Concatenated_dating",
#'   exclude_markers = "phyC"
#' )
#' }
#' @export
run_concatenation_pipeline <- function(input_dir, output_dir, outgroup_pattern = NULL, min_coverage = 0.2,
                                       exclude_markers = NULL) {
  ind_dir     <- file.path(output_dir, "individual_markers")
  concat_dir  <- file.path(output_dir, "concatenated_alignments")
  logs_dir    <- file.path(output_dir, "logs_and_qc")
  tables_dir  <- file.path(output_dir, "final_tables")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ind_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(concat_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  file_super_fasta      <- file.path(concat_dir, "ALIGNMENT_supermatrix.fasta")
  file_dup_groups       <- file.path(logs_dir, "SUPP_TABLE_identical_sequence_groups.csv")
  file_super_phylip     <- file.path(concat_dir, "ALIGNMENT_supermatrix.phy")
  file_super_nexus      <- file.path(concat_dir, "ALIGNMENT_supermatrix.nex")
  
  # Partition exports are named by the format they encode, not by the program expected to read
  # them, except where the file is a configuration script belonging to one specific program.
  # PARTITION_raxml_style.txt holds the single coordinate map ("RAxML-style" is the established
  # name of that format) and is read by RAxML-NG, ModelTest-NG, and IQ-TREE via -q / -p / -Q.
  # PARTITION_nexus_charset.nex holds the same coordinates as a portable NEXUS SETS block.
  # The PartitionFinder2 .cfg and the MrBayes .nex are command scripts, so a program name is
  # the correct identifier for those two.
  file_part_raxml_style <- file.path(concat_dir, "PARTITION_raxml_style.txt")
  file_part_charset     <- file.path(concat_dir, "PARTITION_nexus_charset.nex")
  file_part_pf2         <- file.path(concat_dir, "PARTITION_partitionfinder.cfg")
  file_part_mrbayes     <- file.path(concat_dir, "PARTITION_mrbayes.nex")
  
  file_log              <- file.path(logs_dir, "LOG_concatenation_run.txt")
  file_commands         <- file.path(logs_dir, "LOG_example_commands.txt")
  
  file_name_crosswalk   <- file.path(logs_dir, "SUPP_TABLE_taxon_name_crosswalk.csv")
  
  file_super_stats      <- file.path(tables_dir, "TABLE_supermatrix_statistics.csv")
  file_final_matrix     <- file.path(tables_dir, "TABLE_final_species_alignment_summary.csv")
  file_marker_stats     <- file.path(tables_dir, "TABLE_marker_statistics.csv")
  file_marker_ranges    <- file.path(tables_dir, "TABLE_marker_ranges.tsv")

  log_message <- function(...) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
    cat(msg, "\n")
    write(msg, file = file_log, append = TRUE)
  }

  # A marker whose coverage sits entirely on one side of the root contributes columns that the
  # other side cannot share, and concatenation then assigns branch length where no character is
  # held in common. Reported, never blocking: whether such a marker belongs in the matrix is a
  # decision about the analysis, not something a concatenation routine should settle.
  if (!is.null(outgroup_pattern)) {
    report_marker_group_coverage(
      input_dir       = input_dir,
      outgroup_pattern = outgroup_pattern,
      min_coverage    = min_coverage,
      out_csv         = file.path(logs_dir, "SUPP_TABLE_marker_group_coverage.csv")
    )
  }

  clean_name <- function(x) {
    x <- sapply(strsplit(x, "\\|"), `[`, 1)
    x <- tolower(x); x <- gsub("[^a-z0-9 ]", " ", x); x <- gsub("\\s+", " ", x); x <- trimws(x)
    x <- gsub(" ", "_", x); gsub("^([a-z])", "\\U\\1", x, perl = TRUE)
  }

  extract_sid <- function(x) {
    sapply(strsplit(x, "\\|"), function(p) {
      if (length(p) >= 2) return(trimws(p[length(p)]))
      # fallback: try last piece separated by _ if it contains numbers
      parts <- unlist(strsplit(trimws(p[1]), "_"))
      if (length(parts) > 1 && grepl("[0-9]", parts[length(parts)])) return(parts[length(parts)])
      return("-")
    })
  }

  alignment_to_char <- function(aln, marker_name) {
    aln_char <- as.character(aln)
    split_list <- strsplit(aln_char, "", fixed = TRUE)
    m <- do.call(rbind, split_list)
    if (is.null(dim(m))) m <- matrix(m, nrow = 1)
    rownames(m) <- names(aln)
    storage.mode(m) <- "character"
    toupper(m)
  }

  compute_marker_stats <- function(aln_char, marker_name) {
    cells <- as.vector(aln_char)
    missing_pct <- mean(cells %in% c("N", "n")); gaps_pct <- mean(cells == "-")
    bases <- toupper(cells[cells %in% c("A", "C", "G", "T", "a", "c", "g", "t")])
    gc_pct <- if (length(bases) == 0) NA_real_ else sum(bases %in% c("G", "C")) / length(bases) * 100
    var_sites <- sum(apply(aln_char, 2, function(col) { col2 <- toupper(col[!toupper(col) %in% c("-", "N")]); length(unique(col2)) > 1 }))
    pi_sites <- sum(apply(aln_char, 2, function(col) { col2 <- toupper(col[!toupper(col) %in% c("-", "N")]); if (length(col2) == 0) return(FALSE); sum(table(col2) >= 2) >= 2 }))
    inv_sites <- sum(apply(aln_char, 2, function(col) { col2 <- toupper(col[!toupper(col) %in% c("-", "N")]); if (length(col2) == 0) return(FALSE); length(unique(col2)) == 1 }))
    dplyr::tibble(marker = marker_name, n_taxa = nrow(aln_char), alignment_length = ncol(aln_char), missing_pct = round(missing_pct * 100, 3), gaps_pct = round(gaps_pct * 100, 3), gc_pct = round(gc_pct, 3), variable_sites = var_sites, parsimony_informative = pi_sites, invariant_sites = inv_sites)
  }

  log_message("Starting concatenation pipeline.")
  log_message("Input dir: ", input_dir)
  log_message("Output dir: ", output_dir)

  fasta_files <- sort(list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE))
  marker_raw <- tools::file_path_sans_ext(basename(fasta_files))
  marker_clean <- sub("^Masked_", "", marker_raw)

  # A locus can be worth keeping for one purpose and disqualifying for another. phyC resolves
  # Cactaceae for barcoding and carries no Portulacaceae, so in a dating matrix it adds columns
  # that only one side of the root can occupy. The alignment stays in `input_dir`; the exclusion
  # is declared at the call site, named in the log and written to disk, so the composition of the
  # supermatrix is readable from the run alone rather than from whichever files happened to be
  # in the folder.
  if (!is.null(exclude_markers) && length(exclude_markers) > 0L) {
    exclude_markers <- as.character(exclude_markers)
    drop <- marker_clean %in% exclude_markers
    missing <- setdiff(exclude_markers, marker_clean)
    if (length(missing) > 0L) {
      warning("`exclude_markers` names markers absent from `input_dir`: ",
              paste(missing, collapse = ", "), call. = FALSE)
    }
    if (any(drop)) {
      log_message("Markers excluded by request: ", paste(marker_clean[drop], collapse = ", "))
      readr::write_csv(
        dplyr::tibble(marker = marker_clean[drop], file = basename(fasta_files[drop]),
                      reason = "excluded_by_exclude_markers"),
        file.path(logs_dir, "TABLE_markers_excluded.csv")
      )
      fasta_files  <- fasta_files[!drop]
      marker_raw   <- marker_raw[!drop]
      marker_clean <- marker_clean[!drop]
    }
    if (length(fasta_files) == 0L) {
      stop("`exclude_markers` removed every alignment in `input_dir`; nothing left to concatenate.",
           call. = FALSE)
    }
  }

  all_alignments <- list()
  crosswalk_list <- list()
  for (i in seq_along(fasta_files)) {
    aln <- Biostrings::readDNAStringSet(fasta_files[i])
    original_names <- names(aln)
    cleaned_names <- clean_name(original_names)
    names(aln) <- cleaned_names
    all_alignments[[marker_clean[i]]] <- aln
    crosswalk_list[[i]] <- dplyr::tibble(marker = marker_clean[i], original_name = original_names, cleaned_name = cleaned_names)
  }
  crosswalk_df <- dplyr::bind_rows(crosswalk_list)
  readr::write_csv(crosswalk_df, file_name_crosswalk)


  # Build the comprehensive species-marker alignment summary
  registry_file_tables <- file.path(dirname(dirname(input_dir)), "4_Cleaned", "tables", "TABLE_sequence_registry_with_acceptance.csv")
  registry_file_root   <- file.path(dirname(dirname(input_dir)), "4_Cleaned", "TABLE_sequence_registry_with_acceptance.csv")
  registry_file <- if (file.exists(registry_file_tables)) registry_file_tables else registry_file_root
  if (file.exists(registry_file)) {
    registry_df <- readr::read_csv(registry_file, show_col_types = FALSE)
    status_map <- dplyr::distinct(registry_df, species, source_branch, accepted_or_outgroup) |>
      dplyr::mutate(species_clean = clean_name(species)) |>
      dplyr::distinct(species_clean, .keep_all = TRUE)
      
    sid_map <- registry_df |>
      dplyr::filter(!is.na(species), !is.na(marker_key), !is.na(sid)) |>
      dplyr::mutate(
        sid = as.character(sid), 
        fasta_name = as.character(fasta_name), 
        accepted_rank = dplyr::if_else(accepted_or_outgroup, 1L, 0L), 
        outgroup_rank = dplyr::if_else(is_outgroup_species, 1L, 0L)
      ) |>
      dplyr::arrange(marker_key, species, dplyr::desc(accepted_rank), dplyr::desc(outgroup_rank), sid, fasta_name) |>
      dplyr::distinct(species, marker_key, .keep_all = TRUE) |>
      dplyr::filter(accepted_or_outgroup) |>
      dplyr::mutate(species_clean = clean_name(species)) |>
      dplyr::select(species_clean, marker_key, sid)
      
    matrix_long <- crosswalk_df |>
      dplyr::left_join(sid_map, by = c("cleaned_name" = "species_clean", "marker" = "marker_key")) |>
      dplyr::mutate(
        sid_extracted = extract_sid(original_name),
        sid = dplyr::case_when(
          !is.na(sid) ~ as.character(sid),
          sid_extracted != "-" ~ sid_extracted,
          TRUE ~ "-"
        )
      ) |>
      dplyr::select(-sid_extracted)
  } else {
    status_map <- dplyr::tibble(species = character(), source_branch = character(), accepted_or_outgroup = logical())
    matrix_long <- crosswalk_df |>
      dplyr::mutate(sid = extract_sid(original_name))
  }
  
  species_summary <- matrix_long |>
    dplyr::group_by(cleaned_name) |>
    dplyr::summarise(
      retained_markers = dplyr::n(),
      pct_markers = round((dplyr::n() / length(marker_clean)) * 100, 2),
      .groups = "drop"
    )
  
  matrix_wide <- matrix_long |>
    dplyr::select(cleaned_name, marker, sid) |>
    dplyr::distinct(cleaned_name, marker, .keep_all = TRUE) |>
    tidyr::pivot_wider(names_from = marker, values_from = sid) |>
    dplyr::mutate(dplyr::across(-cleaned_name, ~ tidyr::replace_na(as.character(.), "-")))
  
  final_summary_matrix <- species_summary |>
    dplyr::left_join(matrix_wide, by = "cleaned_name") |>
    dplyr::rename(species = cleaned_name)
    
  if (nrow(status_map) > 0) {
    status_map <- status_map |> dplyr::rename(species_class = source_branch) |> dplyr::select(-species) |> dplyr::rename(species = species_clean)
    final_summary_matrix <- final_summary_matrix |>
      dplyr::left_join(status_map, by = "species") |>
      dplyr::relocate(species_class, accepted_or_outgroup, .after = species)
  }
  
  readr::write_csv(final_summary_matrix, file_final_matrix)
  
  taxa_all <- sort(unique(unlist(lapply(all_alignments, names))))
  concat_list <- list(); stats_list <- list(); part_comp_list <- list(); family_presence_list <- list(); marker_ranges <- dplyr::tibble(marker = character(), start = integer(), end = integer())
  current_pos <- 1L
  
  if (nrow(status_map) > 0) {
    ingroup_taxa_vec <- status_map |> dplyr::filter(species_class == "ingroup") |> dplyr::pull(species)
    outgroup_taxa_vec <- status_map |> dplyr::filter(species_class == "outgroup") |> dplyr::pull(species)
  } else {
    ingroup_taxa_vec <- taxa_all
    outgroup_taxa_vec <- character(0)
  }

  # Family split within the outgroup. Portulaca is the sole Portulacaceae genus mined by
  # assemble_outgroup_phylotar(); Anacampseros/Grahamia/Talinopsis are Anacampserotaceae.
  # These are computed on the REAL final taxon names (post Stage 5 realignment), unlike the
  # equivalent split in integrate_and_clean_markers() (Stage 4), which cannot see which
  # individual outgroup taxa Stage 5's joint MAFFT + DECIPHER realignment later drops.
  # The trailing underscore is required: without it "^Portulaca" also matches Portulacaria
  # (Didiereaceae), which is neither Portulacaceae nor part of the intended rooting sample.
  anacampserotaceae_taxa_vec <- outgroup_taxa_vec[grepl("^(Anacampseros|Grahamia|Talinopsis)_", outgroup_taxa_vec)]
  portulacaceae_taxa_vec <- outgroup_taxa_vec[grepl("^Portulaca_", outgroup_taxa_vec)]
  talinaceae_taxa_vec <- outgroup_taxa_vec[grepl("^(Talinum|Talinella)_", outgroup_taxa_vec)]

  for (marker in marker_clean) {
    aln <- all_alignments[[marker]]
    aln_char <- alignment_to_char(aln, marker)
    aln_len <- ncol(aln_char)
    stats_list[[marker]] <- compute_marker_stats(aln_char, marker)
    
    aln_full <- matrix("-", nrow = length(taxa_all), ncol = aln_len, dimnames = list(taxa_all, NULL))
    common_taxa <- intersect(rownames(aln_char), taxa_all)
    aln_full[common_taxa, ] <- aln_char[common_taxa, , drop = FALSE]
    
    concat_list[[marker]] <- aln_full
    marker_ranges <- dplyr::add_row(marker_ranges, marker = marker, start = current_pos, end = current_pos + aln_len - 1L)
    current_pos <- current_pos + aln_len
    
    # Compute partition comparison
    aln_in <- aln_full[rownames(aln_full) %in% ingroup_taxa_vec, , drop = FALSE]
    aln_out <- aln_full[rownames(aln_full) %in% outgroup_taxa_vec, , drop = FALSE]
    aln_anacamp <- aln_full[rownames(aln_full) %in% anacampserotaceae_taxa_vec, , drop = FALSE]
    aln_portulaca <- aln_full[rownames(aln_full) %in% portulacaceae_taxa_vec, , drop = FALSE]
    aln_talinum <- aln_full[rownames(aln_full) %in% talinaceae_taxa_vec, , drop = FALSE]
    taxa_in_present <- if (nrow(aln_in) > 0) sum(apply(aln_in, 1, function(r) !all(r %in% c("-", "N")))) else 0L
    taxa_out_present <- if (nrow(aln_out) > 0) sum(apply(aln_out, 1, function(r) !all(r %in% c("-", "N")))) else 0L
    taxa_anacamp_present <- if (nrow(aln_anacamp) > 0) sum(apply(aln_anacamp, 1, function(r) !all(r %in% c("-", "N")))) else 0L
    taxa_portulaca_present <- if (nrow(aln_portulaca) > 0) sum(apply(aln_portulaca, 1, function(r) !all(r %in% c("-", "N")))) else 0L
    taxa_talinum_present <- if (nrow(aln_talinum) > 0) sum(apply(aln_talinum, 1, function(r) !all(r %in% c("-", "N")))) else 0L
    
    var_in <- if (nrow(aln_in) > 0) sum(apply(aln_in, 2, function(col) { c_clean <- toupper(col[!toupper(col) %in% c("-", "N")]); length(unique(c_clean)) > 1 })) else 0L
    pi_in  <- if (nrow(aln_in) > 0) sum(apply(aln_in, 2, function(col) { c_clean <- toupper(col[!toupper(col) %in% c("-", "N")]); if (length(c_clean) == 0) return(FALSE); sum(table(c_clean) >= 2) >= 2 })) else 0L
    var_tot <- sum(apply(aln_full, 2, function(col) { c_clean <- toupper(col[!toupper(col) %in% c("-", "N")]); length(unique(c_clean)) > 1 }))
    pi_tot  <- sum(apply(aln_full, 2, function(col) { c_clean <- toupper(col[!toupper(col) %in% c("-", "N")]); if (length(c_clean) == 0) return(FALSE); sum(table(c_clean) >= 2) >= 2 }))
    
    part_comp_list[[marker]] <- dplyr::tibble(
      marker = marker,
      alignment_length = aln_len,
      n_taxa_total = nrow(aln_full),
      n_taxa_ingroup_present = taxa_in_present,
      n_taxa_outgroup_present = taxa_out_present,
      variable_sites_ingroup = var_in,
      variable_sites_total = var_tot,
      parsimony_informative_ingroup = pi_in,
      parsimony_informative_total = pi_tot,
      pi_inflation_ratio = dplyr::if_else(pi_in > 0, round(pi_tot / pi_in, 2), NA_real_)
    )

    # Kept separately (not written to TABLE_partition_informativeness_comparison.csv) -- used
    # only internally below to patch the "Final loci retained" rows in
    # TABLE_dataset_species_summary.csv (Stage 4) with true post-realignment counts.
    family_presence_list[[marker]] <- dplyr::tibble(
      marker = marker,
      anacamp_present = taxa_anacamp_present > 0,
      portulaca_present = taxa_portulaca_present > 0,
      talinum_present = taxa_talinum_present > 0
    )
  }

  readr::write_csv(dplyr::bind_rows(stats_list), file_marker_stats)
  part_comp_df <- dplyr::bind_rows(part_comp_list)
  readr::write_csv(part_comp_df, file.path(tables_dir, "TABLE_partition_informativeness_comparison.csv"))
  utils::write.table(marker_ranges, file = file_marker_ranges, sep = "\t", quote = FALSE, row.names = FALSE)

  log_message("Markers concatenated: ", nrow(marker_ranges))
  log_message("Partitions retaining outgroup: ", sum(part_comp_df$n_taxa_outgroup_present > 0),
              " of ", nrow(part_comp_df))

  # --- Patch the "Final loci retained" rows in TABLE_dataset_species_summary.csv (Stage 4) in place. ---
  # integrate_and_clean_markers() (Stage 4) writes those rows as provisional, pre-realignment
  # counts, because it runs before Stage 5's joint MAFFT + DECIPHER realignment -- which can drop
  # individual taxa from a locus (e.g. a handful of divergent outgroup accessions masked down to
  # nothing) even when the marker file itself survives. Now that the true final, realigned
  # supermatrix exists, overwrite those same rows with the correct counts. A locus counts as
  # retained for a family only if at least one taxon of that family has real (non-gap/non-N)
  # sequence in the final concatenated partition.
  family_presence_df <- dplyr::bind_rows(family_presence_list)
  n_final_loci_cactaceae <- sum(part_comp_df$n_taxa_ingroup_present > 0)
  n_final_loci_anacampserotaceae <- if (nrow(family_presence_df) > 0) sum(family_presence_df$anacamp_present) else 0L
  n_final_loci_portulacaceae <- if (nrow(family_presence_df) > 0) sum(family_presence_df$portulaca_present) else 0L
  n_final_loci_talinaceae <- if (nrow(family_presence_df) > 0) sum(family_presence_df$talinum_present) else 0L
  n_final_loci_total <- nrow(part_comp_df)

  # Species counts must be patched for the same reason as the locus counts: Stage 4 writes them
  # before run_joint_realignment() (Stage 5) drops individual terminals, so its values overcount
  # the outgroup. taxa_all is exactly the set of terminals present in the final supermatrix.
  n_final_sp_anacampserotaceae <- length(intersect(taxa_all, anacampserotaceae_taxa_vec))
  n_final_sp_portulacaceae     <- length(intersect(taxa_all, portulacaceae_taxa_vec))
  n_final_sp_talinaceae        <- length(intersect(taxa_all, talinaceae_taxa_vec))
  n_final_sp_outgroup          <- length(intersect(taxa_all, outgroup_taxa_vec))
  n_final_sp_joint             <- length(taxa_all)

  species_summary_file_tables <- file.path(dirname(dirname(input_dir)), "4_Cleaned", "tables", "TABLE_dataset_species_summary.csv")
  species_summary_file_root   <- file.path(dirname(dirname(input_dir)), "4_Cleaned", "TABLE_dataset_species_summary.csv")
  species_summary_file <- if (file.exists(species_summary_file_tables)) species_summary_file_tables else species_summary_file_root

  if (file.exists(species_summary_file)) {
    species_summary_df <- readr::read_csv(species_summary_file, show_col_types = FALSE)
    final_values <- dplyr::tribble(
      ~metric, ~value, ~details,
      "Final loci retained (Cactaceae Ingroup)", as.character(n_final_loci_cactaceae), "Loci with real (non-gap) sequence for Cactaceae in the final, realigned supermatrix (Stage 6)",
      "Final loci retained (Anacampserotaceae Outgroup)", as.character(n_final_loci_anacampserotaceae), "Loci with real (non-gap) sequence for Anacampserotaceae in the final, realigned supermatrix (Stage 6)",
      "Final loci retained (Portulacaceae Outgroup)", as.character(n_final_loci_portulacaceae), "Loci with real (non-gap) sequence for Portulacaceae in the final, realigned supermatrix (Stage 6)",
      "Final loci retained (Talinaceae Outgroup)", as.character(n_final_loci_talinaceae), "Loci with real (non-gap) sequence for Talinaceae in the final, realigned supermatrix (Stage 6)",
      "Total unique loci in final joint dataset", as.character(n_final_loci_total), "Total partitions in the final concatenated supermatrix (Stage 6)",
      "Unique Anacampserotaceae outgroup species retained", as.character(n_final_sp_anacampserotaceae), "Anacampserotaceae terminals present in the final, realigned supermatrix (Stage 6)",
      "Unique Portulacaceae outgroup species retained", as.character(n_final_sp_portulacaceae), "Portulaca terminals present in the final, realigned supermatrix (Stage 6)",
      "Unique Talinaceae outgroup species retained", as.character(n_final_sp_talinaceae), "Talinaceae terminals present in the final, realigned supermatrix (Stage 6)",
      "Total unique outgroup species retained", as.character(n_final_sp_outgroup), "Total outgroup terminals in the final, realigned supermatrix (Anacampserotaceae + Portulacaceae + Talinaceae, Stage 6)",
      "Total unique species in final dataset (Joint)", as.character(n_final_sp_joint), "Total terminals in the final concatenated supermatrix (Stage 6)"
    )
    for (i in seq_len(nrow(final_values))) {
      row_match <- species_summary_df$metric == final_values$metric[i]
      if (any(row_match)) {
        species_summary_df$value[row_match] <- final_values$value[i]
        species_summary_df$details[row_match] <- final_values$details[i]
      } else {
        species_summary_df <- dplyr::bind_rows(species_summary_df, final_values[i, ])
      }
    }
    readr::write_csv(species_summary_df, species_summary_file)
  } else {
    warning(sprintf("Could not find %s to update with final locus counts; run integrate_and_clean_markers() (Stage 4) first.", species_summary_file), call. = FALSE)
  }

  supermatrix <- do.call(cbind, concat_list)
  cells <- as.vector(supermatrix)
  bases <- toupper(cells[cells %in% c("A", "C", "G", "T", "a", "c", "g", "t")])

  supermatrix_stats <- dplyr::tibble(
    n_taxa = nrow(supermatrix), total_alignment_length = ncol(supermatrix),
    missing_pct = round(mean(toupper(cells) == "N") * 100, 3), gaps_pct = round(mean(cells == "-") * 100, 3),
    gc_pct = round(if (length(bases) == 0) NA_real_ else sum(bases %in% c("G", "C")) / length(bases) * 100, 3),
    variable_sites = sum(apply(supermatrix, 2, function(col) { col2 <- toupper(col[!toupper(col) %in% c("-", "N")]); length(unique(col2)) > 1 })),
    parsimony_informative = sum(apply(supermatrix, 2, function(col) { col2 <- toupper(col[!toupper(col) %in% c("-", "N")]); if (length(col2) == 0) return(FALSE); sum(table(col2) >= 2) >= 2 })),
    invariant_sites = sum(apply(supermatrix, 2, function(col) { col2 <- toupper(col[!toupper(col) %in% c("-", "N")]); length(col2) > 0 && length(unique(col2)) == 1 }))
  )
  readr::write_csv(supermatrix_stats, file_super_stats)

  super_seqs <- apply(supermatrix, 1, paste0, collapse = "")
  super_list <- as.list(super_seqs)
  names(super_list) <- rownames(supermatrix)

  seqinr::write.fasta(sequences = super_list, names = names(super_list), file.out = file_super_fasta)

  # PHYLIP
  taxa <- rownames(supermatrix)
  header <- paste(length(taxa), nchar(super_seqs[1]))
  writeLines(c(header, paste(taxa, super_seqs)), file_super_phylip)

  super_list_nexus <- lapply(super_list, function(s) strsplit(s, "")[[1]])
  ape::write.nexus.data(x = super_list_nexus, file = file_super_nexus, format = "dna", interleave = FALSE)

  # One coordinate map, written once. The datatype token stays "DNA" so that the downstream model
  # selection step decides the substitution model, rather than inheriting a model asserted here.
  part_lines <- dplyr::pull(dplyr::mutate(marker_ranges, line = paste0("DNA, ", marker, " = ", start, "-", end)), line)
  writeLines(part_lines, con = file_part_raxml_style)

  # NEXUS SETS block carrying the same coordinates, without the command blocks that make the
  # MrBayes export specific to MrBayes.
  charset_lines <- c(
    "#NEXUS",
    "",
    "begin sets;",
    sprintf("   charset %s = %d-%d;", marker_ranges$marker, marker_ranges$start, marker_ranges$end),
    sprintf("   charpartition loci = %s;",
            paste(sprintf("%s:%s", marker_ranges$marker, marker_ranges$marker), collapse = ", ")),
    "end;"
  )
  writeLines(charset_lines, con = file_part_charset)

  pf2_text <- c("## PartitionFinder2 configuration file", "", "## ALIGNMENT FILE ##", paste0("alignment = ", basename(file_super_phylip), ";"), "", "## BRANCHLENGTHS: linked | unlinked ##", "branchlengths = linked;", "", "## MODELS OF EVOLUTION ##", "models = all;", "", "## MODEL SELECTION: AIC | AICc | BIC ##", "model_selection = AICc;", "", "## DATA BLOCKS ##", "[data_blocks]", paste0(marker_ranges$marker, " = ", marker_ranges$start, "-", marker_ranges$end, ";"), "", "## SCHEMES ##", "[schemes]", "search = greedy;")
  writeLines(pf2_text, con = file_part_pf2)

  cat("#nexus\nbegin mrbayes;\n", file = file_part_mrbayes)
  for (i in seq_len(nrow(marker_ranges))) cat(sprintf("   charset %s = %d-%d;\n", marker_ranges$marker[i], marker_ranges$start[i], marker_ranges$end[i]), file = file_part_mrbayes, append = TRUE)
  cat(sprintf("\n   partition combined = %d: %s;\n", nrow(marker_ranges), paste(marker_ranges$marker, collapse = ", ")), file = file_part_mrbayes, append = TRUE)
  cat("   set partition = combined;\nend;\n", file = file_part_mrbayes, append = TRUE)

  for (marker in marker_clean) seqinr::write.fasta(sequences = as.list(apply(concat_list[[marker]], 1, paste0, collapse = "")), names = rownames(concat_list[[marker]]), file.out = file.path(ind_dir, paste0(marker, ".fasta")))

  # Terminals sharing an identical concatenated sequence. RAxML-NG's --check / --parse step
  # collapses these into a single tip and writes a *.reduced.phy, so the matrix actually analysed
  # in Stage 7 has fewer terminals than the supermatrix exported here. Exported so the collapsed
  # taxa can be reported or grafted back rather than disappearing silently.
  dup_groups <- split(names(super_seqs), unname(super_seqs))
  dup_groups <- dup_groups[lengths(dup_groups) > 1]
  dup_table <- if (length(dup_groups) == 0) {
    dplyr::tibble(group = integer(), taxon = character())
  } else {
    dplyr::tibble(
      group = rep(seq_along(dup_groups), lengths(dup_groups)),
      taxon = unlist(dup_groups, use.names = FALSE)
    )
  }
  readr::write_csv(dup_table, file_dup_groups)

  n_dup_collapsed <- sum(lengths(dup_groups)) - length(dup_groups)
  log_message("Supermatrix: ", nrow(supermatrix), " taxa x ", ncol(supermatrix), " bp")
  log_message("Identical-sequence groups: ", length(dup_groups),
              "; terminals RAxML-NG will collapse: ", n_dup_collapsed,
              " (effective matrix: ", nrow(supermatrix) - n_dup_collapsed, " taxa)")
  log_message(paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"))

  cmds <- c(
    "# Example commands - adjust paths and threads as needed",
    "",
    "## IQ-TREE2",
    paste0("iqtree2 -s ", basename(file_super_phylip), " -p ", basename(file_part_raxml_style), " -m MFP+MERGE -B 1000 -T 4"),
    "",
    "## RAxML-NG",
    paste0("raxml-ng --all --msa ", basename(file_super_phylip), " --model ", basename(file_part_raxml_style), " --bs-trees 100 --threads 4"),
    "",
    "## ModelTest-NG",
    paste0("modeltest-ng --datatype nt --input ", basename(file_super_phylip), " --partitions ", basename(file_part_raxml_style), " --template raxml"),
    "",
    "## PartitionFinder2",
    paste0("partitionfinder2 ", basename(file_part_pf2)),
    "",
    "## MrBayes",
    paste0("mb ", basename(file_part_mrbayes))
  )
  writeLines(cmds, con = file_commands)

  message("Concatenation pipeline finalized! \U0001f335")
}


#' Report How Each Marker Covers the Ingroup and the Outgroup
#'
#' Counts, per locus, how many ingroup and how many outgroup terminals carry real sequence, and
#' warns when either side is empty.
#'
#' Concatenation assumes that the loci being joined describe the same terminals. A locus sampled
#' almost entirely on one side of the root breaks that assumption without breaking anything
#' visible: it adds columns the other side cannot share, and the branch lengths spanning the
#' bipartition are then estimated from the loci that remain. In the August 2026 supermatrix,
#' `phyC` covered 15 of 24 Anacampserotaceae terminals at 99.3% and no Portulacaceae terminal at
#' all; the branch subtending the outgroup fell from roughly 600 expected substitutions to 0.03,
#' and the two deepest calibrated nodes returned their own bounds rather than an estimate.
#'
#' The function reports and does not block. Whether a locus of that kind belongs in a given matrix
#' depends on which analysis the matrix is for: a nuclear locus with no outgroup coverage is
#' unusable for dating and valuable for species discrimination.
#'
#' @param input_dir Character. Directory of aligned locus FASTA files, one per marker.
#' @param outgroup_pattern Character. Regular expression matched against terminal names, or a named
#'   character vector of regular expressions (e.g. `c(Anacampserotaceae = "...", Portulacaceae = "...", Talinaceae = "...")`)
#'   to additionally compute per-group coverage columns.
#' @param min_coverage Numeric. Fraction of non-gap, non-missing characters at which a terminal
#'   counts as covered by that marker. Defaults to `0.2`.
#' @param out_csv Character or `NULL`. Path to write the table to. Defaults to `NULL`.
#' @return Invisibly, a data frame with one row per marker: alignment length, terminals and
#'   covered terminals on each side, and the median coverage of each side.
#' @examples
#' \dontrun{
#' report_marker_group_coverage(
#'   input_dir = "5_MAFFT_Cleaned/aligned_markers",
#'   outgroup_pattern = c(
#'     Anacampserotaceae = "^(Anacampseros|Grahamia|Talinopsis)_",
#'     Portulacaceae     = "^Portulaca_",
#'     Talinaceae        = "^(Talinum|Talinella)_"
#'   )
#' )
#' }
#' @export
report_marker_group_coverage <- function(input_dir, outgroup_pattern, min_coverage = 0.2,
                                         out_csv = NULL) {
  files <- list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE)
  if (length(files) == 0L) {
    stop("No FASTA files found in '", input_dir, "'.", call. = FALSE)
  }

  combined_pattern <- if (length(outgroup_pattern) > 1L) {
    paste0("(", paste(unname(outgroup_pattern), collapse = ")|("), ")")
  } else {
    outgroup_pattern
  }

  rows <- list()
  for (f in files) {
    marker <- tools::file_path_sans_ext(basename(f))
    aln <- Biostrings::readDNAStringSet(f)
    if (length(aln) == 0L) next
    txt <- as.character(aln)
    len <- nchar(txt[1])
    cov <- vapply(txt, function(x) {
      ch <- strsplit(x, "", fixed = TRUE)[[1]]
      sum(!ch %in% c("-", "?", "N", "n")) / length(ch)
    }, numeric(1))

    is_out <- grepl(combined_pattern, names(txt))
    row_df <- data.frame(
      marker = marker, aln_len = len,
      n_ingroup = sum(!is_out), n_ingroup_covered = sum(cov[!is_out] >= min_coverage),
      n_outgroup = sum(is_out), n_outgroup_covered = sum(cov[is_out] >= min_coverage),
      median_cov_ingroup = if (any(!is_out)) stats::median(cov[!is_out]) else NA_real_,
      median_cov_outgroup = if (any(is_out)) stats::median(cov[is_out]) else NA_real_,
      stringsAsFactors = FALSE
    )
    if (!is.null(names(outgroup_pattern)) && all(nzchar(names(outgroup_pattern)))) {
      for (grp_nm in names(outgroup_pattern)) {
        grp_pat <- outgroup_pattern[[grp_nm]]
        is_grp <- grepl(grp_pat, names(txt))
        row_df[[paste0("n_", grp_nm)]] <- sum(is_grp)
        row_df[[paste0("n_", grp_nm, "_covered")]] <- sum(cov[is_grp] >= min_coverage)
        row_df[[paste0("median_cov_", grp_nm)]] <- if (any(is_grp)) stats::median(cov[is_grp]) else NA_real_
      }
    }
    rows[[marker]] <- row_df
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  empty_out <- out$marker[out$n_outgroup_covered == 0L]
  if (length(empty_out) > 0L) {
    warning("No outgroup terminal reaches ", min_coverage, " coverage in: ",
            paste(empty_out, collapse = ", "),
            ". These loci cannot contribute to the branch separating the outgroup, and any ",
            "calibration addressed by a node spanning it rests on the remaining loci.",
            call. = FALSE)
  }
  empty_in <- out$marker[out$n_ingroup_covered == 0L]
  if (length(empty_in) > 0L) {
    warning("No ingroup terminal reaches ", min_coverage, " coverage in: ",
            paste(empty_in, collapse = ", "), ".", call. = FALSE)
  }

  if (!is.null(out_csv)) {
    utils::write.csv(out, out_csv, row.names = FALSE)
    message("Marker coverage by group written to ", out_csv)
  }
  invisible(out)
}
