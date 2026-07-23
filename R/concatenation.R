#' Concatenate Locus Alignments and Build Partition Coordinate Maps
#'
#' Concatenates individual orthologous locus alignments end-to-end into a unified multilocus supermatrix.
#' Combining independent genomic loci increases statistical power to resolve difficult ancestral nodes while
#' allowing partitioned substitution modeling to account for mutational rate heterogeneity across genomic regions.
#' Exports partition files compatible with `RAxML-NG`, `IQ-TREE`, `PartitionFinder2`, and `MrBayes`.
#'
#' @param input_dir Character. Path to directory containing curated, aligned locus FASTA files.
#' @param output_dir Character. Path to destination root directory for concatenated alignments and partition maps.
#' @return A data frame containing supermatrix dimensions, taxon coverage, and genomic partition bounds.
#' @examples
#' \dontrun{
#' run_concatenation_pipeline(
#'   input_dir = "5_MAFFT_Cleaned/aligned_markers",
#'   output_dir = "6_Concatenated"
#' )
#' }
#' @export
run_concatenation_pipeline <- function(input_dir, output_dir) {
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
  file_super_phylip     <- file.path(concat_dir, "ALIGNMENT_supermatrix.phy")
  file_super_nexus      <- file.path(concat_dir, "ALIGNMENT_supermatrix.nex")
  
  file_part_raxml       <- file.path(concat_dir, "PARTITION_raxml_ng.txt")
  file_part_iqtree      <- file.path(concat_dir, "PARTITION_iqtree.part")
  file_part_pf2         <- file.path(concat_dir, "PARTITION_partitionfinder.cfg")
  file_part_mrbayes     <- file.path(concat_dir, "PARTITION_mrbayes.nex")
  
  file_log              <- file.path(logs_dir, "LOG_concatenation_run.txt")
  file_commands         <- file.path(logs_dir, "LOG_example_commands.txt")
  
  file_name_crosswalk   <- file.path(logs_dir, "SUPP_TABLE_taxon_name_crosswalk.csv")
  
  file_species_marker   <- file.path(tables_dir, "TABLE_species_marker_occupancy.csv")
  file_super_stats      <- file.path(tables_dir, "TABLE_supermatrix_statistics.csv")
  file_final_matrix     <- file.path(tables_dir, "TABLE_final_species_alignment_summary.csv")
  file_marker_stats     <- file.path(tables_dir, "TABLE_marker_statistics.csv")
  file_marker_ranges    <- file.path(tables_dir, "TABLE_marker_ranges.tsv")

  log_message <- function(...) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
    cat(msg, "\n")
    write(msg, file = file_log, append = TRUE)
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

  fasta_files <- sort(list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE))
  marker_raw <- tools::file_path_sans_ext(basename(fasta_files))
  marker_clean <- sub("^Masked_", "", marker_raw)

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
  registry_file <- file.path(dirname(dirname(input_dir)), "4_Cleaned", "TABLE_sequence_registry_with_acceptance.csv")
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
  concat_list <- list(); stats_list <- list(); marker_ranges <- dplyr::tibble(marker = character(), start = integer(), end = integer())
  current_pos <- 1L
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
  }
  
  readr::write_csv(dplyr::bind_rows(stats_list), file_marker_stats)
  write.table(marker_ranges, file = file_marker_ranges, sep = "\t", quote = FALSE, row.names = FALSE)

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

  ape::write.nexus.data(x = super_list, file = file_super_nexus, format = "dna", interleave = FALSE)

  part_lines <- dplyr::pull(dplyr::mutate(marker_ranges, line = paste0("DNA, ", marker, " = ", start, "-", end)), line)
  writeLines(part_lines, con = file_part_raxml)
  writeLines(part_lines, con = file_part_iqtree)

  pf2_text <- c("## PartitionFinder2 configuration file", "", "## ALIGNMENT FILE ##", paste0("alignment = ", basename(file_super_phylip), ";"), "", "## BRANCHLENGTHS: linked | unlinked ##", "branchlengths = linked;", "", "## MODELS OF EVOLUTION ##", "models = all;", "", "## MODEL SELECTION: AIC | AICc | BIC ##", "model_selection = AICc;", "", "## DATA BLOCKS ##", "[data_blocks]", paste0(marker_ranges$marker, " = ", marker_ranges$start, "-", marker_ranges$end, ";"), "", "## SCHEMES ##", "[schemes]", "search = greedy;")
  writeLines(pf2_text, con = file_part_pf2)

  cat("#nexus\nbegin mrbayes;\n", file = file_part_mrbayes)
  for (i in seq_len(nrow(marker_ranges))) cat(sprintf("   charset %s = %d-%d;\n", marker_ranges$marker[i], marker_ranges$start[i], marker_ranges$end[i]), file = file_part_mrbayes, append = TRUE)
  cat(sprintf("\n   partition combined = %d: %s;\n", nrow(marker_ranges), paste(marker_ranges$marker, collapse = ", ")), file = file_part_mrbayes, append = TRUE)
  cat("   set partition = combined;\nend;\n", file = file_part_mrbayes, append = TRUE)

  for (marker in marker_clean) seqinr::write.fasta(sequences = as.list(apply(concat_list[[marker]], 1, paste0, collapse = "")), names = rownames(concat_list[[marker]]), file.out = file.path(ind_dir, paste0(marker, ".fasta")))

  cmds <- c("# Example commands - adjust paths and threads as needed", "", "## IQ-TREE2", paste0("iqtree2 -s ", basename(file_super_phylip), " -spp ", basename(file_part_iqtree), " -m MFP+MERGE -B 1000 -T 4"), "", "## RAxML-NG", paste0("raxml-ng --all --msa ", basename(file_super_phylip), " --model ", basename(file_part_raxml), " --bs-trees 100 --threads 4"), "", "## PartitionFinder2", paste0("partitionfinder2 ", basename(file_part_pf2)), "", "## MrBayes", paste0("mb ", basename(file_part_mrbayes)))
  writeLines(cmds, con = file_commands)

  message("Concatenation pipeline finalized! \U0001f335")
}
