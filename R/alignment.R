#' Infer Positional Homology via MAFFT Alignment
#'
#' Establishes hypotheses of positional homology across unaligned orthologous nucleotide sequence clusters.
#' Positional homology alignment is a crucial prerequisite for maximum-likelihood phylogenetic inference,
#' ensuring that corresponding nucleotide sites derived from common evolutionary ancestry are aligned
#' prior to substitution model evaluation.
#'
#' @param input_fasta Character. Path to unaligned input FASTA file.
#' @param output_fasta Character. Path to destination aligned FASTA output file.
#' @param mafft_exec Character. System command or full path to the executable `MAFFT` binary. Defaults to `"mafft"`.
#' @param mafft_opts Character. Command-line parameters passed directly to `MAFFT`. Defaults to `"--auto"`.
#' @return Invisible numeric exit status code (0 for successful alignment completion).
#' @references
#' Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7:
#' Improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772–780.
#' \doi{10.1093/molbev/mst010}
#' @examples
#' \dontrun{
#' run_mafft(
#'   input_fasta = "raw_cluster.fasta",
#'   output_fasta = "aligned_cluster.fasta",
#'   mafft_opts = "--auto"
#' )
#' }
#' @export
run_mafft <- function(input_fasta, output_fasta, mafft_exec = "mafft", mafft_opts = "--auto") {
  args <- c(strsplit(mafft_opts, "\\s+")[[1]], shQuote(input_fasta))
  stderr_fasta <- paste0(output_fasta, ".stderr.log")
  
  status <- system2(
    mafft_exec,
    args = args,
    stdout = output_fasta,
    stderr = stderr_fasta
  )
  
  if (!file.exists(output_fasta) || file.info(output_fasta)$size == 0) {
    stop("MAFFT did not produce a readable output FASTA. Check logs: ", stderr_fasta, call. = FALSE)
  }
  
  if (!is.null(status) && status != 0) {
    stop(
      sprintf("MAFFT returned non-zero exit status (%s). Check: %s", status, stderr_fasta),
      call. = FALSE
    )
  }
  invisible(status)
}

stop_if_missing_dir <- function(path_dir, label) {
  if (!dir.exists(path_dir)) {
    stop(sprintf("%s does not exist: %s", label, path_dir), call. = FALSE)
  }
}

check_mafft_available <- function(mafft_exec) {
  test <- suppressWarnings(
    system2(mafft_exec, args = "--version", stdout = TRUE, stderr = TRUE)
  )
  status <- attr(test, "status")
  if (!is.null(status) && status != 0) {
    stop("MAFFT is not available or failed to run with '--version'.", call. = FALSE)
  }
  test
}

safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_string_counts <- function(s, pattern) {
  if (length(s) == 0 || is.na(s) || nchar(s) == 0) return(0L)
  sum(strsplit(s, "", fixed = TRUE)[[1]] == pattern)
}

clean_ambiguous <- function(x) {
  seqs <- as.character(x)
  seqs_up <- toupper(seqs)
  seqs_clean <- gsub("[^ACGT-]", "N", seqs_up, perl = TRUE)
  out <- Biostrings::DNAStringSet(seqs_clean)
  names(out) <- names(x)
  out
}

get_alignment_stats <- function(x) {
  seqs <- as.character(x)
  
  if (length(seqs) == 0) {
    return(data.frame(
      n_sequences = 0L,
      min_sequence_length = NA_integer_,
      max_sequence_length = NA_integer_,
      mean_sequence_length = NA_real_,
      mean_fraction_missing = NA_real_,
      mean_fraction_gaps = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  seqlen <- nchar(seqs)
  n_missing <- vapply(seqs, safe_string_counts, integer(1), pattern = "N")
  n_gaps <- vapply(seqs, safe_string_counts, integer(1), pattern = "-")
  
  pct_missing <- ifelse(seqlen > 0, n_missing / seqlen, NA_real_)
  pct_gaps <- ifelse(seqlen > 0, n_gaps / seqlen, NA_real_)
  
  data.frame(
    n_sequences = length(seqs),
    min_sequence_length = min(seqlen),
    max_sequence_length = max(seqlen),
    mean_sequence_length = safe_mean(seqlen),
    mean_fraction_missing = safe_mean(pct_missing),
    mean_fraction_gaps = safe_mean(pct_gaps),
    stringsAsFactors = FALSE
  )
}

write_run_log <- function(log_file, input_folder, output_root, mafft_exec, mafft_opts) {
  mafft_version <- tryCatch(
    paste(check_mafft_available(mafft_exec), collapse = " "),
    error = function(e) paste("Unavailable:", conditionMessage(e))
  )
  
  lines <- c(
    "Alignment run log",
    paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("Input folder:", input_folder),
    paste("Output root:", output_root),
    paste("MAFFT executable:", mafft_exec),
    paste("MAFFT options:", mafft_opts),
    paste("MAFFT version:", mafft_version),
    "",
    "R session info:",
    paste(capture.output(sessionInfo()), collapse = "\n")
  )
  writeLines(lines, con = log_file)
}

build_marker_consolidated_row <- function(marker_id, summary_all, filter_info, final_alignment_written, masking_applied, status) {
  get_stage_row <- function(stage_name) {
    out <- summary_all[summary_all$stage == stage_name, , drop = FALSE]
    if (nrow(out) == 0) {
      return(data.frame(
        stage = stage_name,
        n_sequences = NA_integer_,
        min_sequence_length = NA_integer_,
        max_sequence_length = NA_integer_,
        mean_sequence_length = NA_real_,
        mean_fraction_missing = NA_real_,
        mean_fraction_gaps = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    out[1, , drop = FALSE]
  }
  
  raw_row    <- get_stage_row("raw_input")
  mafft_row  <- get_stage_row("mafft_aligned")
  common_row <- get_stage_row("common_gap_removed")
  masked_row <- get_stage_row("masked_alignment")
  final_row  <- get_stage_row("final_filtered")
  
  n_removed_low_coverage <- if (!is.null(filter_info) && nrow(filter_info) > 0) sum(filter_info$removed_low_coverage, na.rm = TRUE) else 0L
  n_removed_high_missingness <- if (!is.null(filter_info) && nrow(filter_info) > 0) sum(filter_info$removed_high_missingness, na.rm = TRUE) else 0L
  n_removed_total_filter <- if (!is.null(filter_info) && nrow(filter_info) > 0) sum(filter_info$removed_final, na.rm = TRUE) else 0L
  
  raw_n    <- raw_row$n_sequences
  mafft_n  <- mafft_row$n_sequences
  common_n <- common_row$n_sequences
  masked_n <- masked_row$n_sequences
  final_n  <- final_row$n_sequences
  
  raw_len    <- raw_row$max_sequence_length
  mafft_len  <- mafft_row$max_sequence_length
  common_len <- common_row$max_sequence_length
  masked_len <- masked_row$max_sequence_length
  final_len  <- final_row$max_sequence_length
  
  data.frame(
    marker = marker_id,
    
    n_sequences_raw_input = raw_n,
    n_sequences_mafft_aligned = mafft_n,
    n_sequences_common_gap_removed = common_n,
    n_sequences_masked_alignment = masked_n,
    n_sequences_final_filtered = final_n,
    
    alignment_length_raw_input = raw_len,
    alignment_length_mafft_aligned = mafft_len,
    alignment_length_common_gap_removed = common_len,
    alignment_length_masked_alignment = masked_len,
    alignment_length_final_filtered = final_len,
    
    mean_fraction_missing_raw_input = raw_row$mean_fraction_missing,
    mean_fraction_missing_mafft_aligned = mafft_row$mean_fraction_missing,
    mean_fraction_missing_common_gap_removed = common_row$mean_fraction_missing,
    mean_fraction_missing_masked_alignment = masked_row$mean_fraction_missing,
    mean_fraction_missing_final_filtered = final_row$mean_fraction_missing,
    
    mean_fraction_gaps_raw_input = raw_row$mean_fraction_gaps,
    mean_fraction_gaps_mafft_aligned = mafft_row$mean_fraction_gaps,
    mean_fraction_gaps_common_gap_removed = common_row$mean_fraction_gaps,
    mean_fraction_gaps_masked_alignment = masked_row$mean_fraction_gaps,
    mean_fraction_gaps_final_filtered = final_row$mean_fraction_gaps,
    
    n_sequences_removed_total = if (!is.na(raw_n) && !is.na(final_n)) raw_n - final_n else NA_integer_,
    n_sequences_removed_low_coverage = n_removed_low_coverage,
    n_sequences_removed_high_missingness = n_removed_high_missingness,
    n_sequences_removed_by_filtering = n_removed_total_filter,
    
    n_sites_removed_by_common_gap_removal = if (!is.na(mafft_len) && !is.na(common_len)) mafft_len - common_len else NA_integer_,
    n_sites_removed_by_masking = if (!is.na(common_len) && !is.na(masked_len)) common_len - masked_len else NA_integer_,
    n_sites_removed_total_from_mafft = if (!is.na(mafft_len) && !is.na(final_len)) mafft_len - final_len else NA_integer_,
    
    pct_sequences_retained_after_mafft = if (!is.na(raw_n) && raw_n > 0) mafft_n / raw_n else NA_real_,
    pct_sequences_retained_after_common_gap_removal = if (!is.na(raw_n) && raw_n > 0) common_n / raw_n else NA_real_,
    pct_sequences_retained_after_masking = if (!is.na(raw_n) && raw_n > 0) masked_n / raw_n else NA_real_,
    pct_sequences_retained_final = if (!is.na(raw_n) && raw_n > 0) final_n / raw_n else NA_real_,
    
    pct_sites_retained_after_common_gap_removal = if (!is.na(mafft_len) && mafft_len > 0) common_len / mafft_len else NA_real_,
    pct_sites_retained_after_masking = if (!is.na(mafft_len) && mafft_len > 0) masked_len / mafft_len else NA_real_,
    pct_sites_retained_final_vs_mafft = if (!is.na(mafft_len) && mafft_len > 0) final_len / mafft_len else NA_real_,
    
    masking_applied = masking_applied,
    final_alignment_written = final_alignment_written,
    status = status,
    stringsAsFactors = FALSE
  )
}

process_marker_file <- function(
    fasta_file,
    output_aln_dir,
    output_table_dir,
    mafft_exec,
    mafft_opts,
    mask_alignment_regions = TRUE,
    min_non_gap_fraction = 0.30,
    max_missing_fraction = 0.30
) {
  marker_file <- basename(fasta_file)
  marker_id <- tools::file_path_sans_ext(marker_file)
  
  clean_fasta_path <- file.path(output_aln_dir, paste0("ALN_clean_input_", marker_id, ".fasta"))
  mafft_fasta_path <- file.path(output_aln_dir, paste0("ALN_mafft_raw_", marker_id, ".fasta"))
  final_fasta_path <- file.path(output_aln_dir, paste0("ALN_masked_final_", marker_id, ".fasta"))
  
  summary_csv <- file.path(output_table_dir, paste0("TABLE_marker_alignment_summary_", marker_id, ".csv"))
  filter_csv  <- file.path(output_table_dir, paste0("TABLE_marker_sequence_filter_log_", marker_id, ".csv"))
  
  if (file.exists(final_fasta_path) && file.exists(summary_csv)) {
    message("CACHE: Marker already processed: ", marker_id)
    return(invisible(NULL))
  }
  
  raw_in <- Biostrings::readDNAStringSet(fasta_file)
  if (length(raw_in) == 0) stop("Input FASTA contains zero sequences.", call. = FALSE)
  if (anyDuplicated(names(raw_in)) > 0) warning(sprintf("Duplicated sequence names detected in %s", marker_file), call. = FALSE)
  
  stats_before <- get_alignment_stats(raw_in)
  
  raw_in_clean <- clean_ambiguous(raw_in)
  Biostrings::writeXStringSet(raw_in_clean, clean_fasta_path)
  
  run_mafft(input_fasta = clean_fasta_path, output_fasta = mafft_fasta_path, mafft_exec = mafft_exec, mafft_opts = mafft_opts)
  
  aln <- Biostrings::readDNAStringSet(mafft_fasta_path)
  if (length(aln) == 0) stop("Aligned FASTA contains zero sequences after MAFFT.", call. = FALSE)
  
  stats_after_mafft <- get_alignment_stats(aln)
  
  if (!mask_alignment_regions) {
    Biostrings::writeXStringSet(aln, final_fasta_path)
    summary_all <- rbind(
      cbind(stage = "raw_input", stats_before),
      cbind(stage = "mafft_aligned", stats_after_mafft),
      cbind(stage = "common_gap_removed", stats_after_mafft),
      cbind(stage = "masked_alignment", stats_after_mafft),
      cbind(stage = "final_filtered", stats_after_mafft)
    )
    write.csv(summary_all, summary_csv, row.names = FALSE)
    
    filter_info <- data.frame(
      sequence_id = names(aln),
      non_gap_sites_after_masking = NA_integer_,
      fraction_missing_after_masking = NA_real_,
      removed_low_coverage = FALSE,
      removed_high_missingness = FALSE,
      removed_final = FALSE,
      stringsAsFactors = FALSE
    )
    write.csv(filter_info, filter_csv, row.names = FALSE)
    
    return(build_marker_consolidated_row(marker_id, summary_all, filter_info, TRUE, FALSE, "OK_NO_MASK"))
  }
  
  aln_no_gaps <- DECIPHER::RemoveGaps(aln, removeGaps = "common")
  stats_after_remove_gaps <- get_alignment_stats(aln_no_gaps)
  
  aln_masked <- as(DECIPHER::MaskAlignment(aln_no_gaps, correction = (length(aln_no_gaps) < 200)), "DNAStringSet")
  stats_after_mask <- get_alignment_stats(aln_masked)
  
  masked_strings <- as.character(aln_masked)
  masked_widths <- nchar(masked_strings)
  
  if (length(masked_widths) == 0 || max(masked_widths, na.rm = TRUE) == 0) {
    remove_gap_filter <- rep(TRUE, length(masked_strings))
    non_gaps <- rep(0L, length(masked_strings))
    pct_missing_masked <- rep(NA_real_, length(masked_strings))
  } else {
    aln_length_masked <- max(masked_widths, na.rm = TRUE)
    non_gaps <- nchar(gsub("-", "", masked_strings))
    min_non_gap_sites <- aln_length_masked * min_non_gap_fraction
    remove_gap_filter <- non_gaps < min_non_gap_sites
    
    pct_missing_masked <- vapply(masked_strings, function(s) {
      if (is.na(s) || nchar(s) == 0) return(NA_real_)
      safe_string_counts(s, "N") / nchar(s)
    }, numeric(1))
  }
  
  remove_missing_filter <- pct_missing_masked > max_missing_fraction
  remove_missing_filter[is.na(remove_missing_filter)] <- TRUE
  remove_final <- remove_gap_filter | remove_missing_filter
  
  aln_final <- aln_masked[!remove_final]
  stats_final <- get_alignment_stats(aln_final)
  final_written <- FALSE
  if (length(aln_final) > 0) {
    Biostrings::writeXStringSet(aln_final, final_fasta_path)
    final_written <- TRUE
  }
  
  filter_info <- data.frame(
    sequence_id = names(aln_masked),
    non_gap_sites_after_masking = non_gaps,
    fraction_missing_after_masking = pct_missing_masked,
    removed_low_coverage = remove_gap_filter,
    removed_high_missingness = remove_missing_filter,
    removed_final = remove_final,
    stringsAsFactors = FALSE
  )
  write.csv(filter_info, filter_csv, row.names = FALSE)
  
  summary_all <- rbind(
    cbind(stage = "raw_input", stats_before),
    cbind(stage = "mafft_aligned", stats_after_mafft),
    cbind(stage = "common_gap_removed", stats_after_remove_gaps),
    cbind(stage = "masked_alignment", stats_after_mask),
    cbind(stage = "final_filtered", stats_final)
  )
  write.csv(summary_all, summary_csv, row.names = FALSE)
  
  return(build_marker_consolidated_row(marker_id, summary_all, filter_info, final_written, TRUE, if (final_written) "OK" else "ZERO_RETAINED"))
}

#' Execute Complete Alignment and Gap-Masking Pipeline
#'
#' Orchestrates multiple sequence alignment (MSA) and automated quality-control masking across orthologous sequence clusters.
#' Primary alignment hypotheses are inferred using `MAFFT` (Katoh & Standley, 2013). Subsequently, ambiguous sites,
#' poorly aligned terminal fragments, and non-homologous insertions are masked using the `DECIPHER` framework (Wright, 2024),
#' eliminating systematic noise while retaining phylogenetically informative nucleotide positions for downstream supermatrix assembly.
#'
#' @param input_folder Character. Directory path containing raw unaligned orthologous FASTA files.
#' @param output_dir Character. Root destination directory for output subfolders (`alignments/`, `tables/`, `logs/`).
#' @param fasta_pattern Character. Regular expression pattern matching target FASTA files. Defaults to `"\\.fasta$"`.
#' @param mask_alignment_regions Logical. Apply automated alignment masking via `DECIPHER`? Defaults to `TRUE`.
#' @param min_non_gap_fraction Numeric. Minimum allowable proportion of non-gap characters required to retain a site column. Defaults to `0.30`.
#' @param max_missing_fraction Numeric. Maximum allowable proportion of missing or ambiguous characters (`N`) allowed per sequence. Defaults to `0.30`.
#' @param mafft_exec Character. System command or full path to the executable `MAFFT` binary. Defaults to `"mafft"`.
#' @param mafft_opts Character. Command-line parameters passed directly to `MAFFT`. Defaults to `"--auto"`.
#' @return A data frame containing site length, missingness, and sequence retention statistics across processed loci.
#' @references
#' Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7:
#' Improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772–780.
#' \doi{10.1093/molbev/mst010}
#'
#' Wright, E. S. (2024). Fast and Flexible Search for Homologous Biological Sequences with DECIPHER v3.
#' *The R Journal*, 16(2), 191-200. \doi{10.18129/B9.bioc.DECIPHER}
#' @examples
#' \dontrun{
#' run_alignment_pipeline(
#'   input_folder = "1_phylotaR_out_ingroup",
#'   output_dir = "2_MAFFT_Cactaceae",
#'   min_non_gap_fraction = 0.30,
#'   max_missing_fraction = 0.30
#' )
#' }
#' @export
run_alignment_pipeline <- function(
    input_folder,
    output_dir,
    fasta_pattern = "\\.fasta$",
    mask_alignment_regions = TRUE,
    min_non_gap_fraction = 0.30,
    max_missing_fraction = 0.30,
    mafft_exec = "mafft",
    mafft_opts = "--auto"
) {
  stop_if_missing_dir(input_folder, "Input folder")
  
  output_aln_dir <- file.path(output_dir, "alignments")
  output_table_dir <- file.path(output_dir, "tables")
  output_log_dir <- file.path(output_dir, "logs")

  dir.create(output_aln_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_log_dir, recursive = TRUE, showWarnings = FALSE)
  
  check_mafft_available(mafft_exec)
  
  log_file <- file.path(output_log_dir, "LOG_alignment_run_info.txt")
  write_run_log(log_file, input_folder, output_dir, mafft_exec, mafft_opts)
  
  fasta_files <- sort(list.files(input_folder, pattern = fasta_pattern, full.names = TRUE))
  if (length(fasta_files) == 0) stop("No FASTA files found in input folder.", call. = FALSE)
  
  manifest_list <- vector("list", length(fasta_files))
  for (i in seq_along(fasta_files)) {
    f <- fasta_files[i]
    message("\nProcessing: ", basename(f))
    
    manifest_list[[i]] <- tryCatch(
      process_marker_file(f, output_aln_dir, output_table_dir, mafft_exec, mafft_opts, mask_alignment_regions, min_non_gap_fraction, max_missing_fraction),
      error = function(e) {
        data.frame(
          marker = tools::file_path_sans_ext(basename(f)),
          n_sequences_raw_input = NA_integer_, n_sequences_mafft_aligned = NA_integer_, n_sequences_common_gap_removed = NA_integer_,
          n_sequences_masked_alignment = NA_integer_, n_sequences_final_filtered = NA_integer_,
          alignment_length_raw_input = NA_integer_, alignment_length_mafft_aligned = NA_integer_, alignment_length_common_gap_removed = NA_integer_,
          alignment_length_masked_alignment = NA_integer_, alignment_length_final_filtered = NA_integer_,
          mean_fraction_missing_raw_input = NA_real_, mean_fraction_missing_mafft_aligned = NA_real_, mean_fraction_missing_common_gap_removed = NA_real_,
          mean_fraction_missing_masked_alignment = NA_real_, mean_fraction_missing_final_filtered = NA_real_,
          mean_fraction_gaps_raw_input = NA_real_, mean_fraction_gaps_mafft_aligned = NA_real_, mean_fraction_gaps_common_gap_removed = NA_real_,
          mean_fraction_gaps_masked_alignment = NA_real_, mean_fraction_gaps_final_filtered = NA_real_,
          n_sequences_removed_total = NA_integer_, n_sequences_removed_low_coverage = NA_integer_, n_sequences_removed_high_missingness = NA_integer_,
          n_sequences_removed_by_filtering = NA_integer_, n_sites_removed_by_common_gap_removal = NA_integer_, n_sites_removed_by_masking = NA_integer_,
          n_sites_removed_total_from_mafft = NA_integer_, pct_sequences_retained_after_mafft = NA_real_, pct_sequences_retained_after_common_gap_removal = NA_real_,
          pct_sequences_retained_after_masking = NA_real_, pct_sequences_retained_final = NA_real_, pct_sites_retained_after_common_gap_removal = NA_real_,
          pct_sites_retained_after_masking = NA_real_, pct_sites_retained_final_vs_mafft = NA_real_,
          masking_applied = mask_alignment_regions, final_alignment_written = FALSE, status = paste("ERROR:", conditionMessage(e)), stringsAsFactors = FALSE
        )
      }
    )
  }
  
  valid_manifests <- Filter(function(x) !is.null(x), manifest_list)
  manifest <- do.call(rbind, valid_manifests)
  
  manifest_path <- file.path(output_log_dir, "LOG_marker_processing_manifest.csv")
  write.csv(manifest, manifest_path, row.names = FALSE)
  
  final_summary_path <- file.path(output_table_dir, "TABLE_alignment_summary_all_markers.csv")
  write.csv(manifest, final_summary_path, row.names = FALSE)
  
  message("\nRun completed. \U0001f335")
  return(manifest)
}

#' Reconcile and Validate Taxonomic Nomenclature
#'
#' Reconciles sequence tip labels against authoritative botanical checklists (e.g., Caryophyllales.org checklist).
#' Resolves taxonomic synonymies, infraspecific variants, and orthographic errors, guaranteeing nomenclatural stability
#' across public GenBank sequence downloads.
#'
#' @param raw_input_fasta Character. Path to input FASTA file containing raw GenBank sequence accessions.
#' @param checklist_path Character. Path to accepted taxonomic checklist CSV or Excel file.
#' @param output_clean_dir Character. Directory path where standardized FASTA sequence output will be saved.
#' @param force_process Logical. Force reprocessing if output cached file exists? Defaults to `FALSE`.
#' @return Character string path to the generated FASTA file with standardized species binomials.
#' @examples
#' \dontrun{
#' clean_taxonomic_names(
#'   raw_input_fasta = "marker_rbcl.fasta",
#'   checklist_path = "CactaceaeFullList_accepted.csv",
#'   output_clean_dir = "cleaned_names"
#' )
#' }
#' @export
clean_taxonomic_names <- function(raw_input_fasta, checklist_path, output_clean_dir, force_process = FALSE) {
  
  dir.create(output_clean_dir, recursive = TRUE, showWarnings = FALSE)
  marker_id <- tools::file_path_sans_ext(basename(raw_input_fasta))
  target_out <- file.path(output_clean_dir, paste0(marker_id, "_clean_names.fasta"))
  
  if (!force_process && file.exists(target_out)) {
    message("CACHE: Target clean names FASTA already calculated: ", target_out)
    return(target_out)
  }
  
  dna_in <- Biostrings::readDNAStringSet(raw_input_fasta)
  
  if (grepl("\\.xlsx?$", checklist_path)) {
    checklist_df <- readxl::read_excel(checklist_path)
  } else {
    checklist_df <- read.csv(checklist_path, stringsAsFactors = FALSE)
  }
  
  accepted_names <- unique(trimws(checklist_df$pureName))
  
  names_seqs <- names(dna_in)
  binomial_leaves <- sub("^([^ ]+_[^ ]+).*", "\\1", names_seqs)
  keep_idx <- binomial_leaves %in% accepted_names
  dna_clean <- dna_in[keep_idx]
  names(dna_clean) <- sub("_", " ", binomial_leaves[keep_idx])
  
  dna_final <- dna_clean[!duplicated(names(dna_clean))]
  Biostrings::writeXStringSet(dna_final, target_out)
  return(target_out)
}

#' Perform Joint Realignment Across Integrated Ingroup and Outgroup Sequences
#'
#' Re-estimates positional homology alignments (`MAFFT`) after aggregating ingroup and outgroup sequence markers.
#' Joint realignment accommodates positional variations introduced when combining outgroup taxons with the focal ingroup clade.
#'
#' @param input_dir Character. Directory containing merged, taxonomically reconciled locus FASTA files.
#' @param output_fasta_dir Character. Directory path to save output realigned FASTA sequence files.
#' @param output_aln_dir Character. Directory path to store detailed realignment log reports.
#' @return Invisible NULL upon successful execution.
#' @references
#' Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7:
#' Improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772–780.
#' \doi{10.1093/molbev/mst010}
#' @examples
#' \dontrun{
#' run_joint_realignment(
#'   input_dir = "4_Cleaned/cleaned_markers",
#'   output_fasta_dir = "5_MAFFT_Cleaned",
#'   output_aln_dir = "5_MAFFT_Cleaned/aligned_markers"
#' )
#' }
#' @export
run_joint_realignment <- function(input_dir, output_fasta_dir, output_aln_dir) {
  dir.create(output_fasta_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_aln_dir, recursive = TRUE, showWarnings = FALSE)
  
  log_file <- file.path(output_fasta_dir, "LOG_ALN_FINAL_RUN.txt")
  writeLines(paste0("Run started: ", Sys.time()), log_file)
  
  log_write <- function(text) {
    con <- file(log_file, open = "a")
    writeLines(text, con)
    close(con)
  }
  
  files <- sort(list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE))
  log_write(paste("Markers detected:", length(files)))
  
  collapsed_markers <- character()
  
  for (f in files) {
    fn <- basename(f)
    message("Processing joint realignment: ", fn)
    
    raw <- Biostrings::readDNAStringSet(f)
    raw_clean <- clean_ambiguous(raw)
    
    temp_clean <- file.path(output_fasta_dir, paste0("TEMP_", fn))
    Biostrings::writeXStringSet(raw_clean, temp_clean)
    
    raw_out <- file.path(output_fasta_dir, paste0("TEMP_ALN_", fn))
    raw_err <- file.path(output_fasta_dir, paste0("TEMP_ALN_", tools::file_path_sans_ext(fn), ".log"))
    
    status <- system2(
      command = "mafft",
      args = c("--auto", temp_clean),
      stdout = raw_out,
      stderr = raw_err
    )
    
    if (!identical(status, 0L)) {
      err_msg <- if (file.exists(raw_err)) paste(readLines(raw_err, warn = FALSE), collapse = "\n") else "No stderr captured."
      stop("MAFFT failed for: ", fn, "\n", err_msg)
    }
    
    aln <- Biostrings::readDNAStringSet(raw_out)
    aln <- DECIPHER::RemoveGaps(aln, removeGaps = "common")
    aln <- DECIPHER::MaskAlignment(aln, correction = (length(aln) < 200))
    aln <- as(aln, "DNAStringSet")
    
    aln_char <- as.character(aln)
    aln_len <- nchar(aln_char[1])
    
    non_gaps <- nchar(gsub("-", "", aln_char))
    pctMissing <- vapply(aln_char, function(s) sum(strsplit(s, "")[[1]] == "N") / nchar(s), numeric(1))
    
    keep <- non_gaps >= (0.3 * aln_len) & pctMissing <= 0.3
    aln_final <- aln[keep]
    
    if (length(aln_final) == 0) {
      collapsed_markers <- c(collapsed_markers, fn)
    }
    
    if (length(aln_final) > 0) {
      out_fasta <- file.path(output_aln_dir, fn)
      Biostrings::writeXStringSet(aln_final, out_fasta)
    }
    
    # Quality control exports (always export, even if collapsed)
    utils::write.csv(data.frame(
      Seq = names(aln),
      NonGaps = non_gaps,
      pctMissing = pctMissing,
      Retained = keep
    ), file = file.path(output_fasta_dir, paste0("LOG_SEQ_FILTER_", tools::file_path_sans_ext(fn), ".csv")),
    row.names = FALSE)
    
    # Calculate stats
    seqs <- as.character(aln_final)
    if (length(seqs) == 0) {
      stats_df <- data.frame(nSeq = 0, minLen = NA_real_, maxLen = NA_real_, avgLen = NA_real_, pctMissing = NA_real_, pctGaps = NA_real_)
    } else {
      seqlen <- nchar(seqs)
      n_missing <- vapply(strsplit(seqs, "", fixed = TRUE), function(z) sum(z == "N"), numeric(1))
      n_gaps <- vapply(strsplit(seqs, "", fixed = TRUE), function(z) sum(z == "-"), numeric(1))
      stats_df <- data.frame(
        nSeq = length(seqs),
        minLen = min(seqlen),
        maxLen = max(seqlen),
        avgLen = mean(seqlen),
        pctMissing = mean(n_missing / seqlen),
        pctGaps = mean(n_gaps / seqlen)
      )
    }
    utils::write.csv(stats_df,
              file = file.path(output_fasta_dir, paste0("TABLE_ALN_SUMMARY_", tools::file_path_sans_ext(fn), ".csv")),
              row.names = FALSE)
              
    file.remove(temp_clean, raw_out, raw_err)
  }
  
  log_write(paste("Collapsed markers:", paste(collapsed_markers, collapse = ", ")))
  log_write(paste("Run finished:", Sys.time()))
  
  message("\nJoint realignment completed. \U0001f335")
}


