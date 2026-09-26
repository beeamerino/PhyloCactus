# ------------------------------------------------------------------------------
# Molecular diagnostic branch (11_barcoding/), final library step. PhyloCactus 0.5.0.
#
# Decisions of BMM, 2026-09-21 (Phase 2 proposal, sec. 8):
# - Sequences that run_alignment_pipeline() logs as "no_match_either_direction" in
#   LOG_STRAND_<locus>.csv leave the library and are declared in a table.
# - Within each species and locus, accessions with the same aligned sequence are reduced to one, the
#   alphabetically smallest sid; a table maps every accession to its representative. Identical
#   sequences of different species are not touched.
# - A species has replica when it has two or more distinct sequences in the locus. The replication
#   threshold is evaluated again after the collapse.
# The columns of the curated alignment are not modified.
# ------------------------------------------------------------------------------

#' Collapse identical aligned sequences within each species of one locus
#' @param aln Named character vector of aligned sequences, headers Genus_species|sid.
#' @param locus Character. Locus name, copied into the mapping table.
#' @return List with `aln` (representatives only) and `map` (one row per accession: locus, species,
#'   sid, representative, collapsed).
#' @noRd
.bc_collapse_identical <- function(aln, locus) {
  h <- .bc_parse_header(names(aln))
  key <- paste(h$species, toupper(unname(aln)), sep = "\r")
  # C-locale ordering, so the representative does not depend on the session locale
  first_sid <- tapply(h$sid, key, function(x) sort(x, method = "radix")[1])
  representative <- unname(first_sid[key])
  map <- data.frame(locus = rep(locus, length(aln)), species = h$species, sid = h$sid,
                    representative = representative, collapsed = h$sid != representative,
                    stringsAsFactors = FALSE)
  list(aln = aln[!map$collapsed], map = map)
}

#' Build the Final Reference Library of the Molecular Diagnostic Branch
#'
#' Takes the curated alignments written by [curate_barcoding_markers()] and writes the reference
#' library used by the validation. For each locus it (1) removes the sequences that the alignment
#' pipeline logged as matching the locus in neither direction (`no_match_either_direction` in
#' `LOG_STRAND_<locus>.csv`), which are paralogs or other non-homologous copies, not intraspecific
#' replication; (2) collapses accessions with the same aligned sequence within each species to the
#' alphabetically smallest sid; (3) evaluates the replication threshold again, counting a species as
#' replicated when it keeps two or more distinct sequences. Identical sequences of different species
#' are kept. The columns of the curated alignment are not modified, so removing sequences can leave
#' columns made only of gaps.
#'
#' Every input is read from the directory of the step that wrote it, so the library runs on its own
#' whenever those directories exist. The loci are those marked `enters` in the screening table.
#'
#' @param assembly_dir Character. Output directory of [assemble_barcoding_dataset()]; the accession
#'   registry is read from it.
#' @param curated_dir Character. Output directory of [curate_barcoding_markers()]; the curated
#'   alignments and the `LOG_STRAND_<locus>.csv` files are read from it.
#' @param screening_dir Character. Output directory of [screen_barcoding_markers()]; the loci that
#'   entered the screening are read from `TABLE_barcoding_marker_screening.csv`.
#' @param output_dir Character. Destination; refused inside a directory of the phylogeny.
#' @param metadata_file Character. GenBank metadata cache written by [assemble_barcoding_dataset()],
#'   read for the description of each excluded sequence.
#' @param min_species_with_replica Integer. A locus enters the library when at least this number of
#'   species keep two or more distinct sequences. Defaults to `20L` (decision of 2026-09-21).
#' @param paralog_loci Character vector. Loci under paralogy surveillance, flagged in the column
#'   `possible_paralog` of every table; the flag never excludes a locus. Defaults to `"pepC_like"`.
#' @return Invisibly, a list with `summary`, `excluded` and `collapsed`. Writes `LIB_<locus>.fasta`
#'   for every locus that enters, `TABLE_barcoding_library_summary.csv`,
#'   `TABLE_barcoding_excluded_nonhomologous.csv`, `TABLE_barcoding_collapsed_identical.csv` and
#'   `TABLE_barcoding_funnel.csv` (per locus: accessions assembled, curated, excluded, collapsed and
#'   kept, species with replica, and whether the locus entered the screening and the library).
#' @examples
#' \dontrun{
#' lib <- finalize_barcoding_library(
#'   assembly_dir = "11_barcoding/1_assembly",
#'   curated_dir = "11_barcoding/2_curated",
#'   screening_dir = "11_barcoding/3_screening"
#' )
#' }
#' @export
finalize_barcoding_library <- function(assembly_dir = file.path("11_barcoding", "1_assembly"),
                                       curated_dir = file.path("11_barcoding", "2_curated"),
                                       screening_dir = file.path("11_barcoding", "3_screening"),
                                       output_dir = file.path("11_barcoding", "4_library"),
                                       metadata_file = file.path("11_barcoding", "cache",
                                                                 "CACHE_GENBANK_METADATA_BARCODING.csv"),
                                       min_species_with_replica = 20L,
                                       paralog_loci = "pepC_like") {
  .bc_assert_output_dir(output_dir)
  registry <- .bc_read_registry(assembly_dir)
  screening <- .bc_read_step_table(file.path(screening_dir, "TABLE_barcoding_marker_screening.csv"),
                                   "screen_barcoding_markers")
  loci <- screening$locus[screening$enters %in% TRUE]
  if (length(loci) == 0) stop("No locus entered the screening in ", screening_dir, ".", call. = FALSE)
  message("Loci that entered the screening: ", paste(loci, collapse = ", "))
  aln_files <- file.path(curated_dir, "alignments", paste0("ALN_masked_final_", loci, ".fasta"))
  if (any(!file.exists(aln_files))) {
    stop("No curated alignment for locus: ", paste(loci[!file.exists(aln_files)], collapse = ", "),
         ". Run curate_barcoding_markers() first.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  meta <- NULL
  if (!is.null(metadata_file) && file.exists(metadata_file)) {
    meta <- utils::read.csv(metadata_file, stringsAsFactors = FALSE)
  } else {
    message("No GenBank metadata file found; the descriptions of excluded sequences are left empty.")
  }

  excluded <- list()
  collapsed <- list()
  rows <- list()
  n_after <- integer(0)
  for (i in seq_along(loci)) {
    l <- loci[i]
    dna <- Biostrings::readDNAStringSet(aln_files[i])
    aln <- stats::setNames(as.character(dna), names(dna))

    # 1. Non-homologous sequences, as logged by the strand check of the alignment pipeline
    log_file <- file.path(curated_dir, "tables", paste0("LOG_STRAND_", l, ".csv"))
    unmatched <- character(0)
    if (file.exists(log_file)) {
      lg <- utils::read.csv(log_file, stringsAsFactors = FALSE)
      lg <- lg[lg$Action == "no_match_either_direction", , drop = FALSE]
      unmatched <- lg$Seq
      if (nrow(lg) > 0) {
        h <- .bc_parse_header(lg$Seq)
        excluded[[l]] <- data.frame(
          locus = l, sid = h$sid, species = h$species,
          description = if (is.null(meta)) NA_character_ else meta$Description_gb[match(h$sid, meta$sid)],
          reason = "no_match_either_direction",
          share_forward = lg$ShareForward, share_reverse = lg$ShareReverse,
          in_curated_alignment = lg$Seq %in% names(aln),
          stringsAsFactors = FALSE
        )
      }
    } else {
      message("Locus '", l, "': no LOG_STRAND file; no sequence excluded for homology.")
    }
    n_excluded <- sum(names(aln) %in% unmatched)
    aln <- aln[!names(aln) %in% unmatched]

    # 2. Identical sequences within species
    n_rep_before <- sum(table(.bc_parse_header(names(aln))$species) >= 2L)
    col <- .bc_collapse_identical(aln, l)
    collapsed[[l]] <- col$map
    n_after[l] <- length(col$aln)

    # 3. Threshold on distinct sequences
    reps <- .bc_parse_header(names(col$aln))$sid
    reg_l <- registry[registry$locus == l & registry$sid %in% reps, , drop = FALSE]
    ctx <- .bc_locus_summary(reg_l)
    if (is.null(ctx) || nrow(ctx) == 0) {
      ctx <- data.frame(locus = l, total_species = 0L, species_with_replicate = 0L, total_accessions = 0L,
                        total_genera = 0L, genera_with_2plus_species = 0L, stringsAsFactors = FALSE)
    }
    if (nrow(reg_l) != length(reps)) {
      warning(sprintf("Locus '%s': %d representative(s) not found in the registry.", l, length(reps) - nrow(reg_l)),
              call. = FALSE)
    }
    passes <- ctx$species_with_replicate >= min_species_with_replica
    rows[[l]] <- data.frame(ctx,
                            species_with_replicate_before_collapse = n_rep_before,
                            accessions_excluded_nonhomologous = n_excluded,
                            accessions_collapsed = sum(col$map$collapsed),
                            passes_threshold = passes, enters = passes,
                            stringsAsFactors = FALSE)
    if (passes) {
      Biostrings::writeXStringSet(Biostrings::DNAStringSet(col$aln), file.path(output_dir, paste0("LIB_", l, ".fasta")))
    }
    message(sprintf("Locus '%s': %d non-homologous removed, %d identical collapsed, %d species with replica (%s).",
                    l, n_excluded, sum(col$map$collapsed), ctx$species_with_replicate,
                    if (passes) "enters" else "does not enter"))
  }

  summary_tab <- do.call(rbind, rows)
  rownames(summary_tab) <- NULL
  excluded_tab <- if (length(excluded) > 0) do.call(rbind, excluded) else
    data.frame(locus = character(0), sid = character(0), species = character(0), description = character(0),
               reason = character(0), share_forward = numeric(0), share_reverse = numeric(0),
               in_curated_alignment = logical(0), stringsAsFactors = FALSE)
  rownames(excluded_tab) <- NULL
  collapsed_tab <- do.call(rbind, collapsed)
  rownames(collapsed_tab) <- NULL

  # Per-locus funnel (decision of 2026-09-22): every locus assembled or screened, with its count at each
  # stage. A locus that did not enter the screening keeps its row, with NA after that point.
  all_loci <- sort(union(unique(as.character(registry$locus)), as.character(screening$locus)))
  funnel_tab <- do.call(rbind, lapply(all_loci, function(l) {
    f_aln <- file.path(curated_dir, "alignments", paste0("ALN_masked_final_", l, ".fasta"))
    n_cur <- if (file.exists(f_aln)) sum(startsWith(readLines(f_aln), ">")) else NA_integer_
    s <- summary_tab[summary_tab$locus == l, , drop = FALSE]
    done <- nrow(s) == 1L
    data.frame(
      locus = l,
      accessions_assembled = sum(registry$locus == l),
      accessions_curated = as.integer(n_cur),
      enters_screening = l %in% loci,
      accessions_excluded_nonhomologous = if (done) s$accessions_excluded_nonhomologous else NA_integer_,
      accessions_collapsed = if (done) s$accessions_collapsed else NA_integer_,
      accessions_after_collapse = if (done) unname(n_after[l]) else NA_integer_,
      species_with_replicate = if (done) s$species_with_replicate else NA_integer_,
      enters = if (done) s$enters else FALSE,
      stringsAsFactors = FALSE
    )
  }))

  summary_tab <- .bc_flag_provisional(.bc_flag_paralog(summary_tab, paralog_loci))
  excluded_tab <- .bc_flag_provisional(.bc_flag_paralog(excluded_tab, paralog_loci))
  collapsed_tab <- .bc_flag_provisional(.bc_flag_paralog(collapsed_tab, paralog_loci))
  funnel_tab <- .bc_flag_provisional(.bc_flag_paralog(funnel_tab, paralog_loci))
  utils::write.csv(funnel_tab, file.path(output_dir, "TABLE_barcoding_funnel.csv"), row.names = FALSE)

  utils::write.csv(summary_tab, file.path(output_dir, "TABLE_barcoding_library_summary.csv"), row.names = FALSE)
  utils::write.csv(excluded_tab, file.path(output_dir, "TABLE_barcoding_excluded_nonhomologous.csv"), row.names = FALSE)
  utils::write.csv(collapsed_tab, file.path(output_dir, "TABLE_barcoding_collapsed_identical.csv"), row.names = FALSE)

  invisible(list(summary = summary_tab, excluded = excluded_tab, collapsed = collapsed_tab, funnel = funnel_tab))
}
