# ------------------------------------------------------------------------------
# Molecular diagnostic branch (11_barcoding/), curation and screening steps. PhyloCactus 0.5.0.
#
# Curation calls the phylogeny's alignment pipeline, run_alignment_pipeline(), unchanged: it works
# sequence by sequence, does not collapse by species and keeps the Genus_species|sid headers.
# Screening recomputes the per-locus context columns on the curated alignments, applies the branch's
# replication threshold there, flags saturation without excluding, and counts identical sequences.
# ------------------------------------------------------------------------------

#' Read a table written by an earlier step of the branch, naming that step when the file is absent
#' @noRd
.bc_read_step_table <- function(path, step_fn) {
  if (!file.exists(path)) {
    stop("Missing ", path, ". Run ", step_fn, "() first, or point the argument at its output directory.",
         call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

#' Accession registry of step 1, read from its directory
#' @noRd
.bc_read_registry <- function(assembly_dir) {
  .bc_read_step_table(file.path(assembly_dir, "TABLE_barcoding_accession_registry.csv"), "assemble_barcoding_dataset")
}

#' Paralogy surveillance: add `posible_paralogo` to a table with a `locus` column
#'
#' Phase 2, item 5 of the 0.5.0 plan: the loci named in `paralog_loci` (by default `pepC_like`) are
#' flagged in every output of the branch. The flag informs and never excludes.
#' @noRd
.bc_flag_paralog <- function(tab, paralog_loci) {
  if (is.null(tab)) return(tab)
  tab$posible_paralogo <- as.character(tab$locus) %in% paralog_loci
  tab
}

#' Replication threshold, evaluated on the curated accessions of one locus
#' @param assembled Registry rows of one locus (`sid`, `species`).
#' @param curated_sids Sids that survived curation.
#' @noRd
.bc_screen_threshold <- function(assembled, curated_sids, min_species_with_replica = 20L) {
  n_rep <- function(d) sum(table(d$species) >= 2L)
  curated <- assembled[assembled$sid %in% curated_sids, , drop = FALSE]
  data.frame(
    especies_con_replica_ensamblado = n_rep(assembled),
    especies_con_replica_curado = n_rep(curated),
    pasa_umbral = n_rep(curated) >= min_species_with_replica,
    stringsAsFactors = FALSE
  )
}

#' Screening decision: the threshold decides, saturation only flags (decision D3 of 2026-09-21)
#' @noRd
.bc_screen_decision <- function(screen) {
  screen$entra <- screen$pasa_umbral %in% TRUE
  screen$marca_saturacion <- screen$saturado %in% TRUE
  screen
}

#' Accessions and distinct aligned sequences, per locus and species
#' @param aln_list Named list (by locus) of character vectors of aligned sequences named by header.
#' @noRd
.bc_identical_sequences <- function(aln_list) {
  do.call(rbind, lapply(names(aln_list), function(l) {
    a <- aln_list[[l]]
    sp <- .bc_parse_header(names(a))$species
    s <- toupper(unname(a))
    do.call(rbind, lapply(sort(unique(sp)), function(x) {
      data.frame(locus = l, species = x,
                 accesiones = sum(sp == x),
                 secuencias_distintas = length(unique(s[sp == x])),
                 stringsAsFactors = FALSE)
    }))
  }))
}

#' Curate the Loci of the Molecular Diagnostic Branch
#'
#' Aligns and curates the per-locus FASTA files written by [assemble_barcoding_dataset()] with the
#' phylogeny's alignment pipeline, [run_alignment_pipeline()], called unchanged and with the
#' parameters the phylogeny uses for the ingroup. Only the loci named in `loci` are processed. The
#' pipeline works sequence by sequence and keeps the `Genus_species|sid` headers, so every accession
#' remains identified after curation.
#'
#' @param input_dir Character. Directory with one FASTA per locus (`11_barcoding/1_assembly`).
#' @param output_dir Character. Destination; refused inside a directory of the phylogeny.
#' @param loci Character vector or `NULL`. Loci to curate, as FASTA file names without extension.
#'   `NULL` (the default) curates every FASTA in `input_dir`: the loci come out of the data, not out of
#'   a list (decision of 2026-09-22).
#' @inheritParams run_alignment_pipeline
#' @return Invisibly, the summary table returned by [run_alignment_pipeline()].
#' @examples
#' \dontrun{
#' curate_barcoding_markers(
#'   input_dir = "11_barcoding/1_assembly",
#'   output_dir = "11_barcoding/2_curated"
#' )
#' }
#' @export
curate_barcoding_markers <- function(input_dir = file.path("11_barcoding", "1_assembly"),
                                     output_dir = file.path("11_barcoding", "2_curated"),
                                     loci = NULL,
                                     mask_alignment_regions = TRUE,
                                     min_non_gap_fraction = 0.30,
                                     max_missing_fraction = 0.30,
                                     min_masked_alignment_length = 100L,
                                     preserve_iupac = TRUE,
                                     fix_strand = TRUE,
                                     mafft_exec = "mafft",
                                     mafft_opts = "--auto") {
  .bc_assert_output_dir(output_dir)
  if (is.null(loci)) {
    loci <- sub("\\.fasta$", "", list.files(input_dir, pattern = "\\.fasta$"))
    if (length(loci) == 0) {
      stop("No FASTA in ", input_dir, ". Run assemble_barcoding_dataset() first.", call. = FALSE)
    }
    message("Loci to curate (every assembled locus): ", paste(loci, collapse = ", "))
  }
  files <- file.path(input_dir, paste0(loci, ".fasta"))
  if (any(!file.exists(files))) {
    stop("No FASTA for locus: ", paste(loci[!file.exists(files)], collapse = ", "), call. = FALSE)
  }
  esc <- gsub("([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1", loci)
  fasta_pattern <- paste0("^(", paste(esc, collapse = "|"), ")\\.fasta$")

  manifest <- run_alignment_pipeline(
    input_folder = input_dir,
    output_dir = output_dir,
    fasta_pattern = fasta_pattern,
    mask_alignment_regions = mask_alignment_regions,
    min_non_gap_fraction = min_non_gap_fraction,
    max_missing_fraction = max_missing_fraction,
    min_masked_alignment_length = min_masked_alignment_length,
    preserve_iupac = preserve_iupac,
    fix_strand = fix_strand,
    mafft_exec = mafft_exec,
    mafft_opts = mafft_opts
  )
  invisible(manifest)
}

#' Screen the Curated Loci of the Molecular Diagnostic Branch
#'
#' For each curated locus, recomputes the five context columns of the validation plan on the
#' accessions that survived curation, applies the branch's replication threshold to them, flags
#' saturation (which never excludes a locus), and counts identical aligned sequences per species.
#' Inputs are read from the directories of the earlier steps, so the screening runs on its own
#' whenever those directories exist.
#'
#' @param assembly_dir Character. Output directory of [assemble_barcoding_dataset()]; the accession
#'   registry `TABLE_barcoding_accession_registry.csv` is read from it.
#' @param curated_dir Character. Output directory of [curate_barcoding_markers()].
#' @param loci Character vector or `NULL`. Loci to screen. `NULL` (the default) screens every curated
#'   alignment in `curated_dir`; the threshold then decides which loci enter (decision of 2026-09-22).
#' @param output_dir Character. Destination of the screening tables.
#' @param min_species_with_replica Integer. A locus enters the branch when at least this number of
#'   species keep two or more accessions after curation. Defaults to `20L` (decision of 2026-09-21).
#' @param saturation_flag_cutoff Numeric. Slope threshold of the saturation proxy, as in
#'   [run_marker_screening()]. Defaults to `0.3`.
#' @param paralog_loci Character vector. Loci under paralogy surveillance, flagged in the column
#'   `posible_paralogo` of every table; the flag never excludes a locus. Defaults to `"pepC_like"`.
#' @return Invisibly, a list with `screening` and `identical`. Writes
#'   `TABLE_barcoding_marker_screening.csv` and `TABLE_barcoding_identical_sequences.csv`.
#' @export
screen_barcoding_markers <- function(assembly_dir = file.path("11_barcoding", "1_assembly"),
                                     curated_dir = file.path("11_barcoding", "2_curated"),
                                     loci = NULL,
                                     output_dir = file.path("11_barcoding", "3_screening"),
                                     min_species_with_replica = 20L,
                                     saturation_flag_cutoff = 0.3,
                                     paralog_loci = "pepC_like") {
  .bc_assert_output_dir(output_dir)
  registry <- .bc_read_registry(assembly_dir)
  if (is.null(loci)) {
    loci <- sub("^ALN_masked_final_(.*)\\.fasta$", "\\1",
                list.files(file.path(curated_dir, "alignments"), pattern = "^ALN_masked_final_.*\\.fasta$"))
    if (length(loci) == 0) {
      stop("No curated alignment in ", file.path(curated_dir, "alignments"),
           ". Run curate_barcoding_markers() first.", call. = FALSE)
    }
    message("Loci to screen (every curated locus): ", paste(loci, collapse = ", "))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  aln_list <- list()
  rows <- lapply(loci, function(l) {
    f <- file.path(curated_dir, "alignments", paste0("ALN_masked_final_", l, ".fasta"))
    assembled <- registry[registry$locus == l, , drop = FALSE]
    curated_sids <- character(0)
    sat <- list(slope = NA_real_, saturated = NA, reason = "no_alignment")
    if (file.exists(f)) {
      aln <- Biostrings::readDNAStringSet(f)
      curated_sids <- .bc_parse_header(names(aln))$sid
      aln_list[[l]] <<- stats::setNames(as.character(aln), names(aln))
      sat <- .test_saturation_proxy(ape::read.dna(f, format = "fasta"), saturation_method = "corrected",
                                    saturation_flag_cutoff = saturation_flag_cutoff)
    }
    thr <- .bc_screen_threshold(assembled, curated_sids, min_species_with_replica)
    ctx <- .bc_locus_summary(assembled[assembled$sid %in% curated_sids, , drop = FALSE])
    if (is.null(ctx) || nrow(ctx) == 0) {
      ctx <- data.frame(locus = l, especies_totales = 0L, especies_con_replica = 0L, accesiones_totales = 0L,
                        generos_totales = 0L, generos_con_2_o_mas_especies = 0L, stringsAsFactors = FALSE)
    }
    data.frame(ctx, thr[, c("especies_con_replica_ensamblado", "pasa_umbral")],
               accesiones_ensamblado = nrow(assembled),
               saturado = sat$saturated, pendiente_saturacion = sat$slope, motivo_saturacion = sat$reason,
               stringsAsFactors = FALSE)
  })
  screening <- .bc_flag_provisional(.bc_flag_paralog(.bc_screen_decision(do.call(rbind, rows)), paralog_loci))
  identical <- if (length(aln_list) > 0) {
    .bc_flag_provisional(.bc_flag_paralog(.bc_identical_sequences(aln_list), paralog_loci))
  } else {
    NULL
  }

  utils::write.csv(screening, file.path(output_dir, "TABLE_barcoding_marker_screening.csv"), row.names = FALSE)
  if (!is.null(identical)) {
    utils::write.csv(identical, file.path(output_dir, "TABLE_barcoding_identical_sequences.csv"), row.names = FALSE)
  }
  invisible(list(screening = screening, identical = identical))
}
