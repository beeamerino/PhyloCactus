#' Own records of step 1, read from a folder of own data (decisions O1 to O3 of BMM, 29-09)
#'
#' The folder (`0_own/` of the tutorial) holds `SAMPLES.csv`, one row per sample (`sample`,
#' `species`, `voucher`, `data_type`, `published`, `notes`), and one FASTA per library locus, named
#' as the locus (`matK.fasta`), with headers `>sample` (text after the first space is ignored). A
#' record enters with the sid `own:<sample>:<locus>`, `source = "own"` and the sample as its
#' `genbank_name`; the species is resolved with the checklist (Korotkova *et al*. 2021) as GenBank
#' names are. A record whose sample is not in the sheet, whose species is not in the checklist or
#' whose locus is not a locus of the registry is left out and listed.
#' @param own_dir Character. The folder.
#' @param checklist_path As in [assemble_barcoding_dataset()].
#' @param loci Character. The loci of the registry.
#' @return A list with `registry` (rows in the columns of the registry, plus `source`), `sequences`
#'   (named by sid) and `left_out` (`sample`, `locus`, `reason`).
#' @noRd
.bc_own_records <- function(own_dir, checklist_path, loci) {
  sheet_file <- file.path(own_dir, "SAMPLES.csv")
  if (!file.exists(sheet_file)) stop("No SAMPLES.csv in ", own_dir, ".", call. = FALSE)
  sheet <- utils::read.csv(sheet_file, stringsAsFactors = FALSE, colClasses = "character")
  req <- c("sample", "species", "voucher", "data_type", "published", "notes")
  miss <- setdiff(req, names(sheet))
  if (length(miss)) stop("SAMPLES.csv lacks columns: ", paste(miss, collapse = ", "), call. = FALSE)
  sheet$sample <- trimws(sheet$sample)
  if (anyDuplicated(sheet$sample)) {
    stop("A sample appears more than once in SAMPLES.csv: ",
         paste(unique(sheet$sample[duplicated(sheet$sample)]), collapse = ", "), call. = FALSE)
  }
  sheet$resolved <- .bc_resolve_names(sheet$species, checklist_path)

  rows <- list()
  seqs <- character(0)
  left <- list()
  for (f in sort(list.files(own_dir, pattern = "\\.fasta$", full.names = TRUE))) {
    locus <- sub("\\.fasta$", "", basename(f))
    fa <- Biostrings::readDNAStringSet(f)
    smp <- sub("[[:space:]].*$", "", trimws(names(fa)))
    if (anyDuplicated(smp)) {
      stop("A sample appears more than once in ", basename(f), ": ",
           paste(unique(smp[duplicated(smp)]), collapse = ", "), call. = FALSE)
    }
    i <- match(smp, sheet$sample)
    reason <- ifelse(is.na(i), "sample_not_in_sheet",
              ifelse(rep(!locus %in% loci, length(smp)), "locus_not_in_library",
              ifelse(is.na(sheet$resolved[i]), "species_not_in_checklist", NA_character_)))
    if (any(!is.na(reason))) {
      left[[f]] <- data.frame(sample = smp[!is.na(reason)], locus = locus, reason = reason[!is.na(reason)],
                              stringsAsFactors = FALSE)
    }
    ok <- is.na(reason)
    if (!any(ok)) next
    sid <- paste0("own:", smp[ok], ":", locus)
    rows[[f]] <- data.frame(sid = sid, species = sheet$resolved[i[ok]], genus = .bc_genus(sheet$resolved[i[ok]]),
                            locus = locus, cluster_id = NA_integer_, genbank_name = smp[ok], source = "own",
                            stringsAsFactors = FALSE)
    seqs <- c(seqs, stats::setNames(as.character(fa)[ok], sid))
  }
  registry <- if (length(rows)) do.call(rbind, rows) else
    data.frame(sid = character(0), species = character(0), genus = character(0), locus = character(0),
               cluster_id = integer(0), genbank_name = character(0), source = character(0), stringsAsFactors = FALSE)
  rownames(registry) <- NULL
  left_out <- if (length(left)) do.call(rbind, left) else
    data.frame(sample = character(0), locus = character(0), reason = character(0), stringsAsFactors = FALSE)
  rownames(left_out) <- NULL
  list(registry = registry, sequences = seqs[registry$sid], left_out = left_out)
}
