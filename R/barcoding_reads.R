# ------------------------------------------------------------------------------
# Molecular diagnostic branch (11_barcoding/), identification of paired reads (route A, T4).
# PhyloCactus 0.5.0.
#
# fastp, GetOrganelle (plastome and nrDNA) and step 10 on the assembled paths, with the parameters of
# the route probe of Phase 8. The same function builds 0_genomic_raw/reads/<run>/ for F3 (M3).
#
# Governing documents: 01_minutes_and_phases/04_barcoding_0.5.0/ of the audit repository, route A
# proposal of 2026-09-28 (decisions RA1 to RA6 of BMM).
# ------------------------------------------------------------------------------

#' Default runner of the external steps
#' @noRd
.bc_reads_runner <- function(cmd, args, stdout = "", stderr = "") {
  system2(cmd, args, stdout = stdout, stderr = stderr)
}

#' Median read length over the first reads of a FASTQ file (plain or gzipped)
#' @noRd
.bc_fastq_median_length <- function(path, n = 10000L) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "r") else file(path, "r")
  on.exit(close(con), add = TRUE)
  x <- readLines(con, n = 4L * n, warn = FALSE)
  s <- x[seq(2L, length(x), by = 4L)]
  if (length(s) == 0L) return(NA_real_)
  stats::median(nchar(s))
}

#' Identify paired Illumina reads through route A
#'
#' Step 10 for genome-skim or WGS reads (input type T4, route A of the route probe of Phase 8).
#' The reads are cleaned with `fastp`; the plastome (`embplant_pt`) and the nrDNA (`embplant_nr`) are
#' assembled with GetOrganelle (`-R 15` and `-R 10`, `-k 21,45,65,85,105`); the path sequences of each
#' (the `graph1.1` path first) are identified with [identify_barcoding_query()], and for each locus the
#' row with the longest region is kept.
#'
#' Accepted: paired Illumina reads, two FASTQ files, plain or gzipped. Single-end reads stop (not
#' tested in 0.5.0), and so do long reads (median length above 600 bases over the first 10 000 reads;
#' not admitted in 0.5.0). Aligned reads (BAM, CRAM, SAM) stop with the command that converts them to
#' FASTQ. Target capture, RAD-seq and GBS, RNA-seq and mixed amplicons are rejected in
#' 0.5.0 (matrix of accepted inputs of 28-09) but cannot be told from the files; they should not be
#' given to this function.
#'
#' A GetOrganelle step that fails, or exits 0 without writing a path, is reported with a message and
#' its stderr is kept; the loci not found then have reason `assembly_failed`, not `no_overlap`, since
#' their absence was not observed. Before any step runs, the executables and the GetOrganelle seed and
#' label databases (`embplant_pt`, `embplant_mt`, `embplant_nr`) are checked.
#'
#' Writes to `output_dir/<run_name>/`: `fastp.json`, `plastome_path.fasta`, `nrdna_path.fasta`, the
#' `STDERR_<step>.txt` of each step, `MANIFEST.csv` (md5 of inputs and outputs, versions of the tools,
#' commands) and `TABLE_barcoding_identify_<run_name>.csv`. With `keep_intermediate = FALSE` the
#' cleaned reads, the fastp HTML report and the GetOrganelle work folders are deleted; the input reads
#' are never touched. The identification models are cached in `output_dir/models`.
#'
#' @param reads1,reads2 Character. The two FASTQ files of the run; `reads2 = NULL` (single-end) stops.
#' @param library_dir,metrics_dir,threshold_dir,output_dir,threshold,min_overlap,seed As in
#'   [identify_barcoding_query()].
#' @param run_name Character. Name of the run folder and of the table.
#' @param threads Integer. Threads of `fastp` and GetOrganelle.
#' @param keep_intermediate Logical. Keep the cleaned reads and the GetOrganelle work folders.
#' @param fastp,getorganelle Character. The executables.
#' @param getorg_path Character. The GetOrganelle folder of databases (`GETORG_PATH`, by default
#'   `~/.GetOrganelle`).
#' @param report Logical. Write the identification report ([report_barcoding_identification()]).
#'   The run record `RUN_<run_name>.csv` is written in any case.
#' @param runner Function `(cmd, args, stdout, stderr)` that runs an external step, as [system2()];
#'   the default runs the tools. The tests pass a fake one.
#' @return Invisibly, the table of [identify_barcoding_query()], one row per locus, with the columns
#'   `contig` and `organelle` (`plastome` or `nrdna`).
#' @seealso [identify_barcoding_query()], [identify_barcoding_assembly()].
#' @examples
#' \dontrun{
#' identify_barcoding_reads("SRR_1.fastq.gz", "SRR_2.fastq.gz", run_name = "SRR", threads = 8)
#' }
#' @export
identify_barcoding_reads <- function(reads1,
                                     reads2,
                                     library_dir = file.path("11_barcoding", "4_library"),
                                     metrics_dir = file.path("11_barcoding", "11_metrics"),
                                     threshold_dir = file.path("11_barcoding", "9_threshold"),
                                     output_dir = file.path("11_barcoding", "10_identify"),
                                     run_name = "reads",
                                     threads = 4L,
                                     keep_intermediate = FALSE,
                                     fastp = "fastp",
                                     getorganelle = "get_organelle_from_reads.py",
                                     getorg_path = Sys.getenv("GETORG_PATH", "~/.GetOrganelle"),
                                     runner = .bc_reads_runner,
                                     threshold = 60,
                                     min_overlap = 100L,
                                     seed = 1L,
                                     report = TRUE) {
  # RA1: input
  if (is.null(reads2)) {
    stop("Single-end reads are not tested in 0.5.0 (matrix of accepted inputs of 28-09): give the two ",
         "FASTQ files of a paired Illumina run.", call. = FALSE)
  }
  for (f in c(reads1, reads2)) .bc_stop_aligned_reads(f)
  for (f in c(reads1, reads2)) if (!file.exists(f)) stop("Reads not found: ", f, ".", call. = FALSE)
  ml <- .bc_fastq_median_length(reads1)
  if (!is.na(ml) && ml > 600) {
    stop("Median read length ", ml, " bases: long reads are not admitted in 0.5.0 (matrix of accepted ",
         "inputs of 28-09).", call. = FALSE)
  }
  # RA3: tools and databases before anything runs
  if (identical(runner, .bc_reads_runner)) {
    exe <- Sys.which(c(fastp, getorganelle))
    if (any(!nzchar(exe))) {
      stop("Not found on the PATH: ", paste(c(fastp, getorganelle)[!nzchar(exe)], collapse = ", "),
           ". Route A needs fastp and GetOrganelle (bioconda).", call. = FALSE)
    }
    fastp <- exe[[1]]
    getorganelle <- exe[[2]]
  }
  dbs <- c("embplant_pt", "embplant_mt", "embplant_nr")
  need <- file.path(path.expand(getorg_path), rep(c("SeedDatabase", "LabelDatabase"), each = 3L), paste0(dbs, ".fasta"))
  if (!all(file.exists(need))) {
    stop("GetOrganelle databases missing in ", getorg_path, ": ", paste(basename(dirname(need[!file.exists(need)])),
         basename(need[!file.exists(need)]), sep = "/", collapse = ", "),
         ". Install them with: get_organelle_config.py -a embplant_pt,embplant_mt,embplant_nr", call. = FALSE)
  }
  .bc_assert_output_dir(output_dir)
  lib <- .bc_library_from_dir(library_dir)
  loci <- sort(unique(lib$locus), method = "radix")
  restore_rng <- .bc_rng_state()
  on.exit(restore_rng(), add = TRUE)
  rdir <- file.path(output_dir, run_name)
  dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
  commands <- character(0)
  run_step <- function(step, cmd, args) {
    commands <<- c(commands, paste(basename(cmd), paste(args, collapse = " ")))
    st <- runner(cmd, args, stdout = FALSE, stderr = file.path(rdir, paste0("STDERR_", step, ".txt")))
    if (is.null(st)) 0L else as.integer(st)
  }

  # RA2: fastp
  cl <- file.path(rdir, c("clean_1.fq.gz", "clean_2.fq.gz"))
  st <- run_step("fastp", fastp, c("-i", reads1, "-I", reads2, "-o", cl[1], "-O", cl[2], "-w", as.integer(threads),
                                   "-j", file.path(rdir, "fastp.json"), "-h", file.path(rdir, "fastp.html")))
  if (st != 0L || !all(file.exists(cl))) stop("fastp failed (exit status ", st, "); see STDERR_fastp.txt.", call. = FALSE)

  # RA2 and RA3: GetOrganelle, plastome and nrDNA
  org <- c(embplant_pt = "plastome", embplant_nr = "nrdna")
  failed <- character(0)
  contigs <- character(0)
  for (db in names(org)) {
    go <- file.path(rdir, paste0("getorganelle_", db))
    st <- run_step(paste0("getorganelle_", db), getorganelle,
                   c("-1", cl[1], "-2", cl[2], "-F", db, "-o", go, "-t", as.integer(threads),
                     "-R", if (db == "embplant_pt") 15L else 10L, "-k", "21,45,65,85,105", "--overwrite"))
    pf <- list.files(go, pattern = "path_sequence\\.fasta$", full.names = TRUE)
    pf <- pf[order(!grepl("graph1\\.1\\.", pf), pf, method = "radix")]
    if (st != 0L || length(pf) == 0L) {
      failed <- c(failed, db)
      message("GetOrganelle ", db, ": ", if (st != 0L) paste0("exit status ", st) else "exit 0 but no path written",
              "; its loci not found are reported as assembly_failed (see STDERR_getorganelle_", db, ".txt).")
      next
    }
    x <- Biostrings::readDNAStringSet(pf[1])
    names(x) <- paste0(org[[db]], "_", seq_along(x))
    Biostrings::writeXStringSet(x, file.path(rdir, paste0(org[[db]], "_path.fasta")))
    contigs <- c(contigs, stats::setNames(as.character(x), names(x)))
  }

  # Step 10 on the contigs; per locus, the row with the longest region
  if (length(contigs)) {
    tab <- .bc_identify_table(contigs, loci, lib, library_dir, metrics_dir, output_dir, threshold, min_overlap, seed,
                              threshold_dir = threshold_dir)
  } else {
    tab <- .bc_identify_table(stats::setNames("", paste0(run_name, "__none")), loci, lib, library_dir, metrics_dir,
                              output_dir, threshold, min_overlap, seed, threshold_dir = threshold_dir)
  }
  tab$contig <- ifelse(tab$reason %in% "no_overlap", NA_character_, tab$query)
  tab$organelle <- ifelse(is.na(tab$contig), NA_character_, sub("_[0-9]+$", "", tab$contig))
  len <- ifelse(is.na(tab$region_length), -1L, tab$region_length)
  tab <- tab[order(tab$locus, is.na(tab$contig), -len, tab$query, method = "radix"), , drop = FALSE]
  tab <- tab[!duplicated(tab$locus), , drop = FALSE]
  if (length(failed)) {
    miss <- is.na(tab$contig)
    tab$reason[miss] <- "assembly_failed"
  }
  tab$query[is.na(tab$contig)] <- run_name
  rownames(tab) <- NULL

  # RA4: manifest and retention
  outs <- file.path(rdir, c("fastp.json", "plastome_path.fasta", "nrdna_path.fasta"))
  outs <- outs[file.exists(outs)]
  versions <- vapply(c(fastp, getorganelle), function(cmd) {
    v <- tryCatch(runner(cmd, "--version", stdout = TRUE, stderr = TRUE), error = function(e) NA_character_)
    paste(trimws(v[nzchar(trimws(v))]), collapse = " ")
  }, "")
  manifest <- rbind(
    data.frame(kind = "input", name = basename(c(reads1, reads2)), value = unname(tools::md5sum(c(reads1, reads2)))),
    data.frame(kind = "output", name = basename(outs), value = unname(tools::md5sum(outs))),
    data.frame(kind = "tool", name = basename(c(fastp, getorganelle)), value = unname(versions)),
    data.frame(kind = "command", name = seq_along(commands), value = commands))
  utils::write.csv(manifest, file.path(rdir, "MANIFEST.csv"), row.names = FALSE)
  if (!isTRUE(keep_intermediate)) {
    unlink(c(cl, file.path(rdir, "fastp.html"), file.path(rdir, paste0("getorganelle_", names(org)))), recursive = TRUE)
  }
  out_file <- file.path(rdir, paste0("TABLE_barcoding_identify_", run_name, ".csv"))
  utils::write.csv(tab, out_file, row.names = FALSE)
  .bc_write_run_record(rdir, run_name, input_path = paste(c(reads1, reads2), collapse = "; "),
                       input_md5 = paste(unname(tools::md5sum(c(reads1, reads2))), collapse = "; "), input_type = "reads",
                       route = c("fastp", "GetOrganelle embplant_pt", "GetOrganelle embplant_nr", "crop", "cut to locus core (J1)",
                                 "IdTaxa", "cut to genus core (J3b)", "longest region per locus"),
                       library_dir = library_dir, loci = loci, models_dir = file.path(output_dir, "models"),
                       threshold = threshold, seed = seed, min_overlap = min_overlap,
                       extra = manifest[manifest$kind != "input", , drop = FALSE])
  .bc_identify_banner(tab, out_file)
  if (isTRUE(report)) {
    report_barcoding_identification(run_name, results_dir = rdir, library_dir = library_dir, metrics_dir = metrics_dir)
  }
  invisible(tab)
}
