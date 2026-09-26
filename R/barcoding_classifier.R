# Classifier of the molecular diagnostic branch (11_barcoding/7_classifier/).
# Phase 5A of PhyloCactus 0.5.0. This file implements the classifier and computes no accuracy:
# the negative controls are Phase 5B and the real figures are Phase 6, in that order and for a
# reason (a classifier written while its accuracy is on screen ends up fitted to it).

#' Empty prediction row, so every path returns the same columns
#'
#' The two methods write the same columns, which is what makes them comparable column by column.
#' `nn_distance` and `margin` are the scores of the nearest neighbour; `confidence` is the score
#' of IdTaxa. Each method fills its own and leaves the other empty. The threshold is not applied
#' here: it is swept in Phase 6 over these scores, without running the classifier again.
#' @noRd
.bc_prediction_row <- function(state, species = NA_character_, genus = NA_character_,
                               candidates = NA_character_, distance = NA_real_, margin = NA_real_,
                               confidence = NA_real_, reason = NA_character_) {
  data.frame(state = as.integer(state), predicted_species = species, predicted_genus = genus,
             candidates = candidates, nn_distance = distance, margin = margin,
             confidence = confidence, reason = reason, stringsAsFactors = FALSE)
}

#' Three states from the structure of the tie at the minimum distance
#'
#' One species at the minimum distance is state 1. Several species of a single genus is state 2,
#' with the candidate set declared: the data do not tell those species apart, and saying which one
#' it is would be inventing it. A tie across genera is state 3. No threshold takes part in this
#' decision, which is why the threshold can be swept afterwards.
#'
#' The function refuses a query that sits in its own training set. The folds of Phase 3 never put it
#' there, and the checker proves it, but the classifier does not rely on that: the defect that
#' retired the section in v0.4.2 was exactly this one.
#'
#' `allow_resubstitution = TRUE` lifts that guard, and only CN3 uses it (validation plan, sec. 4).
#' It is an argument and not a silent case on purpose: resubstitution is the defect itself, measured
#' on purpose to report its size, and whoever runs it has to ask for it by name.
#' @noRd
.bc_classify_nn <- function(dmat, train_ids, test_id, species, allow_resubstitution = FALSE) {
  if (test_id %in% train_ids && !isTRUE(allow_resubstitution)) {
    stop("Query '", test_id, "' is in its own training set. ",
         "The classifier refuses to evaluate a sequence it was trained on.", call. = FALSE)
  }
  train_ids <- train_ids[train_ids %in% rownames(dmat)]
  if (length(train_ids) == 0L) {
    return(.bc_prediction_row(3L, reason = "no_training"))
  }
  d <- dmat[test_id, train_ids]
  names(d) <- train_ids
  d <- d[!is.na(d)]
  if (length(d) == 0L) {
    return(.bc_prediction_row(3L, reason = "no_comparable_positions"))
  }
  d_min <- min(d)
  tied <- sort(unique(unname(species[names(d)[d == d_min]])), method = "radix")
  candidates <- paste(tied, collapse = "|")
  genera <- unique(.bc_genus(tied))

  if (length(tied) == 1L) {
    others <- d[unname(species[names(d)]) != tied]
    margin <- if (length(others) == 0L) NA_real_ else min(others) - d_min
    .bc_prediction_row(1L, species = tied, genus = genera, candidates = candidates,
                       distance = d_min, margin = margin)
  } else if (length(genera) == 1L) {
    .bc_prediction_row(2L, genus = genera, candidates = candidates,
                       distance = d_min, margin = 0)
  } else {
    .bc_prediction_row(3L, candidates = candidates, distance = d_min, margin = 0,
                       reason = "tie_across_genera")
  }
}

#' Seed of one query, a pure function of the run seed and of where the query sits
#'
#' IdTaxa draws k-mers at random for its confidence, and the draw depends on the state of the
#' generator when `IdTaxa()` is called. Measured on 2026-09-26 (DECIPHER 3.9.4): without a seed two
#' calls disagree; with the same seed they agree; and several queries in one call do not give what
#' one call per query gives. So each query gets its own seed, derived from the seed of the run, the
#' locus, the scheme, the fold and the `sid`, and from nothing else: not from the other folds, not
#' from the chunk, not from the order of the run (decision K1 of BMM).
#'
#' The four fields are joined by a carriage return and hashed base 31 modulo the Mersenne prime
#' 2^31 - 1, in doubles, which are exact below 2^53; the result is combined with the run seed. It
#' touches no random state.
#' @noRd
.bc_query_seed <- function(seed, locus, scheme, fold, sid) {
  m <- 2147483647
  key <- paste(locus, scheme, as.character(fold), sid, sep = "\r")
  h <- 0
  for (b in as.numeric(utf8ToInt(enc2utf8(key)))) h <- (h * 31 + b) %% m
  s <- ((as.numeric(seed) %% m) * 1000003 + h) %% m
  as.integer(if (s == 0) 1 else s)
}

#' The session's random state, as a function that puts it back
#'
#' Used as `restore <- .bc_rng_state(); on.exit(restore(), add = TRUE)` by the steps that train or
#' classify with IdTaxa: LearnTaxa() and IdTaxa() draw from the generator, and a step of the package
#' leaves the session's generator as it found it.
#' @noRd
.bc_rng_state <- function() {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  old_kind <- RNGkind()
  function() {
    suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
    if (had) assign(".Random.seed", old_seed, envir = globalenv())
    else suppressWarnings(rm(list = ".Random.seed", envir = globalenv()))
    invisible(TRUE)
  }
}

#' Evaluate an expression with a declared generator and seed, and give the session its state back
#' @noRd
.bc_with_seed <- function(seed, expr) {
  restore <- .bc_rng_state()
  on.exit(restore(), add = TRUE)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
  set.seed(seed)
  force(expr)
}

#' Training of IdTaxa for one fold
#'
#' `DECIPHER::LearnTaxa()` with the genus and the species as the two ranks below `Root`. A fold is
#' learnt once and every query of the fold is classified against the same training (decision K2 of
#' BMM).
#'
#' LearnTaxa() draws from the generator: it classifies its own training sequences in rounds with
#' `sample()`, and where labels conflict, as GenBank's do, the training it returns depends on the
#' state of the generator. The first measurements of 2026-09-26 missed it because their labels were
#' clean; the pilot on real data found it (two identical runs, different tables). So the training
#' has its own seed, `train_seed`, from `.bc_query_seed(seed, locus, scheme, fold, "<training>")`.
#' @noRd
.bc_idtaxa_train <- function(train, train_labels, train_seed = NULL) {
  if (length(train) != length(train_labels)) {
    stop("train and train_labels have different lengths.", call. = FALSE)
  }
  genera <- .bc_genus(unname(train_labels))
  taxonomy_str <- paste0("Root;", genera, ";", unname(train_labels))
  dna <- Biostrings::DNAStringSet(toupper(unname(train)))
  names(dna) <- names(train)
  learn <- function() DECIPHER::LearnTaxa(train = dna, taxonomy = taxonomy_str, verbose = FALSE)
  if (is.null(train_seed)) learn() else .bc_with_seed(train_seed, learn())
}

#' One row of IdTaxa from the confidences of each rank, cut at a threshold
#'
#' IdTaxa runs at threshold 0, so the genus and the species it reaches and the confidence of each
#' are always written (`genus_idtaxa`, `species_idtaxa`, `genus_confidence`, `species_confidence`).
#' The three states then come from `threshold`, as IdTaxa itself would have cut them: species at or
#' above it is state 1; genus alone at or above it is state 2, with the species of that genus in the
#' training set as candidates; neither is state 3. Cutting afterwards gives exactly what a run at
#' that threshold gives (measured on 2026-09-26), which is why the threshold can be swept in step 9.
#' @noRd
.bc_idtaxa_row <- function(genus_idtaxa, species_idtaxa, genus_confidence, species_confidence,
                           threshold, train_labels) {
  if (!is.na(species_confidence) && species_confidence >= threshold) {
    r <- .bc_prediction_row(1L, species = species_idtaxa, genus = genus_idtaxa,
                            candidates = species_idtaxa, confidence = species_confidence)
  } else if (!is.na(genus_confidence) && genus_confidence >= threshold) {
    labels <- unique(unname(train_labels))
    candidates <- sort(labels[.bc_genus(labels) == genus_idtaxa], method = "radix")
    r <- .bc_prediction_row(2L, genus = genus_idtaxa, candidates = paste(candidates, collapse = "|"),
                            confidence = genus_confidence)
  } else {
    r <- .bc_prediction_row(3L, reason = "low_confidence")
  }
  cbind(r, data.frame(genus_idtaxa = genus_idtaxa, species_idtaxa = species_idtaxa,
                      genus_confidence = genus_confidence, species_confidence = species_confidence,
                      stringsAsFactors = FALSE))
}

#' The same three states, from IdTaxa, the classifier of real use
#'
#' `train` and `train_labels` are the training sequences and their species; `trained`, when given,
#' is their training from `.bc_idtaxa_train()`, learnt once per fold; without it the training is
#' learnt here, with `train_seed`. `query_seed` is the seed of
#' this query from `.bc_query_seed()`: it is set just before the call to `IdTaxa()`, and the
#' session's random state is given back afterwards. Without it the answer is not reproducible.
#' @noRd
.bc_classify_idtaxa <- function(train, train_labels, query, threshold = 60, processors = 1L,
                                trained = NULL, query_seed = NULL, train_seed = NULL) {
  if (is.null(trained)) trained <- .bc_idtaxa_train(train, train_labels, train_seed = train_seed)
  q <- Biostrings::DNAStringSet(toupper(unname(query)))
  call_idtaxa <- function() {
    DECIPHER::IdTaxa(q, trained, strand = "top", threshold = 0, processors = processors,
                     verbose = FALSE)
  }
  ids <- if (is.null(query_seed)) call_idtaxa() else .bc_with_seed(query_seed, call_idtaxa())
  taxon <- ids[[1]]$taxon
  conf <- ids[[1]]$confidence
  ok <- !grepl("^unclassified", taxon)
  taxon <- taxon[ok]
  conf <- conf[ok]
  .bc_idtaxa_row(genus_idtaxa = if (length(taxon) >= 2L) taxon[2] else NA_character_,
                 species_idtaxa = if (length(taxon) >= 3L) taxon[3] else NA_character_,
                 genus_confidence = if (length(conf) >= 2L) conf[2] else NA_real_,
                 species_confidence = if (length(conf) >= 3L) conf[3] else NA_real_,
                 threshold = threshold, train_labels = train_labels)
}

#' Training alignment of one fold: the training rows, and no column left with gaps only
#'
#' A column that only the queries of the fold had opened carries no base once they are out. It is
#' dropped because `MAFFT --add --keeplength` would drop it too, and the check of
#' `.bc_align_to_library()` that the library keeps its width would then stop the run.
#' @noRd
.bc_training_alignment <- function(dna, train_ids) {
  m <- as.matrix(dna)
  missing <- setdiff(train_ids, rownames(m))
  if (length(missing) > 0L) {
    stop("Training sequences not in the alignment: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  m <- m[train_ids, , drop = FALSE]
  m[, colSums(as.character(m) != "-") > 0L, drop = FALSE]
}

#' Nearest neighbour of one query measured the way the identification will measure it
#'
#' The query arrives unaligned. It is oriented against the training set with the rule of
#' `.normalise_strand()`, added to the training alignment with `MAFFT --add --keeplength`, one call
#' for this query alone, and classified with the three states of `.bc_classify_nn()` under E13.
#' CN2 and the add path of step 7 both go through here, so an alien query and a legitimate one are
#' measured by the same operation.
#'
#' A query that matches the training set in neither direction is not aligned and not classified: it
#' is state 3 with the reason `no_match`, no species and no distance (decision D4 of BMM,
#' 2026-09-26). A sequence with no homology to the locus has no nearest neighbour worth publishing.
#'
#' @param train_dna Aligned training set, `sid` as row names.
#' @param species Character vector of species named by `sid`; must cover every training row.
#' @param query Character. The query sequence; gaps are removed.
#' @param pool Strand reference pool of the training set, or `NULL` to build it here.
#' @return One-row data frame: `orientation` followed by the columns of `.bc_prediction_row()`.
#' @noRd
.bc_classify_by_add <- function(train_dna, species, query, model, min_comparable, pool = NULL,
                                mafft_exec = "mafft", mafft_opts = "--auto") {
  train_ids <- rownames(train_dna)
  if (is.null(pool)) pool <- .bc_strand_pool(train_dna)
  q <- .bc_dnabin_row(gsub("-", "", query, fixed = TRUE), "query")
  o <- .bc_orient_to_library(q, pool)
  if (o$orientation == "no_match") {
    r <- .bc_prediction_row(3L, reason = "no_match")
  } else {
    joined <- .bc_align_to_library(train_dna, o$query, mafft_exec = mafft_exec, mafft_opts = mafft_opts)
    dm <- .bc_classifier_matrix(joined, model, min_comparable)
    r <- .bc_classify_nn(dm, train_ids, "query", c(species[train_ids], query = "query"))
  }
  cbind(data.frame(orientation = o$orientation, stringsAsFactors = FALSE), r)
}

#' Distance matrix of a locus with the short pairs left without value
#' @noRd
.bc_classifier_matrix <- function(dna, model, min_comparable) {
  dmat <- .bc_distance_matrix(dna, model)
  comp <- .bc_comparable_sites(dna)
  dmat[comp < min_comparable] <- NA_real_
  diag(dmat) <- 0
  dmat
}

#' Folds of one scheme, rebuilt from the table of step 5
#'
#' The training set is the rest of the locus: in scheme E everything but the query, in scheme G
#' everything but the evaluated species, which is exactly what the queries of the fold are. Rebuilt
#' this way the hard rule cannot be broken by a mistake here.
#' @noRd
.bc_folds_from_table <- function(tab, sids_locus) {
  key <- paste(tab$locus, tab$fold, sep = "\r")
  lapply(sort(unique(key), method = "radix"), function(k) {
    d <- tab[key == k, , drop = FALSE]
    list(locus = d$locus[1], fold = d$fold[1], stratum = d$stratum[1],
         test_ids = d$sid, train_ids = setdiff(sids_locus, d$sid))
  })
}

#' Classify the Validation Folds of the Molecular Diagnostic Branch
#'
#' Step 7 of the branch. Reads the library of step 4 and the folds of step 5, classifies every query
#' of every fold and writes the raw predictions, one row per query. It aggregates nothing: no
#' accuracy, no rate, no mean. The negative controls (Phase 5B) and the real metrics (Phase 6) read
#' these tables afterwards.
#'
#' The output has three states (validation plan, sec. 1): species, genus with the specific ambiguity
#' declared, and not assignable. The state comes from the structure of the tie at the minimum
#' distance, not from a threshold, so the threshold can be swept later over `nn_distance`,
#' `margin` and `confidence` without classifying again.
#'
#' @param library_dir Character. Output directory of [finalize_barcoding_library()].
#' @param folds_dir Character. Output directory of [build_barcoding_folds()].
#' @param output_dir Character. Destination; refused inside a directory of the phylogeny.
#' @param method Character. `"nn"`, the branch classifier, or `"idtaxa"`, the contrast.
#' @param model Character. Distance model of [ape::dist.dna()], for `"nn"`.
#' @param min_comparable Integer. Minimum number of comparable positions for a pair to have a value
#'   (E13 of the validation plan). Pairs below it are left without value, never imputed.
#' @param schemes Character vector. `"species"` (scheme E) and `"genus"` (scheme G).
#' @param loci Character vector or `NULL` for every locus of the library.
#' @param max_folds Integer or `NULL`. Folds per locus and scheme; a declared random subset, for the
#'   contrast with IdTaxa, which retrains once per fold. `NULL` means all of them.
#' @param seed Integer. Seed of that subset, so it is reproducible, and, for `"idtaxa"`, the seed
#'   from which the seed of every query is derived with its locus, scheme, fold and `sid`, so that
#'   each answer is reproducible and independent of the rest of the run.
#' @param threshold Numeric. Confidence threshold applied to the output of [DECIPHER::IdTaxa()], for
#'   `"idtaxa"`. IdTaxa itself always runs at threshold 0 and the confidence of each rank is written
#'   (`genus_confidence`, `species_confidence`), so any other threshold can be applied afterwards
#'   with the same result as a run at that threshold.
#' @param chunk,n_chunks Integers, or `NULL` for a run without chunks. With both, only chunk
#'   `chunk` of `n_chunks` is classified: inside each locus and scheme the folds, in the order of
#'   step 5, are dealt round-robin to the chunks. The chunk writes its tables in `chunks/` of
#'   `output_dir`, with the suffix `_chunk<k>of<n>`, together with its session; the tables of a whole
#'   run are then written by [merge_barcoding_chunks()]. For a job array on a cluster, see
#'   [generate_barcoding_job_scripts()].
#' @param alignment Character. How the distance of a query is measured, for `"nn"`. `"library"`,
#'   the default, reads it from the joint alignment of the library, in which the query took part.
#'   `"add"` removes the query's gaps, orients it against the training set of its fold and adds it
#'   to that training alignment with `MAFFT --add --keeplength`, one call per query: the path of CN2
#'   and of a query brought by a user. The probe of 2026-09-25 found the two paths agreeing on 699 of
#'   720 queries and 4 answers changed. `"add"` writes its own tables, with the suffix `_add`, and
#'   never touches those of `"library"`.
#' @param mafft_exec,mafft_opts Command of the `MAFFT` binary and its options, for `alignment = "add"`.
#' @param notify Logical. Send an email when the run ends, whether it finished or failed. Intended
#'   for `method = "idtaxa"`, which retrains once per fold and is measured in hours. A notification
#'   that cannot be sent is reported and ignored: it never fails the run. See
#'   [send_run_notification()] for what has to be configured. Defaults to `FALSE`.
#' @param notify_to,notify_credentials Passed to [send_run_notification()] as `to` and
#'   `credentials`. Both default to `NULL`, which reads the `MY_EMAIL` and `PHYLOCACTUS_SMTP_CREDS`
#'   environment variables.
#' @return Invisibly, a list of prediction tables by scheme. Writes
#'   `TABLE_barcoding_predictions_<scheme>_<method>.csv` and `TABLE_barcoding_timing_<method>.csv`,
#'   one row per locus and scheme with the folds, the queries and the seconds they took; with
#'   `alignment = "add"` both names carry the suffix `_add` and the prediction tables add the columns
#'   `alignment` and `orientation`.
#' @examples
#' \dontrun{
#' classify_barcoding_folds(
#'   library_dir = "11_barcoding/4_library",
#'   folds_dir = "11_barcoding/5_folds",
#'   output_dir = "11_barcoding/7_classifier"
#' )
#' }
#' @export
classify_barcoding_folds <- function(library_dir = file.path("11_barcoding", "4_library"),
                                     folds_dir = file.path("11_barcoding", "5_folds"),
                                     output_dir = file.path("11_barcoding", "7_classifier"),
                                     method = c("nn", "idtaxa"),
                                     model = "raw",
                                     min_comparable = 100L,
                                     schemes = c("species", "genus"),
                                     loci = NULL,
                                     max_folds = NULL,
                                     seed = 1L,
                                     threshold = 60,
                                     chunk = NULL,
                                     n_chunks = NULL,
                                     alignment = c("library", "add"),
                                     mafft_exec = "mafft",
                                     mafft_opts = "--auto",
                                     notify = FALSE,
                                     notify_to = NULL,
                                     notify_credentials = NULL) {
  method <- match.arg(method)
  alignment <- match.arg(alignment)
  .bc_check_chunk(chunk, n_chunks)
  if (alignment == "add" && method != "nn") {
    stop("alignment = \"add\" measures the distance of the nearest neighbour and applies only to ",
         "method = \"nn\". IdTaxa does not align.", call. = FALSE)
  }
  # MAFFT is checked before anything is read or written, so a missing binary costs nothing
  if (alignment == "add") .bc_assert_mafft(mafft_exec)
  suffix <- paste0(method, if (alignment == "add") "_add" else "")
  chunk_tag <- if (is.null(chunk)) "" else paste0("_chunk", as.integer(chunk), "of", as.integer(n_chunks))
  started <- Sys.time()
  call_run <- function() {
    # Training and classification touch the generator; the session gets its state back
    restore_rng <- .bc_rng_state()
    on.exit(restore_rng(), add = TRUE)
    .bc_assert_output_dir(output_dir)
    lib <- .bc_library_from_dir(library_dir)
    for (sc in schemes) {
      f <- file.path(folds_dir, paste0("TABLE_barcoding_folds_", sc, ".csv"))
      if (!file.exists(f)) {
        stop("No ", basename(f), " in ", folds_dir, ". Run build_barcoding_folds() first.", call. = FALSE)
      }
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write_dir <- if (is.null(chunk)) output_dir else file.path(output_dir, "chunks")
    dir.create(write_dir, recursive = TRUE, showWarnings = FALSE)
    loci_todos <- sort(unique(lib$locus), method = "radix")
    loci <- if (is.null(loci)) loci_todos else intersect(loci_todos, loci)

    out <- list()
    timings <- list()
    for (sc in schemes) {
      tab <- utils::read.csv(file.path(folds_dir, paste0("TABLE_barcoding_folds_", sc, ".csv")),
                             stringsAsFactors = FALSE)
      rows <- list()
      for (l in loci) {
        d <- lib[lib$locus == l, , drop = FALSE]
        species <- stats::setNames(d$species, d$sid)
        folds <- .bc_folds_from_table(tab[tab$locus == l, , drop = FALSE], d$sid)
        if (length(folds) == 0L) next
        folds <- .bc_sample_folds(folds, max_folds, seed)
        if (!is.null(chunk)) folds <- folds[.bc_chunk_assign(length(folds), n_chunks) == chunk]
        if (length(folds) == 0L) next

        dna <- ape::read.dna(file.path(library_dir, paste0("LIB_", l, ".fasta")), format = "fasta",
                             as.matrix = TRUE)
        rownames(dna) <- .bc_parse_header(labels(dna))$sid
        dmat <- if (method == "nn" && alignment == "library") {
          .bc_classifier_matrix(dna, model, min_comparable)
        } else NULL
        seqs <- if (method == "idtaxa") {
          s <- stats::setNames(toupper(apply(as.character(as.matrix(dna)), 1, paste, collapse = "")),
                               rownames(dna))
          gsub("-", "", s, fixed = TRUE)
        } else NULL

        t0 <- Sys.time()
        for (p in folds) {
          if (alignment == "add") {
            train_aln <- .bc_training_alignment(dna, p$train_ids)
            pool <- .bc_strand_pool(train_aln)
          }
          # IdTaxa learns the fold once (decision K2); LearnTaxa() draws from the generator, so the
          # training has a seed of its own
          trained <- if (method == "idtaxa") {
            .bc_idtaxa_train(seqs[p$train_ids], species[p$train_ids],
                             train_seed = .bc_query_seed(seed, l, sc, p$fold, "<training>"))
          } else NULL
          for (q in p$test_ids) {
            if (alignment == "add") {
              s <- paste(as.character(dna[q, ]), collapse = "")
              r <- .bc_classify_by_add(train_aln, species, s, model, min_comparable, pool = pool,
                                       mafft_exec = mafft_exec, mafft_opts = mafft_opts)
              rows[[length(rows) + 1L]] <- cbind(
                data.frame(locus = l, scheme = sc, fold = p$fold, stratum = p$stratum, sid = q,
                           true_species = unname(species[q]),
                           true_genus = .bc_genus(unname(species[q])),
                           method = method, alignment = "add", stringsAsFactors = FALSE),
                r)
              next
            }
            r <- if (method == "nn") {
              .bc_classify_nn(dmat, p$train_ids, q, species)
            } else {
              .bc_classify_idtaxa(seqs[p$train_ids], species[p$train_ids], seqs[[q]],
                                  threshold = threshold, trained = trained,
                                  query_seed = .bc_query_seed(seed, l, sc, p$fold, q))
            }
            rows[[length(rows) + 1L]] <- cbind(
              data.frame(locus = l, scheme = sc, fold = p$fold, stratum = p$stratum, sid = q,
                         true_species = unname(species[q]),
                         true_genus = .bc_genus(unname(species[q])),
                         method = method, stringsAsFactors = FALSE),
              r)
          }
        }
        seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        n_queries <- sum(vapply(folds, function(p) length(p$test_ids), integer(1)))
        timings[[length(timings) + 1L]] <- data.frame(
          locus = l, scheme = sc, method = method, alignment = alignment,
          folds = length(folds), queries = n_queries, seconds = seconds,
          stringsAsFactors = FALSE)
        message(sprintf("Locus '%s', scheme %s, method %s, alignment %s: %d folds, %d predictions written in %.1f s.",
                        l, sc, method, alignment, length(folds), n_queries, seconds))
      }
      pred <- if (length(rows) > 0L) do.call(rbind, rows) else NULL
      if (!is.null(pred)) rownames(pred) <- NULL
      out[[sc]] <- pred
      utils::write.csv(if (is.null(pred)) data.frame() else pred,
                       file.path(write_dir, paste0("TABLE_barcoding_predictions_", sc, "_", suffix,
                                                   chunk_tag, ".csv")),
                       row.names = FALSE)
    }
    # The running time is written, not left in prose: the projection of the complete IdTaxa run
    # (about 34 h) and of the add path (8.5 h) rested on rates noted by hand.
    utils::write.csv(if (length(timings) > 0L) do.call(rbind, timings) else data.frame(),
                     file.path(write_dir, paste0("TABLE_barcoding_timing_", suffix, chunk_tag, ".csv")),
                     row.names = FALSE)
    if (!is.null(chunk)) {
      # The session of every chunk is kept next to its tables: on a cluster each task may run on
      # another node, and the versions it ran with are part of the result
      writeLines(c(paste0("Chunk ", chunk, " of ", n_chunks, ", method ", suffix, ", seed ", seed,
                          ", finished ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                          " on ", Sys.info()[["nodename"]]),
                   utils::capture.output(utils::sessionInfo())),
                 file.path(write_dir, paste0("SESSION_barcoding_", suffix, chunk_tag, ".txt")))
    }
    .bc_classifier_banner(out, write_dir, paste0(suffix, chunk_tag))
    invisible(out)
  }

  if (!isTRUE(notify)) return(call_run())

  result <- tryCatch(call_run(), error = function(e) e)

  if (inherits(result, "error")) {
    note <- .compose_run_notification(
      analysis = paste0("barcoding classification of the folds, method ", method),
      status = "failed",
      started = started,
      outputs = c("Classifier output directory" = output_dir),
      error_message = conditionMessage(result)
    )
    send_run_notification(note$subject, note$body,
                          to = notify_to, credentials = notify_credentials)
    # Re-raised as the original condition, so the caller sees the same error it would have seen
    # with notify = FALSE, with its class and its call intact.
    stop(result)
  }

  note <- .compose_run_notification(
    analysis = paste0("barcoding classification of the folds, method ", method),
    status = "finished",
    started = started,
    outputs = c(stats::setNames(
      file.path(output_dir, paste0("TABLE_barcoding_predictions_", names(result), "_", suffix, ".csv")),
      paste0("Predictions, scheme ", names(result))),
      "Classifier output directory" = output_dir)
  )
  send_run_notification(note$subject, note$body,
                        to = notify_to, credentials = notify_credentials)
  invisible(result)
}

#' Closing banner of step 7, in the format the phylogeny already uses
#'
#' It reports what was written and where, and nothing else. No accuracy, no rate: Phase 5A does not
#' measure any, and a banner is exactly where a number nobody asked for would slip in. `method` is the
#' suffix of the tables, `nn`, `nn_add` or `idtaxa`.
#' @noRd
.bc_classifier_banner <- function(out, output_dir, method) {
  cat("\n====================================================\n")
  cat("  Barcoding Classification Complete \U0001f335\n")
  cat("====================================================\n")
  cat("  Method:                ", method, "\n")
  for (sc in names(out)) {
    cat("  Scheme ", sc, ": ", NROW(out[[sc]]), " predictions in ",
        file.path(output_dir, paste0("TABLE_barcoding_predictions_", sc, "_", method, ".csv")),
        "\n", sep = "")
  }
  cat("  Output directory:      ", output_dir, "\n")
  cat("====================================================\n\n")
  invisible(TRUE)
}

#' A declared subset of folds, reproducible and without touching the session's random state
#' @noRd
.bc_sample_folds <- function(folds, max_folds, seed) {
  if (is.null(max_folds) || length(folds) <= max_folds) return(folds)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(list = ".Random.seed", envir = globalenv())), add = TRUE)
  }
  set.seed(seed)
  folds[sort(sample.int(length(folds), as.integer(max_folds)))]
}
