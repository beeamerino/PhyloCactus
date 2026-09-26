# Classifier of the molecular diagnostic branch (11_barcoding/7_classifier/).
# Phase 5A of PhyloCactus 0.5.0. This file implements the classifier and computes no accuracy:
# the negative controls are Phase 5B and the real figures are Phase 6, in that order and for a
# reason (a classifier written while its accuracy is on screen ends up fitted to it).

#' Empty prediction row, so every path returns the same columns
#'
#' The two methods write the same columns, which is what makes them comparable column by column.
#' `distancia_vecino` and `margen` are the scores of the nearest neighbour; `confianza` is the score
#' of IdTaxa. Each method fills its own and leaves the other empty. The threshold is not applied
#' here: it is swept in Phase 6 over these scores, without running the classifier again.
#' @noRd
.bc_prediction_row <- function(estado, especie = NA_character_, genero = NA_character_,
                               candidatas = NA_character_, distancia = NA_real_, margen = NA_real_,
                               confianza = NA_real_, motivo = NA_character_) {
  data.frame(estado = as.integer(estado), especie_predicha = especie, genero_predicho = genero,
             candidatas = candidatas, distancia_vecino = distancia, margen = margen,
             confianza = confianza, motivo = motivo, stringsAsFactors = FALSE)
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
    return(.bc_prediction_row(3L, motivo = "sin_entrenamiento"))
  }
  d <- dmat[test_id, train_ids]
  names(d) <- train_ids
  d <- d[!is.na(d)]
  if (length(d) == 0L) {
    return(.bc_prediction_row(3L, motivo = "sin_posiciones_comparables"))
  }
  d_min <- min(d)
  empatadas <- sort(unique(unname(species[names(d)[d == d_min]])), method = "radix")
  candidatas <- paste(empatadas, collapse = "|")
  generos <- unique(.bc_genus(empatadas))

  if (length(empatadas) == 1L) {
    otras <- d[unname(species[names(d)]) != empatadas]
    margen <- if (length(otras) == 0L) NA_real_ else min(otras) - d_min
    .bc_prediction_row(1L, especie = empatadas, genero = generos, candidatas = candidatas,
                       distancia = d_min, margen = margen)
  } else if (length(generos) == 1L) {
    .bc_prediction_row(2L, genero = generos, candidatas = candidatas,
                       distancia = d_min, margen = 0)
  } else {
    .bc_prediction_row(3L, candidatas = candidatas, distancia = d_min, margen = 0,
                       motivo = "empate_entre_generos")
  }
}

#' The same three states, from IdTaxa, as a contrast for the nearest neighbour
#'
#' `DECIPHER::LearnTaxa()` and `IdTaxa()` with the genus and the species as the two ranks below
#' `Root`. IdTaxa returns only the ranks above its own confidence threshold, so the depth it reaches
#' gives the state directly: species is state 1, genus alone is state 2 with the species of that
#' genus in the training set as candidates, and neither is state 3.
#'
#' It reports `confianza` and leaves the distances empty, and it is a contrast: the branch
#' classifier is the nearest neighbour (decision of BMM, 2026-09-22).
#' @noRd
.bc_classify_idtaxa <- function(train, train_labels, query, threshold = 60, processors = 1L) {
  if (length(train) != length(train_labels)) {
    stop("train and train_labels have different lengths.", call. = FALSE)
  }
  generos <- .bc_genus(unname(train_labels))
  taxonomia <- paste0("Root;", generos, ";", unname(train_labels))
  dna <- Biostrings::DNAStringSet(toupper(unname(train)))
  names(dna) <- names(train)
  entrenado <- DECIPHER::LearnTaxa(train = dna, taxonomy = taxonomia, verbose = FALSE)
  q <- Biostrings::DNAStringSet(toupper(unname(query)))
  ids <- DECIPHER::IdTaxa(q, entrenado, strand = "top", threshold = threshold,
                          processors = processors, verbose = FALSE)
  taxon <- ids[[1]]$taxon
  conf <- ids[[1]]$confidence
  ok <- !grepl("^unclassified", taxon)
  taxon <- taxon[ok]
  conf <- conf[ok]

  if (length(taxon) >= 3L) {
    .bc_prediction_row(1L, especie = taxon[3], genero = taxon[2], candidatas = taxon[3],
                       confianza = conf[3])
  } else if (length(taxon) == 2L) {
    candidatas <- sort(unique(unname(train_labels)[generos == taxon[2]]), method = "radix")
    .bc_prediction_row(2L, genero = taxon[2], candidatas = paste(candidatas, collapse = "|"),
                       confianza = conf[2])
  } else {
    .bc_prediction_row(3L, motivo = "confianza_insuficiente")
  }
}

#' Training alignment of one fold: the training rows, and no column left with gaps only
#'
#' A column that only the queries of the fold had opened carries no base once they are out. It is
#' dropped because `MAFFT --add --keeplength` would drop it too, and the check of
#' [.bc_align_to_library()] that the library keeps its width would then stop the run.
#' @noRd
.bc_training_alignment <- function(dna, train_ids) {
  m <- as.matrix(dna)
  faltan <- setdiff(train_ids, rownames(m))
  if (length(faltan) > 0L) {
    stop("Training sequences not in the alignment: ", paste(faltan, collapse = ", "), call. = FALSE)
  }
  m <- m[train_ids, , drop = FALSE]
  m[, colSums(as.character(m) != "-") > 0L, drop = FALSE]
}

#' Nearest neighbour of one query measured the way the identification will measure it
#'
#' The query arrives unaligned. It is oriented against the training set with the rule of
#' `.normalise_strand()`, added to the training alignment with `MAFFT --add --keeplength`, one call
#' for this query alone, and classified with the three states of [.bc_classify_nn()] under E13.
#' CN2 and the add path of step 7 both go through here, so an alien query and a legitimate one are
#' measured by the same operation.
#'
#' A query that matches the training set in neither direction is not aligned and not classified: it
#' is state 3 with the reason `sin_coincidencia`, no species and no distance (decision D4 of BMM,
#' 2026-09-26). A sequence with no homology to the locus has no nearest neighbour worth publishing.
#'
#' @param train_dna Aligned training set, `sid` as row names.
#' @param species Character vector of species named by `sid`; must cover every training row.
#' @param query Character. The query sequence; gaps are removed.
#' @param pool Strand reference pool of the training set, or `NULL` to build it here.
#' @return One-row data frame: `orientacion` followed by the columns of [.bc_prediction_row()].
#' @noRd
.bc_classify_by_add <- function(train_dna, species, query, model, min_comparable, pool = NULL,
                                mafft_exec = "mafft", mafft_opts = "--auto") {
  train_ids <- rownames(train_dna)
  if (is.null(pool)) pool <- .bc_strand_pool(train_dna)
  q <- .bc_dnabin_row(gsub("-", "", query, fixed = TRUE), "consulta")
  o <- .bc_orient_to_library(q, pool)
  if (o$orientacion == "sin_coincidencia") {
    r <- .bc_prediction_row(3L, motivo = "sin_coincidencia")
  } else {
    junto <- .bc_align_to_library(train_dna, o$query, mafft_exec = mafft_exec, mafft_opts = mafft_opts)
    dm <- .bc_classifier_matrix(junto, model, min_comparable)
    r <- .bc_classify_nn(dm, train_ids, "consulta", c(species[train_ids], consulta = "consulta"))
  }
  cbind(data.frame(orientacion = o$orientacion, stringsAsFactors = FALSE), r)
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
  key <- paste(tab$locus, tab$pliegue, sep = "\r")
  lapply(sort(unique(key), method = "radix"), function(k) {
    d <- tab[key == k, , drop = FALSE]
    list(locus = d$locus[1], pliegue = d$pliegue[1], estrato = d$estrato[1],
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
#' distance, not from a threshold, so the threshold can be swept later over `distancia_vecino`,
#' `margen` and `confianza` without classifying again.
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
#' @param seed Integer. Seed of that subset, so it is reproducible.
#' @param threshold Numeric. Confidence threshold of [DECIPHER::IdTaxa()], for `"idtaxa"`.
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
#'   `alineamiento` and `orientacion`.
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
                                     alignment = c("library", "add"),
                                     mafft_exec = "mafft",
                                     mafft_opts = "--auto",
                                     notify = FALSE,
                                     notify_to = NULL,
                                     notify_credentials = NULL) {
  method <- match.arg(method)
  alignment <- match.arg(alignment)
  if (alignment == "add" && method != "nn") {
    stop("alignment = \"add\" measures the distance of the nearest neighbour and applies only to ",
         "method = \"nn\". IdTaxa does not align.", call. = FALSE)
  }
  # MAFFT is checked before anything is read or written, so a missing binary costs nothing
  if (alignment == "add") .bc_assert_mafft(mafft_exec)
  sufijo <- paste0(method, if (alignment == "add") "_add" else "")
  started <- Sys.time()
  call_run <- function() {
    .bc_assert_output_dir(output_dir)
    lib <- .bc_library_from_dir(library_dir)
    for (sc in schemes) {
      f <- file.path(folds_dir, paste0("TABLE_barcoding_folds_", sc, ".csv"))
      if (!file.exists(f)) {
        stop("No ", basename(f), " in ", folds_dir, ". Run build_barcoding_folds() first.", call. = FALSE)
      }
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    loci_todos <- sort(unique(lib$locus), method = "radix")
    loci <- if (is.null(loci)) loci_todos else intersect(loci_todos, loci)

    out <- list()
    tiempos <- list()
    for (sc in schemes) {
      tab <- utils::read.csv(file.path(folds_dir, paste0("TABLE_barcoding_folds_", sc, ".csv")),
                             stringsAsFactors = FALSE)
      filas <- list()
      for (l in loci) {
        d <- lib[lib$locus == l, , drop = FALSE]
        species <- stats::setNames(d$species, d$sid)
        pliegues <- .bc_folds_from_table(tab[tab$locus == l, , drop = FALSE], d$sid)
        if (length(pliegues) == 0L) next
        pliegues <- .bc_sample_folds(pliegues, max_folds, seed)

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
        for (p in pliegues) {
          if (alignment == "add") {
            ent <- .bc_training_alignment(dna, p$train_ids)
            pool <- .bc_strand_pool(ent)
          }
          for (q in p$test_ids) {
            if (alignment == "add") {
              s <- paste(as.character(dna[q, ]), collapse = "")
              r <- .bc_classify_by_add(ent, species, s, model, min_comparable, pool = pool,
                                       mafft_exec = mafft_exec, mafft_opts = mafft_opts)
              filas[[length(filas) + 1L]] <- cbind(
                data.frame(locus = l, esquema = sc, pliegue = p$pliegue, estrato = p$estrato, sid = q,
                           especie_verdadera = unname(species[q]),
                           genero_verdadero = .bc_genus(unname(species[q])),
                           metodo = method, alineamiento = "add", stringsAsFactors = FALSE),
                r)
              next
            }
            r <- if (method == "nn") {
              .bc_classify_nn(dmat, p$train_ids, q, species)
            } else {
              .bc_classify_idtaxa(seqs[p$train_ids], species[p$train_ids], seqs[[q]],
                                  threshold = threshold)
            }
            filas[[length(filas) + 1L]] <- cbind(
              data.frame(locus = l, esquema = sc, pliegue = p$pliegue, estrato = p$estrato, sid = q,
                         especie_verdadera = unname(species[q]),
                         genero_verdadero = .bc_genus(unname(species[q])),
                         metodo = method, stringsAsFactors = FALSE),
              r)
          }
        }
        segundos <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        n_consultas <- sum(vapply(pliegues, function(p) length(p$test_ids), integer(1)))
        tiempos[[length(tiempos) + 1L]] <- data.frame(
          locus = l, esquema = sc, metodo = method, alineamiento = alignment,
          pliegues = length(pliegues), consultas = n_consultas, segundos = segundos,
          stringsAsFactors = FALSE)
        message(sprintf("Locus '%s', scheme %s, method %s, alignment %s: %d folds, %d predictions written in %.1f s.",
                        l, sc, method, alignment, length(pliegues), n_consultas, segundos))
      }
      pred <- do.call(rbind, filas)
      rownames(pred) <- NULL
      out[[sc]] <- pred
      utils::write.csv(pred, file.path(output_dir,
                                       paste0("TABLE_barcoding_predictions_", sc, "_", sufijo, ".csv")),
                       row.names = FALSE)
    }
    # The running time is written, not left in prose: the projection of the complete IdTaxa run
    # (about 34 h) and of the add path (8.5 h) rested on rates noted by hand.
    utils::write.csv(do.call(rbind, tiempos),
                     file.path(output_dir, paste0("TABLE_barcoding_timing_", sufijo, ".csv")),
                     row.names = FALSE)
    .bc_classifier_banner(out, output_dir, sufijo)
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
      file.path(output_dir, paste0("TABLE_barcoding_predictions_", names(result), "_", sufijo, ".csv")),
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
    cat("  Scheme ", sc, ": ", nrow(out[[sc]]), " predictions in ",
        file.path(output_dir, paste0("TABLE_barcoding_predictions_", sc, "_", method, ".csv")),
        "\n", sep = "")
  }
  cat("  Output directory:      ", output_dir, "\n")
  cat("====================================================\n\n")
  invisible(TRUE)
}

#' A declared subset of folds, reproducible and without touching the session's random state
#' @noRd
.bc_sample_folds <- function(pliegues, max_folds, seed) {
  if (is.null(max_folds) || length(pliegues) <= max_folds) return(pliegues)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(list = ".Random.seed", envir = globalenv())), add = TRUE)
  }
  set.seed(seed)
  pliegues[sort(sample.int(length(pliegues), as.integer(max_folds)))]
}
