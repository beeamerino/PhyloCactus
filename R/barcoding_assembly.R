# ------------------------------------------------------------------------------
# Molecular diagnostic branch (11_barcoding/), assembly step. PhyloCactus 0.5.0.
#
# Reads the same phylotaR workspace as the phylogeny (0_phylotaR_raw_Ingroup/), through the same
# function, and builds a library with every accession of every species: none of the three
# reductions to one sequence per species of the phylogeny is applied. Each sequence keeps its
# GenBank accession (sid) in the FASTA header, as Genus_species|sid.
#
# Governing documents: personal/review_log/phylocactus/04_barcoding_0.5.0/, validation plan of
# 2026-09-10 with its amendments of 2026-09-19 and 2026-09-21, and the Phase 2 design of 2026-09-21.
# ------------------------------------------------------------------------------

#' Refuse an output directory inside the phylogeny's directories
#' @noRd
.bc_assert_output_dir <- function(output_dir) {
  parts <- strsplit(normalizePath(output_dir, mustWork = FALSE), "[/\\\\]")[[1]]
  hit <- parts[grepl("^0_phylotaR_raw_|^1_phylotaR_out_|^4_Cleaned$", parts)]
  if (length(hit) > 0) {
    stop("output_dir is inside a directory of the phylogeny branch (", hit[1], "): ", output_dir,
         ". The molecular diagnostic branch writes only under its own directory.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Comparison key of clean_taxonomic_names(): spaces, hyphens and underscores collapsed
#' @noRd
.bc_checklist_key <- function(x) {
  x <- sub("\\|.*$", "", as.character(x))
  x <- gsub("[[:space:]_-]+", "_", trimws(x))
  gsub("^_+|_+$", "", x)
}

#' Checklist names, read as clean_taxonomic_names() reads them
#'
#' Every sheet of an Excel checklist except the `Facts` sheets, every rank, column `pureName`.
#' @noRd
.bc_read_checklist_names <- function(checklist_path) {
  if (grepl("\\.xlsx?$", checklist_path, ignore.case = TRUE)) {
    sheets <- readxl::excel_sheets(checklist_path)
    sheets <- sheets[!grepl("^facts", sheets, ignore.case = TRUE)]
    pure <- unlist(lapply(sheets, function(sh) {
      as.character(readxl::read_excel(checklist_path, sheet = sh)$pureName)
    }), use.names = FALSE)
  } else {
    pure <- as.character(utils::read.csv(checklist_path, stringsAsFactors = FALSE)$pureName)
  }
  unique(trimws(pure[!is.na(pure)]))
}

#' Accepted name for each GenBank name, with the rule of clean_taxonomic_names()
#'
#' Taxonomic IDs of the Outgroup of the Molecular Diagnostic Branch
#'
#' The six genus-level NCBI Taxonomy IDs that CN2 queries the library with (validation plan, sec. 4):
#' the three families that surround Cactaceae.
#'
#' - **Anacampserotaceae:** *Talinopsis* (107598), *Grahamia* (107617), *Anacampseros* (107583).
#' - **Portulacaceae:** *Portulaca* (3582).
#' - **Talinaceae:** *Talinum* (107600), *Talinella* (108056).
#'
#' *Talinopsis* is Anacampserotaceae, not Talinaceae, despite the name and the adjacent id.
#' *Amphipetalum* (1835425) is not requested: GenBank holds no nucleotide records for it, so the
#' request would only produce an empty download.
#'
#' **This is the branch's own copy** (decision of BMM, 2026-09-23). The phylogeny keeps its list in
#' `assemble_outgroup_phylotar()` and is not touched, so the branch runs without it. A test compares
#' the two and fails if they ever drift apart.
#'
#' @return Character vector of six taxonomic IDs.
#' @examples
#' barcoding_outgroup_taxids()
#' @export
barcoding_outgroup_taxids <- function() {
  c("107598", "107617", "107583", "3582", "107600", "108056")
}

#' What to mine, which is not the same as which parent resolves duplicates
#' @noRd
.bc_mining_taxids <- function(taxids, preferred_parent) {
  if (is.null(taxids)) preferred_parent else as.character(taxids)
}

#' Resolve the Names of a Set of Accessions, With or Without a Checklist
#'
#' Two behaviours, and the difference matters outside Cactaceae. With a checklist, the one the
#' branch has used since Phase 2: the checklist spelling for what matches and `NA` for what does
#' not. With `checklist_path = NULL`, no checklist at all: the GenBank names come back cleaned to
#' the `Genus_species` convention and nothing is dropped.
#'
#' The second is what CN2 needs (validation plan, sec. 4). An outgroup query is, by construction, a
#' name that a checklist of Cactaceae does not hold, so resolving it against that checklist would
#' discard the whole control. Added in Phase 5B; the assembly of the ingroup keeps its checklist and
#' its behaviour untouched.
#'
#' @param genbank_names Character vector of names as GenBank writes them.
#' @param checklist_path Path to a checklist, a character vector of accepted names, or `NULL` for
#'   no checklist.
#' @return Character vector, one entry per input name.
#' @noRd
.bc_resolve_names <- function(genbank_names, checklist_path = NULL) {
  if (is.null(checklist_path) || (length(checklist_path) == 1L && is.na(checklist_path))) {
    limpio <- gsub("[[:space:]]+", " ", trimws(as.character(genbank_names)))
    return(gsub(" ", "_", limpio, fixed = TRUE))
  }
  nombres <- if (length(checklist_path) == 1L && file.exists(checklist_path)) {
    .bc_read_checklist_names(checklist_path)
  } else {
    as.character(checklist_path)
  }
  .bc_match_names(genbank_names, nombres)
}

#' Returns the checklist spelling with underscores for a name that matches, and `NA` for one that
#' does not. Synonyms are not resolved: the checklist holds accepted names only.
#' @noRd
.bc_match_names <- function(genbank_names, checklist_names) {
  checklist_names <- unique(trimws(checklist_names[!is.na(checklist_names)]))
  lookup <- stats::setNames(checklist_names, .bc_checklist_key(checklist_names))
  lookup <- lookup[!duplicated(names(lookup))]
  hit <- unname(lookup[.bc_checklist_key(genbank_names)])
  gsub("[[:space:]]+", "_", hit)
}

#' FASTA header of the branch, Genus_species|sid
#' @noRd
.bc_header <- function(species, sid) paste0(species, "|", sid)

#' Species and sid from a branch header. Splits on "|" only, so hyphenated epithets survive.
#' @noRd
.bc_parse_header <- function(header) {
  data.frame(species = sub("\\|.*$", "", header),
             sid = sub("^[^|]*\\|", "", header),
             stringsAsFactors = FALSE)
}

#' Genus of an underscored binomial, ignoring a hybrid sign
#' @noRd
.bc_genus <- function(species) {
  tok <- strsplit(species, "_", fixed = TRUE)
  vapply(tok, function(t) {
    t <- t[!t %in% c("x", "X", "\u00d7")]
    if (length(t)) t[1] else NA_character_
  }, character(1))
}

#' Resolve sids present in more than one cluster, with the rule of assemble_ingroup_phylotar()
#'
#' A sid in some cluster whose parent is `preferred_parent` is kept only there; otherwise it is kept
#' in the cluster with the lowest identifier.
#' @param records Data frame with `cluster_id` and `sid`, one row per cluster membership.
#' @param parent_of Named character vector: parent taxon of each cluster, named by cluster id.
#' @return List with `kept` and `removed` memberships.
#' @noRd
.bc_resolve_cluster_duplicates <- function(records, parent_of, preferred_parent = "3593") {
  prnt <- unname(parent_of[as.character(records$cluster_id)])
  dup <- records$sid %in% records$sid[duplicated(records$sid)]
  keep <- rep(TRUE, nrow(records))
  if (any(dup)) {
    groups <- split(which(dup), records$sid[dup])
    for (ix in groups) {
      if (any(prnt[ix] == preferred_parent, na.rm = TRUE)) {
        keep[ix] <- prnt[ix] %in% preferred_parent
      } else {
        keep[ix] <- records$cluster_id[ix] == min(records$cluster_id[ix])
      }
    }
  }
  list(kept = records[keep, , drop = FALSE], removed = records[!keep, , drop = FALSE])
}

#' Apply the curated accession exclusions, by (cluster, sid) pair
#'
#' A row whose pair does not exist in `records` removes nothing and is reported in `missing_pairs`.
#' @param removed_auto Memberships already removed as cluster duplicates, to count the effective
#'   exclusions.
#' @noRd
.bc_apply_manual_exclusions <- function(records, exclusions, removed_auto = NULL) {
  rec_key <- paste(records$cluster_id, records$sid)
  exc_key <- unique(paste(as.integer(exclusions$cluster_id), exclusions$sid))
  existing <- exc_key[exc_key %in% rec_key]
  auto_key <- if (is.null(removed_auto)) character(0) else paste(removed_auto$cluster_id, removed_auto$sid)
  exc_cluster <- as.integer(sub(" .*$", "", exc_key))
  list(
    kept = records[!rec_key %in% exc_key, , drop = FALSE],
    excluded_keys = existing,
    missing_pairs = setdiff(exc_key, rec_key),
    counts = list(
      n_file = length(exc_key),
      n_in_clusters = sum(exc_cluster %in% records$cluster_id),
      n_existing = length(existing),
      n_effective = length(setdiff(existing, auto_key))
    )
  )
}

#' Accession registry: one row per (locus, sid), no sid in two loci
#' @noRd
.bc_build_registry <- function(records) {
  req <- c("cluster_id", "sid", "locus", "species", "nombre_genbank")
  miss <- setdiff(req, names(records))
  if (length(miss) > 0) stop("Registry input lacks columns: ", paste(miss, collapse = ", "), call. = FALSE)
  r <- records[!duplicated(paste(records$locus, records$sid)), req, drop = FALSE]
  multi <- unique(r$sid[duplicated(r$sid)])
  if (length(multi) > 0) {
    stop("A sid is assigned to more than one locus: ", paste(utils::head(multi, 10), collapse = ", "), call. = FALSE)
  }
  out <- data.frame(
    sid = r$sid,
    species = r$species,
    genus = .bc_genus(r$species),
    locus = r$locus,
    cluster_id = as.integer(r$cluster_id),
    nombre_genbank = r$nombre_genbank,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$locus, out$species, out$sid), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Locus name of every cluster, keeping the clusters whose gene is not recognised
#'
#' A cluster with no standardized marker name is kept under the provisional name `cluster_<ID>` and
#' listed with the GenBank description of its seed (decision of 2026-09-22).
#' @param clusters Data frame with `ID`, `Marker_std` and `Description` (seed description).
#' @return List with `map` (`ID`, `Marker_std`, `locus_provisional`) and `unnamed` (`ID`, `locus`,
#'   `descripcion_semilla`).
#' @noRd
.bc_name_clusters <- function(clusters) {
  id <- as.character(clusters$ID)
  marker <- trimws(as.character(clusters$Marker_std))
  provisional <- is.na(marker) | !nzchar(marker)
  marker[provisional] <- paste0("cluster_", id[provisional])
  desc <- if ("Description" %in% names(clusters)) as.character(clusters$Description) else rep(NA_character_, length(id))
  list(
    map = data.frame(ID = id, Marker_std = marker, locus_provisional = provisional, stringsAsFactors = FALSE),
    unnamed = data.frame(ID = id[provisional], locus = marker[provisional], descripcion_semilla = desc[provisional],
                         stringsAsFactors = FALSE)
  )
}

#' Flag the loci named provisionally after their cluster (`cluster_<ID>`)
#' @noRd
.bc_flag_provisional <- function(tab) {
  if (is.null(tab)) return(tab)
  tab$locus_provisional <- grepl("^cluster_[0-9]+$", as.character(tab$locus))
  tab
}

#' Cluster funnel of step 1: how many clusters survive each stage
#' @noRd
.bc_cluster_funnel <- function(n_workspace, n_selected, n_with_sequences, marker_map) {
  data.frame(
    etapa = c("clusteres_en_workspace", "clusteres_con_mas_de_min_species", "clusteres_con_secuencias_tras_filtros",
              "clusteres_con_nombre", "clusteres_sin_nombre_conservados", "loci_ensamblados"),
    n = as.integer(c(n_workspace, n_selected, n_with_sequences, sum(!marker_map$locus_provisional),
                     sum(marker_map$locus_provisional), length(unique(marker_map$Marker_std)))),
    stringsAsFactors = FALSE
  )
}

#' The five context columns of the validation plan (sec. 3), per locus
#'
#' A species has replication in a locus when it has two or more distinct accessions in that locus.
#' Nothing is counted across loci.
#' @noRd
.bc_locus_summary <- function(registry) {
  loci <- sort(unique(registry$locus))
  do.call(rbind, lapply(loci, function(l) {
    d <- registry[registry$locus == l, , drop = FALSE]
    n_by_sp <- table(d$species)
    sp_by_gen <- tapply(d$species, d$genus, function(x) length(unique(x)))
    data.frame(
      locus = l,
      especies_totales = length(n_by_sp),
      especies_con_replica = sum(n_by_sp >= 2L),
      accesiones_totales = nrow(d),
      generos_totales = length(sp_by_gen),
      generos_con_2_o_mas_especies = sum(sp_by_gen >= 2L),
      stringsAsFactors = FALSE
    )
  }))
}

#' Write one FASTA per locus, headers Genus_species|sid
#' @param sequences Character vector of sequences named by sid.
#' @noRd
.bc_write_locus_fastas <- function(registry, sequences, out_dir) {
  .bc_assert_output_dir(out_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  miss <- setdiff(registry$sid, names(sequences))
  if (length(miss) > 0) stop("No sequence for sid: ", paste(utils::head(miss, 10), collapse = ", "), call. = FALSE)
  paths <- character(0)
  for (l in sort(unique(registry$locus))) {
    d <- registry[registry$locus == l, , drop = FALSE]
    seqs <- Biostrings::DNAStringSet(toupper(unname(sequences[d$sid])))
    names(seqs) <- .bc_header(d$species, d$sid)
    f <- file.path(out_dir, paste0(l, ".fasta"))
    Biostrings::writeXStringSet(seqs, f)
    paths <- c(paths, f)
  }
  invisible(paths)
}

#' Seed metadata with a cache of the branch: only seeds missing from the cache are fetched
#'
#' A seed whose fetch failed is not cached, so the next run tries again.
#' @noRd
.bc_seed_metadata_cached <- function(seeds, cache_csv, fetcher = cp_fetch_seed_metadata, ...) {
  seeds <- as.character(seeds)
  cached <- if (file.exists(cache_csv)) {
    utils::read.csv(cache_csv, stringsAsFactors = FALSE, colClasses = "character")
  } else {
    data.frame(Seed = character(0), Species = character(0), Description = character(0), stringsAsFactors = FALSE)
  }
  missing <- setdiff(unique(seeds), cached$Seed)
  if (length(missing) > 0) {
    fetched <- as.data.frame(fetcher(missing, ...), stringsAsFactors = FALSE)[, c("Seed", "Species", "Description")]
    fetched <- fetched[!(is.na(fetched$Species) & is.na(fetched$Description)), , drop = FALSE]
    cached <- rbind(cached, fetched)
    dir.create(dirname(cache_csv), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(cached, cache_csv, row.names = FALSE)
  }
  idx <- match(seeds, cached$Seed)
  tibble::tibble(Seed = seeds, Species = cached$Species[idx], Description = cached$Description[idx])
}

#' Compare the branch's cluster-to-locus assignment with the phylogeny's
#' @noRd
.bc_compare_marker_maps <- function(branch, phylogeny) {
  b <- data.frame(ID = as.character(branch$ID), marker_rama = trimws(as.character(branch$Marker_std)),
                  stringsAsFactors = FALSE)
  p <- data.frame(ID = as.character(phylogeny$ID), marker_filogenia = trimws(as.character(phylogeny$Marker_std)),
                  stringsAsFactors = FALSE)
  m <- merge(b, p, by = "ID", all = TRUE)
  m$estado <- ifelse(is.na(m$marker_filogenia), "solo_rama",
              ifelse(is.na(m$marker_rama), "solo_filogenia",
              ifelse(m$marker_rama == m$marker_filogenia, "igual", "distinto")))
  m <- m[order(suppressWarnings(as.integer(m$ID)), m$ID), , drop = FALSE]
  rownames(m) <- NULL
  m
}

#' Assemble the Reference Library of the Molecular Diagnostic Branch
#'
#' Builds the multi-accession reference library of the molecular diagnostic section from the
#' `phylotaR` workspace shared with the phylogeny. The workspace is read, or created with the same
#' `phylotaR` call as [assemble_ingroup_phylotar()] when it does not exist, so either branch can run
#' first. Clusters are selected with the phylogeny's criterion (more than `min_species` species).
#' Sequences present in several clusters are resolved with the phylogeny's rule and the curated
#' accession exclusions are applied, but none of the reductions to one sequence per species is: every
#' accession of every species is kept. Each cluster is assigned to a locus with the same code as the
#' phylogeny, from GenBank metadata cached under `output_dir/cache/`. Names are reconciled against the
#' checklist with the rule of [clean_taxonomic_names()]; names that do not match are dropped and
#' reported. Synonyms are not resolved.
#'
#' @param wd_path Character. The `phylotaR` workspace, shared with the phylogeny (e.g.
#'   `0_phylotaR_raw_Ingroup`).
#' @param output_dir Character. Root of the branch. Refused if it lies inside a directory of the
#'   phylogeny. Defaults to `"11_barcoding"`.
#' @param target_genes_file,genes_map_file Character. Marker definitions; default to the files
#'   distributed with the package, as in [assemble_ingroup_phylotar()].
#' @param manual_exclusions_file Character. Curated accession exclusions; defaults to
#'   `inst/extdata/manual_exclusions_ingroup.csv`.
#' @param apply_manual_exclusions Logical. Apply the curated exclusions? Defaults to `TRUE`.
#' @param checklist_path Character. Accepted checklist (Korotkova *et al*. 2021); `NULL`, the default,
#'   uses the checklist distributed with the package. `NA` applies no checklist and keeps every
#'   GenBank name, cleaned to `Genus_species`: that is what the outgroup of CN2 needs, since a
#'   checklist of Cactaceae would discard every outgroup name.
#' @param min_species Integer. A cluster is retained when it holds more than this number of species.
#'   Defaults to `50`, as in the phylogeny.
#' @param preferred_parent Character. NCBI Taxonomy ID of the focal ingroup, used to decide which
#'   cluster keeps a sid that sits in two. Defaults to `"3593"`.
#' @param taxids Character vector or `NULL`. What to mine from GenBank when the workspace of
#'   `wd_path` does not exist yet, which is what lets the branch run without the phylogeny having
#'   run first. `NULL` mines `preferred_parent`, the behaviour since Phase 2. Several IDs are mined
#'   together, as [barcoding_outgroup_taxids()] returns them for CN2.
#' @param ncbi_dr Character or `NULL`. Directory of the BLAST+ binaries, used only when the workspace
#'   has to be created.
#' @param force_download Logical. Mine GenBank again even if the workspace exists. Defaults to `FALSE`.
#' @param phylogeny_map_file Character or `NULL`. Cluster-to-locus table of the phylogeny, compared
#'   with the branch's own when it exists. Defaults to
#'   `1_phylotaR_out_Ingroup/TABLE_CLUSTER_MARKER_ASSIGNMENT_INGROUP.csv` next to `wd_path`.
#' @return Invisibly, a list with `registry`, `summary`, `marker_map`, `comparison` and `counts`.
#'   Files written under `output_dir/1_assembly/`: one FASTA per locus,
#'   `TABLE_barcoding_accession_registry.csv`, `TABLE_barcoding_cluster_marker_assignment.csv`,
#'   `TABLE_barcoding_marker_assignment_comparison.csv` (when the phylogeny's table exists),
#'   `TABLE_barcoding_discarded_names.csv` and `TABLE_barcoding_replication_summary_assembly.csv`.
#' @references
#' Korotkova, N., Aquino, D., Arias, S., Eggli, U., Franck, A., Gómez-Hinostrosa, C., Guerrero, P. C.,
#' Hernández, H. M., Kohlbecker, A., Köhler, M., Luther, K., Majure, L. C., Müller, A., Metzing, D.,
#' Nyffeler, R., Sánchez, D., Schlumpberger, B. O., & Berendsohn, W. G. (2021). Cactaceae at
#' Caryophyllales.org - a dynamic online species-level taxonomic backbone for the family.
#' *Willdenowia*, 51(2), 251-270. \doi{10.3372/wi.51.51208}
#' @examples
#' \dontrun{
#' assemble_barcoding_dataset(
#'   wd_path = "0_phylotaR_raw_Ingroup",
#'   output_dir = "11_barcoding"
#' )
#' }
#' @export
assemble_barcoding_dataset <- function(wd_path,
                                       output_dir = "11_barcoding",
                                       target_genes_file = NULL,
                                       genes_map_file = NULL,
                                       manual_exclusions_file = NULL,
                                       apply_manual_exclusions = TRUE,
                                       checklist_path = NULL,
                                       min_species = 50,
                                       preferred_parent = "3593",
                                       taxids = NULL,
                                       ncbi_dr = NULL,
                                       force_download = FALSE,
                                       phylogeny_map_file = NULL) {
  .bc_assert_output_dir(output_dir)

  if (is.null(target_genes_file)) target_genes_file <- system.file("extdata", "target_genes.txt", package = "PhyloCactus")
  if (is.null(genes_map_file)) genes_map_file <- system.file("extdata", "genes_map.csv", package = "PhyloCactus")
  if (is.null(manual_exclusions_file)) manual_exclusions_file <- system.file("extdata", "manual_exclusions_ingroup.csv", package = "PhyloCactus")
  # NULL keeps the checklist of Cactaceae, which is the behaviour the closed phases depend on.
  # NA is the way to ask for no checklist at all, which is what the outgroup of CN2 needs.
  if (is.null(checklist_path)) checklist_path <- system.file("extdata", "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx", package = "PhyloCactus")
  if (is.null(phylogeny_map_file)) {
    phylogeny_map_file <- file.path(dirname(normalizePath(wd_path, mustWork = FALSE)), "1_phylotaR_out_Ingroup",
                                    "TABLE_CLUSTER_MARKER_ASSIGNMENT_INGROUP.csv")
  }
  ficheros <- c(target_genes_file, genes_map_file, checklist_path)
  ficheros <- ficheros[!is.na(ficheros)]
  for (f in ficheros) {
    if (!nzchar(f) || !file.exists(f)) stop("Input file not found: ", f, call. = FALSE)
  }

  dir_asm <- file.path(output_dir, "1_assembly")
  dir_cache <- file.path(output_dir, "cache")
  dir_logs <- file.path(output_dir, "logs")
  for (d in c(dir_asm, dir_cache, dir_logs)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  path_log <- file.path(dir_logs, "LOG_barcoding_assembly.txt")
  if (file.exists(path_log)) file.remove(path_log)
  log_message <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
    cat(msg, "\n")
    write(msg, file = path_log, append = TRUE)
  }
  log_message("Assembling the molecular diagnostic library (no reduction to one sequence per species).")

  # Marker definitions, read as in assemble_ingroup_phylotar()
  markers_df <- utils::read.table(target_genes_file, header = TRUE, stringsAsFactors = FALSE)
  markers <- unique(cp_normalize_marker(as.character(markers_df$markers)))
  pattern <- cp_build_pattern(markers)
  genes_map_df <- utils::read.csv(genes_map_file, stringsAsFactors = FALSE) |>
    dplyr::mutate(search = cp_normalize_marker(search), replace = trimws(replace)) |>
    dplyr::distinct(search, .keep_all = TRUE)
  gene_lookup <- genes_map_df |> dplyr::select(search_gene = search, Gene_std = replace)
  marker_lookup <- genes_map_df |> dplyr::select(search_marker = search, Marker_std = replace)

  # 1. Workspace, shared with the phylogeny
  phylota <- .phylotar_load_or_mine(wd_path, preferred_parent = preferred_parent,
                                    txid = .bc_mining_taxids(taxids, preferred_parent), ncbi_dr = ncbi_dr,
                                    force_download = force_download, log_message = log_message)
  log_message("Workspace: ", length(phylota@cids), " clusters, ", length(phylota@sids), " sequences in clusters.")

  # 2. Clusters, with the phylogeny's criterion. The species count is taken on the reduced object, as
  # in assemble_ingroup_phylotar(); the records are then taken from the unreduced clusters.
  species_reduced <- phylotaR::drop_by_rank(phylota, rnk = "species", n = 1)
  cluster_ids <- species_reduced@cids
  ntaxa <- phylotaR::get_ntaxa(species_reduced, cid = cluster_ids, rnk = "species")
  keep_clusters <- cluster_ids[ntaxa > min_species]
  selected <- phylotaR::drop_clstrs(phylota, cid = keep_clusters)
  log_message("Clusters with more than ", min_species, " species: ", length(selected@cids))

  records <- do.call(rbind, lapply(selected@cids, function(cid) {
    data.frame(cluster_id = as.integer(cid), sid = as.character(selected@clstrs[[cid]]@sids),
               stringsAsFactors = FALSE)
  }))
  log_message("Cluster memberships: ", nrow(records), "; distinct sids: ", length(unique(records$sid)))

  # 3. Accession-level filters of the phylogeny, unchanged
  parent_of <- vapply(selected@cids, function(cid) cp_get_cluster_parent(selected, cid), character(1))
  names(parent_of) <- selected@cids
  dd <- .bc_resolve_cluster_duplicates(records, parent_of, preferred_parent = preferred_parent)
  log_message("Redundant cluster memberships removed (not accessions): ", nrow(dd$removed))
  kept <- dd$kept
  counts_manual <- NULL
  if (isTRUE(apply_manual_exclusions) && nzchar(manual_exclusions_file) && file.exists(manual_exclusions_file)) {
    exclusions <- utils::read.csv(manual_exclusions_file, stringsAsFactors = FALSE)
    me <- .bc_apply_manual_exclusions(records, exclusions, removed_auto = dd$removed)
    kept <- kept[!paste(kept$cluster_id, kept$sid) %in% me$excluded_keys, , drop = FALSE]
    counts_manual <- me$counts
    log_message("Manual exclusions: ", me$counts$n_file, " in the file; ", me$counts$n_in_clusters,
                " in the selected clusters; ", me$counts$n_existing, " whose (cluster, sid) pair exists; ",
                me$counts$n_effective, " effective.")
    if (length(me$missing_pairs) > 0) {
      log_message("Exclusion rows whose pair does not exist: ", paste(me$missing_pairs, collapse = "; "))
    }
  } else if (!isTRUE(apply_manual_exclusions)) {
    log_message("Manual exclusions DISABLED (apply_manual_exclusions = FALSE).")
  }

  txids <- vapply(kept$sid, function(s) as.character(selected@sqs[[s]]@txid), character(1), USE.NAMES = FALSE)
  kept$nombre_genbank <- cp_clean_species_name(phylotaR::get_tx_slot(selected, txid = txids, slt_nm = "scnm"))

  # 4. Cluster-to-locus assignment, with the phylogeny's code and the branch's own caches
  metadata <- cp_download_all_metadata(unique(kept$sid),
                                       file.path(dir_cache, "CACHE_GENBANK_METADATA_BARCODING.csv"),
                                       batch_size = 200, sleep_time = 0.5, max_retries = 5,
                                       log_message = log_message)
  df_meta <- data.frame(cluster_id = kept$cluster_id, species = kept$nombre_genbank, sid = kept$sid,
                        stringsAsFactors = FALSE) |>
    dplyr::left_join(metadata, by = "sid") |>
    cp_annotate_marker_text(pattern = pattern, genes_map_df = genes_map_df)

  cleaned_clstrs <- selected@clstrs
  for (i in seq_along(cleaned_clstrs@clstrs)) {
    cid <- names(cleaned_clstrs@clstrs)[i]
    cleaned_clstrs@clstrs[[i]]@sids <- intersect(cleaned_clstrs@clstrs[[i]]@sids,
                                                 kept$sid[kept$cluster_id == as.integer(cid)])
  }
  selected_clean <- selected
  selected_clean@clstrs <- cleaned_clstrs
  nonempty <- as.character(sort(unique(kept$cluster_id)))
  phylota_final <- phylotaR::drop_clstrs(phylota = selected_clean, cid = nonempty)
  log_message("Clusters with sequences after the accession-level filters: ", length(phylota_final@cids))

  smmry <- phylotaR::summary(phylota_final) |> dplyr::mutate(ID = as.integer(ID))
  cluster_gene_summary <- cp_summarise_cluster_markers(df_meta) |> dplyr::mutate(cluster_id = as.integer(cluster_id))
  smmry <- smmry |> dplyr::left_join(cluster_gene_summary, by = c("ID" = "cluster_id"))
  seed_meta <- .bc_seed_metadata_cached(smmry$Seed, file.path(dir_cache, "CACHE_SEED_METADATA_BARCODING.csv"),
                                        log_message = log_message)
  smmry <- .cp_enrich_cluster_summary(smmry, seed_meta, pattern, gene_lookup, marker_lookup)

  # Clusters whose gene is not recognised in target_genes.txt stay in the funnel under a provisional
  # name (decision of 2026-09-22), so a locus that is new in GenBank is not lost at this point.
  named <- .bc_name_clusters(data.frame(ID = smmry$ID, Marker_std = smmry$Marker_std,
                                        Description = smmry$Description, stringsAsFactors = FALSE))
  marker_map <- named$map
  if (nrow(named$unnamed) > 0) {
    log_message("Clusters with no recognised gene, kept under a provisional name: ",
                paste(named$unnamed$locus, collapse = ", "))
  }
  utils::write.csv(named$unnamed, file.path(dir_asm, "TABLE_barcoding_unnamed_clusters.csv"), row.names = FALSE)
  utils::write.csv(marker_map, file.path(dir_asm, "TABLE_barcoding_cluster_marker_assignment.csv"), row.names = FALSE)
  funnel <- .bc_cluster_funnel(n_workspace = length(phylota@cids), n_selected = length(selected@cids),
                               n_with_sequences = length(phylota_final@cids), marker_map = marker_map)
  utils::write.csv(funnel, file.path(dir_asm, "TABLE_barcoding_cluster_funnel.csv"), row.names = FALSE)

  comparison <- NULL
  if (file.exists(phylogeny_map_file)) {
    comparison <- .bc_compare_marker_maps(marker_map, utils::read.csv(phylogeny_map_file, stringsAsFactors = FALSE))
    utils::write.csv(comparison, file.path(dir_asm, "TABLE_barcoding_marker_assignment_comparison.csv"), row.names = FALSE)
    tab <- table(comparison$estado)
    log_message("Comparison with the phylogeny's assignment: ",
                paste(names(tab), as.integer(tab), sep = " ", collapse = "; "))
  } else {
    log_message("No cluster-to-locus table of the phylogeny found; no comparison written.")
  }

  kept$locus <- marker_map$Marker_std[match(as.character(kept$cluster_id), marker_map$ID)]
  kept <- kept[!is.na(kept$locus), , drop = FALSE]
  kept <- kept[!duplicated(paste(kept$locus, kept$sid)), , drop = FALSE]

  # 5. Names, with the rule of clean_taxonomic_names()
  kept$species <- .bc_resolve_names(kept$nombre_genbank, checklist_path)
  discarded <- kept[is.na(kept$species), , drop = FALSE]
  discarded_tab <- as.data.frame(table(nombre_genbank = discarded$nombre_genbank, locus = discarded$locus),
                                 stringsAsFactors = FALSE)
  discarded_tab <- discarded_tab[discarded_tab$Freq > 0, , drop = FALSE]
  names(discarded_tab)[names(discarded_tab) == "Freq"] <- "accesiones"
  utils::write.csv(discarded_tab, file.path(dir_asm, "TABLE_barcoding_discarded_names.csv"), row.names = FALSE)
  log_message("Accessions dropped because their name is not in the checklist: ", nrow(discarded),
              " in ", length(unique(discarded$nombre_genbank)), " names.")

  # 6. Registry, FASTA files and per-locus summary
  registry <- .bc_build_registry(kept[!is.na(kept$species), , drop = FALSE])
  utils::write.csv(registry, file.path(dir_asm, "TABLE_barcoding_accession_registry.csv"), row.names = FALSE)

  sequences <- vapply(registry$sid, function(s) rawToChar(selected@sqs[[s]]@sq), character(1))
  names(sequences) <- registry$sid
  .bc_write_locus_fastas(registry, sequences, dir_asm)

  summary_tab <- .bc_locus_summary(registry)
  utils::write.csv(summary_tab, file.path(dir_asm, "TABLE_barcoding_replication_summary_assembly.csv"), row.names = FALSE)
  log_message("Registry: ", nrow(registry), " accessions of ", length(unique(registry$species)), " species in ",
              length(unique(registry$locus)), " loci.")

  invisible(list(
    registry = registry,
    summary = summary_tab,
    marker_map = marker_map,
    comparison = comparison,
    funnel = funnel,
    counts = list(clusters_selected = length(selected@cids), memberships = nrow(records),
                  sids = length(unique(records$sid)), redundant_removed = nrow(dd$removed),
                  manual = counts_manual, discarded_accessions = nrow(discarded))
  ))
}
