#' Assemble Ingroup Sequence Clusters via phylotaR Similarity Mining
#'
#' Retrieves orthologous sequence clusters for a focal taxonomic ingroup (e.g., family **Cactaceae**, NCBI Taxonomy ID: 3593)
#' directly from GenBank using similarity clustering via `phylotaR` (Bennett *et al.*, 2018).
#' Clusters are formed by sequence similarity, not by GenBank locus annotations, so sequences are
#' not lost to gene-name synonyms or mislabeled records.
#'
#' @param wd_path Character. Path to the `phylotaR` workspace directory storing local parameters and database caches.
#' @param target_genes_file Character. Path to the target locus list text file. If `NULL`, defaults to package `inst/extdata/target_genes.txt`.
#' @param genes_map_file Character. Path to the gene synonymy mapping CSV file. If `NULL`, defaults to package `inst/extdata/genes_map.csv`.
#' @param manual_exclusions_file Character. Path to the accession exclusion CSV file. If `NULL`, defaults to package `inst/extdata/manual_exclusions_ingroup.csv`.
#' @param apply_manual_exclusions Logical. Apply the curated accession exclusion list? **Defaults to `TRUE`.** The lists shipped with the package are the product of manual curation by the expert team supporting `PhyloCactus`: each excluded accession was inspected and removed on taxonomic or sequence-quality grounds that are not recoverable from GenBank metadata alone. The lists are keyed by locus and accession (`locus`, `sid`, `reason`, `source`), so that they survive a new mining that renumbers the clusters: an accession is removed from every cluster of the locus named (decision X2 of 29-09). Across both lists they hold 81 accessions, in 5 ingroup loci and 13 outgroup loci, and they live in `inst/extdata/manual_exclusions_ingroup.csv` and `inst/extdata/manual_exclusions_outgroup.csv`. Setting this to `FALSE` reproduces the uncurated cluster set, which is the way to quantify what the curation actually removes; the applied list is always written to `TABLE_MANUAL_EXCLUSIONS_*.csv` in the output directory, so any run documents its own curation state.
#' @param min_species Integer. Species-richness threshold for cluster retention. Clusters are retained when they contain **strictly more than** `min_species` distinct species, so `min_species = 50` keeps clusters with 51 species or more. Defaults to `50`.
#' @param preferred_parent Character. NCBI Taxonomy ID of the focal ingroup parent node. Defaults to `"3593"` (Cactaceae).
#' @param ncbi_dr Character. Path to local `BLAST+` binaries directory. If `NULL`, attempts system environment auto-detection.
#' @param force_download Logical. Force fresh database retrieval instead of using local cache? Defaults to `FALSE`.
#' @param out_dir Character. Output directory path to save cluster summaries, FASTA sequence matrices, and log reports.
#' @param notify Logical. Send an email through [send_run_notification()] when a mining ends or fails;
#'   reading an existing workspace sends nothing. Defaults to `FALSE`.
#' @param notify_to,notify_credentials Passed to [send_run_notification()] as `to` and `credentials`.
#' @return A data frame summarizing sequence occupancy, taxon representation, and cluster characteristics across retained loci.
#' @references
#' Bennett, D. J., Hettling, H., Silvestro, D., Zizka, A., Bacon, C. D., Faurby, S., ... & Antonelli, A. (2018).
#' phylotaR: An automated pipeline for retrieving orthologous DNA sequences from GenBank in R.
#' *Life*, 8(2), 20. \doi{10.3390/life8020020}
#' @examples
#' \dontrun{
#' assemble_ingroup_phylotar(
#'   wd_path = "0_phylotaR_raw_Ingroup",
#'   min_species = 50,
#'   preferred_parent = "3593"
#' )
#' }
#' @export
assemble_ingroup_phylotar <- function(wd_path, target_genes_file = NULL, genes_map_file = NULL, manual_exclusions_file = NULL, apply_manual_exclusions = TRUE, min_species = 50, preferred_parent = "3593", ncbi_dr = NULL, force_download = FALSE, out_dir = "1_phylotaR_out_Ingroup", notify = FALSE, notify_to = NULL, notify_credentials = NULL) {
  
  
  # Resolve inputs
  if (is.null(target_genes_file)) target_genes_file <- system.file("extdata", "target_genes.txt", package = "PhyloCactus")
  if (is.null(genes_map_file)) genes_map_file <- system.file("extdata", "genes_map.csv", package = "PhyloCactus")
  if (is.null(manual_exclusions_file)) manual_exclusions_file <- system.file("extdata", "manual_exclusions_ingroup.csv", package = "PhyloCactus")
  
  if (!file.exists(target_genes_file)) stop("Target genes file not found at: ", target_genes_file)
  if (!file.exists(genes_map_file)) stop("Genes map file not found at: ", genes_map_file)
  
  dir_out_base <- out_dir
  dir_out_cluster_fasta <- file.path(dir_out_base, "Cluster_raw")
  dir_out_logs <- file.path(dir_out_base, "logs")
  dir_out_cache <- file.path(dir_out_base, "cache")
  
  dir.create(wd_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_out_base, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_out_cluster_fasta, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_out_logs, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_out_cache, recursive = TRUE, showWarnings = FALSE)
  
  path_log <- file.path(dir_out_logs, "LOG_INGROUP_PHYLOTAR_ASSEMBLY.txt")
  path_table_cluster_summary_raw <- file.path(dir_out_base, "TABLE_CLUSTER_SUMMARY_INGROUP_RAW.csv")
  path_table_cluster_summary_clean <- file.path(dir_out_base, "TABLE_CLUSTER_SUMMARY_INGROUP_CLEAN.csv")
  path_table_species_cluster_map_clean <- file.path(dir_out_base, "TABLE_SPECIES_CLUSTER_MAP_INGROUP_CLEAN.csv")
  path_table_accession_occupancy_clean <- file.path(dir_out_base, "TABLE_ACCESSION_OCCUPANCY_INGROUP_CLEAN.csv")
  path_table_marker_summary <- file.path(dir_out_base, "TABLE_MARKER_SUMMARY_INGROUP.csv")
  path_table_duplicate_conflicts <- file.path(dir_out_base, "TABLE_DUPLICATE_SID_CONFLICTS_INGROUP.csv")
  path_table_manual_exclusions <- file.path(dir_out_base, "TABLE_MANUAL_EXCLUSIONS_INGROUP.csv")
  path_table_cluster_marker_assignment <- file.path(dir_out_base, "TABLE_CLUSTER_MARKER_ASSIGNMENT_INGROUP.csv")
  
  path_metadata_cache_csv <- file.path(dir_out_cache, "CACHE_GENBANK_METADATA_INGROUP.csv")
  path_phylota_clean_rdata <- file.path(dir_out_base, "RDATA_PHYLOTAR_INGROUP_CLEANED.RData")
  
  log_message <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
    cat(msg, "\n")
    write(msg, file = path_log, append = TRUE)
  }
  
  if (file.exists(path_log)) file.remove(path_log)
  log_message("Starting ingroup phylotaR assembly refactor. \U0001f335")
  
  # Read markers
  markers_df <- utils::read.table(target_genes_file, header = TRUE, stringsAsFactors = FALSE)
  markers <- unique(cp_normalize_marker(as.character(markers_df$markers)))
  pattern <- cp_build_pattern(markers)
  
  genes_map_df <- utils::read.csv(genes_map_file, stringsAsFactors = FALSE) |>
    dplyr::mutate(search = cp_normalize_marker(search), replace = trimws(replace)) |>
    dplyr::distinct(search, .keep_all = TRUE)
    
  gene_lookup <- genes_map_df |> dplyr::select(search_gene = search, Gene_std = replace)
  marker_lookup <- genes_map_df |> dplyr::select(search_marker = search, Marker_std = replace)
  
  log_message("Inputs loaded successfully.")
  
  # Setup / Load phylotaR. Shared with assemble_barcoding_dataset(), so that both branches make the
  # same phylotaR call and share the same raw workspace, whichever of them runs first.
  phylota <- .phylotar_load_or_mine(wd_path, preferred_parent = preferred_parent, ncbi_dr = ncbi_dr,
                                    force_download = force_download, log_message = log_message,
                                    notify = notify, notify_to = notify_to, notify_credentials = notify_credentials)
  
  # 8. SPECIES REDUCTION AND CLUSTER FILTERING
  species_reduced <- phylotaR::drop_by_rank(phylota, rnk = "species", n = 1)
  cluster_ids <- species_reduced@cids
  # Y1 (BMM, 29-09): the clusters of one locus are pooled before the cut, so that a locus split into
  # small clusters by the mining is not lost; a cluster with no locus name is judged alone
  cluster_locus <- .cp_cluster_loci(phylota, pattern, genes_map_df, marker_lookup)
  cluster_species <- .cp_cluster_species(species_reduced)
  # Y2: only the clusters that overlap the main cluster of their locus are pooled
  overlap <- .cp_cluster_overlap(phylota, cluster_locus, cluster_species)
  by_locus <- .cp_select_clusters_by_locus(cluster_species, cluster_locus, min_species, overlap = overlap)
  if (nrow(by_locus$fragments) > 0) {
    log_message("Clusters not pooled with their locus (no overlap with its main cluster): ",
                paste(by_locus$fragments$cluster_id, by_locus$fragments$locus, sep = " ", collapse = "; "))
  }
  keep_clusters <- intersect(cluster_ids, as.character(by_locus$keep))
  selected <- phylotaR::drop_clstrs(species_reduced, cid = keep_clusters)
  
  log_message("Clusters before >", min_species, " filter: ", length(cluster_ids))
  log_message("Clusters retained after >", min_species, " filter (clusters of one locus pooled): ", length(selected@cids),
              " in ", sum(by_locus$table$kept), " loci or unnamed clusters")
  
  # 9. RAW CLUSTER TABLES AND METADATA
  df_species_clusters <- cp_extract_cluster_species_sid(selected)
  metadata_raw <- cp_download_all_metadata(unique(df_species_clusters$sid), path_metadata_cache_csv, batch_size = 200, sleep_time = 0.5, max_retries = 5, log_message = log_message)
  
  df_species_clusters_metadata <- df_species_clusters |>
    dplyr::left_join(metadata_raw, by = "sid") |>
    cp_annotate_marker_text(pattern = pattern, genes_map_df = genes_map_df)
    
  smmry_sel <- phylotaR::summary(selected) |> dplyr::mutate(ID = as.integer(ID))
  cluster_gene_summary <- cp_summarise_cluster_markers(df_species_clusters_metadata) |> dplyr::mutate(cluster_id = as.integer(cluster_id))
  smmry_sel_enriched <- smmry_sel |> dplyr::left_join(cluster_gene_summary, by = c("ID" = "cluster_id"))
  
  seed_meta <- cp_fetch_seed_metadata(smmry_sel_enriched$Seed, log_message = log_message)

  smmry_sel_enriched <- .cp_enrich_cluster_summary(smmry_sel_enriched, seed_meta, pattern, gene_lookup, marker_lookup)
  
  readr::write_csv(tibble::as_tibble(smmry_sel_enriched), path_table_cluster_summary_raw)
  
  # 10. RAW SPECIES / ACCESSION MAPS
  species_cluster_map <- df_species_clusters_metadata |>
    dplyr::group_by(species) |>
    dplyr::summarise(clusters = paste(unique(cluster_id), collapse = ","), Genes = paste(unique(Genes_text), collapse = "; "), n_clusters = dplyr::n_distinct(cluster_id), .groups = "drop")
    
  species_sid_matrix <- df_species_clusters_metadata |>
    dplyr::select(species, cluster_id, sid) |>
    tidyr::pivot_wider(names_from = cluster_id, names_prefix = "cid_", values_from = sid, values_fill = "")
    
  species_table <- df_species_clusters_metadata |> dplyr::left_join(species_sid_matrix, by = "species")
  
  sid_conflicts <- df_species_clusters_metadata |>
    dplyr::group_by(sid) |>
    dplyr::summarise(n_clusters = dplyr::n(), clusters = paste(unique(cluster_id), collapse = ","), species = paste(unique(species), collapse = ","), .groups = "drop") |>
    dplyr::filter(n_clusters > 1)
    
  readr::write_csv(tibble::as_tibble(sid_conflicts), path_table_duplicate_conflicts)
  
  species_with_sid_dup <- df_species_clusters_metadata |> dplyr::filter(sid %in% sid_conflicts$sid) |> dplyr::distinct(species) |> dplyr::pull(species)
  df_species_clusters_metadata <- df_species_clusters_metadata |>
    dplyr::left_join(sid_conflicts |> dplyr::select(sid, duplicated = n_clusters), by = "sid") |>
    dplyr::mutate(duplicated = ifelse(is.na(duplicated), 0L, as.integer(duplicated)))
    
  species_table <- species_table |>
    dplyr::left_join(sid_conflicts |> dplyr::select(sid, duplicated = n_clusters), by = "sid") |>
    dplyr::mutate(duplicated = ifelse(is.na(duplicated), 0L, as.integer(duplicated)))
    
  species_cluster_map <- species_cluster_map |>
    dplyr::mutate(duplicated = ifelse(species %in% species_with_sid_dup, 1L, 0L))
    
  # 11. DUPLICATE-RESOLUTION CLEANING
  dup_sids <- sid_conflicts$sid
  dup_records <- df_species_clusters_metadata |>
    dplyr::filter(sid %in% dup_sids) |>
    dplyr::mutate(
      cluster_id = as.integer(cluster_id),
      prnt = vapply(cluster_id, function(cid) cp_get_cluster_parent(selected, cid), FUN.VALUE = character(1))
    )
    
  keepers <- dup_records |>
    dplyr::group_by(sid) |>
    dplyr::mutate(
      keep = dplyr::case_when(
        any(prnt == preferred_parent, na.rm = TRUE) ~ (prnt == preferred_parent),
        TRUE ~ (cluster_id == min(cluster_id, na.rm = TRUE))
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(keep) |>
    dplyr::select(cluster_id, sid) |>
    dplyr::distinct()
    
  rows_to_delete_auto <- dup_records |>
    dplyr::anti_join(keepers, by = c("cluster_id", "sid")) |>
    dplyr::select(cluster_id, sid) |>
    dplyr::distinct()
    
  manual_rows <- .cp_read_manual_exclusions(apply_manual_exclusions, manual_exclusions_file, path_table_manual_exclusions,
                                            df_species_clusters_metadata, cluster_locus, log_message)
  
  rows_to_delete <- dplyr::bind_rows(
    rows_to_delete_auto |> dplyr::mutate(reason = "automatic duplicate resolution"),
    manual_rows
  ) |> dplyr::distinct(cluster_id, sid, .keep_all = TRUE)
  
  log_message("Automatic duplicate removals: ", nrow(rows_to_delete_auto))
  log_message("Manual duplicate removals: ", nrow(manual_rows))
  log_message("Total duplicate removals applied: ", nrow(rows_to_delete))
  
  # 12. CLEAN TABLES
  df_species_clusters_clean <- df_species_clusters_metadata |>
    dplyr::anti_join(rows_to_delete |> dplyr::select(cluster_id, sid), by = c("cluster_id", "sid"))
  
  species_table_clean <- species_table |>
    dplyr::mutate(cluster_id = as.integer(cluster_id)) |>
    dplyr::anti_join(rows_to_delete |> dplyr::select(cluster_id, sid), by = c("cluster_id", "sid"))
    
  species_cluster_map_clean <- df_species_clusters_clean |>
    dplyr::group_by(species) |>
    dplyr::summarise(clusters = paste(unique(cluster_id), collapse = ","), n_clusters = dplyr::n_distinct(cluster_id), Genes = paste(unique(Genes_text), collapse = "; "), .groups = "drop") |>
    dplyr::left_join(species_sid_matrix, by = "species")
    
  # 13. UPDATED SUMMARY COUNTS
  smmry_sel_enriched_fixed <- smmry_sel_enriched |> dplyr::mutate(cluster_id = as.integer(ID))
  
  n_original <- species_sid_matrix |>
    tidyr::pivot_longer(cols = dplyr::starts_with("cid_"), names_to = "cid_col", values_to = "sid", values_drop_na = FALSE) |>
    dplyr::mutate(cluster_id = as.integer(sub("^cid_", "", cid_col))) |>
    dplyr::filter(sid != "") |> dplyr::group_by(cluster_id) |> dplyr::summarise(n_original = dplyr::n(), .groups = "drop")
    
  n_removed <- rows_to_delete |> dplyr::group_by(cluster_id) |> dplyr::summarise(n_removed = dplyr::n(), .groups = "drop")
  n_final <- df_species_clusters_clean |> dplyr::group_by(cluster_id) |> dplyr::summarise(n_final = dplyr::n(), .groups = "drop")
  
  smmry_sel_enriched_updated <- smmry_sel_enriched_fixed |>
    dplyr::left_join(n_original, by = "cluster_id") |> dplyr::left_join(n_removed, by = "cluster_id") |> dplyr::left_join(n_final, by = "cluster_id") |>
    dplyr::mutate(n_original = ifelse(is.na(n_original), 0L, as.integer(n_original)), n_removed = ifelse(is.na(n_removed), 0L, as.integer(n_removed)), n_final = ifelse(is.na(n_final), 0L, as.integer(n_final)))
    
  # 14. CLEAN PHYLOTA OBJECT
  cleaned_clstrs <- selected@clstrs
  for (i in seq_along(cleaned_clstrs@clstrs)) {
    cid <- names(cleaned_clstrs@clstrs)[i]
    sids_now <- cleaned_clstrs@clstrs[[i]]@sids
    sids_remove <- rows_to_delete |> dplyr::filter(cluster_id == as.integer(cid)) |> dplyr::pull(sid)
    cleaned_clstrs@clstrs[[i]]@sids <- setdiff(sids_now, sids_remove)
  }
  selected_clean <- selected
  selected_clean@clstrs <- cleaned_clstrs
  
  clusters_to_keep_after_cleaning <- smmry_sel_enriched_updated |> dplyr::filter(n_final > 0) |> dplyr::pull(ID) |> as.character()
  phylota_final <- phylotaR::drop_clstrs(phylota = selected_clean, cid = clusters_to_keep_after_cleaning)
  
  log_message("Clusters retained after duplicate cleaning: ", length(phylota_final@cids))
  
  # 15. FINAL TABLES AFTER CLEANING
  final_cids <- phylota_final@cids
  ntaxa_final <- phylotaR::get_ntaxa(phylota = phylota_final, cid = final_cids, rnk = "species")
  ntaxa_final_df <- tibble::tibble(cluster_id = final_cids, n_taxa = ntaxa_final)
  
  df_species_clusters_final <- cp_extract_cluster_species_sid(phylota_final)
  metadata_final <- metadata_raw
  
  df_species_clusters_metadata_final <- df_species_clusters_final |>
    dplyr::left_join(metadata_final, by = "sid") |> cp_annotate_marker_text(pattern = pattern, genes_map_df = genes_map_df)
    
  smmry_final <- phylotaR::summary(phylota_final) |> dplyr::mutate(ID = as.integer(ID))
  cluster_gene_summary_final <- cp_summarise_cluster_markers(df_species_clusters_metadata_final) |> dplyr::mutate(cluster_id = as.integer(cluster_id))
  smmry_final_enriched <- smmry_final |> dplyr::left_join(cluster_gene_summary_final, by = c("ID" = "cluster_id"))
  
  seed_meta_final <- cp_fetch_seed_metadata(smmry_final_enriched$Seed, log_message = log_message)
  smmry_final_enriched <- .cp_enrich_cluster_summary(smmry_final_enriched, seed_meta_final, pattern, gene_lookup, marker_lookup)
  
  species_sid_matrix_final <- df_species_clusters_metadata_final |> dplyr::select(species, cluster_id, sid) |> tidyr::pivot_wider(names_from = cluster_id, names_prefix = "cid_", values_from = sid, values_fill = "")
  species_genes_final <- df_species_clusters_metadata_final |> dplyr::group_by(species) |> dplyr::summarise(Genes_text = paste(unique(Genes_text), collapse = "; "), .groups = "drop")
  species_table_final <- species_sid_matrix_final |> dplyr::left_join(species_genes_final, by = "species") |> dplyr::relocate(Genes_text, .after = species)
  
  # 16. EXPORT CLEAN CLUSTER FASTAS
  cp_write_cluster_fastas(phylota_final, dir_out_cluster_fasta)
  log_message("Cluster FASTA export complete.")
  
  # 17. MERGE CLUSTERS BY STANDARDIZED MARKER
  df_map <- smmry_final_enriched |> dplyr::mutate(ID = as.integer(ID)) |> dplyr::select(ID, Marker_std) |> dplyr::distinct()
  readr::write_csv(df_map, path_table_cluster_marker_assignment)
  
  cluster_fasta_files <- list.files(dir_out_cluster_fasta, full.names = TRUE, pattern = "^CLUSTER_.*\\.fasta$") |> sort()
  df_files <- tibble::tibble(file = cluster_fasta_files, ID = as.integer(sub("^CLUSTER_([0-9]+)\\.fasta$", "\\1", basename(cluster_fasta_files))))
  df_joined <- df_files |> dplyr::left_join(df_map, by = "ID")
  
  # A cluster whose description matches no rule in genes_map.csv comes back as an empty string, not
  # as NA, so testing for NA alone let it through: paste0("", ".fasta") then wrote every unmapped
  # cluster into a single file named ".fasta", which no later stage reads and no later report
  # mentions. Found on 2026-09-04 in the outgroup run, where rpoB and psbB of Talinum paniculatum
  # and Portulaca oleracea disappeared that way. Harmless there, because neither locus has an
  # ingroup counterpart, but a variant description of a real marker would vanish identically.
  df_joined <- df_joined |> dplyr::mutate(Marker_std = trimws(as.character(Marker_std)))
  unmapped <- df_joined |> dplyr::filter(is.na(Marker_std) | !nzchar(Marker_std))
  if (nrow(unmapped) > 0) {
    warning("Clusters with no standardized marker name, dropped from the marker FASTAs: IDs ",
            paste(sort(unique(unmapped$ID)), collapse = ", "),
            ". Their descriptions match no rule in genes_map.csv. Check the Description column of ",
            basename(path_table_cluster_summary_clean),
            " and add a rule if any of them is a locus this analysis should use.", call. = FALSE)
    log_message("Clusters without Marker_std, dropped: ", paste(sort(unique(unmapped$ID)), collapse = ", "))
  }

  df_joined <- df_joined |> dplyr::filter(!is.na(Marker_std), nzchar(Marker_std))
  markers_final <- sort(unique(df_joined$Marker_std))
  
  for (mk in markers_final) {
    log_message("Processing marker FASTA merge: ", mk)
    files_mk <- df_joined |> dplyr::filter(Marker_std == mk) |> dplyr::arrange(file) |> dplyr::pull(file)
    seqs_list <- lapply(files_mk, Biostrings::readDNAStringSet)
    seqs_all <- do.call(c, seqs_list)
    sp_names <- sub("^([^ ]+_[^ ]+).*", "\\1", names(seqs_all))
    seqs_clean <- seqs_all[!duplicated(sp_names)]
    out_file <- file.path(dir_out_base, paste0(mk, ".fasta"))
    Biostrings::writeXStringSet(seqs_clean, out_file)
  }
  log_message("Marker FASTA export complete.")
  
  # 18. FINAL OCCUPANCY AND MARKER SUMMARIES
  df_species_clusters_metadata_final <- df_species_clusters_metadata_final |> dplyr::mutate(cluster_id = as.integer(cluster_id))
  df_joined <- df_joined |> dplyr::mutate(ID = as.integer(ID))
  smmry_final_enriched <- smmry_final_enriched |> dplyr::mutate(ID = as.integer(as.character(ID)))
  
  occ_table <- df_species_clusters_metadata_final |> dplyr::left_join(df_joined |> dplyr::select(ID, Marker_std), by = c("cluster_id" = "ID"))
  cluster_level_stats <- smmry_final_enriched |> dplyr::mutate(cluster_id = as.integer(ID)) |> dplyr::select(cluster_id, Marker_std, n_sequences)
  
  marker_summary <- cluster_level_stats |>
    dplyr::group_by(Marker_std) |>
    dplyr::summarise(num_clusters = dplyr::n(), num_sequences = sum(n_sequences, na.rm = TRUE), mean_seq_per_cluster = mean(n_sequences), median_seq_per_cluster = stats::median(n_sequences), min_seq_in_cluster = min(n_sequences), max_seq_in_cluster = max(n_sequences), .groups = "drop") |>
    dplyr::left_join(occ_table |> dplyr::group_by(Marker_std) |> dplyr::summarise(num_species = dplyr::n_distinct(Species_gb), .groups = "drop"), by = "Marker_std") |>
    dplyr::select(Marker_std, num_clusters, num_species, num_sequences, mean_seq_per_cluster, median_seq_per_cluster, min_seq_in_cluster, max_seq_in_cluster) |>
    dplyr::arrange(Marker_std)
    
  # 19. EXPORT TABLES
  readr::write_csv(tibble::as_tibble(smmry_final_enriched), path_table_cluster_summary_clean)
  readr::write_csv(tibble::as_tibble(species_table_final), path_table_species_cluster_map_clean)
  readr::write_csv(tibble::as_tibble(occ_table), path_table_accession_occupancy_clean)
  readr::write_csv(tibble::as_tibble(marker_summary), path_table_marker_summary)
  
  # 20. SAVE OBJECTS
  save(phylota_final, smmry_final_enriched, species_table_final, occ_table, marker_summary, file = path_phylota_clean_rdata)
  log_message("Assembly and export fully completed. \U0001f335")
  
  res <- list(
    selected_clusters = phylota_final,
    retained_cids = phylota_final@cids
  )
  message("\nIngroup assembly and export fully completed. \U0001f335")
  return(res)
}

#' Assemble Outgroup Sequence Clusters via phylotaR
#'
#' Retrieves orthologous sequence clusters for specified outgroup lineages (e.g., *Portulaca*, *Anacampseros*, *Talinopsis*, *Grahamia*, *Talinum*, *Talinella*)
#' matching the locus target constraints defined for the focal ingroup. Outer reference sampling provides phylogenetically
#' informative root positions necessary for maximum-likelihood tree search and divergence time estimation.
#'
#' @param wd_path Character. Path to the `phylotaR` workspace directory.
#' @param target_genes_file Character. Path to the target locus list text file. If `NULL`, defaults to package `inst/extdata/target_genes.txt`.
#' @param genes_map_file Character. Path to the gene synonymy mapping CSV file. If `NULL`, defaults to package `inst/extdata/genes_map.csv`.
#' @param manual_exclusions_file Character. Path to outgroup exclusions CSV file. If `NULL`, defaults to package `inst/extdata/manual_exclusions_outgroup.csv`.
#' @param apply_manual_exclusions Logical. Apply the curated accession exclusion list? Defaults to `TRUE`; see [assemble_ingroup_phylotar()] for the provenance of these lists and for what setting it to `FALSE` is useful for.
#' @param outgroups Character vector of NCBI Taxonomy IDs for outgroup lineages. Defaults to `c("107598", "107617", "107583", "3582", "107600", "108056")`.
#' @param force_download Logical. Force fresh database retrieval instead of using local cache? Defaults to `FALSE`.
#' @param out_dir Character. Output directory path to save outgroup cluster tables and FASTA sequence files.
#' @param notify Logical. Send an email through [send_run_notification()] when a mining ends or fails;
#'   reading an existing workspace sends nothing. Defaults to `FALSE`.
#' @param notify_to,notify_credentials Passed to [send_run_notification()] as `to` and `credentials`.
#' @return A list containing the processed outgroup cluster objects and retained cluster IDs.
#' @references
#' Bennett, D. J., Hettling, H., Silvestro, D., Zizka, A., Bacon, C. D., Faurby, S., ... & Antonelli, A. (2018).
#' phylotaR: An automated pipeline for retrieving orthologous DNA sequences from GenBank in R.
#' *Life*, 8(2), 20. \doi{10.3390/life8020020}
#' @examples
#' \dontrun{
#' assemble_outgroup_phylotar(
#'   wd_path = "0_phylotaR_raw_Outgroup",
#'   outgroups = c("3582", "107583")
#' )
#' }
#' @export
assemble_outgroup_phylotar <- function(wd_path, target_genes_file = NULL, genes_map_file = NULL, manual_exclusions_file = NULL, apply_manual_exclusions = TRUE, outgroups = c("107598", "107617", "107583", "3582", "107600", "108056"), force_download = FALSE, out_dir = "1_phylotaR_out_Outgroup", notify = FALSE, notify_to = NULL, notify_credentials = NULL) {
  
  
  if (is.null(target_genes_file)) target_genes_file <- system.file("extdata", "target_genes.txt", package = "PhyloCactus")
  if (is.null(genes_map_file)) genes_map_file <- system.file("extdata", "genes_map.csv", package = "PhyloCactus")
  if (is.null(manual_exclusions_file)) manual_exclusions_file <- system.file("extdata", "manual_exclusions_outgroup.csv", package = "PhyloCactus")
  
  dir_out_base <- out_dir
  dir_out_cluster_fasta <- file.path(dir_out_base, "Cluster_raw_outgroup")
  dir_out_logs <- file.path(dir_out_base, "logs")
  dir_out_cache <- file.path(dir_out_base, "cache")
  
  dir.create(wd_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_out_base, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_out_cluster_fasta, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_out_logs, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_out_cache, recursive = TRUE, showWarnings = FALSE)
  
  path_log <- file.path(dir_out_logs, "LOG_OUTGROUP_PHYLOTAR_ASSEMBLY.txt")
  path_table_cluster_summary_raw <- file.path(dir_out_base, "TABLE_CLUSTER_SUMMARY_OUTGROUP_RAW.csv")
  path_table_cluster_summary_clean <- file.path(dir_out_base, "TABLE_CLUSTER_SUMMARY_OUTGROUP_CLEAN.csv")
  path_table_species_cluster_map_clean <- file.path(dir_out_base, "TABLE_SPECIES_CLUSTER_MAP_OUTGROUP_CLEAN.csv")
  path_table_accession_occupancy_clean <- file.path(dir_out_base, "TABLE_ACCESSION_OCCUPANCY_OUTGROUP_CLEAN.csv")
  path_table_marker_summary <- file.path(dir_out_base, "TABLE_MARKER_SUMMARY_OUTGROUP.csv")
  path_table_duplicate_conflicts <- file.path(dir_out_base, "TABLE_DUPLICATE_SID_CONFLICTS_OUTGROUP.csv")
  path_table_manual_exclusions <- file.path(dir_out_base, "TABLE_MANUAL_EXCLUSIONS_OUTGROUP.csv")
  path_table_cluster_marker_assignment <- file.path(dir_out_base, "TABLE_CLUSTER_MARKER_ASSIGNMENT_OUTGROUP.csv")
  
  path_metadata_cache_csv <- file.path(dir_out_cache, "CACHE_GENBANK_METADATA_OUTGROUP.csv")
  path_phylota_clean_rdata <- file.path(dir_out_base, "RDATA_PHYLOTAR_OUTGROUP_CLEANED.RData")
  
  log_message <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
    cat(msg, "\n")
    write(msg, file = path_log, append = TRUE)
  }
  
  if (file.exists(path_log)) file.remove(path_log)
  log_message("Starting outgroup phylotaR assembly refactor. \U0001f335")
  
  markers_df <- utils::read.table(target_genes_file, header = TRUE, stringsAsFactors = FALSE)
  markers <- unique(cp_normalize_marker(as.character(markers_df$markers)))
  pattern <- cp_build_pattern(markers)
  genes_map_df <- utils::read.csv(genes_map_file, stringsAsFactors = FALSE) |>
    dplyr::mutate(search = cp_normalize_marker(search), replace = trimws(replace)) |>
    dplyr::distinct(search, .keep_all = TRUE)
    
  gene_lookup <- genes_map_df |> dplyr::select(search_gene = search, Gene_std = replace)
  marker_lookup <- genes_map_df |> dplyr::select(search_marker = search, Marker_std = replace)
  
  # Workspace through the function shared with the ingroup and the barcoding branch: one set of
  # parameters for both halves of folder 0 (E1, BMM, 28-09)
  phylota <- .phylotar_load_or_mine(wd_path, preferred_parent = outgroups[1], txid = outgroups,
                                    force_download = force_download, log_message = log_message,
                                    notify = notify, notify_to = notify_to, notify_credentials = notify_credentials)
  
  # 8. SPECIES REDUCTION AND CLUSTER FILTERING
  species_reduced <- phylotaR::drop_by_rank(phylota, rnk = "species", n = 1)
  selected <- phylotaR::drop_clstrs(species_reduced, cid = species_reduced@cids)
  
  log_message("Clusters before >0 filter: ", length(species_reduced@cids))
  log_message("Clusters retained after >0 filter: ", length(selected@cids))
  
  # 9. RAW CLUSTER TABLES AND METADATA
  df_species_clusters <- cp_extract_cluster_species_sid(selected)
  metadata_raw <- cp_download_all_metadata(unique(df_species_clusters$sid), path_metadata_cache_csv, log_message = log_message)
  
  df_species_clusters_metadata <- df_species_clusters |>
    dplyr::left_join(metadata_raw, by = "sid") |>
    cp_annotate_marker_text(pattern = pattern, genes_map_df = genes_map_df)
    
  smmry_sel <- phylotaR::summary(selected) |> dplyr::mutate(ID = as.integer(ID))
  cluster_gene_summary <- cp_summarise_cluster_markers(df_species_clusters_metadata) |> dplyr::mutate(cluster_id = as.integer(cluster_id))
  smmry_sel_enriched <- smmry_sel |> dplyr::left_join(cluster_gene_summary, by = c("ID" = "cluster_id"))
  
  seed_meta <- cp_fetch_seed_metadata(smmry_sel_enriched$Seed, log_message = log_message)
  smmry_sel_enriched <- smmry_sel_enriched |>
    dplyr::mutate(Species = seed_meta$Species, Description = cp_normalize_marker(seed_meta$Description)) |>
    dplyr::mutate(
      Genes_raw = stringr::str_extract_all(Description, pattern),
      Genes_text = dplyr::if_else(
        is.na(Description), NA_character_,
        purrr::map_chr(Genes_raw, ~ paste(unique(cp_normalize_marker(.x)), collapse = ", "))
      ),
      top_marker = cp_normalize_marker(top_marker)
    ) |>
    dplyr::select(-Genes_raw) |>
    dplyr::left_join(gene_lookup, by = c("Genes_text" = "search_gene")) |> dplyr::left_join(marker_lookup, by = c("top_marker" = "search_marker")) |>
    dplyr::mutate(Gene_std = ifelse(is.na(Gene_std), Genes_text, Gene_std), Marker_std = ifelse(is.na(Marker_std), top_marker, Marker_std))
  
  readr::write_csv(tibble::as_tibble(smmry_sel_enriched), path_table_cluster_summary_raw)
  
  # 10. RAW SPECIES / ACCESSION MAPS
  species_cluster_map <- df_species_clusters_metadata |> dplyr::group_by(species) |> dplyr::summarise(clusters = paste(unique(cluster_id), collapse = ","), Genes = paste(unique(Genes_text), collapse = "; "), n_clusters = dplyr::n_distinct(cluster_id), .groups = "drop")
  species_sid_matrix <- df_species_clusters_metadata |> dplyr::select(species, cluster_id, sid) |> tidyr::pivot_wider(names_from = cluster_id, names_prefix = "cid_", values_from = sid, values_fill = "")
  species_table <- df_species_clusters_metadata |> dplyr::left_join(species_sid_matrix, by = "species")
  sid_conflicts <- df_species_clusters_metadata |> dplyr::group_by(sid) |> dplyr::summarise(n_clusters = dplyr::n(), clusters = paste(unique(cluster_id), collapse = ","), species = paste(unique(species), collapse = ","), .groups = "drop") |> dplyr::filter(n_clusters > 1)
  
  readr::write_csv(tibble::as_tibble(sid_conflicts), path_table_duplicate_conflicts)
  
  species_with_sid_dup <- df_species_clusters_metadata |> dplyr::filter(sid %in% sid_conflicts$sid) |> dplyr::distinct(species) |> dplyr::pull(species)
  df_species_clusters_metadata <- df_species_clusters_metadata |> dplyr::left_join(sid_conflicts |> dplyr::select(sid, duplicated = n_clusters), by = "sid") |> dplyr::mutate(duplicated = ifelse(is.na(duplicated), 0L, as.integer(duplicated)))
  species_table <- species_table |> dplyr::left_join(sid_conflicts |> dplyr::select(sid, duplicated = n_clusters), by = "sid") |> dplyr::mutate(duplicated = ifelse(is.na(duplicated), 0L, as.integer(duplicated)))
  species_cluster_map <- species_cluster_map |> dplyr::mutate(duplicated = ifelse(species %in% species_with_sid_dup, 1L, 0L))
  
  # 11. DUPLICATE-RESOLUTION CLEANING
  dup_sids <- sid_conflicts$sid
  dup_records <- df_species_clusters_metadata |>
    dplyr::filter(sid %in% dup_sids) |>
    dplyr::mutate(
      cluster_id = as.integer(cluster_id),
      prnt = vapply(cluster_id, function(cid) cp_get_cluster_parent(selected, cid), FUN.VALUE = character(1))
    )
  keepers <- dup_records |>
    dplyr::group_by(sid) |>
    dplyr::arrange(cluster_id, .by_group = TRUE) |>
    dplyr::mutate(
      has_preferred = any(prnt == "866800", na.rm = TRUE),
      keep_rank = dplyr::case_when(
        has_preferred & prnt == "866800" ~ 1L,
        has_preferred & prnt != "866800" ~ 2L,
        !has_preferred ~ 1L,
        TRUE ~ 3L
      )
    ) |>
    dplyr::arrange(keep_rank, cluster_id, .by_group = TRUE) |>
    dplyr::mutate(keep = dplyr::row_number() == 1L) |>
    dplyr::ungroup() |>
    dplyr::filter(keep) |>
    dplyr::select(cluster_id, sid) |>
    dplyr::distinct()
  rows_to_delete_auto <- dup_records |> dplyr::anti_join(keepers, by = c("cluster_id", "sid")) |> dplyr::select(cluster_id, sid) |> dplyr::distinct()
  
  manual_rows <- .cp_read_manual_exclusions(apply_manual_exclusions, manual_exclusions_file, path_table_manual_exclusions,
                                            df_species_clusters_metadata, .cp_cluster_loci(phylota, pattern, genes_map_df, marker_lookup),
                                            log_message)
  
  rows_to_delete <- dplyr::bind_rows(rows_to_delete_auto |> dplyr::mutate(reason = "auto"), manual_rows) |> dplyr::distinct(cluster_id, sid, .keep_all = TRUE)
  log_message("Automatic duplicate removals: ", nrow(rows_to_delete_auto))
  log_message("Manual duplicate removals: ", nrow(manual_rows))
  log_message("Total duplicate removals applied: ", nrow(rows_to_delete))
  
  # 12. CLEAN TABLES
  df_species_clusters_clean <- df_species_clusters_metadata |> dplyr::anti_join(rows_to_delete |> dplyr::select(cluster_id, sid), by = c("cluster_id", "sid"))
  species_table_clean <- species_table |> dplyr::mutate(cluster_id = as.integer(cluster_id)) |> dplyr::anti_join(rows_to_delete |> dplyr::select(cluster_id, sid), by = c("cluster_id", "sid"))
  species_cluster_map_clean <- df_species_clusters_clean |> dplyr::group_by(species) |> dplyr::summarise(clusters = paste(unique(cluster_id), collapse = ","), n_clusters = dplyr::n_distinct(cluster_id), Genes = paste(unique(Genes_text), collapse = "; "), .groups = "drop") |> dplyr::left_join(species_sid_matrix, by = "species")
  
  # 13. UPDATED SUMMARY COUNTS
  smmry_sel_enriched_fixed <- smmry_sel_enriched |> dplyr::mutate(cluster_id = as.integer(ID))
  n_original <- species_sid_matrix |> tidyr::pivot_longer(cols = dplyr::starts_with("cid_"), names_to = "cid_col", values_to = "sid", values_drop_na = FALSE) |> dplyr::mutate(cluster_id = as.integer(sub("^cid_", "", cid_col))) |> dplyr::filter(sid != "") |> dplyr::group_by(cluster_id) |> dplyr::summarise(n_original = dplyr::n(), .groups = "drop")
  n_removed <- rows_to_delete |> dplyr::group_by(cluster_id) |> dplyr::summarise(n_removed = dplyr::n(), .groups = "drop")
  n_final <- df_species_clusters_clean |> dplyr::group_by(cluster_id) |> dplyr::summarise(n_final = dplyr::n(), .groups = "drop")
  smmry_sel_enriched_updated <- smmry_sel_enriched_fixed |> dplyr::left_join(n_original, by = "cluster_id") |> dplyr::left_join(n_removed, by = "cluster_id") |> dplyr::left_join(n_final, by = "cluster_id") |> dplyr::mutate(n_original = ifelse(is.na(n_original), 0L, as.integer(n_original)), n_removed = ifelse(is.na(n_removed), 0L, as.integer(n_removed)), n_final = ifelse(is.na(n_final), 0L, as.integer(n_final)))
  
  # 14. CLEAN PHYLOTA OBJECT
  cleaned_clstrs <- selected@clstrs
  for (i in seq_along(cleaned_clstrs@clstrs)) {
    cid <- names(cleaned_clstrs@clstrs)[i]
    sids_remove <- rows_to_delete |> dplyr::filter(cluster_id == as.integer(cid)) |> dplyr::pull(sid)
    cleaned_clstrs@clstrs[[i]]@sids <- setdiff(cleaned_clstrs@clstrs[[i]]@sids, sids_remove)
  }
  selected_clean <- selected
  selected_clean@clstrs <- cleaned_clstrs
  
  clusters_to_keep_after_cleaning <- smmry_sel_enriched_updated |> dplyr::filter(n_final > 0) |> dplyr::pull(ID) |> as.character()
  phylota_final <- phylotaR::drop_clstrs(phylota = selected_clean, cid = clusters_to_keep_after_cleaning)
  log_message("Clusters retained after duplicate cleaning: ", length(phylota_final@cids))
  
  # 15. FINAL TABLES AFTER CLEANING
  final_cids <- phylota_final@cids
  ntaxa_final <- phylotaR::get_ntaxa(phylota = phylota_final, cid = final_cids, rnk = "species")
  ntaxa_final_df <- tibble::tibble(cluster_id = final_cids, n_taxa = ntaxa_final)
  df_species_clusters_final <- cp_extract_cluster_species_sid(phylota_final)
  metadata_final <- metadata_raw
  df_species_clusters_metadata_final <- df_species_clusters_final |> dplyr::left_join(metadata_final, by = "sid") |> cp_annotate_marker_text(pattern = pattern, genes_map_df = genes_map_df)
  smmry_final <- phylotaR::summary(phylota_final) |> dplyr::mutate(ID = as.integer(ID))
  cluster_gene_summary_final <- cp_summarise_cluster_markers(df_species_clusters_metadata_final) |> dplyr::mutate(cluster_id = as.integer(cluster_id))
  smmry_final_enriched <- smmry_final |> dplyr::left_join(cluster_gene_summary_final, by = c("ID" = "cluster_id"))
  seed_meta_final <- cp_fetch_seed_metadata(smmry_final_enriched$Seed, log_message = log_message)

  smmry_final_enriched <- smmry_final_enriched |>
    dplyr::mutate(Species = seed_meta_final$Species, Description = cp_normalize_marker(seed_meta_final$Description)) |>
    dplyr::mutate(
      Genes_raw = stringr::str_extract_all(Description, pattern),
      Genes_text = dplyr::if_else(
        is.na(Description), NA_character_,
        purrr::map_chr(Genes_raw, ~ paste(unique(cp_normalize_marker(.x)), collapse = ", "))
      ),
      top_marker = cp_normalize_marker(top_marker)
    ) |>
    dplyr::select(-Genes_raw) |>
    dplyr::left_join(gene_lookup, by = c("Genes_text" = "search_gene")) |> dplyr::left_join(marker_lookup, by = c("top_marker" = "search_marker")) |>
    dplyr::mutate(Gene_std = ifelse(is.na(Gene_std), Genes_text, Gene_std), Marker_std = ifelse(is.na(Marker_std), top_marker, Marker_std))
  
  species_sid_matrix_final <- df_species_clusters_metadata_final |> dplyr::select(species, cluster_id, sid) |> tidyr::pivot_wider(names_from = cluster_id, names_prefix = "cid_", values_from = sid, values_fill = "")
  species_genes_final <- df_species_clusters_metadata_final |> dplyr::group_by(species) |> dplyr::summarise(Genes_text = paste(unique(Genes_text), collapse = "; "), .groups = "drop")
  species_table_final <- species_sid_matrix_final |> dplyr::left_join(species_genes_final, by = "species") |> dplyr::relocate(Genes_text, .after = species)
  
  # 16. EXPORT CLEAN CLUSTER FASTAS
  cp_write_cluster_fastas(phylota_final, dir_out_cluster_fasta)
  
  # 17. MERGE CLUSTERS BY STANDARDIZED MARKER
  df_map <- smmry_final_enriched |> dplyr::mutate(ID = as.integer(ID)) |> dplyr::select(ID, Marker_std) |> dplyr::distinct()
  readr::write_csv(df_map, path_table_cluster_marker_assignment)
  
  cluster_fasta_files <- list.files(dir_out_cluster_fasta, full.names = TRUE, pattern = "^CLUSTER_.*\\.fasta$") |> sort()
  df_files <- tibble::tibble(file = cluster_fasta_files, ID = as.integer(sub("^CLUSTER_([0-9]+)\\.fasta$", "\\1", basename(cluster_fasta_files))))
  df_joined <- df_files |> dplyr::left_join(df_map, by = "ID")

  # Same guard as the ingroup assembly, and this branch previously had none at all: an unmapped
  # cluster was dropped in silence. See the note there for what that cost on 2026-09-04.
  df_joined <- df_joined |> dplyr::mutate(Marker_std = trimws(as.character(Marker_std)))
  unmapped <- df_joined |> dplyr::filter(is.na(Marker_std) | !nzchar(Marker_std))
  if (nrow(unmapped) > 0) {
    warning("Clusters with no standardized marker name, dropped from the marker FASTAs: IDs ",
            paste(sort(unique(unmapped$ID)), collapse = ", "),
            ". Their descriptions match no rule in genes_map.csv. Check the Description column of ",
            basename(path_table_cluster_summary_clean),
            " and add a rule if any of them is a locus this analysis should use.", call. = FALSE)
    log_message("Clusters without Marker_std, dropped: ", paste(sort(unique(unmapped$ID)), collapse = ", "))
  }

  df_joined <- df_joined |> dplyr::filter(!is.na(Marker_std), nzchar(Marker_std))
  markers_final <- sort(unique(df_joined$Marker_std))
  
  for (mk in markers_final) {
    log_message("Processing marker FASTA merge: ", mk)
    files_mk <- df_joined |> dplyr::filter(Marker_std == mk) |> dplyr::arrange(file) |> dplyr::pull(file)
    seqs_list <- lapply(files_mk, Biostrings::readDNAStringSet)
    seqs_all <- do.call(c, seqs_list)
    sp_names <- sub("^([^ ]+_[^ ]+).*", "\\1", names(seqs_all))
    seqs_clean <- seqs_all[!duplicated(sp_names)]
    out_file <- file.path(dir_out_base, paste0(mk, ".fasta"))
    Biostrings::writeXStringSet(seqs_clean, out_file)
  }
  
  # 18. FINAL OCCUPANCY AND MARKER SUMMARIES
  df_species_clusters_metadata_final <- df_species_clusters_metadata_final |> dplyr::mutate(cluster_id = as.integer(cluster_id))
  df_joined <- df_joined |> dplyr::mutate(ID = as.integer(ID))
  smmry_final_enriched <- smmry_final_enriched |> dplyr::mutate(ID = as.integer(as.character(ID)))
  occ_table <- df_species_clusters_metadata_final |> dplyr::left_join(df_joined |> dplyr::select(ID, Marker_std), by = c("cluster_id" = "ID"))
  cluster_level_stats <- smmry_final_enriched |> dplyr::mutate(cluster_id = as.integer(ID)) |> dplyr::select(cluster_id, Marker_std, n_sequences)
  
  marker_summary <- cluster_level_stats |>
    dplyr::group_by(Marker_std) |>
    dplyr::summarise(num_clusters = dplyr::n(), num_sequences = sum(n_sequences, na.rm = TRUE), mean_seq_per_cluster = mean(n_sequences), median_seq_per_cluster = stats::median(n_sequences), min_seq_in_cluster = min(n_sequences), max_seq_in_cluster = max(n_sequences), .groups = "drop") |>
    dplyr::left_join(occ_table |> dplyr::group_by(Marker_std) |> dplyr::summarise(num_species = dplyr::n_distinct(Species_gb), .groups = "drop"), by = "Marker_std") |>
    dplyr::select(Marker_std, num_clusters, num_species, num_sequences, mean_seq_per_cluster, median_seq_per_cluster, min_seq_in_cluster, max_seq_in_cluster) |>
    dplyr::arrange(Marker_std)
    
  # 19. EXPORT TABLES
  readr::write_csv(tibble::as_tibble(smmry_final_enriched), path_table_cluster_summary_clean)
  readr::write_csv(tibble::as_tibble(species_table_final), path_table_species_cluster_map_clean)
  readr::write_csv(tibble::as_tibble(occ_table), path_table_accession_occupancy_clean)
  readr::write_csv(tibble::as_tibble(marker_summary), path_table_marker_summary)
  
  # 20. SAVE OBJECTS
  save(phylota_final, smmry_final_enriched, species_table_final, occ_table, marker_summary, file = path_phylota_clean_rdata)
  log_message("Assembly and export fully completed. \U0001f335")
  
  res <- list(
    selected_clusters = phylota_final,
    retained_cids = phylota_final@cids
  )
  message("\nOutgroup assembly and export fully completed. \U0001f335")
  return(res)
}

#' Fetch GenBank Sequence Metadata via NCBI Entrez Utilities
#'
#' Retrieves the organism name and the definition line of each GenBank sequence identifier (SID)
#' with `ape::read.GenBank()`, in batches with retries, and caches the result.
#'
#' @param sids Character vector of GenBank Sequence Identifiers (SIDs) to query.
#' @param cache_file Character. File path to store and load cached metadata tables.
#' @param batch_size Integer. Number of sequence IDs requested per HTTP batch query. Defaults to `200`.
#' @param sleep_time Numeric. Pause duration in seconds between consecutive batch requests to respect NCBI rate limits. Defaults to `0.5`.
#' @param max_retries Integer. Maximum retry attempts permitted per batch before failing. Defaults to `5`.
#' @param force_download Logical. Not used in this version: records already in `cache_file` are always reused. Defaults to `FALSE`.
#' @return A data frame with columns `sid`, `Species_gb` and `Description_gb`.
#' @examples
#' \dontrun{
#' fetch_genbank_metadata(
#'   sids = c("AY123456", "AY123457"),
#'   cache_file = "cache_metadata.csv"
#' )
#' }
#' @export
fetch_genbank_metadata <- function(sids, cache_file, batch_size = 200, sleep_time = 0.5, max_retries = 5, force_download = FALSE) {
  cp_download_all_metadata(sids, cache_file, batch_size, sleep_time, max_retries, log_message = message)
}

# --- Shared Utilities (Unexported) ---

cp_clean_species_name <- function(x) {
  trimws(gsub("\\s+", "_", gsub("\\.", "", x)))
}

cp_normalize_marker <- function(x) {
  stringr::str_replace_all(stringr::str_replace_all(trimws(x), "\\s+", " "), "-|-", "-")
}

cp_build_pattern <- function(markers) {
  paste(stringr::str_replace_all(markers[order(nchar(markers), decreasing = TRUE)], "([.|()\\^{}+$*?\\[\\]\\\\-])", "\\\\\\1"), collapse = "|")
}

cp_get_cluster_parent <- function(phylota_obj, cid) {
  cid_chr <- as.character(cid)
  if (!cid_chr %in% names(phylota_obj@clstrs@clstrs)) return(NA_character_)
  as.character(phylota_obj@clstrs@clstrs[[cid_chr]]@prnt)
}

cp_download_all_metadata <- function(sid_vector, cache_csv, batch_size = 200L, sleep_time = 0.5, max_retries = 5L, log_message = function(...) {}) {
  sid_vector <- sort(unique(as.character(sid_vector)))
  if (length(sid_vector) == 0) return(tibble::tibble(sid = character(), Species_gb = character(), Description_gb = character()))
  
  cached <- tibble::tibble(sid = character(), Species_gb = character(), Description_gb = character())
  if (file.exists(cache_csv)) {
    cached <- suppressMessages(readr::read_csv(cache_csv, show_col_types = FALSE)) |> dplyr::distinct(sid, .keep_all = TRUE)
  }
  
  missing_sids <- setdiff(sid_vector, cached$sid)
  new_metadata <- list()
  if (length(missing_sids) > 0) {
    blocks <- split(missing_sids, ceiling(seq_along(missing_sids) / batch_size))
    for (i in seq_along(blocks)) {
      Sys.sleep(sleep_time)
      pause <- sleep_time
      gb <- NULL
      for (att in seq_len(max_retries)) {
        gb <- tryCatch({ suppressWarnings(ape::read.GenBank(blocks[[i]])) }, error = function(e) NULL)
        if (!is.null(gb) && !is.null(attr(gb, "species"))) break
        Sys.sleep(pause)
        pause <- pause * 2
      }
      if (!is.null(gb)) {
        new_metadata[[i]] <- tibble::tibble(sid = names(gb), Species_gb = attr(gb, "species"), Description_gb = attr(gb, "description"))
      }
    }
  }
  
  all_metadata <- dplyr::bind_rows(cached, dplyr::bind_rows(new_metadata)) |> dplyr::distinct(sid, .keep_all = TRUE)
  if(nrow(all_metadata) > 0) readr::write_csv(all_metadata, cache_csv)
  
  all_metadata |> dplyr::filter(sid %in% sid_vector)
}

cp_fetch_seed_metadata <- function(sids, batch_size = 200L, sleep_time = 0.5, max_retries = 5L, log_message = function(...) {}) {
  sids <- as.character(sids)
  out <- tibble::tibble(Seed = sids, Species = NA_character_, Description = NA_character_)
  if (length(sids) == 0) return(out)

  blocks <- split(sids, ceiling(seq_along(sids) / batch_size))
  fetched <- list()
  for (i in seq_along(blocks)) {
    if (i > 1) Sys.sleep(sleep_time)
    pause <- sleep_time
    gb <- NULL
    for (att in seq_len(max_retries)) {
      gb <- tryCatch({ suppressWarnings(ape::read.GenBank(blocks[[i]])) }, error = function(e) NULL)
      if (!is.null(gb) && !is.null(attr(gb, "species"))) break
      log_message("Seed metadata batch ", i, "/", length(blocks), " attempt ", att, " failed, retrying...")
      Sys.sleep(pause)
      pause <- pause * 2
    }
    if (!is.null(gb)) {
      fetched[[i]] <- tibble::tibble(Seed = names(gb), Species = attr(gb, "species"), Description = attr(gb, "description"))
    } else {
      log_message("Seed metadata batch ", i, "/", length(blocks), " unrecovered after ", max_retries, " attempts (", length(blocks[[i]]), " accessions).")
    }
  }

  fetched_df <- dplyr::bind_rows(fetched)
  if (nrow(fetched_df) > 0) {
    out <- out |>
      dplyr::select(Seed) |>
      dplyr::left_join(dplyr::distinct(fetched_df, Seed, .keep_all = TRUE), by = "Seed")
  }
  out
}

cp_summarise_cluster_markers <- function(df_cluster_metadata) {
  df_cluster_metadata |>
    dplyr::group_by(cluster_id) |>
    dplyr::summarise(
      all_markers = paste(sort(unique(Genes_text)), collapse = "; "),
      n_markers = dplyr::n_distinct(Genes_text),
      top_marker = {
        valid_genes <- Genes_text[!is.na(Genes_text) & Genes_text != ""]
        if (length(valid_genes) > 0) {
          tab <- table(valid_genes)
          names(which.max(tab))
        } else {
          ""
        }
      },
      top_marker_freq = {
        valid_genes <- Genes_text[!is.na(Genes_text) & Genes_text != ""]
        if (length(valid_genes) > 0) {
          tab <- table(valid_genes)
          max(tab) / length(Genes_text)
        } else {
          NA_real_
        }
      },
      n_sequences = dplyr::n(),
      n_species = dplyr::n_distinct(species),
      .groups = "drop"
    )
}

cp_extract_cluster_species_sid <- function(phylota_obj) {
  purrr::map_dfr(phylota_obj@cids, function(cid) {
    sids <- phylota_obj@clstrs[[cid]]@sids
    txids <- sapply(sids, function(sid) phylota_obj@sqs[[sid]]@txid)
    species <- cp_clean_species_name(phylotaR::get_tx_slot(phylota_obj, txid = txids, slt_nm = "scnm"))
    tibble::tibble(cluster_id = as.integer(cid), species = species, sid = as.character(sids))
  })
}

cp_annotate_marker_text <- function(df, pattern, genes_map_df) {
  out <- df |>
    dplyr::mutate(
      Description_gb = cp_normalize_marker(ifelse(is.na(Description_gb), "", Description_gb)),
      Genes_raw = stringr::str_extract_all(Description_gb, pattern),
      Genes_text = purrr::map_chr(Genes_raw, ~ paste(unique(cp_normalize_marker(.x)), collapse = ", "))
    ) |> dplyr::select(-Genes_raw)
  out |>
    dplyr::left_join(dplyr::select(genes_map_df, search, replace), by = c("Genes_text" = "search")) |>
    dplyr::mutate(Genes_text = ifelse(is.na(replace), Genes_text, replace)) |>
    dplyr::select(-replace)
}

cp_write_cluster_fastas <- function(phylota_obj, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (cid in sort(as.integer(phylota_obj@cids))) {
    sids <- phylota_obj@clstrs[[as.character(cid)]]@sids
    txids <- sapply(sids, function(sid) phylota_obj@sqs[[sid]]@txid)
    scientific_names <- cp_clean_species_name(phylotaR::get_tx_slot(phylota_obj, txid = txids, slt_nm = "scnm"))
    if (length(scientific_names) != length(sids)) {
      scientific_names <- scientific_names[seq_along(sids)]
    }
    phylotaR::write_sqs(phylota_obj, sid = sids, sq_nm = scientific_names, outfile = file.path(outdir, paste0("CLUSTER_", cid, ".fasta")))
  }
}

# ------------------------------------------------------------------------------
# Shared by assemble_ingroup_phylotar() and assemble_barcoding_dataset(). Moved here in 0.5.0
# without changing their behaviour, so that the phylogeny and the molecular diagnostic branch make
# the same phylotaR call on the same workspace and assign clusters to loci with the same code.
# ------------------------------------------------------------------------------

#' phylotaR search term shared by both branches
#'
#' Excludes predicted, unverified, whole-genome shotgun, synthetic, RefSeq and TSA records, and
#' records named at open nomenclature (sp., aff., cf.), below species (var., subsp.) or as hybrids.
#' @return Character string passed to `phylotaR::setup(srch_trm = )`.
#' @noRd
.phylotar_search_terms <- function() {
  paste0(
    "NOT predicted[TI] ",
    "NOT \"whole genome shotgun\"[TI] ",
    "NOT unverified[TI] ",
    "NOT \"synthetic construct\"[Organism] ",
    "NOT refseq[filter] ",
    "NOT TSA[Keyword] ",
    "NOT \"sp.\"[TI] ",
    "NOT \"sp.\"[Organism] ",
    "NOT \"sp\"[Organism] ",
    "NOT \"sp\"[Organism] ",
    "NOT \"aff.\"[TI] ",
    "NOT \"aff\"[Organism] ",
    "NOT \"cf.\"[TI] ",
    "NOT \"cf\"[Organism] ",
    "NOT \"var.\"[TI] ",
    "NOT \"var\"[TI] ",
    "NOT \"var\"[Organism] ",
    "NOT \"var.\"[Organism] ",
    "NOT \"variety\"[TI] ",
    "NOT \"subsp.\"[TI] ",
    "NOT \"subsp\"[TI] ",
    "NOT \"subsp.\"[Organism] ",
    "NOT \"subsp\"[Organism] ",
    "NOT \"subspecies\"[Organism] ",
    "NOT \"x\"[Organism] ",
    "NOT \" x \"[Organism]"
  )
}

#' Read a phylotaR workspace, or mine GenBank into it when it does not exist
#'
#' Reads the workspace with `phylotaR::read_phylota()`. When that fails, or when
#' `force_download = TRUE`, runs `phylotaR::setup()` and `phylotaR::run()` with the parameters of
#' the phylogeny and reads the result. The reader and the two pipeline steps are arguments only so
#' that the tests can stand in for them.
#' @noRd
.phylotar_load_or_mine <- function(wd_path, preferred_parent = "3593", txid = NULL, ncbi_dr = NULL,
                                   force_download = FALSE, log_message = function(...) {},
                                   mnsql = 100L, mxsql = 5000L,
                                   notify = FALSE, notify_to = NULL, notify_credentials = NULL,
                                   notify_fn = send_run_notification,
                                   reader = phylotaR::read_phylota,
                                   setup_fn = phylotaR::setup,
                                   run_fn = phylotaR::run) {
  # What is mined and which parent wins when a sid sits in two clusters are two different things.
  # They coincided while only the ingroup used this, with its single focal clade. The outgroup of
  # CN2 is several clades and has no preferred parent, so `txid` carries what to mine and defaults
  # to `preferred_parent`, which leaves every existing call doing exactly what it did.
  if (is.null(txid)) txid <- preferred_parent
  mine <- function() {
    log_message("No valid phylota object found or force_download=TRUE. Running setup and phylotaR pipeline...")
    if (is.null(ncbi_dr)) {
      env_blast <- Sys.getenv("BLAST_PATH")
      if (env_blast != "") {
        ncbi_dr <- env_blast
      } else {
        blastn_path <- Sys.which("blastn")
        if (blastn_path != "") ncbi_dr <- dirname(blastn_path)
      }
    }

    setup_fn(
      wd = wd_path, txid = txid, ncbi_dr = ncbi_dr, v = TRUE, ncps = 1, mncvrg = 80,
      mnsql = mnsql, mxsql = mxsql, srch_trm = .phylotar_search_terms()
    )

    log_message("Executing phylotaR::run()...")
    tryCatch({
      run_fn(wd = wd_path)
    }, error = function(erun) {
      log_message("Error in phylotaR::run(): ", erun$message)
    })

    log_message("Attempting to load phylota object after run()...")
    tryCatch({
      reader(wd_path)
    }, error = function(e2) {
      log_message("read_phylota error: ", e2$message)
      stop("Could not load phylota object from ", wd_path)
    })
  }
  tryCatch({
    if (force_download) stop("Force download enabled")
    reader(wd_path)
  }, error = function(e) {
    if (!isTRUE(notify)) return(mine())
    # L4 (BMM, 28-09): the mining takes hours on the Mac; it says when it ends or fails
    started <- Sys.time()
    analysis <- paste0("phylotaR mining of ", paste(txid, collapse = ", "))
    res <- tryCatch(mine(), error = function(err) {
      note <- .compose_run_notification(analysis, status = "failed", started = started,
                                        outputs = c(Workspace = wd_path), error_message = conditionMessage(err))
      notify_fn(note$subject, note$body, to = notify_to, credentials = notify_credentials)
      stop(err)
    })
    note <- .compose_run_notification(analysis, status = "finished", started = started, outputs = c(Workspace = wd_path))
    notify_fn(note$subject, note$body, to = notify_to, credentials = notify_credentials)
    res
  })
}

#' Assign each cluster its standardised gene and marker names
#'
#' Adds the species and description of each cluster seed, extracts the gene names the description
#' mentions, and maps the dominant marker of the cluster (`top_marker`) through `genes_map.csv`.
#' A marker that `genes_map.csv` does not map keeps its own name.
#' @param smmry Cluster summary with `Seed` and `top_marker`, one row per cluster.
#' @param seed_meta Seed metadata (`Species`, `Description`), one row per row of `smmry`, same order.
#' @noRd
.cp_enrich_cluster_summary <- function(smmry, seed_meta, pattern, gene_lookup, marker_lookup) {
  smmry |>
    dplyr::mutate(
      Species = seed_meta$Species,
      Description = cp_normalize_marker(seed_meta$Description)
    ) |>
    dplyr::mutate(
      Genes_raw = stringr::str_extract_all(Description, pattern),
      Genes_text = dplyr::if_else(
        is.na(Description), NA_character_,
        purrr::map_chr(Genes_raw, ~ paste(unique(cp_normalize_marker(.x)), collapse = ", "))
      ),
      top_marker = cp_normalize_marker(top_marker)
    ) |>
    dplyr::select(-Genes_raw) |>
    dplyr::left_join(gene_lookup, by = c("Genes_text" = "search_gene")) |>
    dplyr::left_join(marker_lookup, by = c("top_marker" = "search_marker")) |>
    dplyr::mutate(
      Gene_std = ifelse(is.na(Gene_std), Genes_text, Gene_std),
      Marker_std = ifelse(is.na(Marker_std), top_marker, Marker_std)
    )
}

#' The locus of every cluster, read from the definition lines of its sequences
#'
#' The same reading as the naming of the kept clusters (`cp_annotate_marker_text()`, the dominant
#' marker of `cp_summarise_cluster_markers()`, `marker_lookup`), on the definition lines the
#' workspace already holds, so that every cluster can be named before the cut of min_species
#' (decision Y1 of BMM, 29-09) without downloading metadata.
#' @return A data frame with `cluster_id` and `locus` (`NA` when no gene is recognised).
#' @noRd
.cp_cluster_loci <- function(phylota, pattern, genes_map_df, marker_lookup, variants = .cp_locus_variants()) {
  cids <- phylota@cids
  if (length(cids) == 0) return(data.frame(cluster_id = integer(0), locus = character(0)))
  df <- do.call(rbind, lapply(cids, function(cid) {
    sids <- phylota@clstrs[[cid]]@sids
    data.frame(cluster_id = as.integer(cid), sid = sids, species = "",
               Description_gb = vapply(sids, function(s) phylota@sqs[[s]]@dfln, "", USE.NAMES = FALSE),
               stringsAsFactors = FALSE)
  }))
  s <- cp_summarise_cluster_markers(cp_annotate_marker_text(df, pattern = pattern, genes_map_df = genes_map_df))
  top <- cp_normalize_marker(s$top_marker)
  locus <- marker_lookup$Marker_std[match(top, marker_lookup$search_marker)]
  locus <- ifelse(is.na(locus), top, locus)
  locus[!nzchar(locus)] <- NA_character_
  .cp_apply_locus_variants(data.frame(cluster_id = as.integer(s$cluster_id), locus = locus, stringsAsFactors = FALSE),
                           variants)
}

#' Variant names of a locus and the locus they are pooled with (decision Y3 of BMM, 29-09)
#'
#' The table fixed before the probe X3 ran (`inst/extdata/locus_name_variants.csv`): names that
#' `.cp_cluster_loci()` gives to clusters of a library locus (`trnK` for `matK`, `trnL-rpl32` for
#' `rpl32-trnL`...). Applied after `genes_map.csv`, which is not changed.
#' @noRd
.cp_locus_variants <- function() {
  utils::read.csv(system.file("extdata", "locus_name_variants.csv", package = "PhyloCactus"),
                  stringsAsFactors = FALSE, colClasses = "character")
}

#' A variant name takes the name of its locus; other names and `NA` are kept
#' @noRd
.cp_apply_locus_variants <- function(cluster_locus, variants) {
  if (is.null(variants) || !nrow(variants)) return(cluster_locus)
  hit <- match(cluster_locus$locus, variants$name)
  cluster_locus$locus[!is.na(hit)] <- variants$locus[hit[!is.na(hit)]]
  cluster_locus
}

#' Clusters over the cut of min_species, the clusters of one locus pooled (decision Y1 of BMM, 29-09)
#'
#' A locus is kept when its clusters together hold more than `min_species` species; a cluster with
#' no locus name is judged alone, as before.
#' @param cluster_species Data frame with `cluster_id` and `species` (one row per species of a
#'   cluster, counted on the workspace reduced to one sequence per species).
#' @param cluster_locus Data frame with `cluster_id` and `locus`.
#' @return A list with `keep` (cluster ids) and `table` (per locus or unnamed cluster: clusters and
#'   species).
#' @noRd
.cp_select_clusters_by_locus <- function(cluster_species, cluster_locus, min_species, overlap = NULL) {
  fragments <- data.frame(cluster_id = integer(0), locus = character(0))
  if (!is.null(overlap)) {
    # Y2 (BMM, 29-09): a named cluster that does not overlap the main cluster of its locus is not pooled
    out <- overlap$cluster_id[!overlap$overlaps_main]
    loc_all <- cluster_locus$locus[match(out, cluster_locus$cluster_id)]
    out <- out[!is.na(loc_all) & nzchar(loc_all)]
    fragments <- data.frame(cluster_id = as.integer(out),
                            locus = cluster_locus$locus[match(out, cluster_locus$cluster_id)], stringsAsFactors = FALSE)
    cluster_species <- cluster_species[!cluster_species$cluster_id %in% out, , drop = FALSE]
  }
  loc <- cluster_locus$locus[match(cluster_species$cluster_id, cluster_locus$cluster_id)]
  named <- !is.na(loc) & nzchar(loc)
  unit <- ifelse(named, loc, paste0("cluster_", cluster_species$cluster_id))
  sp <- tapply(cluster_species$species, unit, function(x) length(unique(x)))
  cl <- tapply(cluster_species$cluster_id, unit, function(x) length(unique(x)))
  table <- data.frame(locus = names(sp), named = !startsWith(names(sp), "cluster_") | names(sp) %in% loc[named],
                      clusters = as.integer(cl[names(sp)]), species = as.integer(sp), stringsAsFactors = FALSE)
  table$kept <- table$species > min_species
  keep <- sort(unique(cluster_species$cluster_id[unit %in% table$locus[table$kept]]))
  list(keep = as.integer(keep), table = table, fragments = fragments)
}

#' Whether a cluster overlaps the main cluster of its locus (Y2)
#'
#' TRUE when at least `min_fraction` of `seqs` share `min_shared` or more 20-mers with the sequences
#' of the main cluster, on either strand.
#' @noRd
.cp_overlaps_main <- function(seqs, main_seqs, k = 20L, min_shared = 5L, min_fraction = 0.5) {
  kmers <- function(x) { x <- toupper(x); n <- nchar(x); if (n < k) character(0) else unique(substring(x, 1:(n - k + 1L), k:n)) }
  rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
  pool <- unique(unlist(lapply(c(main_seqs, rc(main_seqs)), kmers), use.names = FALSE))
  shared <- vapply(seqs, function(x) sum(kmers(x) %in% pool), 0L)
  length(seqs) > 0 && mean(shared >= min_shared) >= min_fraction
}

#' For every cluster, whether it overlaps the main cluster (most species) of its locus (Y2)
#'
#' Unnamed clusters are judged alone and marked TRUE. At most `max_seqs` sequences of a cluster are
#' read, the first by sid, so that the check stays fast on large clusters.
#' @noRd
.cp_cluster_overlap <- function(phylota, cluster_locus, cluster_species, max_seqs = 100L) {
  seqs_of <- function(cid) {
    sids <- sort(phylota@clstrs[[as.character(cid)]]@sids, method = "radix")
    sids <- utils::head(sids, max_seqs)
    vapply(sids, function(s) rawToChar(phylota@sqs[[s]]@sq), "", USE.NAMES = FALSE)
  }
  n_sp <- tapply(cluster_species$species, cluster_species$cluster_id, function(x) length(unique(x)))
  out <- data.frame(cluster_id = as.integer(cluster_locus$cluster_id), overlaps_main = TRUE)
  loci <- unique(stats::na.omit(cluster_locus$locus[nzchar(cluster_locus$locus)]))
  for (l in loci) {
    cids <- cluster_locus$cluster_id[cluster_locus$locus %in% l]
    if (length(cids) < 2L) next
    sp <- as.integer(n_sp[as.character(cids)]); sp[is.na(sp)] <- 0L
    main <- cids[order(-sp, cids)][1]
    main_seqs <- seqs_of(main)
    for (cid in setdiff(cids, main)) out$overlaps_main[out$cluster_id == cid] <- .cp_overlaps_main(seqs_of(cid), main_seqs)
  }
  out
}

#' The curated exclusions keyed by (locus, sid) as the (cluster, sid) pairs of this workspace (X2)
#'
#' An accession is removed from every cluster of the locus named in the list; its clusters of other
#' loci are left alone.
#' @return A list with `pairs` (cluster_id, sid, reason), `unmatched` (rows of the list with no
#'   pair) and `counts` (`n_file`, `n_matched` rows of the list, `n_pairs`).
#' @noRd
.cp_exclusion_pairs <- function(exclusions, records, cluster_locus) {
  rec <- data.frame(cluster_id = as.integer(records$cluster_id), sid = as.character(records$sid), stringsAsFactors = FALSE)
  rec$locus <- cluster_locus$locus[match(rec$cluster_id, cluster_locus$cluster_id)]
  key_exc <- paste(exclusions$locus, exclusions$sid)
  hit <- paste(rec$locus, rec$sid) %in% key_exc
  pairs <- rec[hit, c("cluster_id", "sid"), drop = FALSE]
  pairs$reason <- exclusions$reason[match(paste(rec$locus[hit], rec$sid[hit]), key_exc)]
  pairs <- pairs[!duplicated(paste(pairs$cluster_id, pairs$sid)), , drop = FALSE]
  rownames(pairs) <- NULL
  matched <- key_exc %in% paste(rec$locus, rec$sid)
  list(pairs = pairs, unmatched = exclusions[!matched, , drop = FALSE],
       counts = list(n_file = nrow(exclusions), n_matched = sum(matched), n_pairs = nrow(pairs)))
}

#' Species of every cluster on the workspace reduced to one sequence per species
#' @noRd
.cp_cluster_species <- function(species_reduced) {
  do.call(rbind, lapply(species_reduced@cids, function(cid) {
    tx <- unique(as.character(phylotaR::get_txids(species_reduced, cid = cid, rnk = "species")))
    data.frame(cluster_id = rep(as.integer(cid), length(tx)), species = tx, stringsAsFactors = FALSE)
  }))
}

#' Read the curated exclusions, keyed by (locus, sid), as the (cluster, sid) pairs to remove (X2)
#'
#' The list is written to `path_table` as read; the log says how many of its rows are found in this
#' workspace, so that a run reports the curation it actually applied.
#' @noRd
.cp_read_manual_exclusions <- function(apply, file, path_table, records, cluster_locus, log_message) {
  empty <- tibble::tibble(cluster_id = integer(), sid = character(), reason = character())
  if (!isTRUE(apply)) {
    log_message("Manual exclusions DISABLED (apply_manual_exclusions = FALSE): the expert-curated accession list is not applied.")
    return(empty)
  }
  if (is.null(file) || !file.exists(file)) return(empty)
  exc <- utils::read.csv(file, stringsAsFactors = FALSE)
  utils::write.csv(exc, path_table, row.names = FALSE)
  mp <- .cp_exclusion_pairs(exc, records, cluster_locus)
  log_message("Manual exclusions: ", mp$counts$n_file, " (locus, accession) rows in the list; ", mp$counts$n_matched,
              " found in this workspace, as ", mp$counts$n_pairs, " cluster memberships.")
  tibble::as_tibble(mp$pairs) |> dplyr::mutate(cluster_id = as.integer(cluster_id))
}
