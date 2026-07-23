#' Assemble Ingroup Sequence Clusters via phylotaR Similarity Mining
#'
#' Retrieves orthologous sequence clusters for a focal taxonomic ingroup (e.g., family **Cactaceae**, NCBI Taxonomy ID: 3593)
#' directly from GenBank using similarity clustering via `phylotaR` (Bennett *et al.*, 2018).
#' Relying on sequence similarity rather than inconsistent locus annotations prevents missing orthologous sequence data
#' caused by gene synonymy or mislabeling in public sequence repositories.
#'
#' @param wd_path Character. Path to the `phylotaR` workspace directory storing local parameters and database caches.
#' @param target_genes_file Character. Path to the target locus list text file. If `NULL`, defaults to package `inst/extdata/target_genes.txt`.
#' @param genes_map_file Character. Path to the gene synonymy mapping CSV file. If `NULL`, defaults to package `inst/extdata/genes_map.csv`.
#' @param manual_exclusions_file Character. Path to the accession exclusion CSV file. If `NULL`, defaults to package `inst/extdata/manual_exclusions_ingroup.csv`.
#' @param min_species Integer. Minimum number of distinct species required to retain an orthologous sequence cluster. Defaults to `50`.
#' @param preferred_parent Character. NCBI Taxonomy ID of the focal ingroup parent node. Defaults to `"3593"` (Cactaceae).
#' @param ncbi_dr Character. Path to local `BLAST+` binaries directory. If `NULL`, attempts system environment auto-detection.
#' @param force_download Logical. Force fresh database retrieval instead of using local cache? Defaults to `FALSE`.
#' @param out_dir Character. Output directory path to save cluster summaries, FASTA sequence matrices, and log reports.
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
assemble_ingroup_phylotar <- function(wd_path, target_genes_file = NULL, genes_map_file = NULL, manual_exclusions_file = NULL, min_species = 50, preferred_parent = "3593", ncbi_dr = NULL, force_download = FALSE, out_dir = "1_phylotaR_out_Ingroup") {
  
  
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
  markers_df <- read.table(target_genes_file, header = TRUE, stringsAsFactors = FALSE)
  markers <- unique(cp_normalize_marker(as.character(markers_df$markers)))
  pattern <- cp_build_pattern(markers)
  
  genes_map_df <- read.csv(genes_map_file, stringsAsFactors = FALSE) |>
    dplyr::mutate(search = cp_normalize_marker(search), replace = trimws(replace)) |>
    dplyr::distinct(search, .keep_all = TRUE)
    
  gene_lookup <- genes_map_df |> dplyr::select(search_gene = search, Gene_std = replace)
  marker_lookup <- genes_map_df |> dplyr::select(search_marker = search, Marker_std = replace)
  
  log_message("Inputs loaded successfully.")
  
  # Setup / Load phylotaR
  phylota <- tryCatch({
    if (force_download) stop("Force download enabled")
    phylotaR::read_phylota(wd_path)
  }, error = function(e) {
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
    
    phylotaR::setup(
      wd = wd_path, txid = preferred_parent, ncbi_dr = ncbi_dr, v = TRUE, ncps = 1, mncvrg = 80,
      srch_trm = paste0(
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
    )
    
    log_message("Executing phylotaR::run()...")
    tryCatch({
      phylotaR::run(wd = wd_path)
    }, error = function(erun) {
      log_message("Error in phylotaR::run(): ", erun$message)
    })
    
    log_message("Attempting to load phylota object after run()...")
    tryCatch({
      phylotaR::read_phylota(wd_path)
    }, error = function(e2) {
      log_message("read_phylota error: ", e2$message)
      stop("Could not load phylota object from ", wd_path)
    })
  })
  
  # 8. SPECIES REDUCTION AND CLUSTER FILTERING
  species_reduced <- phylotaR::drop_by_rank(phylota, rnk = "species", n = 1)
  cluster_ids <- species_reduced@cids
  ntaxa <- phylotaR::get_ntaxa(species_reduced, cid = cluster_ids, rnk = "species")
  keep_clusters <- cluster_ids[ntaxa > min_species]
  selected <- phylotaR::drop_clstrs(species_reduced, cid = keep_clusters)
  
  log_message("Clusters before >", min_species, " filter: ", length(cluster_ids))
  log_message("Clusters retained after >", min_species, " filter: ", length(selected@cids))
  
  # 9. RAW CLUSTER TABLES AND METADATA
  df_species_clusters <- cp_extract_cluster_species_sid(selected)
  metadata_raw <- cp_download_all_metadata(unique(df_species_clusters$sid), path_metadata_cache_csv, batch_size = 200, sleep_time = 0.5, max_retries = 5, log_message = log_message)
  
  df_species_clusters_metadata <- df_species_clusters |>
    dplyr::left_join(metadata_raw, by = "sid") |>
    cp_annotate_marker_text(pattern = pattern, genes_map_df = genes_map_df)
    
  smmry_sel <- phylotaR::summary(selected) |> dplyr::mutate(ID = as.integer(ID))
  cluster_gene_summary <- cp_summarise_cluster_markers(df_species_clusters_metadata) |> dplyr::mutate(cluster_id = as.integer(cluster_id))
  smmry_sel_enriched <- smmry_sel |> dplyr::left_join(cluster_gene_summary, by = c("ID" = "cluster_id"))
  
  seed_data <- tryCatch({ ape::read.GenBank(smmry_sel_enriched$Seed) }, error = function(e) NULL)
  
  if (!is.null(seed_data)) {
    smmry_sel_enriched <- smmry_sel_enriched |>
      dplyr::mutate(
        Species = attr(seed_data, "species"),
        Description = attr(seed_data, "description"),
        Description = cp_normalize_marker(Description)
      ) |>
      dplyr::mutate(
        Genes_raw = stringr::str_extract_all(Description, pattern),
        Genes_text = purrr::map_chr(Genes_raw, ~ paste(unique(cp_normalize_marker(.x)), collapse = ", ")),
        top_marker = cp_normalize_marker(top_marker)
      ) |>
      dplyr::select(-Genes_raw) |>
      dplyr::left_join(gene_lookup, by = c("Genes_text" = "search_gene")) |>
      dplyr::left_join(marker_lookup, by = c("top_marker" = "search_marker")) |>
      dplyr::mutate(
        Gene_std = ifelse(is.na(Gene_std), Genes_text, Gene_std),
        Marker_std = ifelse(is.na(Marker_std), top_marker, Marker_std)
      )
  } else {
    smmry_sel_enriched$Species <- NA
    smmry_sel_enriched$Description <- NA
    smmry_sel_enriched$Genes_text <- NA
    smmry_sel_enriched$Gene_std <- NA
    smmry_sel_enriched$Marker_std <- NA
  }
  
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
    
  manual_rows <- tibble::tibble(cluster_id = integer(), sid = character(), reason = character())
  if (!is.null(manual_exclusions_file) && file.exists(manual_exclusions_file)) {
    manual_rows <- readr::read_csv(manual_exclusions_file, show_col_types = FALSE) |> dplyr::mutate(cluster_id = as.integer(cluster_id))
    readr::write_csv(manual_rows, path_table_manual_exclusions)
  }
  
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
  ntaxa_final <- get_ntaxa(phylota = phylota_final, cid = final_cids, rnk = "species")
  ntaxa_final_df <- tibble::tibble(cluster_id = final_cids, n_taxa = ntaxa_final)
  
  df_species_clusters_final <- cp_extract_cluster_species_sid(phylota_final)
  metadata_final <- metadata_raw
  
  df_species_clusters_metadata_final <- df_species_clusters_final |>
    dplyr::left_join(metadata_final, by = "sid") |> cp_annotate_marker_text(pattern = pattern, genes_map_df = genes_map_df)
    
  smmry_final <- phylotaR::summary(phylota_final) |> dplyr::mutate(ID = as.integer(ID))
  cluster_gene_summary_final <- cp_summarise_cluster_markers(df_species_clusters_metadata_final) |> dplyr::mutate(cluster_id = as.integer(cluster_id))
  smmry_final_enriched <- smmry_final |> dplyr::left_join(cluster_gene_summary_final, by = c("ID" = "cluster_id"))
  
  seed_data_final <- tryCatch({ ape::read.GenBank(smmry_final_enriched$Seed) }, error = function(e) NULL)
  if (!is.null(seed_data_final)) {
    smmry_final_enriched <- smmry_final_enriched |>
      dplyr::mutate(Species = attr(seed_data_final, "species"), Description = attr(seed_data_final, "description"), Description = cp_normalize_marker(Description)) |>
      dplyr::mutate(Genes_raw = stringr::str_extract_all(Description, pattern), Genes_text = purrr::map_chr(Genes_raw, ~ paste(unique(cp_normalize_marker(.x)), collapse = ", ")), top_marker = cp_normalize_marker(top_marker)) |>
      dplyr::select(-Genes_raw) |>
      dplyr::left_join(gene_lookup, by = c("Genes_text" = "search_gene")) |>
      dplyr::left_join(marker_lookup, by = c("top_marker" = "search_marker")) |>
      dplyr::mutate(Gene_std = ifelse(is.na(Gene_std), Genes_text, Gene_std), Marker_std = ifelse(is.na(Marker_std), top_marker, Marker_std))
  }
  
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
  
  missing_marker_assignment <- df_joined |> dplyr::filter(is.na(Marker_std))
  if (nrow(missing_marker_assignment) > 0) log_message("Warning: clusters without Marker_std: ", paste(missing_marker_assignment$ID, collapse = ", "))
  
  df_joined <- df_joined |> dplyr::filter(!is.na(Marker_std))
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
    dplyr::summarise(num_clusters = dplyr::n(), num_sequences = sum(n_sequences, na.rm = TRUE), mean_seq_per_cluster = mean(n_sequences), median_seq_per_cluster = median(n_sequences), min_seq_in_cluster = min(n_sequences), max_seq_in_cluster = max(n_sequences), .groups = "drop") |>
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
#' Retrieves orthologous sequence clusters for specified outgroup lineages (e.g., *Portulaca*, *Anacampseros*, *Talinopsis*, *Grahamia*)
#' matching the locus target constraints defined for the focal ingroup. Outer reference sampling provides phylogenetically
#' informative root positions necessary for maximum-likelihood tree search and divergence time estimation.
#'
#' @param wd_path Character. Path to the `phylotaR` workspace directory.
#' @param target_genes_file Character. Path to the target locus list text file. If `NULL`, defaults to package `inst/extdata/target_genes.txt`.
#' @param genes_map_file Character. Path to the gene synonymy mapping CSV file. If `NULL`, defaults to package `inst/extdata/genes_map.csv`.
#' @param manual_exclusions_file Character. Path to outgroup exclusions CSV file. If `NULL`, defaults to package `inst/extdata/manual_exclusions_outgroup.csv`.
#' @param outgroups Character vector of NCBI Taxonomy IDs for outgroup lineages. Defaults to `c("107598", "107617", "107583", "3582")`.
#' @param force_download Logical. Force fresh database retrieval instead of using local cache? Defaults to `FALSE`.
#' @param out_dir Character. Output directory path to save outgroup cluster tables and FASTA sequence files.
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
assemble_outgroup_phylotar <- function(wd_path, target_genes_file = NULL, genes_map_file = NULL, manual_exclusions_file = NULL, outgroups = c("107598", "107617", "107583", "3582"), force_download = FALSE, out_dir = "1_phylotaR_out_Outgroup") {
  
  
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
  
  markers_df <- read.table(target_genes_file, header = TRUE, stringsAsFactors = FALSE)
  markers <- unique(cp_normalize_marker(as.character(markers_df$markers)))
  pattern <- cp_build_pattern(markers)
  genes_map_df <- read.csv(genes_map_file, stringsAsFactors = FALSE) |>
    dplyr::mutate(search = cp_normalize_marker(search), replace = trimws(replace)) |>
    dplyr::distinct(search, .keep_all = TRUE)
    
  gene_lookup <- genes_map_df |> dplyr::select(search_gene = search, Gene_std = replace)
  marker_lookup <- genes_map_df |> dplyr::select(search_marker = search, Marker_std = replace)
  
  phylota <- tryCatch({
    if (force_download) stop("Force download enabled")
    phylotaR::read_phylota(wd_path)
  }, error = function(e) {
    log_message("No valid phylota object found or force_download=TRUE. Running setup and phylotaR pipeline...")
    ncbi_dr_local <- NULL
    env_blast <- Sys.getenv("BLAST_PATH")
    if (env_blast != "") {
      ncbi_dr_local <- env_blast
    } else {
      blastn_path <- Sys.which("blastn")
      if (blastn_path != "") ncbi_dr_local <- dirname(blastn_path)
    }
    
    phylotaR::setup(
      wd = wd_path, txid = outgroups, ncbi_dr = ncbi_dr_local, v = TRUE, ncps = 1, mncvrg = 80,
      srch_trm = paste0(
        "NOT predicted[TI] NOT \"whole genome shotgun\"[TI] NOT unverified[TI] ",
        "NOT \"synthetic construct\"[Organism] NOT refseq[filter] NOT TSA[Keyword] ",
        "NOT \"sp.\"[TI] NOT \"sp.\"[Organism] NOT \"sp\"[Organism] NOT \"aff.\"[TI] ",
        "NOT \"aff\"[Organism] NOT \"cf.\"[TI] NOT \"cf\"[Organism] NOT \"var.\"[TI] ",
        "NOT \"var\"[TI] NOT \"var\"[Organism] NOT \"var.\"[Organism] NOT \"variety\"[TI] ",
        "NOT \"subsp.\"[TI] NOT \"subsp\"[TI] NOT \"subsp.\"[Organism] ",
        "NOT \"subsp\"[Organism] NOT \"subspecies\"[Organism] NOT \"x\"[Organism] NOT \" x \"[Organism]"
      )
    )
    
    log_message("Executing phylotaR::run() for outgroups...")
    tryCatch({
      phylotaR::run(wd = wd_path)
    }, error = function(erun) {
      log_message("Error in outgroup phylotaR::run(): ", erun$message)
    })
    
    log_message("Attempting to load outgroup phylota object after run()...")
    tryCatch({
      phylotaR::read_phylota(wd_path)
    }, error = function(e2) {
      log_message("read_phylota error: ", e2$message)
      stop("Could not load outgroup phylota object from ", wd_path)
    })
  })
  
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
  
  seed_data <- tryCatch({ ape::read.GenBank(smmry_sel_enriched$Seed) }, error = function(e) NULL)
  if(!is.null(seed_data)) {
    smmry_sel_enriched <- smmry_sel_enriched |>
      dplyr::mutate(Species = attr(seed_data, "species"), Description = attr(seed_data, "description"), Description = cp_normalize_marker(Description)) |>
      dplyr::mutate(Genes_raw = stringr::str_extract_all(Description, pattern), Genes_text = purrr::map_chr(Genes_raw, ~ paste(unique(cp_normalize_marker(.x)), collapse = ", ")), top_marker = cp_normalize_marker(top_marker)) |>
      dplyr::select(-Genes_raw) |>
      dplyr::left_join(gene_lookup, by = c("Genes_text" = "search_gene")) |> dplyr::left_join(marker_lookup, by = c("top_marker" = "search_marker")) |>
      dplyr::mutate(Gene_std = ifelse(is.na(Gene_std), Genes_text, Gene_std), Marker_std = ifelse(is.na(Marker_std), top_marker, Marker_std))
  }
  
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
  dup_records <- df_species_clusters_metadata |> dplyr::filter(sid %in% dup_sids) |> dplyr::mutate(cluster_id = as.integer(cluster_id), prnt = vapply(cluster_id, function(cid) cp_get_cluster_parent(selected, cid), FUN.VALUE = character(1)))
  keepers <- dup_records |> dplyr::group_by(sid) |> dplyr::arrange(cluster_id, .by_group = TRUE) |> dplyr::mutate(has_preferred = any(prnt == "866800", na.rm = TRUE), keep_rank = dplyr::case_when(has_preferred & prnt == "866800" ~ 1L, has_preferred & prnt != "866800" ~ 2L, !has_preferred ~ 1L, TRUE ~ 3L)) |> dplyr::arrange(keep_rank, cluster_id, .by_group = TRUE) |> dplyr::mutate(keep = row_number() == 1L) |> dplyr::ungroup() |> dplyr::filter(keep) |> dplyr::select(cluster_id, sid) |> dplyr::distinct()
  rows_to_delete_auto <- dup_records |> dplyr::anti_join(keepers, by = c("cluster_id", "sid")) |> dplyr::select(cluster_id, sid) |> dplyr::distinct()
  
  manual_rows <- tibble::tibble(cluster_id = integer(), sid = character(), reason = character())
  if (!is.null(manual_exclusions_file) && file.exists(manual_exclusions_file)) {
    manual_rows <- readr::read_csv(manual_exclusions_file, show_col_types = FALSE) |> dplyr::mutate(cluster_id = as.integer(cluster_id))
    readr::write_csv(manual_rows, path_table_manual_exclusions)
  }
  
  rows_to_delete <- dplyr::bind_rows(rows_to_delete_auto |> dplyr::mutate(reason = "auto"), manual_rows) |> dplyr::distinct(cluster_id, sid, .keep_all = TRUE)
  log_message("Automatic duplicate removals: ", nrow(rows_to_delete_auto))
  log_message("Manual duplicate removals: ", nrow(manual_rows))
  
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
  ntaxa_final <- get_ntaxa(phylota = phylota_final, cid = final_cids, rnk = "species")
  ntaxa_final_df <- tibble::tibble(cluster_id = final_cids, n_taxa = ntaxa_final)
  df_species_clusters_final <- cp_extract_cluster_species_sid(phylota_final)
  metadata_final <- metadata_raw
  df_species_clusters_metadata_final <- df_species_clusters_final |> dplyr::left_join(metadata_final, by = "sid") |> cp_annotate_marker_text(pattern = pattern, genes_map_df = genes_map_df)
  smmry_final <- phylotaR::summary(phylota_final) |> dplyr::mutate(ID = as.integer(ID))
  cluster_gene_summary_final <- cp_summarise_cluster_markers(df_species_clusters_metadata_final) |> dplyr::mutate(cluster_id = as.integer(cluster_id))
  smmry_final_enriched <- smmry_final |> dplyr::left_join(cluster_gene_summary_final, by = c("ID" = "cluster_id"))
  seed_data_final <- tryCatch({ ape::read.GenBank(smmry_final_enriched$Seed) }, error = function(e) NULL)
  
  if(!is.null(seed_data_final)) {
    smmry_final_enriched <- smmry_final_enriched |>
      dplyr::mutate(Species = attr(seed_data_final, "species"), Description = attr(seed_data_final, "description"), Description = cp_normalize_marker(Description)) |>
      dplyr::mutate(Genes_raw = stringr::str_extract_all(Description, pattern), Genes_text = purrr::map_chr(Genes_raw, ~ paste(unique(cp_normalize_marker(.x)), collapse = ", ")), top_marker = cp_normalize_marker(top_marker)) |>
      dplyr::select(-Genes_raw) |>
      dplyr::left_join(gene_lookup, by = c("Genes_text" = "search_gene")) |> dplyr::left_join(marker_lookup, by = c("top_marker" = "search_marker")) |>
      dplyr::mutate(Gene_std = ifelse(is.na(Gene_std), Genes_text, Gene_std), Marker_std = ifelse(is.na(Marker_std), top_marker, Marker_std))
  }
  
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
  df_joined <- df_joined |> dplyr::filter(!is.na(Marker_std))
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
    dplyr::summarise(num_clusters = dplyr::n(), num_sequences = sum(n_sequences, na.rm = TRUE), mean_seq_per_cluster = mean(n_sequences), median_seq_per_cluster = median(n_sequences), min_seq_in_cluster = min(n_sequences), max_seq_in_cluster = max(n_sequences), .groups = "drop") |>
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
#' Queries NCBI Entrez Utilities to retrieve sequence lengths, organism taxonomy, publication titles, and accession IDs
#' for a collection of GenBank sequence identifiers (SIDs). Metadata retrieval enriches raw sequence clusters with verifiable audit data.
#'
#' @param sids Character vector of GenBank Sequence Identifiers (SIDs) to query.
#' @param cache_file Character. File path to store and load cached metadata tables.
#' @param batch_size Integer. Number of sequence IDs requested per HTTP batch query. Defaults to `200`.
#' @param sleep_time Numeric. Pause duration in seconds between consecutive batch requests to respect NCBI rate limits. Defaults to `0.5`.
#' @param max_retries Integer. Maximum retry attempts permitted per batch before failing. Defaults to `5`.
#' @param force_download Logical. Force fresh Entrez queries instead of loading local cache? Defaults to `FALSE`.
#' @return A data frame containing fetched GenBank metadata fields for each requested sequence ID.
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
        gb <- tryCatch({ ape::read.GenBank(blocks[[i]]) }, error = function(e) NULL)
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
