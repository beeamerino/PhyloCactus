# -------------------------------------------------------------
# PhyloCactus: Tutorial 3 - Phylogenetics Pipeline: Data Visualization and IUCN Summaries
# -------------------------------------------------------------
# This script covers Stages 11 and 12:
# Data Visualization, IUCN Enrichment, and Manuscript Figures.
# It assumes you have already run Tutorial 1 and Tutorial 2 (up to Stage 10) and
# your working directory is set to the tutorial folder.
# -------------------------------------------------------------

# ============================================================
# LIBRARIES
# ============================================================
library(PhyloCactus)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(forcats)
library(purrr)
library(scales)
library(ggtree)
library(ape)
library(treeio)

tutorial_dir <- "~/Desktop/PhyloCactus_Tutorial"
setwd(tutorial_dir)

# ============================================================
# Stage 11: Integrating Biological Metadata and Evaluating Hypothesis Visualization
# ============================================================
cat("\n=======================================================\n")
cat("Stage 11: Integrating Biological Metadata and Evaluating Hypothesis Visualization\n")
cat("=======================================================\n")

# Create Output Directory
dir.create("9_Visualization/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("9_Visualization/tables", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# INPUT FILES (Referenced from tutorial root paths)
# ============================================================
species_file <- "6_Concatenated/final_tables/TABLE_final_species_alignment_summary.csv"
output_file <- "9_Visualization/tables/TABLE_species_all_data.csv"
prioritization_file <- "9_Visualization/tables/TABLE_conservation_prioritization_gaps.csv"
marker_stats_file <- "6_Concatenated/final_tables/TABLE_marker_statistics.csv"
marker_ranges_file <- "6_Concatenated/final_tables/TABLE_marker_ranges.tsv"

# Canonical inputs that do not change across runs
accepted_list_file <- system.file("extdata", "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx", package = "PhyloCactus")
constraints_csv_path <- system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus")

# Tree outputs from Stages 9 and 10
supp_tree_file <- if (file.exists("7_Phylogenetics/cactus_support.raxml.support")) {
  "7_Phylogenetics/cactus_support.raxml.support"
} else {
  "7_Phylogenetics/cactus.raxml.bestTree"
}
dated_tree_file <- "8_Dating/BestTree_treePL.tree"
chrono_file <- "8_Dating/dated_summary_hpd.tree"

# ============================================================
# LOAD SPECIES TABLE & ACCEPTED CHECKLIST (Korotkova et al. 2021)
# ============================================================
message("Reading species genetics table from supermatrix outputs...")
species_summary <- read_csv(species_file, show_col_types = FALSE) %>%
  mutate(species_clean = str_replace_all(species, "_", " "))

message("Reading accepted checklist (Korotkova et al. 2021, Caryophyllales.org)...")
if (grepl("\\.xlsx?$", accepted_list_file, ignore.case = TRUE)) {
  sheets <- readxl::excel_sheets(accepted_list_file)
  checklist_sheets <- sheets[!grepl("^facts", sheets, ignore.case = TRUE)]
  list_df <- lapply(checklist_sheets, function(sh) {
    df <- readxl::read_excel(accepted_list_file, sheet = sh)
    df$family_sheet <- sh
    df
  })
  accepted_list <- dplyr::bind_rows(list_df)
} else {
  accepted_list <- read.csv(accepted_list_file, stringsAsFactors = FALSE)
  if (!"family_sheet" %in% names(accepted_list)) accepted_list$family_sheet <- "Cactaceae"
}

# Normalization helper ensuring perfect match with molecular alignment keys
normalize_species <- function(x) {
  x |>
    stringr::str_replace_all("\u00d7", "") |>
    stringr::str_replace_all("\\bx[[:space:]]+", "") |>
    stringr::str_trim() |>
    stringr::str_replace_all("[ -]+", "_")
}

accepted_species_table <- accepted_list %>%
  filter(tolower(family_sheet) == "cactaceae", RANK == "Species") %>%
  mutate(pureName = str_squish(as.character(pureName)),
         fullName = str_squish(as.character(fullName))) %>%
  filter(!str_detect(pureName, regex("\\b(subsp|ssp|var|forma|f\\.|subg|sect|ser|cf\\.|aff\\.|sp\\.|spp\\.|nr\\.)\\b", ignore_case = TRUE))) %>%
  mutate(species = normalize_species(pureName)) %>%
  filter(nzchar(species)) %>%
  select(species, checklist_pureName = pureName, checklist_fullName = fullName,
         checklist_author = author, checklist_rank = RANK, checklist_taxon = taxon, checklist_uuid = uuid) %>%
  distinct(species, .keep_all = TRUE)

message("Reading phylogeny tree tips...")
ml_tree_tips <- ape::read.tree(supp_tree_file)$tip.label
chrono_tips <- if (file.exists(dated_tree_file)) {
  ape::read.tree(dated_tree_file)$tip.label
} else if (file.exists(chrono_file)) {
  as.phylo(treeio::read.beast(chrono_file))$tip.label
} else {
  ml_tree_tips
}

# Identify identical sequence groups collapsed during ML tree inference
raxml_log <- "7_Phylogenetics/cactus.raxml.log"
rx_lines <- if (file.exists(raxml_log)) readLines(raxml_log) else character(0)
warn_lines <- rx_lines[grepl("are exactly identical", rx_lines)]
collapsed_taxa <- if (length(warn_lines) > 0) {
  unique(stringr::str_match(warn_lines, "and ([A-Za-z0-9_]+) are")[, 2])
} else if (file.exists("6_Concatenated/logs_and_qc/SUPP_TABLE_identical_sequence_groups.csv")) {
  read_csv("6_Concatenated/logs_and_qc/SUPP_TABLE_identical_sequence_groups.csv", show_col_types = FALSE)$taxon
} else {
  setdiff(species_summary$species, ml_tree_tips)
}

message("Reconciling taxonomy, sequences, and phylogenetic tips...")
species_summary <- species_summary %>%
  full_join(accepted_species_table, by = "species") %>%
  mutate(
    species_clean = if_else(!is.na(species_clean), species_clean, str_replace_all(species, "_", " ")),
    in_checklist = !is.na(checklist_pureName),
    in_alignment = !is.na(species_class),
    in_ml_tree = species %in% ml_tree_tips,
    in_chrono_tree = species %in% chrono_tips,
    is_outgroup = if_else(in_alignment, species_class == "outgroup", FALSE),
    is_ingroup = if_else(in_alignment, species_class == "ingroup", !is_outgroup),
    # Backward-compatible flags
    has_sequences = in_alignment,
    in_phylogeny = in_ml_tree,
    accepted_in_checklist = in_checklist,
    accepted_or_outgroup = in_checklist | is_outgroup,
    # Categorical pipeline status
    status_pipeline = case_when(
      is_outgroup ~ "outgroup_reference",
      in_checklist & in_alignment & in_chrono_tree ~ "sampled_and_dated",
      in_alignment & !in_ml_tree & (species %in% collapsed_taxa) ~ "sequenced_collapsed_duplicate",
      in_alignment & !in_ml_tree ~ "sequenced_unplaced",
      in_checklist & !in_alignment ~ "unsampled_checklist_gap",
      TRUE ~ "other"
    ),
    record_status = case_when(
      is_outgroup ~ "accepted_outgroup",
      in_ml_tree & in_checklist ~ "sampled_ingroup",
      in_ml_tree & !in_checklist ~ "rejected_ingroup_in_tree",
      in_alignment & in_checklist ~ "sequenced_not_in_tree",
      in_alignment & !in_checklist ~ "rejected_ingroup_no_tree",
      !in_alignment & in_checklist ~ "unsampled_ingroup",
      TRUE ~ "other"
    )
  ) %>%
  distinct(species, .keep_all = TRUE) %>%
  arrange(species_clean)

# ============================================================
# RUN IUCN EXTRACTION OR LOAD INTEGRATED SPECIES TABLE
# ============================================================
if (file.exists(output_file)) {
  message("Reading integrated species table: ", output_file)
  species_summary_iucn <- read_csv(output_file, show_col_types = FALSE)
} else {
  message("Extracting IUCN data...")
  iucn_data <- purrr::map_dfr(species_summary$species_clean, PhyloCactus::get_iucn_data)
  species_summary_iucn <- species_summary %>%
    left_join(iucn_data, by = "species_clean")
  write_csv(species_summary_iucn, output_file)
}

# ============================================================
# CONSERVATION PRIORITIZATION: THREATENED UNSAMPLED GAPS
# ============================================================
message("Generating conservation prioritization table for unsampled gaps...")
TABLE_conservation_prioritization_gaps <- species_summary_iucn %>%
  filter(status_pipeline == "unsampled_checklist_gap", iucn_category %in% c("CR", "EN", "VU")) %>%
  mutate(genus = stringr::str_extract(species, "^[A-Za-z]+")) %>%
  select(species, species_clean, genus, iucn_category, iucn_category_name,
         iucn_criteria, iucn_population_trend, iucn_locations, iucn_endemic_locations,
         iucn_is_endemic, checklist_author, checklist_uuid) %>%
  arrange(factor(iucn_category, levels = c("CR", "EN", "VU")), genus, species)

write_csv(TABLE_conservation_prioritization_gaps, prioritization_file)
message("Conservation prioritization table written with ", nrow(TABLE_conservation_prioritization_gaps), " priority taxa.")

# ============================================================
# LOAD DATA FOR FIGURES
# ============================================================
marker_stats <- read_csv(marker_stats_file, show_col_types = FALSE)
marker_ranges <- read_tsv(marker_ranges_file, show_col_types = FALSE)
marker_cols <- intersect(marker_stats$marker, names(species_summary_iucn))

# ============================================================
# DATASET SUMMARY TABLE
# ============================================================
TABLE_dataset_summary <- tibble(
  metric = c(
    "Total ingroup species in checklist",
    "Total species with sequences (supermatrix)",
    "Ingroup species with sequences (supermatrix)",
    "Outgroup species with sequences (supermatrix)",
    "Ingroup species sampled and dated in tree",
    "Outgroup reference species in tree",
    "Sequenced ingroup species collapsed duplicates",
    "Unsampled ingroup species (checklist gaps)",
    "Total species in final phylogeny",
    "Ingroup species mapped to IUCN assessment",
    "Ingroup species without IUCN assessment",
    "Threatened species in phylogeny (VU/EN/CR) (Ingroup)",
    "Threatened unsequenced species (VU/EN/CR) (Priority gaps)",
    "Endemic species in phylogeny (Ingroup)"
  ),
  value = c(
    sum(species_summary_iucn$in_checklist & species_summary_iucn$is_ingroup, na.rm = TRUE),
    sum(species_summary_iucn$in_alignment, na.rm = TRUE),
    sum(species_summary_iucn$is_ingroup & species_summary_iucn$in_alignment, na.rm = TRUE),
    sum(species_summary_iucn$is_outgroup & species_summary_iucn$in_alignment, na.rm = TRUE),
    sum(species_summary_iucn$status_pipeline == "sampled_and_dated", na.rm = TRUE),
    sum(species_summary_iucn$status_pipeline == "outgroup_reference", na.rm = TRUE),
    sum(species_summary_iucn$status_pipeline == "sequenced_collapsed_duplicate", na.rm = TRUE),
    sum(species_summary_iucn$status_pipeline == "unsampled_checklist_gap", na.rm = TRUE),
    sum(species_summary_iucn$in_ml_tree, na.rm = TRUE),
    sum(species_summary_iucn$is_ingroup & species_summary_iucn$iucn_found == TRUE, na.rm = TRUE),
    sum(species_summary_iucn$is_ingroup & (is.na(species_summary_iucn$iucn_found) | species_summary_iucn$iucn_found == FALSE), na.rm = TRUE),
    sum(species_summary_iucn$status_pipeline == "sampled_and_dated" & species_summary_iucn$iucn_category %in% c("VU", "EN", "CR"), na.rm = TRUE),
    sum(species_summary_iucn$status_pipeline == "unsampled_checklist_gap" & species_summary_iucn$iucn_category %in% c("VU", "EN", "CR"), na.rm = TRUE),
    sum(species_summary_iucn$in_ml_tree & species_summary_iucn$is_ingroup & species_summary_iucn$iucn_is_endemic == TRUE, na.rm = TRUE)
  )
)
write_csv(TABLE_dataset_summary, "9_Visualization/tables/TABLE_dataset_summary.csv")

# ============================================================
# HARMONIZED PLOT THEME
# ============================================================
theme_phylocactus <- function(base_size = 12, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.background = element_rect(fill = "#F5F5F5", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "white", linewidth = 0.6),
      panel.grid.minor = element_line(color = "white", linewidth = 0.3),
      axis.ticks = element_line(color = "grey70", linewidth = 0.4),
      axis.title = element_text(face = "italic", color = "grey20"),
      axis.title.x = element_text(face = "italic", margin = margin(t = 8)),
      axis.title.y = element_text(face = "italic", margin = margin(r = 8)),
      axis.text = element_text(color = "grey30"),
      plot.title = element_text(face = "bold", color = "grey15", hjust = 0),
      plot.subtitle = element_text(color = "grey40", hjust = 0, margin = margin(b = 8)),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.title = element_text(face = "bold", size = rel(0.9)),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold", color = "grey20", size = rel(0.95))
    )
}

# ============================================================
# FIGURES: SEQUENTIAL SCIENTIFIC PRESENTATION (Figures 1 to 17)
# ============================================================

# -------------------------------------------------------------
# Part 1: Molecular Dataset and Alignment Properties (Figures 1 - 3)
# -------------------------------------------------------------
# FIGURE 1: Marker Coverage Across Species
message("Rendering Figure 1 (Marker Coverage)...")
p1 <- ggplot(marker_stats, aes(x = reorder(marker, n_taxa), y = n_taxa)) +
  geom_col(fill = "grey35", width = 0.7) +
  coord_flip() +
  theme_phylocactus(base_size = 12) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(title = "Marker Coverage Across Species", x = "Marker", y = "Number of species")
ggsave("9_Visualization/figures/Figure_1_marker_coverage.pdf", p1, width = 7, height = 5)

# FIGURE 2: Distribution of Marker Completeness
message("Rendering Figure 2 (Marker Completeness)...")
if ("retained_markers" %in% names(species_summary_iucn)) {
  p2 <- ggplot(species_summary_iucn %>% filter(in_ml_tree == TRUE), aes(x = retained_markers)) +
    geom_histogram(binwidth = 1, fill = "grey35", color = "white") +
    scale_x_continuous(breaks = seq(1, 11, 1)) +
    theme_phylocactus(base_size = 12) +
    labs(title = "Distribution of Marker Completeness", x = "Number of retained markers", y = "Number of species")
  ggsave("9_Visualization/figures/Figure_2_marker_completeness.pdf", p2, width = 7, height = 5)
}

# FIGURE 3: Marker Presence Heatmap
message("Rendering Figure 3 (Marker Presence Heatmap)...")
heatmap_data <- species_summary_iucn %>%
  filter(in_ml_tree == TRUE) %>%
  select(species, all_of(marker_cols)) %>%
  pivot_longer(cols = all_of(marker_cols), names_to = "marker", values_to = "accession") %>%
  mutate(present = accession != "-")
p3 <- ggplot(heatmap_data, aes(x = marker, y = fct_rev(species), fill = present)) +
  scale_fill_grey(start = 0.95, end = 0.25, name = "Present", labels = c("No", "Yes")) +
  geom_tile() +
  theme_phylocactus(base_size = 10) +
  labs(title = "Marker Presence Heatmap", x = "Marker", y = "Species") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave("9_Visualization/figures/Figure_3_marker_heatmap.pdf", p3, width = 8, height = 14)

# -------------------------------------------------------------
# Part 2: Supermatrix Concatenation Structure (Figures 4 - 5)
# -------------------------------------------------------------
# FIGURE 4: Supermatrix Structure
message("Rendering Figure 4 (Supermatrix Structure)...")
if (exists("marker_ranges") && nrow(marker_ranges) > 0) {
  marker_ranges <- marker_ranges %>% mutate(length = end - start + 1)
  p4 <- ggplot(marker_ranges) +
    geom_segment(aes(x = start, xend = end, y = marker, yend = marker), linewidth = 6, color = "grey35") +
    theme_phylocactus(base_size = 12) +
    theme(axis.text.y = element_text(face = "italic")) +
    labs(title = "Supermatrix Structure", x = "Alignment position (bp)", y = "Marker")
  ggsave("9_Visualization/figures/Figure_4_supermatrix_structure.pdf", p4, width = 10, height = 5)
}

# FIGURE 5: Marker Informativeness
message("Rendering Figure 5 (Marker Informativeness)...")
if ("parsimony_informative" %in% names(marker_stats)) {
  p5 <- ggplot(marker_stats, aes(x = alignment_length, y = parsimony_informative, label = marker)) +
    geom_point(size = 3.5, color = "grey35") +
    geom_text(nudge_y = 20, size = 3.2, fontface = "italic", color = "grey20") +
    theme_phylocactus(base_size = 12) +
    labs(title = "Marker Informativeness", x = "Alignment length (bp)", y = "Parsimony informative sites")
  ggsave("9_Visualization/figures/Figure_5_marker_informativeness.pdf", p5, width = 7, height = 5)
}

# -------------------------------------------------------------
# Part 3: Inferred Phylogeny and Divergence Timescale (Figures 6 - 9)
# -------------------------------------------------------------
# FIGURE 6: Maximum-Likelihood Tree with Bootstrap Support (FBP)
message("Rendering Figure 6 (maximum-likelihood tree with FBP bootstrap support)...")
if (file.exists(supp_tree_file)) {
  ml_tree <- ape::read.tree(supp_tree_file)
  ml_tree <- ape::ladderize(ml_tree)
  
  if (file.exists(constraints_csv_path) && nzchar(constraints_csv_path)) {
    constraints_tbl <- read_csv(constraints_csv_path, show_col_types = FALSE)
    if (is.null(ml_tree$node.label)) {
      ml_tree$node.label <- rep(NA_character_, ml_tree$Nnode)
    }

    MAIN_LEVEL4_VALUES <- c("Anacampseros", "Grahamia", "Talinopsis", "Leuenbergeria", "Pereskia", "Maihuenia", "Blossfeldia", "Opuntieae", "Cylindropuntieae", "Tephrocacteae", "Cacteae", "Core I", "Core II", "Copiapoa", "Calymmanthium", "Rhipsalis", "Portulaca", "Notocacteae", "BCT", "Rhipsalideae")
    ML_MAIN_LAYER_SPECS <- tibble::tribble(
      ~reg_name, ~fontsize, ~barsize, ~offset, ~offset_text, ~fontface, ~sort_desc, ~angle, ~align,
      "level_4_main", 2.5, 0.34, 0.010, 0.0038, "plain", TRUE, 0, TRUE,
      "level_3", 3.5, 0.50, 0.170, 0.0055, "bold", TRUE, 270, TRUE,
      "level_2", 4.0, 0.75, 0.220, 0.0100, "bold", TRUE, 270, TRUE,
      "level_1", 4.5, 0.85, 0.270, 0.0100, "bold", TRUE, 270, TRUE
    )
    
    registry <- PhyloCactus::build_annotation_registry(
      tree = ml_tree,
      constraints_tbl = constraints_tbl,
      main_level4_values = MAIN_LEVEL4_VALUES,
      supp_level4_min_tips = 3L
    )

    max_raw <- max(suppressWarnings(as.numeric(ml_tree$node.label)), na.rm = TRUE)
    is_proportion <- !is.na(max_raw) && max_raw <= 1.0

    # Identify topologically constrained nodes (cactus_constraints.tree)
    # to avoid presenting algorithmic constraints as empirical bootstrap support
    constraint_tree_file <- file.path(tutorial_dir, "7_Phylogenetics/cactus_constraints.tree")
    constrained_nodes_vec <- integer(0)
    if (file.exists(constraint_tree_file) && nzchar(constraint_tree_file)) {
      constrained_df <- tryCatch(
        PhyloCactus::classify_constrained_nodes(ml_tree, constraint_tree_file),
        error = function(e) NULL
      )
      if (!is.null(constrained_df)) {
        constrained_nodes_vec <- constrained_df$node[constrained_df$constrained]
      }
    }

    ml_tree_data <- tidytree::as_tibble(ml_tree) %>%
      mutate(
        raw_val = suppressWarnings(as.numeric(label)),
        pct_val = if (is_proportion) raw_val * 100 else raw_val,
        is_constrained = node %in% constrained_nodes_vec,
        support_label = ifelse(is_constrained | is.na(pct_val), NA_character_, sprintf("%.0f", pct_val)),
        support_class = case_when(
          is_constrained ~ "Constrained",
          pct_val >= 90  ~ ">= 90",
          pct_val >= 70  ~ "70-89",
          TRUE           ~ "< 70"
        ),
        support_class = factor(support_class, levels = c(">= 90", "70-89", "< 70", "Constrained"))
      )
    ml_tree_plot <- tidytree::as.treedata(ml_tree_data)

    p6 <- ggtree(ml_tree_plot, linewidth = 0.3) +
      theme_tree() +
      geom_nodepoint(aes(fill = support_class), shape = 21, size = 2, stroke = 0.2, na.rm = TRUE) +
      scale_fill_manual(
        values = c(">= 90" = "white", "70-89" = "grey", "< 70" = "black", "Constrained" = "steelblue"),
        name = "Support (FBP)",
        na.translate = FALSE
      ) +
      theme(legend.position = "bottom")

    p6 <- PhyloCactus::apply_clade_label_layers(p6, registry, ML_MAIN_LAYER_SPECS)
    xmax <- max(p6$data$x, na.rm = TRUE)
    p6 <- p6 + coord_cartesian(xlim = c(0, xmax + xmax * 0.35), clip = "off") + theme(plot.margin = margin(12, 12, 12, 12))

    ggsave("9_Visualization/figures/Figure_6_ML_Tree.pdf", p6, width = 8.5, height = 11)

    # FIGURE 7: Extended Constrained ML Phylogeny (A0 Poster)
    ML_SUPP_LAYER_SPECS <- tibble::tribble(
      ~reg_name, ~fontsize, ~barsize, ~offset, ~offset_text, ~fontface, ~sort_desc, ~angle, ~align,
      "level_4_supp", 6, 0.60, 0.030, 0.0038, "plain", TRUE, 0, TRUE,
      "level_3", 8, 0.65, 0.150, 0.0038, "bold", TRUE, 270, TRUE,
      "level_2", 9, 0.70, 0.200, 0.0100, "bold", TRUE, 270, TRUE,
      "level_1", 10, 0.80, 0.250, 0.0100, "bold", TRUE, 270, TRUE
    )
    p7_supp <- ggtree(ml_tree_plot, linewidth = 0.2) +
      theme_tree() +
      geom_tiplab(size = 0.70, align = FALSE, linetype = "solid", linesize = 0.15, colour = "grey20") +
      geom_nodelab(aes(label = support_label), size = 0.70, hjust = 0.7, colour = "black", na.rm = TRUE) +
      theme(legend.position = "none")
    p7_supp <- PhyloCactus::apply_clade_label_layers(p7_supp, registry, ML_SUPP_LAYER_SPECS)
    xmax_supp <- max(p7_supp$data$x, na.rm = TRUE)
    p7_supp <- p7_supp + coord_cartesian(xlim = c(0, xmax_supp + xmax_supp * 0.30), clip = "off") + theme(plot.margin = margin(12, 12, 12, 12))
    ggsave("9_Visualization/figures/Figure_7_ML_Tree_Extended.pdf", p7_supp, width = 36, height = 48)

    annotated_treedata <- PhyloCactus::augment_treedata_with_registry(ml_tree_plot, registry)
    treeio::write.beast(annotated_treedata, file = "9_Visualization/figures/annotated_ml_tree.tree")
  }
}

# FIGURE 8: Time-Calibrated Chronogram with HPD intervals
message("Rendering Figure 8 (Chronogram)...")
if (file.exists(chrono_file)) {
  chrono_beast <- treeio::read.beast(chrono_file)
  
  if (file.exists(constraints_csv_path) && nzchar(constraints_csv_path)) {
    constraints_tbl <- read_csv(constraints_csv_path, show_col_types = FALSE)
    tree_phy <- as.phylo(chrono_beast)
    
    if (is.null(tree_phy$node.label)) {
      tree_phy$node.label <- rep(NA_character_, tree_phy$Nnode)
    }

    MAIN_LEVEL4_VALUES <- c("Anacampseros", "Grahamia", "Talinopsis", "Leuenbergeria", "Pereskia", "Maihuenia", "Blossfeldia", "Opuntieae", "Cylindropuntieae", "Tephrocacteae", "Cacteae", "Core I", "Core II", "Copiapoa", "Calymmanthium", "Rhipsalis", "Portulaca", "Notocacteae", "BCT", "Rhipsalideae")
    CHRONO_MAIN_LAYER_SPECS <- tibble::tribble(
      ~reg_name, ~fontsize, ~barsize, ~offset, ~offset_text, ~fontface, ~sort_desc, ~angle, ~align,
      "level_4_main", 2.5, 0.34, 0.010, 0.0038, "plain", TRUE, 0, TRUE,
      "level_3", 3.5, 0.50, 0.170, 0.0055, "bold", TRUE, 270, TRUE,
      "level_2", 4.0, 0.75, 0.220, 0.0100, "bold", TRUE, 270, TRUE,
      "level_1", 4.5, 0.85, 0.270, 0.0100, "bold", TRUE, 270, TRUE
    )
    
    registry_chrono <- PhyloCactus::build_annotation_registry(
      tree = tree_phy,
      constraints_tbl = constraints_tbl,
      main_level4_values = MAIN_LEVEL4_VALUES,
      supp_level4_min_tips = 3L
    )

    p8 <- ggtree(chrono_beast, linewidth = 0.3) +
      theme_tree2() + 
      geom_range(range = 'height_0.95_HPD', color = 'gray60', alpha = 0.4, linewidth = 1.5) +
      theme_tree()
    
    p8 <- PhyloCactus::apply_clade_label_layers(p8, registry_chrono, CHRONO_MAIN_LAYER_SPECS)
    p8 <- PhyloCactus::add_chronogram_axis(p8, tree_phy, by = 5, digits = 0L, segment_size = 3, title_margin_top = 30)
    
    xmax <- max(p8$data$x, na.rm = TRUE)
    p8 <- p8 + coord_cartesian(xlim = c(0, xmax + xmax * 0.35), clip = "off") + theme(plot.margin = margin(12, 12, 12, 12))
      
    ggsave("9_Visualization/figures/Figure_8_Chronogram.pdf", p8, width = 8.5, height = 11)

    # FIGURE 9: Extended Time-Calibrated Chronogram (A0 Poster)
    CHRONO_SUPP_LAYER_SPECS <- tibble::tribble(
      ~reg_name, ~fontsize, ~barsize, ~offset, ~offset_text, ~fontface, ~sort_desc, ~angle, ~align,
      "level_4_supp", 6, 0.60, 0.030, 0.0038, "plain", TRUE, 0, TRUE,
      "level_3", 8, 0.65, 0.150, 0.0038, "bold", TRUE, 270, TRUE,
      "level_2", 9, 0.70, 0.200, 0.0100, "bold", TRUE, 270, TRUE,
      "level_1", 10, 0.80, 0.250, 0.0100, "bold", TRUE, 270, TRUE
    )
    p9_supp <- ggtree(chrono_beast, linewidth = 0.2) +
      theme_tree2() +
      geom_tiplab(size = 0.85, align = FALSE, linetype = "solid", linesize = 0.15, colour = "grey20") +
      geom_range(range = 'height_0.95_HPD', color = 'gray60', alpha = 0.4, linewidth = 0.8) +
      theme_tree()
    p9_supp <- PhyloCactus::apply_clade_label_layers(p9_supp, registry_chrono, CHRONO_SUPP_LAYER_SPECS)
    p9_supp <- PhyloCactus::add_chronogram_axis(p9_supp, tree_phy, by = 5, digits = 0L, segment_size = 10, title_margin_top = 28, bar_size = 30)
    xmax_chrono_supp <- max(p9_supp$data$x, na.rm = TRUE)
    p9_supp <- p9_supp + coord_cartesian(xlim = c(0, xmax_chrono_supp + xmax_chrono_supp * 0.30), clip = "off") + theme(plot.margin = margin(12, 12, 32, 12))
    ggsave("9_Visualization/figures/Figure_9_Chronogram_Extended.pdf", p9_supp, width = 36, height = 48)

    annotated_chrono <- PhyloCactus::augment_treedata_with_registry(chrono_beast, registry_chrono)
    treeio::write.beast(annotated_chrono, file = "9_Visualization/figures/annotated_chronogram.tree")
  }
}

# -------------------------------------------------------------
# Part 4: Taxonomic Reconciliation and Species Cadastre (Figure 10)
# -------------------------------------------------------------
# FIGURE 10: Species Composition by Pipeline Status
message("Rendering Figure 10 (Species Composition across Pipeline)...")
p10 <- ggplot(species_summary_iucn %>% count(status_pipeline), aes(x = reorder(status_pipeline, n), y = n)) +
  geom_col(fill = "grey35", width = 0.7) +
  coord_flip() +
  geom_text(aes(label = scales::comma(n)), hjust = -0.15, size = 3.8, color = "grey20") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_phylocactus(base_size = 12) +
  labs(
    title = "Species Composition",
    subtitle = "Methodological status across checklist, supermatrix, and tree",
    x = NULL,
    y = "Number of species"
  )
ggsave("9_Visualization/figures/Figure_10_species_composition.pdf", p10, width = 7, height = 4.5)

# -------------------------------------------------------------
# Part 5: IUCN Conservation and Macroevolutionary Syntheses (Figures 11 - 17)
# -------------------------------------------------------------
if (any(species_summary_iucn$iucn_found, na.rm = TRUE)) {
  TABLE_iucn_categories <- species_summary_iucn %>%
    filter(iucn_found == TRUE, in_ml_tree == TRUE, is_ingroup == TRUE) %>%
    count(iucn_category, iucn_category_name, sort = TRUE)
  write_csv(TABLE_iucn_categories, "9_Visualization/tables/TABLE_iucn_categories.csv")
    
  # FIGURE 11: IUCN Categories
  message("Rendering Figure 11 (IUCN Categories)...")
  p11 <- ggplot(TABLE_iucn_categories, aes(x = reorder(iucn_category, n), y = n, fill = iucn_category)) +
    PhyloCactus::scale_fill_iucn(name = "IUCN Category") +
    geom_col(width = 0.7) +
    theme_phylocactus(base_size = 12) +
    labs(title = "IUCN Categories", x = "Category", y = "Number of species")
  ggsave("9_Visualization/figures/Figure_11_iucn_categories.pdf", p11, width = 6, height = 4)
    
  # FIGURE 12: Completeness vs IUCN
  message("Rendering Figure 12 (Completeness vs IUCN)...")
  p12 <- ggplot(
    species_summary_iucn %>%
      filter(status_pipeline == "sampled_and_dated",
             iucn_found == TRUE,
             !is.na(iucn_category)),
    aes(x = iucn_category, y = pct_markers, fill = iucn_category)
  ) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.2, outlier.size = 0.6, alpha = 0.8) +
    PhyloCactus::scale_fill_iucn(name = "IUCN Category") +
    theme_phylocactus(base_size = 12) +
    labs(
      title = "Molecular Completeness vs IUCN Category",
      subtitle = "Marker occupancy across IUCN threat categories for sampled ingroup taxa (N = 985)",
      x = "IUCN Category",
      y = "Marker completeness (%)"
    )
  suppressWarnings(ggsave("9_Visualization/figures/Figure_12_completeness_vs_iucn.pdf", p12, width = 7, height = 5))

  # FIGURE 13: Conservation Prioritization of Unsampled Gaps by Genus
  message("Rendering Figure 13 (Conservation Prioritization Gaps)...")
  plot_gaps_genus <- TABLE_conservation_prioritization_gaps %>%
    count(genus, iucn_category) %>%
    mutate(iucn_category = factor(iucn_category, levels = c("CR", "EN", "VU")))

  top_gap_genera <- TABLE_conservation_prioritization_gaps %>%
    count(genus, sort = TRUE) %>%
    slice_max(n, n = 12) %>%
    pull(genus)

  p13 <- ggplot(plot_gaps_genus %>% filter(genus %in% top_gap_genera),
                aes(x = factor(genus, levels = rev(top_gap_genera)), y = n, fill = iucn_category)) +
    geom_col(width = 0.7) +
    coord_flip() +
    PhyloCactus::scale_fill_iucn(name = "Threat Category") +
    theme_phylocactus(base_size = 12) +
    theme(axis.text.y = element_text(face = "italic")) +
    labs(
      title = "Threatened Unsampled Species by Genus",
      subtitle = "Priority Cactaceae taxa (CR, EN, VU) lacking molecular sequences",
      x = "Genus",
      y = "Number of threatened unsequenced species"
    )
  ggsave("9_Visualization/figures/Figure_13_conservation_prioritization_gaps.pdf", p13, width = 7, height = 5)
    
  # FIGURE 14: Endemism by Country
  message("Rendering Figure 14 (Endemism by Country)...")
  TABLE_endemic_species <- species_summary_iucn %>%
    filter(iucn_is_endemic == TRUE, in_ml_tree == TRUE, is_ingroup == TRUE) %>%
    select(species, iucn_category, iucn_endemic_locations)
  write_csv(TABLE_endemic_species, "9_Visualization/tables/TABLE_endemic_species.csv")

  plot_data_endimism <- species_summary_iucn %>%
    filter(iucn_is_endemic == TRUE, in_ml_tree == TRUE, is_ingroup == TRUE) %>%
    separate_rows(iucn_endemic_locations, sep = "; ") %>%
    count(iucn_endemic_locations, sort = TRUE) %>%
    slice_max(n, n = 20)

  p14 <- ggplot(plot_data_endimism, aes(x = reorder(iucn_endemic_locations, n), y = n)) +
    geom_col(fill = "grey35", width = 0.7) +
    coord_flip() +
    theme_phylocactus(base_size = 12) +
    labs(title = "Endemic Species by Country", x = "Country", y = "Number of endemic species")
  ggsave("9_Visualization/figures/Figure_14_endemism_by_country.pdf", p14, width = 8, height = 6)
  
  # FIGURE 15: IUCN Categories by Country (Endemic)
  message("Rendering Figure 15 (IUCN Categories by Country)...")
  plot_iucn_country <- species_summary_iucn %>%
    filter(iucn_is_endemic == TRUE, in_ml_tree == TRUE, is_ingroup == TRUE, !is.na(iucn_endemic_locations)) %>%
    separate_rows(iucn_endemic_locations, sep = "; ") %>%
    count(iucn_endemic_locations, iucn_category)
  top_countries <- plot_data_endimism$iucn_endemic_locations
  plot_iucn_country <- plot_iucn_country %>% filter(iucn_endemic_locations %in% top_countries)
  if (nrow(plot_iucn_country) > 0) {
    p15 <- ggplot(plot_iucn_country, aes(x = factor(iucn_endemic_locations, levels = rev(top_countries)), y = n, fill = iucn_category)) +
      geom_col(width = 0.7) +
      coord_flip() +
      PhyloCactus::scale_fill_iucn(name = "IUCN Category") +
      theme_phylocactus(base_size = 12) +
      labs(title = "Endemic Species by Country and IUCN Category", x = "Country", y = "Number of endemic species")
    ggsave("9_Visualization/figures/Figure_15_iucn_categories_by_country.pdf", p15, width = 8, height = 6)
  }

  # FIGURE 16: Primary Ecological Threats
  message("Rendering Figure 16 (Primary Ecological Threats)...")
  plot_threats <- species_summary_iucn %>%
    filter(!is.na(iucn_threat_names), in_ml_tree == TRUE, is_ingroup == TRUE) %>%
    separate_rows(iucn_threat_names, sep = "; ") %>%
    count(iucn_threat_names, sort = TRUE) %>%
    slice_max(n, n = 15)
  if (nrow(plot_threats) > 0) {
    p16 <- ggplot(plot_threats, aes(x = reorder(iucn_threat_names, n), y = n)) +
      geom_col(fill = "grey35", width = 0.7) +
      coord_flip() +
      theme_phylocactus(base_size = 12) +
      labs(title = "Most Common Threats", x = "Threat", y = "Number of species")
    ggsave("9_Visualization/figures/Figure_16_common_threats.pdf", p16, width = 9, height = 6)
  }
    
  # FIGURE 17: Primary Ecological Habitats
  message("Rendering Figure 17 (Primary Ecological Habitats)...")
  plot_habitats <- species_summary_iucn %>%
    filter(!is.na(iucn_habitat_names), in_ml_tree == TRUE, is_ingroup == TRUE) %>%
    separate_rows(iucn_habitat_names, sep = "; ") %>%
    count(iucn_habitat_names, sort = TRUE) %>%
    slice_max(n, n = 15)
  if (nrow(plot_habitats) > 0) {
    p17 <- ggplot(plot_habitats, aes(x = reorder(iucn_habitat_names, n), y = n)) +
      geom_col(fill = "grey35", width = 0.7) +
      coord_flip() +
      theme_phylocactus(base_size = 12) +
      labs(title = "Most Common Habitats", x = "Habitat", y = "Number of species")
    ggsave("9_Visualization/figures/Figure_17_common_habitats.pdf", p17, width = 9, height = 6)
  }
}

message("Visualization pipeline completed successfully (Figures 1 to 17 generated).")
