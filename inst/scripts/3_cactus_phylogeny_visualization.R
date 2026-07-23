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
# INPUT FILES (Ensure these exist from your Tutorial 1 run)
# ============================================================
species_file <- "6_Concatenated/final_tables/TABLE_final_species_alignment_summary.csv"
output_file <- "9_Visualization/tables/TABLE_species_all_data.csv"
marker_stats_file <- "6_Concatenated/final_tables/TABLE_marker_statistics.csv"
marker_ranges_file <- "6_Concatenated/final_tables/TABLE_marker_ranges.tsv"
accepted_list_file <- system.file("extdata", "CactaceaeFullList_accepted.csv", package = "PhyloCactus")
supp_tree_file <- "7_Phylogenetics/cactus_support.raxml.support"
dated_tree_file <- "8_Dating/BestTree_treePL.tree"
chrono_file <- "8_Dating/dated_summary_hpd.tree"

# ============================================================
# LOAD SPECIES TABLE & ACCEPTED CHECKLIST
# ============================================================
message("Reading species genetics table...")
species_summary <- read_csv(species_file, show_col_types = FALSE) %>%
  mutate(species_clean = str_replace_all(species, "_", " "))

message("Reading accepted checklist...")
accepted_list <- read.csv(accepted_list_file, stringsAsFactors = FALSE)
accepted_species_table <- accepted_list %>%
  filter(RANK == "Species") %>%
  mutate(pureName = str_squish(as.character(pureName)),
         fullName = str_squish(as.character(fullName))) %>%
  filter(!str_detect(pureName, regex("\\b(subsp|ssp|var|forma|f\\.|subg|sect|ser|cf\\.|aff\\.|sp\\.|spp\\.|nr\\.)\\b|×|\\bx\\b", ignore_case = TRUE))) %>%
  filter(str_detect(pureName, "^[A-Z][a-zA-Z-]+\\s[a-z-]+$")) %>%
  mutate(species = str_replace_all(pureName, " ", "_")) %>%
  select(species, checklist_pureName = pureName, checklist_fullName = fullName,
         checklist_author = author, checklist_rank = RANK, checklist_taxon = taxon, checklist_uuid = uuid) %>%
  distinct()

message("Reading phylogeny tree to get tips...")
ml_tree_tips <- ape::read.tree(supp_tree_file)$tip.label

message("Merging checklist metadata...")
species_summary <- species_summary %>%
  full_join(accepted_species_table, by = "species") %>%
  mutate(
    species_clean = if_else(is.na(species_clean), str_replace_all(species, "_", " "), species_clean),
    has_sequences = !is.na(species_class),
    in_phylogeny = species %in% ml_tree_tips,
    is_outgroup = if_else(has_sequences, species_class == "outgroup", FALSE),
    is_ingroup = if_else(has_sequences, species_class == "ingroup", !is_outgroup),
    accepted_in_checklist = !is.na(checklist_pureName),
    accepted_or_outgroup = accepted_in_checklist | is_outgroup,
    record_status = case_when(
      is_outgroup ~ "accepted_outgroup",
      in_phylogeny & accepted_in_checklist ~ "sampled_ingroup",
      in_phylogeny & !accepted_in_checklist ~ "rejected_ingroup_in_tree",
      has_sequences & accepted_in_checklist ~ "sequenced_not_in_tree",
      has_sequences & !accepted_in_checklist ~ "rejected_ingroup_no_tree",
      !has_sequences & accepted_in_checklist ~ "unsampled_ingroup",
      TRUE ~ "other"
    )
  )

# ============================================================
# RUN IUCN EXTRACTION
# ============================================================
# Note: By default, it processes the first 20 species to avoid long waits in the tutorial.
# Change [1:20] to use all species: species_summary$species_clean
message("Extracting IUCN data...")
iucn_data <- map_dfr(species_summary$species_clean[1:20], PhyloCactus::get_iucn_data)

message("Building final table...")
species_summary_iucn <- species_summary %>%
  left_join(iucn_data, by = "species_clean")

message("Exporting final table...")
write_csv(species_summary_iucn, output_file)

message("IUCN enrichment completed successfully.")

# ============================================================
# LOAD DATA FOR FIGURES
# ============================================================
# Reload the file just to ensure it works exactly as a standalone module
species_summary_iucn <- read_csv(output_file, show_col_types = FALSE)
marker_stats <- read_csv(marker_stats_file, show_col_types = FALSE)
marker_ranges <- read_tsv(marker_ranges_file, show_col_types = FALSE)
marker_cols <- marker_stats$marker

# ============================================================
# DATASET SUMMARY TABLE
# ============================================================
TABLE_dataset_summary <- tibble(
  metric = c("Total ingroup species in checklist", "Total species with sequences", 
             "Ingroup species with sequences", "Outgroup species with sequences", 
             "Ingroup species in final phylogeny", "Outgroup species in final phylogeny",
             "Ingroup species mapped to IUCN assessment", "Ingroup species without IUCN assessment",
             "Ingroup species in phylogeny with IUCN assessment",
             "Threatened species in phylogeny (VU/EN/CR) (Ingroup)", "Endemic species in phylogeny (Ingroup)"),
  value = c(sum(species_summary_iucn$accepted_in_checklist, na.rm = TRUE),
            sum(species_summary_iucn$has_sequences, na.rm = TRUE),
            sum(species_summary_iucn$is_ingroup & species_summary_iucn$has_sequences, na.rm = TRUE),
            sum(species_summary_iucn$is_outgroup & species_summary_iucn$has_sequences, na.rm = TRUE),
            sum(species_summary_iucn$in_phylogeny & species_summary_iucn$is_ingroup, na.rm = TRUE),
            sum(species_summary_iucn$in_phylogeny & species_summary_iucn$is_outgroup, na.rm = TRUE),
            sum(species_summary_iucn$is_ingroup & species_summary_iucn$iucn_found == TRUE, na.rm = TRUE),
            sum(species_summary_iucn$is_ingroup & (is.na(species_summary_iucn$iucn_found) | species_summary_iucn$iucn_found == FALSE), na.rm = TRUE),
            sum(species_summary_iucn$in_phylogeny & species_summary_iucn$is_ingroup & species_summary_iucn$iucn_found == TRUE, na.rm = TRUE),
            sum(species_summary_iucn$in_phylogeny & species_summary_iucn$is_ingroup & species_summary_iucn$iucn_category %in% c("VU", "EN", "CR"), na.rm = TRUE),
            sum(species_summary_iucn$in_phylogeny & species_summary_iucn$is_ingroup & species_summary_iucn$iucn_is_endemic == TRUE, na.rm = TRUE))
)
write_csv(TABLE_dataset_summary, "9_Visualization/tables/TABLE_dataset_summary.csv")

# ============================================================
# FIGURES
# ============================================================

# FIGURE 1: Species Composition
p1 <- ggplot(species_summary_iucn %>% filter(accepted_or_outgroup) %>% count(record_status), aes(x = record_status, y = n)) +
  geom_col(width = 0.7) + coord_flip() + theme_minimal(base_size = 12) +
  labs(title = "Species Composition (Phylogeny vs Checklist)", x = NULL, y = "Number of species") + theme(legend.position = "none")

print(p1)
ggsave("9_Visualization/figures/Figure_1_species_composition.pdf", p1, width = 5, height = 4)

# FIGURE 2: Marker Coverage
p2 <- ggplot(marker_stats, aes(x = reorder(marker, n_taxa), y = n_taxa)) +
  geom_col() + coord_flip() + theme_minimal(base_size = 12) +
  labs(title = "Marker Coverage Across Species", x = "Marker", y = "Number of species")
print(p2)
ggsave("9_Visualization/figures/Figure_2_marker_coverage.pdf", p2, width = 7, height = 5)

# FIGURE 3: Marker Completeness
# Use 'retained_markers' or proxy
if("retained_markers" %in% names(species_summary_iucn)){
  p3 <- ggplot(species_summary_iucn %>% filter(in_phylogeny == TRUE), aes(x = retained_markers)) +
      geom_histogram(binwidth = 1, color = "white") +
      scale_x_continuous(breaks = seq(1, 11, 1)) + theme_minimal(base_size = 12) +
      labs(title = "Distribution of Marker Completeness", x = "Number of retained markers", y = "Number of species")
  print(p3)
  ggsave("9_Visualization/figures/Figure_3_marker_completeness.pdf", p3, width = 7, height = 5)
}


# FIGURE 4: Marker Heatmap
heatmap_data <- species_summary_iucn %>% filter(in_phylogeny == TRUE) %>% select(species, all_of(marker_cols)) %>%
  pivot_longer(cols = all_of(marker_cols), names_to = "marker", values_to = "accession") %>% mutate(present = accession != "-")
p4 <- ggplot(heatmap_data, aes(x = marker, y = fct_rev(species), fill = present)) +
  scale_fill_grey(start = 0.95, end = 0.25) + geom_tile() + theme_minimal(base_size = 10) +
  labs(title = "Marker Presence Heatmap", x = "Marker", y = "Species") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
print(p4)
ggsave("9_Visualization/figures/Figure_4_marker_heatmap.pdf", p4, width = 8, height = 14)

# FIGURE 5: Supermatrix Structure
  if (exists("marker_ranges")) {
    marker_ranges <- marker_ranges %>% mutate(length = end - start + 1)
    p5 <- ggplot(marker_ranges) + geom_segment(aes(x = start, xend = end, y = marker, yend = marker), linewidth = 6) +
      theme_minimal(base_size = 12) + labs(title = "Supermatrix Structure", x = "Alignment position", y = "Marker")
    print(p5)
    ggsave("9_Visualization/figures/Figure_5_supermatrix_structure.pdf", p5, width = 10, height = 5)
  }
  
  # FIGURE 6: Marker Informativeness
  if("parsimony_informative" %in% names(marker_stats)) {
    p6 <- ggplot(marker_stats, aes(x = alignment_length, y = parsimony_informative, label = marker)) +
        geom_point(size = 3) + geom_text(nudge_y = 20, size = 3) + theme_minimal(base_size = 12) +
        labs(title = "Marker Informativeness", x = "Alignment length", y = "Parsimony informative sites")
    print(p6)
    ggsave("9_Visualization/figures/Figure_6_marker_informativeness.pdf", p6, width = 7, height = 5)
  }

# IUCN Categories list
if(any(species_summary_iucn$iucn_found, na.rm=TRUE)) {
  TABLE_iucn_categories <- species_summary_iucn %>% filter(iucn_found == TRUE, in_phylogeny == TRUE, is_ingroup == TRUE) %>% count(iucn_category, iucn_category_name, sort = TRUE)
  write_csv(TABLE_iucn_categories, "9_Visualization/tables/TABLE_iucn_categories.csv")
    
  # FIGURE 7: IUCN Categories
  p7 <- ggplot(TABLE_iucn_categories, aes(x = reorder(iucn_category, n), y = n, fill = iucn_category)) +
      scale_fill_iucn(name = "IUCN Category") + geom_col() + theme_minimal(base_size = 12) +
      labs(title = "IUCN Categories (Ingroup Phylogeny)", x = "Category", y = "Number of species")
  print(p7)
  ggsave("9_Visualization/figures/Figure_7_iucn_categories.pdf", p7, width = 6, height = 4)
    
  # FIGURE 8: Completeness VS IUCN
  p8 <- ggplot(species_summary_iucn %>% filter(iucn_found == TRUE, in_phylogeny == TRUE, is_ingroup == TRUE), aes(x = iucn_category, y = pct_markers, fill = iucn_category)) +
      geom_violin(trim = FALSE) + scale_fill_iucn(name = "IUCN Category") + geom_boxplot(width = 0.15, outlier.size = 0.5) +
      theme_minimal(base_size = 12) + labs(title = "Molecular Completeness vs IUCN Category", x = "IUCN category", y = "Marker completeness (%)")
  print(p8)
  ggsave("9_Visualization/figures/Figure_8_completeness_vs_iucn.pdf", p8, width = 7, height = 5)
    
  # FIGURE 9: Endemism
  TABLE_endemic_species <- species_summary_iucn %>% filter(iucn_is_endemic == TRUE, in_phylogeny == TRUE, is_ingroup == TRUE) %>% select(species, iucn_category, iucn_endemic_locations)
  write_csv(TABLE_endemic_species, "9_Visualization/tables/TABLE_endemic_species.csv")
  plot_data_endimism <- species_summary_iucn %>% filter(iucn_is_endemic == TRUE, in_phylogeny == TRUE, is_ingroup == TRUE) %>% separate_rows(iucn_endemic_locations, sep = "; ") %>% count(iucn_endemic_locations, sort = TRUE) %>% slice_max(n, n = 20)
  p9 <- ggplot(plot_data_endimism, aes(x = reorder(iucn_endemic_locations, n), y = n)) + geom_col() + coord_flip() + theme_minimal(base_size = 12) +
      labs(title = "Endemic Species by Country (Ingroup Phylogeny)", x = "Country", y = "Number of endemic species")
  print(p9)
  ggsave("9_Visualization/figures/Figure_9_endemism_by_country.pdf", p9, width = 8, height = 6)
  
  # FIGURE 10: IUCN Categories by Country (Endemic)
  plot_iucn_country <- species_summary_iucn %>% filter(iucn_is_endemic == TRUE, in_phylogeny == TRUE, is_ingroup == TRUE, !is.na(iucn_endemic_locations)) %>% separate_rows(iucn_endemic_locations, sep = "; ") %>% count(iucn_endemic_locations, iucn_category)
  top_countries <- plot_data_endimism$iucn_endemic_locations
  plot_iucn_country <- plot_iucn_country %>% filter(iucn_endemic_locations %in% top_countries)
  if(nrow(plot_iucn_country) > 0) {
    p10 <- ggplot(plot_iucn_country, aes(x = factor(iucn_endemic_locations, levels = rev(top_countries)), y = n, fill = iucn_category)) +
      geom_col() + coord_flip() + scale_fill_iucn(name = "IUCN Category") + theme_minimal(base_size = 12) +
      labs(title = "Endemic Species by Country and IUCN Category", x = "Country", y = "Number of endemic species")
    print(p10)
    ggsave("9_Visualization/figures/Figure_10_iucn_categories_by_country.pdf", p10, width = 8, height = 6)
  }

  # FIGURE 11 & 12: Threats and Habitats
  plot_threats <- species_summary_iucn %>% filter(!is.na(iucn_threat_names), in_phylogeny == TRUE, is_ingroup == TRUE) %>% separate_rows(iucn_threat_names, sep = "; ") %>% count(iucn_threat_names, sort = TRUE) %>% slice_max(n, n = 15)
  if(nrow(plot_threats) > 0) {
    p11 <- ggplot(plot_threats, aes(x = reorder(iucn_threat_names, n), y = n)) + geom_col() + coord_flip() + theme_minimal(base_size = 12) +
        labs(title = "Most Common Threats (Ingroup Phylogeny)", x = "Threat", y = "Number of species")
    print(p11)
    ggsave("9_Visualization/figures/Figure_11_common_threats.pdf", p11, width = 9, height = 6)
  }
    
  plot_habitats <- species_summary_iucn %>% filter(!is.na(iucn_habitat_names), in_phylogeny == TRUE, is_ingroup == TRUE) %>% separate_rows(iucn_habitat_names, sep = "; ") %>% count(iucn_habitat_names, sort = TRUE) %>% slice_max(n, n = 15)
  if(nrow(plot_habitats) > 0) {
    p12 <- ggplot(plot_habitats, aes(x = reorder(iucn_habitat_names, n), y = n)) + geom_col() + coord_flip() + theme_minimal(base_size = 12) +
        labs(title = "Most Common Habitats (Ingroup Phylogeny)", x = "Habitat", y = "Number of species")
    print(p12)
    ggsave("9_Visualization/figures/Figure_12_common_habitats.pdf", p12, width = 9, height = 6)
  }
}

write_csv(species_summary_iucn, "9_Visualization/tables/TABLE_species_iucn_complete.csv")
message("Visualization module completed successfully.")



# ============================================================
# Stage 12: Synthesizing Macroevolutionary Cadastres (Manuscript Figures)
# ============================================================
cat("\n=======================================================\n")
cat("Stage 12: Synthesizing Macroevolutionary Cadastres (Manuscript Figures)\n")
cat("=======================================================\n")


# -------------------------------------------------------------
# Figure 13: Maximum-Likelihood Tree with TBE Support
# -------------------------------------------------------------
message("Rendering Figure 13 (maximum-likelihood tree)...")
# Load the support tree evaluated in Stage 9
supp_tree_file <- "7_Phylogenetics/cactus_support.raxml.support"

if (file.exists(supp_tree_file)) {
  ml_tree <- read.tree(supp_tree_file)
  ml_tree <- ape::ladderize(ml_tree)
  
  # Load constraints table required for annotation
  constraints_csv_path <- system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus")
  if (file.exists(constraints_csv_path) && nzchar(constraints_csv_path)) {
    constraints_tbl <- read_csv(constraints_csv_path, show_col_types = FALSE)
    
    # Ensure tree node labels match Nnode
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

    # Parse numerical support and map it to visual classes
    ml_tree_data <- tidytree::as_tibble(ml_tree) %>%
      mutate(
        support = suppressWarnings(as.numeric(label)),
        support_label = ifelse(!is.na(support), sprintf("%.3f", support), NA_character_),
        support_class = case_when(
          support >= 0.90 ~ ">= 90",
          support >= 0.70 ~ "70-89",
          TRUE ~ "< 70"
        ),
        support_class = factor(support_class, levels = c(">= 90", "70-89", "< 70"))
      )
    ml_tree_plot <- tidytree::as.treedata(ml_tree_data)

    p_ml <- ggtree(ml_tree_plot, size = 0.3) +
      theme_tree() +
      geom_nodepoint(aes(fill = support_class), shape = 21, size = 2, stroke = 0.2, na.rm = TRUE) +
      scale_fill_manual(
        values = c(">= 90" = "white", "70-89" = "grey", "< 70" = "black"),
        name = "TBE Support",
        na.translate = FALSE
      ) +
      theme(legend.position = "bottom")

    # Apply complex formatting using original script logic
    p_ml <- PhyloCactus::apply_clade_label_layers(p_ml, registry, ML_MAIN_LAYER_SPECS)
    # Add expansion for large trees
    xmax <- max(p_ml$data$x, na.rm = TRUE)
    p_ml <- p_ml + coord_cartesian(xlim = c(0, xmax + xmax * 0.35), clip = "off") + theme(plot.margin = margin(12, 12, 12, 12))

    ggsave("9_Visualization/figures/Figure_13_ML_Tree.pdf", p_ml, width = 8.5, height = 11)
    print(p_ml)

    # Export Figure 14 in A0 format
    ML_SUPP_LAYER_SPECS <- tibble::tribble(
      ~reg_name, ~fontsize, ~barsize, ~offset, ~offset_text, ~fontface, ~sort_desc, ~angle, ~align,
      "level_4_supp", 6, 0.60, 0.030, 0.0038, "plain", TRUE, 0, TRUE,
      "level_3", 8, 0.65, 0.150, 0.0038, "bold", TRUE, 270, TRUE,
      "level_2", 9, 0.70, 0.200, 0.0100, "bold", TRUE, 270, TRUE,
      "level_1", 10, 0.80, 0.250, 0.0100, "bold", TRUE, 270, TRUE
    )
    p_ml_supp <- ggtree(ml_tree_plot, size = 0.2) +
      theme_tree() +
      geom_tiplab(size = 0.70, align = FALSE, linetype = "solid", linesize = 0.15, colour = "grey20") +
      geom_nodelab(aes(label = support_label), size = 0.70, hjust = 0.7, colour = "black", na.rm = TRUE) +
      theme(legend.position = "none")
    p_ml_supp <- PhyloCactus::apply_clade_label_layers(p_ml_supp, registry, ML_SUPP_LAYER_SPECS)
    xmax_supp <- max(p_ml_supp$data$x, na.rm = TRUE)
    p_ml_supp <- p_ml_supp + coord_cartesian(xlim = c(0, xmax_supp + xmax_supp * 0.30), clip = "off") + theme(plot.margin = margin(12, 12, 12, 12))
    ggsave("9_Visualization/figures/Figure_14_ML_Tree_Extended.pdf", p_ml_supp, width = 36, height = 48)

    annotated_treedata <- PhyloCactus::augment_treedata_with_registry(ml_tree_plot, registry)
    treeio::write.beast(annotated_treedata, file = "9_Visualization/figures/annotated_ml_tree.tree")

  } else {
    message("Constraints CSV not found. Skipping annotation overlay on Figure 13.")
  }
} else {
  message("Skipping Figure 13: Support tree file not found (", supp_tree_file, ")")
}

# -------------------------------------------------------------
# Figure 15: Time-Calibrated Chronogram with HPD intervals
# -------------------------------------------------------------
message("Rendering Figure 15 (Chronogram)...")
# Load the TreeAnnotator summary chronogram from Stage 10
chrono_file <- "8_Dating/dated_summary_hpd.tree"

if (file.exists(chrono_file)) {
  chrono_beast <- treeio::read.beast(chrono_file)
  
  constraints_csv_path <- system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus")
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

    p_chrono <- ggtree(chrono_beast, size = 0.3) +
      theme_tree2() + 
      geom_range(range = 'height_0.95_HPD', color = 'gray60', alpha = 0.4, linewidth = 1.5) +
      theme_tree()
    
    p_chrono <- PhyloCactus::apply_clade_label_layers(p_chrono, registry_chrono, CHRONO_MAIN_LAYER_SPECS)
    p_chrono <- PhyloCactus::add_chronogram_axis(p_chrono, tree_phy, by = 5, digits = 0L, segment_size = 3, title_margin_top = 30)
    
    xmax <- max(p_chrono$data$x, na.rm = TRUE)
    p_chrono <- p_chrono + coord_cartesian(xlim = c(0, xmax + xmax * 0.35), clip = "off") + theme(plot.margin = margin(12, 12, 12, 12))
      
    ggsave("9_Visualization/figures/Figure_15_Chronogram.pdf", p_chrono, width = 8.5, height = 11)
    print(p_chrono)

    # Export Figure 16 in A0 format
    CHRONO_SUPP_LAYER_SPECS <- tibble::tribble(
      ~reg_name, ~fontsize, ~barsize, ~offset, ~offset_text, ~fontface, ~sort_desc, ~angle, ~align,
      "level_4_supp", 6, 0.60, 0.030, 0.0038, "plain", TRUE, 0, TRUE,
      "level_3", 8, 0.65, 0.150, 0.0038, "bold", TRUE, 270, TRUE,
      "level_2", 9, 0.70, 0.200, 0.0100, "bold", TRUE, 270, TRUE,
      "level_1", 10, 0.80, 0.250, 0.0100, "bold", TRUE, 270, TRUE
    )
    p_chrono_supp <- ggtree(chrono_beast, size = 0.2) +
      theme_tree2() +
      geom_tiplab(size = 0.85, align = FALSE, linetype = "solid", linesize = 0.15, colour = "grey20") +
      geom_range(range = 'height_0.95_HPD', color = 'gray60', alpha = 0.4, linewidth = 0.8) +
      theme_tree()
    p_chrono_supp <- PhyloCactus::apply_clade_label_layers(p_chrono_supp, registry_chrono, CHRONO_SUPP_LAYER_SPECS)
    p_chrono_supp <- PhyloCactus::add_chronogram_axis(p_chrono_supp, tree_phy, by = 5, digits = 0L, segment_size = 10, title_margin_top = 28, bar_size = 30)
    xmax_chrono_supp <- max(p_chrono_supp$data$x, na.rm = TRUE)
    p_chrono_supp <- p_chrono_supp + coord_cartesian(xlim = c(0, xmax_chrono_supp + xmax_chrono_supp * 0.30), clip = "off") + theme(plot.margin = margin(12, 12, 32, 12))
    ggsave("9_Visualization/figures/Figure_16_Chronogram_Extended.pdf", p_chrono_supp, width = 36, height = 48)

    annotated_chrono <- PhyloCactus::augment_treedata_with_registry(chrono_beast, registry_chrono)
    treeio::write.beast(annotated_chrono, file = "9_Visualization/figures/annotated_chronogram.tree")

  }
} else {
  message("Skipping Figure 15: Chronogram file not found (", chrono_file, ")")
}

message("Manuscript figures rendered successfully in '9_Visualization/figures' folder.")