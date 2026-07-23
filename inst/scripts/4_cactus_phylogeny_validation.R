# -------------------------------------------------------------
# PhyloCactus: Tutorial 4 - Phylogenetic Validation and Sub-tree Comparisons
# -------------------------------------------------------------
# This script covers Stage 13 of the phylogenetic pipeline:
# Sub-trees, ASTRAL-III integration, and Phylogenetic Validation.
# -------------------------------------------------------------
library(PhyloCactus)
library(ape)
library(stringr)
library(readr)
library(dplyr)

# -------------------------------------------------------------
# Setup: Creating a Clean Workspace
# -------------------------------------------------------------
tutorial_dir <- "~/Desktop/PhyloCactus_Tutorial"
dir.create(tutorial_dir, showWarnings = FALSE)
setwd(tutorial_dir)

# Create validation directories
dir.create("10_Validation/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("10_Validation/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("10_Validation/trees", recursive = TRUE, showWarnings = FALSE)
dir.create("10_Validation/logs", recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------
# Downloading External Resources
# -------------------------------------------------------------
# The reference trees and ASTRAL-III are not bundled due to size limits.
# You must download them manually and provide their paths.
# See vignette 4 for the DOIs and download instructions:
# - Amaral et al. (2022) DOI: 10.1016/j.biocon.2022.109677 (Tree_80MD.tre)
# - Thompson et al. (2024) DOI: 10.1038/s41467-024-51666-2 (ultra_cacti_JT.tre)
# - de Vos et al. (2025) DOI: 10.1007/s00606-025-01948-z (QC.bestTreeCollapsed.trees & metadata)

cat("\n=======================================================\n")
cat("Stage 13: Quantifying Topological Congruence and Validating Evolutionary Hypotheses\n")
cat("=======================================================\n")

# To keep your local paths secure, it is recommended to set these in ~/.Renviron
# e.g., ASTRAL_PATH="/path/to/astral.5.7.8.jar"
# Here, we will assume you have placed the files in the package's extdata folder 
# for ease of execution during this tutorial.

# -------------------------------------------------------------
# ASTRAL-III PIPELINE FROM R
# -------------------------------------------------------------
cat("\n--- Running ASTRAL-III Pipeline ---\n")

# Context: Why use ASTRAL-III?
# de Vos et al. (2025) is a phylogenomic study using Angiosperms353 to reconstruct 
# the family phylogeny at the genus level, making it a great element for comparison. 
# We use ASTRAL-III here simply to extract their gene trees and reconstruct the 
# species tree because the authors did not share the final coalescent species tree.

astral_jar <- Sys.getenv("ASTRAL_PATH") # Provide path if installed
devos_gene_trees <- Sys.getenv("DEVOS_GENETREES_PATH") # Path to downloaded QC.bestTreeCollapsed.trees
devos_metadata <- Sys.getenv("DEVOS_METADATA_PATH") # Path to downloaded 606_2025_1948_MOESM1_ESM.csv

if(nchar(astral_jar) > 0 && nchar(devos_gene_trees) > 0 && file.exists(astral_jar) && file.exists(devos_gene_trees)) {
  
  output_tree <- "10_Validation/QC.Species_tree_astral.tree"
  renamed_gene_trees_file <- "10_Validation/QC.Species_tree_astral_rename.tree"
  
  meta <- read.csv(devos_metadata, stringsAsFactors = FALSE, check.names = FALSE)
  meta$species_name <- gsub(" ", "_", meta$`Scientific name`)
  meta$Sample_from_tree <- str_extract(meta$`Name in tree`, "P[0-9]+")
  
  map_sample <- setNames(meta$species_name, meta$Sample_from_tree)
  map_srr    <- setNames(meta$species_name, meta$`ENA run acc.`)
  
  gene_trees <- read.tree(devos_gene_trees)
  
  rename_tips <- function(tree){
    tips <- tree$tip.label
    paftol_to_sample <- function(x){
      num <- str_extract(x, "[0-9]+")
      paste0("P", str_sub(num, -5))
    }
    sample_ids <- ifelse(str_detect(tips, "PAFTOL"), paftol_to_sample(tips), NA)
    species_vec <- ifelse(str_detect(tips, "PAFTOL"), map_sample[sample_ids], map_srr[tips])
    species_vec[is.na(species_vec)] <- tips[is.na(species_vec)]
    tree$tip.label <- species_vec
    return(tree)
  }
  
  gene_trees <- lapply(gene_trees, rename_tips)
  class(gene_trees) <- "multiPhylo"
  write.tree(gene_trees, file = renamed_gene_trees_file)
  
  cmd <- paste("java -jar", astral_jar, "-i", renamed_gene_trees_file, "-o", output_tree)
  system(cmd)
  cat("ASTRAL-III execution complete.\n")
} else {
  cat("Skipping ASTRAL-III execution. ASTRAL_PATH or DEVOS_GENETREES_PATH is not set, or files were not found.\n")
}

# -------------------------------------------------------------
# PHYLOGENETIC VALIDATION PIPELINE
# -------------------------------------------------------------
cat("\n--- Running Phylogenetic Validation Pipeline ---\n")

# Fetch the absolute paths to reference trees configured in your ~/.Renviron
tree_paths <- list(
  FocalTree   = "8_Dating/dated_summary_hpd.tree",
  Thompson    = Sys.getenv("THOMPSON_TREE_PATH"),
  Amaral      = Sys.getenv("AMARAL_TREE_PATH"),
  deVos       = "10_Validation/QC.Species_tree_astral.tree"
)

# Filter missing trees
existing_trees <- tree_paths[sapply(tree_paths, function(x) file.exists(x) && nzchar(x))]

if(length(existing_trees) > 1) {
  
  checklist_file   <- system.file("extdata", "CactaceaeFullList_accepted.csv", package = "PhyloCactus")
  constraints_file <- system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus")
  
  validation_outputs <- validate_phylogenies(
    trees_mapping_list = existing_trees,
    checklist_csv = checklist_file,
    constraints_map = constraints_file,
    out_dir = "10_Validation"
  )
  
  if (!is.null(validation_outputs$plot)) {
    print(validation_outputs$plot)
  }
  
  cat("Validation Complete! Results exported to 10_Validations/\n")
} else {
  cat("Insufficient trees found at configured .Renviron paths to run full validation pipeline.\n")
  cat("Please ensure THOMPSON_TREE_PATH, AMARAL_TREE_PATH, and DEVOS_SPECIES_TREE_PATH are correctly set.\n")
}

