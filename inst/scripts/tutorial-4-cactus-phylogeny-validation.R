# -------------------------------------------------------------
# PhyloCactus: Tutorial 4. Phylogenetic Validation and Comparative Analyses
# -------------------------------------------------------------
# Module 13: comparison of the focal tree with published phylogenies of Cactaceae.
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
# Reference Phylogenies and Optional ASTRAL-III Re-inference
# -------------------------------------------------------------
# The reference trees (Amaral et al. 2022, Thompson et al. 2024, Zuntini et al. 2024,
# and de Vos et al. 2025) are not distributed with the package. Obtain each from its
# source (table in Tutorial 4), save the four in one folder and set its path in ~/.Renviron:
# - PHYLOCACTUS_REFERENCE_TREES="/path/to/reference_trees", holding
#   reference_amaral_2022.tree      Amaral et al. 2022: dataset+tree_80MD/Tree_80MD.tre of
#                                   https://doi.org/10.6084/m9.figshare.17057816
#                                   (md5 185926b62a0a01b3187750e664e0fcf7)
#   reference_thompson_2024.tree    Thompson et al. 2024: Alignment, accessions and tree/treepl_cacti.tre of
#                                   https://github.com/jamie-thompson/cactaceae
#                                   (md5 df479c182a14028aa710b508aebaa458)
#   reference_zuntini_2024.tree     Zuntini et al. 2024: treeoflife.all_support_values.4.0.tree of the Kew
#                                   Tree of Life, release 4.0 (sftp.kew.org, pub/treeoflife), pruned to
#                                   its Cactaceae as Tutorial 4 shows (md5 of the reference run's file
#                                   319e6985225858b836ab2ff37054b3e5)
#   reference_devos_2025_astral.tree  de Vos et al. 2025: QC.Species_tree_astral.tree of Online Resource 4,
#                                   doi:10.1007/s00606-025-01948-z (md5 5ceb901a0555521da1a14f0815787362)
# - PHYLOCACTUS_CHECKLIST, the checklist of Korotkova et al. (2021) (Tutorial 1).
# Re-estimating the de Vos et al. (2025) species tree from its gene trees requires ASTRAL-III and Java:
# - Set ASTRAL_PATH in ~/.Renviron (e.g., ASTRAL_PATH="/path/to/astral.5.7.8.jar")
# - Set DEVOS_GENETREES_PATH and DEVOS_METADATA_PATH for the raw de Vos inputs

cat("\n=======================================================\n")
cat("Module 13: Comparison with published phylogenies\n")
cat("=======================================================\n")

# Set these paths in ~/.Renviron so that the script contains no local paths
# e.g., ASTRAL_PATH="/path/to/astral.5.7.8.jar"

# -------------------------------------------------------------
# Module 13a: Species tree of de Vos et al. (2025) with ASTRAL-III
# -------------------------------------------------------------
cat("\n--- Running ASTRAL-III Pipeline ---\n")

# de Vos et al. (2025) sequenced the Angiosperms353 loci. ASTRAL-III estimates the
# species tree from their gene trees under the multispecies coalescent, which
# accounts for gene tree discordance caused by incomplete lineage sorting.

astral_jar <- Sys.getenv("ASTRAL_PATH")

devos_gene_trees <- Sys.getenv("DEVOS_GENETREES_PATH")

devos_metadata <- Sys.getenv("DEVOS_METADATA_PATH")

java_bin <- Sys.which("java")
if (!nzchar(java_bin) && file.exists("/opt/homebrew/opt/openjdk/bin/java")) {
  java_bin <- "/opt/homebrew/opt/openjdk/bin/java"
}

# Set to TRUE to run ASTRAL-III; FALSE (default) copies the precomputed tree.
run_astral_locally <- FALSE

output_tree <- "10_Validation/QC.Species_tree_astral.tree"
renamed_gene_trees_file <- "10_Validation/QC.Species_tree_astral_rename.tree"

if (run_astral_locally && nzchar(astral_jar) && file.exists(astral_jar) && 
    nzchar(devos_gene_trees) && file.exists(devos_gene_trees) && 
    nzchar(devos_metadata) && file.exists(devos_metadata) && 
    nzchar(java_bin)) {
  
  cat("Reconstructing de Vos et al. (2025) species tree using ASTRAL-III...\n")
  meta <- read.csv(devos_metadata, stringsAsFactors = FALSE, check.names = FALSE)
  meta$species_name <- gsub(" ", "_", meta$`Scientific name`)
  meta$Sample_from_tree <- str_extract(meta$`Name in tree`, "P[0-9]+")
  
  map_sample <- setNames(meta$species_name, meta$Sample_from_tree)
  map_srr    <- setNames(meta$species_name, meta$`ENA run acc.`)
  
  gene_trees <- read.tree(devos_gene_trees)
  
  rename_tips <- function(tree) {
    tips <- tree$tip.label
    paftol_to_sample <- function(x) {
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
  
  cmd <- paste(shQuote(java_bin), "-jar", shQuote(astral_jar), 
               "-i", shQuote(renamed_gene_trees_file), 
               "-o", shQuote(output_tree), "2>/dev/null")
  system(cmd)
  
  # Also deposit in trees/ subdirectory
  file.copy(output_tree, "10_Validation/trees/QC.Species_tree_astral.tree", overwrite = TRUE)
  cat("ASTRAL-III execution complete: tree saved to", output_tree, "\n")
} else {
  cat("External ASTRAL-III binary or raw gene trees not detected.\n")
  cat("Copying the de Vos et al. (2025) species tree from the folder of reference trees...\n")
  ref_tree_local <- file.path(Sys.getenv("PHYLOCACTUS_REFERENCE_TREES"), "reference_devos_2025_astral.tree")
  file.copy(ref_tree_local, output_tree, overwrite = TRUE)
  file.copy(ref_tree_local, "10_Validation/trees/QC.Species_tree_astral.tree", overwrite = TRUE)
  cat("De Vos et al. (2025) tree created in tutorial directory:", output_tree, "\n")
}

# -------------------------------------------------------------
# Module 13b: Topological comparison
# -------------------------------------------------------------
cat("\n--- Running Phylogenetic Validation Pipeline ---\n")

# The focal tree is the chronogram of this run (Module 10); the reference trees are read from the
# folder in PHYLOCACTUS_REFERENCE_TREES. When 8_Dating/ is absent, the reference chronogram
# distributed with the package is used as the focal tree.
focal_tree <- file.path("8_Dating", "dated_summary_hpd.tree")
if (!file.exists(focal_tree)) {
  focal_tree <- system.file("extdata", "phylocactus_chronogram_hpd.tree", package = "PhyloCactus")
}
reference_dir <- Sys.getenv("PHYLOCACTUS_REFERENCE_TREES")
tree_paths <- list(
  FocalTree   = focal_tree,
  Zuntini     = file.path(reference_dir, "reference_zuntini_2024.tree"),
  Thompson    = file.path(reference_dir, "reference_thompson_2024.tree"),
  Amaral      = file.path(reference_dir, "reference_amaral_2022.tree"),
  deVos       = file.path(reference_dir, "reference_devos_2025_astral.tree")
)

# Filter missing trees
existing_trees <- tree_paths[sapply(tree_paths, function(x) nzchar(x) && file.exists(x))]

if(length(existing_trees) > 1) {
  
  checklist_file   <- phylocactus_checklist()
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
  
  cat("Validation Complete! Results exported to 10_Validation/\n")
} else {
  cat("Fewer than two trees found: set PHYLOCACTUS_REFERENCE_TREES to the folder of the reference trees.\n")
}
