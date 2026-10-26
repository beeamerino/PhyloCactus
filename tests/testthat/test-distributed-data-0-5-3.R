# Closure of 0.5.2, phase 2 (audit of 26-10-2026, C2 and m10; plan, sections 5.2 and 5.4): the package
# distributes the IUCN category of each species and counts, not the assessments (IUCN Red List Terms
# of Use v3.1, sections 3 and 4); the checklist of Caryophyllales.org (Korotkova et al. 2021) and the
# four published reference trees are read from local paths, not from inst/extdata. Written before
# the code.

.dd_extdata <- function(name = NULL) {
  d <- system.file("extdata", package = "PhyloCactus")
  skip_if(!nzchar(d), "extdata not installed")
  if (is.null(name)) d else file.path(d, name)
}

.dd_root <- function(...) testthat::test_path("..", "..", ...)

.dd_read <- function(...) {
  f <- .dd_root(...)
  skip_if(!file.exists(f), paste(basename(f), "is not in this tree"))
  paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

.dd_synthetic_checklist <- function() testthat::test_path("testdata", "checklist_synthetic.xlsx")

.dd_kept_iucn <- c("iucn_found", "iucn_category", "iucn_category_name")

# ---- IUCN (plan 5.2) ---------------------------------------------------------------------------

test_that("the distributed registry keeps only the IUCN category among the iucn_ columns", {
  sp <- utils::read.csv(.dd_extdata("phylocactus_table_species.csv"), nrows = 5, check.names = FALSE)
  iucn <- grep("^iucn_", names(sp), value = TRUE)
  expect_setequal(iucn, .dd_kept_iucn)
})

test_that("no distributed file carries the text, places or threats of an assessment", {
  hits <- character(0)
  for (f in list.files(.dd_extdata(), pattern = "\\.(csv|tsv|txt)$", full.names = TRUE)) {
    x <- readLines(f, warn = FALSE)
    if (any(grepl("iucn_doc|iucn_locations|iucn_threat_names|iucn_habitat_names|iucn_endemic_locations", x)))
      hits <- c(hits, basename(f))
  }
  expect_identical(hits, character(0))
})

test_that("the two count tables are distributed with the version and the access date of the Red List", {
  by_country <- utils::read.csv(.dd_extdata("phylocactus_table_iucn_endemic_by_country.csv"), stringsAsFactors = FALSE)
  classes <- utils::read.csv(.dd_extdata("phylocactus_table_iucn_classes.csv"), stringsAsFactors = FALSE)
  expect_identical(names(by_country), c("country", "iucn_category", "n", "redlist_version", "access_date"))
  expect_identical(names(classes), c("field", "class", "n", "redlist_version", "access_date"))
  expect_setequal(unique(classes$field), c("threat", "habitat", "endemic_species_in_tree"))
  for (d in list(by_country, classes)) {
    expect_identical(unique(d$redlist_version), "2026-1")
    expect_identical(unique(d$access_date), "2026-10-15")
    expect_true(all(d$n > 0))
  }
})

test_that("the count tables give the counts of the reference run of 0.5.2", {
  by_country <- utils::read.csv(.dd_extdata("phylocactus_table_iucn_endemic_by_country.csv"), stringsAsFactors = FALSE)
  classes <- utils::read.csv(.dd_extdata("phylocactus_table_iucn_classes.csv"), stringsAsFactors = FALSE)
  # 633 endemic ingroup species in the tree (Table 1 of Tutorial 3 in 0.5.2)
  expect_identical(classes$n[classes$field == "endemic_species_in_tree"], 633L)
  per_country <- tapply(by_country$n, by_country$country, sum)
  expect_identical(unname(per_country[c("Mexico", "Brazil", "Chile", "Argentina", "Peru")]),
                   c(310L, 95L, 59L, 54L, 51L))
})

test_that("summarise_iucn_registry() gives the known answer on a synthetic table", {
  sp <- data.frame(
    species = c("A_a", "A_b", "B_c", "B_d", "C_e", "O_f"),
    status_pipeline = c("sampled_and_dated", "sampled_and_dated", "sampled_and_dated",
                        "unsampled_checklist_gap", "sampled_and_dated", "outgroup_reference"),
    in_ml_tree = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
    is_ingroup = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
    iucn_found = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
    iucn_category = c("CR", "LC", "EN", "VU", NA, "LC"),
    iucn_category_name = c("Critically Endangered", "Least Concern", "Endangered", "Vulnerable", NA, "Least Concern"),
    iucn_is_endemic = c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE),
    iucn_endemic_locations = c("Mexico", "Mexico; Chile", NA, "Peru", NA, "Chile"),
    iucn_threat_names = c("Mining; Roads", "Roads", NA, "Mining", NA, "Roads"),
    iucn_habitat_names = c("Desert", NA, "Desert; Shrubland", "Desert", NA, "Desert"),
    iucn_doc_range = c("text", "text", "text", "text", NA, "text"),
    stringsAsFactors = FALSE
  )
  out <- summarise_iucn_registry(sp, redlist_version = "2026-1", access_date = "2026-10-15")

  # The registry loses every iucn_ column but the category
  expect_identical(names(out$registry),
                   c("species", "status_pipeline", "in_ml_tree", "is_ingroup", .dd_kept_iucn))
  expect_identical(nrow(out$registry), 6L)

  # Endemic ingroup species in the tree, by country and category: A_a and A_b only
  bc <- as.data.frame(out$endemic_by_country)
  expect_identical(bc$country, c("Chile", "Mexico", "Mexico"))
  expect_identical(bc$iucn_category, c("LC", "CR", "LC"))
  expect_identical(bc$n, c(1L, 1L, 1L))
  expect_identical(unique(bc$redlist_version), "2026-1")
  expect_identical(unique(bc$access_date), "2026-10-15")

  # Classes over the ingroup species in the tree (A_a, A_b, B_c, C_e)
  cl <- as.data.frame(out$classes)
  get <- function(field, class) cl$n[cl$field == field & cl$class == class]
  expect_identical(get("threat", "Roads"), 2L)
  expect_identical(get("threat", "Mining"), 1L)
  expect_identical(get("habitat", "Desert"), 2L)
  expect_identical(get("habitat", "Shrubland"), 1L)
  expect_identical(get("endemic_species_in_tree", "all"), 2L)
  expect_false(any(cl$class %in% NA))

  # Written to a folder when asked
  d <- withr::local_tempdir()
  summarise_iucn_registry(sp, redlist_version = "2026-1", access_date = "2026-10-15", output_dir = d)
  expect_setequal(list.files(d), c("TABLE_species_registry.csv", "TABLE_iucn_endemic_by_country.csv",
                                   "TABLE_iucn_classes.csv"))
})

test_that("summarise_iucn_registry() asks for the version and the access date of the Red List", {
  sp <- data.frame(species = "A_a", in_ml_tree = TRUE, is_ingroup = TRUE, iucn_found = TRUE,
                   iucn_category = "LC", iucn_category_name = "Least Concern", iucn_is_endemic = FALSE,
                   iucn_endemic_locations = NA, iucn_threat_names = NA, iucn_habitat_names = NA)
  expect_error(summarise_iucn_registry(sp, access_date = "2026-10-15"), "redlist_version")
  expect_error(summarise_iucn_registry(sp, redlist_version = "2026-1"), "access_date")
})

test_that("the help page of get_iucn_data() says what it returns and that it is not redistributed (m10)", {
  rd <- .dd_read("man", "get_iucn_data.Rd")
  expect_false(grepl("verified macroevolutionary", rd, fixed = TRUE))
  expect_match(rd, "summarise_iucn_registry", fixed = TRUE)
})

.dd_assessment_columns <- "iucn_endemic_locations|iucn_threat_names|iucn_habitat_names|iucn_is_endemic"

# The evaluated chunks of a vignette: the code the website runs on the distributed files
.dd_evaluated_chunks <- function(v) {
  lines <- strsplit(v, "\n", fixed = TRUE)[[1]]
  out <- character(0)
  for (h in grep("^```\\{r", lines)) {
    end <- h + match(TRUE, grepl("^```$", lines[(h + 1):length(lines)]))
    if (!grepl("eval *= *FALSE", lines[h])) out <- c(out, lines[h:end])
  }
  paste(out, collapse = "\n")
}

test_that("Tutorial 3 draws Figures 14 to 17 and the endemism row of Table 1 from the count tables", {
  v <- .dd_read("vignettes", "tutorial-3-cactus-phylogeny-visualization.Rmd")
  ev <- .dd_evaluated_chunks(v)
  expect_false(grepl(.dd_assessment_columns, ev))
  expect_match(ev, "phylocactus_table_iucn_endemic_by_country.csv", fixed = TRUE)
  expect_match(ev, "phylocactus_table_iucn_classes.csv", fixed = TRUE)
  # The script of a user's run derives the two tables after Module 11 and draws the figures from them
  s <- strsplit(.dd_read("inst", "scripts", "tutorial-3-cactus-phylogeny-visualization.R"), "\n", fixed = TRUE)[[1]]
  call <- grep("summarise_iucn_registry(", s, fixed = TRUE)
  expect_length(call, 1L)
  after <- paste(s[(call + 1):length(s)], collapse = "\n")
  expect_false(grepl(.dd_assessment_columns, after))
  expect_match(after, "endemic_by_country", fixed = TRUE)
  expect_match(after, "classes", fixed = TRUE)
  expect_match(v, "Version 2026-1", fixed = TRUE)
  expect_match(v, "rredlist", fixed = TRUE)
})

# ---- Checklist (plan 5.4) ----------------------------------------------------------------------

test_that("no workbook and no reference tree is distributed in inst/extdata or declared in the manifest", {
  present <- list.files(.dd_extdata())
  expect_identical(grep("\\.xlsx?$", present, value = TRUE, ignore.case = TRUE), character(0))
  expect_identical(grep("^reference_.*\\.tree$", present, value = TRUE), character(0))
  m <- utils::read.csv(.dd_extdata("MANIFEST.csv"), stringsAsFactors = FALSE)
  expect_false(any(grepl("\\.xlsx?$|^reference_", m$canonical)))
})

test_that("phylocactus_checklist() stops without a path and names where to obtain the checklist", {
  withr::local_envvar(PHYLOCACTUS_CHECKLIST = "")
  err <- tryCatch(phylocactus_checklist(), error = function(e) conditionMessage(e))
  expect_type(err, "character")
  expect_match(err, "PHYLOCACTUS_CHECKLIST", fixed = TRUE)
  expect_match(err, "Caryophyllales.org", fixed = TRUE)
  expect_match(err, "Korotkova", fixed = TRUE)
  withr::local_envvar(PHYLOCACTUS_CHECKLIST = file.path(tempdir(), "no_such_checklist.xlsx"))
  expect_error(phylocactus_checklist(), "does not exist")
})

test_that("phylocactus_checklist() returns the path of the variable or of its argument, and warns on another md5", {
  f <- normalizePath(.dd_synthetic_checklist())
  withr::local_envvar(PHYLOCACTUS_CHECKLIST = f)
  w <- NULL
  p <- withCallingHandlers(suppressMessages(phylocactus_checklist()),
                           warning = function(cnd) { w <<- conditionMessage(cnd); invokeRestart("muffleWarning") })
  expect_identical(p, f)
  # The synthetic workbook is not the checklist of the reference run
  expect_match(w, "0ab7855cef6fe7db0e0105646d5fee08", fixed = TRUE)
  expect_match(w, unname(tools::md5sum(f)), fixed = TRUE)
  # The argument comes before the variable
  withr::local_envvar(PHYLOCACTUS_CHECKLIST = "")
  expect_identical(suppressWarnings(suppressMessages(phylocactus_checklist(f))), f)
})

test_that("phylocactus_checklist() reports the md5 of the file it returns", {
  f <- normalizePath(.dd_synthetic_checklist())
  expect_message(suppressWarnings(phylocactus_checklist(f)), unname(tools::md5sum(f)), fixed = TRUE)
})

test_that("no code, script or vignette reads the checklist or a reference tree from extdata", {
  files <- c(list.files(.dd_root("R"), "\\.R$", full.names = TRUE),
             list.files(.dd_root("inst", "scripts"), "\\.R$", full.names = TRUE),
             list.files(.dd_root("vignettes"), "\\.Rmd$", full.names = TRUE))
  skip_if(length(files) == 0L, "the sources are not in this tree")
  bad <- Filter(function(f) any(grepl('"extdata", *"(CactaceaeFullList|reference_)',
                                      readLines(f, warn = FALSE))), files)
  expect_identical(basename(bad), character(0))
})

test_that("the defaults of the barcoding branch take the checklist from phylocactus_checklist()", {
  expect_match(deparse(body(assemble_barcoding_dataset)), "phylocactus_checklist()", fixed = TRUE, all = FALSE)
  expect_identical(deparse(formals(extract_barcoding_plastome_loci)$checklist_path), "phylocactus_checklist()")
})

test_that("the readers of the checklist read the synthetic workbook", {
  skip_if_not_installed("readxl")
  f <- .dd_synthetic_checklist()
  cac <- .checklist_family_names(f, family = "Cactaceae")
  expect_true(all(c("Leuenbergeria_bleo", "Parodia_graessneri", "Chamaecereus_silvestrii") %in% cac))
  expect_false(any(c("Anacampseros_kurtzii", "Talinopsis_frutescens") %in% cac))
  nm <- .bc_read_checklist_names(f)
  expect_true(all(c("Leuenbergeria bleo", "Anacampseros kurtzii") %in% nm))
  # The Facts sheets are not read as names
  expect_identical(length(nm), 11L)
})

test_that("README and Tutorial 1 say how to obtain the checklist, with its version, md5 and variable", {
  for (f in list(c("README.md"), c("vignettes", "tutorial-1-cactus-phylogeny-prep.Rmd"))) {
    x <- do.call(.dd_read, as.list(f))
    for (k in c("Caryophyllales.org", "Korotkova", "Willdenowia", "PHYLOCACTUS_CHECKLIST",
                "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx", "0ab7855cef6fe7db0e0105646d5fee08"))
      expect_match(x, k, fixed = TRUE, info = paste(f[length(f)], k))
  }
})

# ---- Reference trees (plan 5.4) ----------------------------------------------------------------

test_that("Tutorial 4 reads the reference trees from local paths and draws Figure 1 from the run", {
  for (f in list(c("vignettes", "tutorial-4-cactus-phylogeny-validation.Rmd"),
                 c("inst", "scripts", "tutorial-4-cactus-phylogeny-validation.R"))) {
    x <- do.call(.dd_read, as.list(f))
    expect_match(x, "PHYLOCACTUS_REFERENCE_TREES", fixed = TRUE, info = f[length(f)])
    expect_match(x, "phylocactus_checklist(", fixed = TRUE, info = f[length(f)])
  }
  v <- .dd_read("vignettes", "tutorial-4-cactus-phylogeny-validation.Rmd")
  # Every chunk that needs a reference tree or the checklist is shown, not run
  expect_false(grepl("PHYLOCACTUS_REFERENCE_TREES|phylocactus_checklist\\(|validate_phylogenies\\(|ASTRAL_PATH",
                     .dd_evaluated_chunks(v)))
  expect_match(v, "phylocactus_table_tree_space.csv", fixed = TRUE)
  expect_match(v, "plot_tree_space(", fixed = TRUE)
  for (doi in c("10.1016/j.biocon.2022.109677", "10.1038/s41467-024-51666-2",
                "10.1038/s41586-024-07324-0", "10.1007/s00606-025-01948-z"))
    expect_match(v, doi, fixed = TRUE)
})

test_that("the tree-space coordinates of the reference run are distributed and redraw Figure 1", {
  co <- utils::read.csv(.dd_extdata("phylocactus_table_tree_space.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("analysis_id", "analysis_label", "tree", "Dimension_1", "Dimension_2", "metric") %in% names(co)))
  expect_setequal(unique(co$analysis_id), c("A1_species_all4", "A2_species_3",
                                            "A3_all_species_in_shared_genera", "A4_representative_shared_genera"))
  expect_setequal(unique(co$metric), c("RF", "DPI"))
  p <- plot_tree_space(co, metric = "RF")
  expect_s3_class(p, "ggplot")
})
