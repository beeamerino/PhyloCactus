# ============================================================
# IUCN Helper Functions
# ============================================================

#' Collapse unique values
#' @noRd
collapse_unique <- function(x) {
  x <- unique(x)
  x <- x[!is.na(x)]
  if(length(x) == 0){ return(NA_character_) }
  paste(x, collapse = "; ")
}

#' Extract description and code from IUCN data frame
#' @noRd
extract_desc_code <- function(df) {
  if(is.null(df) || !is.data.frame(df) || nrow(df) == 0){
    return(list(names = NA_character_, codes = NA_character_))
  }
  names_out <- NA_character_
  codes_out <- NA_character_
  if("description" %in% names(df)){ names_out <- collapse_unique(df$description$en) }
  if("code" %in% names(df)){ codes_out <- collapse_unique(df$code) }
  list(names = names_out, codes = codes_out)
}

#' Extract locations
#' @noRd
extract_locations <- function(locs) {
  empty <- list(names = NA_character_, codes = NA_character_, n = NA_integer_,
                endemic_names = NA_character_, endemic_codes = NA_character_,
                endemic_n = NA_integer_, is_endemic = FALSE)
  if(is.null(locs) || !is.data.frame(locs) || nrow(locs) == 0){ return(empty) }
  
  all_names <- unique(locs$description$en)
  all_codes <- unique(locs$code)
  all_names <- all_names[!is.na(all_names)]
  all_codes <- all_codes[!is.na(all_codes)]
  
  endemic_rows <- locs$is_endemic %in% TRUE
  endemic_names <- unique(locs$description$en[endemic_rows])
  endemic_codes <- unique(locs$code[endemic_rows])
  endemic_names <- endemic_names[!is.na(endemic_names)]
  endemic_codes <- endemic_codes[!is.na(endemic_codes)]
  
  list(names = collapse_unique(all_names),
       codes = collapse_unique(all_codes),
       n = length(all_names),
       endemic_names = collapse_unique(endemic_names),
       endemic_codes = collapse_unique(endemic_codes),
       endemic_n = length(endemic_names),
       is_endemic = length(endemic_names) > 0)
}

#' Fetch IUCN Red List Metadata via rredlist API
#'
#' Queries the IUCN Red List API through `rredlist` and returns the latest assessment of one species: the
#' Red List category (e.g., CR, EN, VU) and its criteria, the population trend, the countries of occurrence
#' and of endemism, habitats, threats, and the texts of the assessment.
#' Enriching phylogenetic datasets with IUCN conservation statuses enables evolutionary distinctness and extinction vulnerability analyses.
#'
#' The table returned is Red List Data under the IUCN Red List Terms of Use: it is for the user's own
#' analysis, queried with the user's own token, and is not redistributed (sections 1 and 4 of the terms).
#' `summarise_iucn_registry()` reduces it to the category of each species and to counts, which is what
#' the package distributes.
#'
#' @param sp_name Character. Standardized binomial scientific name of the target species (e.g., `"Astrophytum asterias"`).
#' @return A single-row `tibble` containing extracted IUCN conservation metadata fields.
#' @examples
#' \dontrun{
#' get_iucn_data("Astrophytum asterias")
#' }
#' @export
get_iucn_data <- function(sp_name) {
  message("Processing: ", sp_name)
  out <- tibble::tibble(species_clean = sp_name, iucn_found = FALSE, iucn_category = NA_character_,
                iucn_category_name = NA_character_, iucn_year = NA_character_, iucn_latest = NA,
                iucn_criteria = NA_character_, iucn_url = NA_character_, iucn_family = NA_character_,
                iucn_authority = NA_character_, iucn_common_name = NA_character_, iucn_synonyms = NA_character_,
                iucn_population_trend = NA_character_, iucn_population_trend_code = NA_character_,
                iucn_locations = NA_character_, iucn_location_codes = NA_character_, iucn_num_locations = NA_integer_,
                iucn_endemic_locations = NA_character_, iucn_endemic_location_codes = NA_character_, iucn_num_endemic_locations = NA_integer_,
                iucn_is_endemic = FALSE, iucn_habitat_names = NA_character_, iucn_habitat_codes = NA_character_,
                iucn_threat_names = NA_character_, iucn_threat_codes = NA_character_, iucn_system_names = NA_character_,
                iucn_system_codes = NA_character_, iucn_realm_names = NA_character_, iucn_realm_codes = NA_character_,
                iucn_growth_form_names = NA_character_, iucn_growth_form_codes = NA_character_, iucn_use_trade_names = NA_character_,
                iucn_use_trade_codes = NA_character_, iucn_doc_range = NA_character_, iucn_doc_habitats = NA_character_,
                iucn_doc_threats = NA_character_, iucn_doc_rationale = NA_character_)
  if (!requireNamespace("rredlist", quietly = TRUE)) {
    stop("Package 'rredlist' is required for IUCN Red List queries. Please install it with install.packages('rredlist').", call. = FALSE)
  }
  
  out$query_status <- "success"
  out$error_message <- NA_character_
  tryCatch({
    parts <- strsplit(sp_name, " ")[[1]]
    if(length(parts) < 2){ return(out) }
    res <- rredlist::rl_species(parts[1], parts[2])
    if(is.null(res$assessments) || nrow(res$assessments) == 0){ return(out) }
    
    out$iucn_found <- TRUE
    out$iucn_family <- res$taxon$family_name
    out$iucn_authority <- res$taxon$authority
    out$iucn_common_name <- collapse_unique(res$taxon$common_names$name)
    out$iucn_synonyms <- collapse_unique(res$taxon$synonyms$name)
    
    ass <- res$assessments
    ass_latest <- ass[ass$latest == TRUE, ]
    if(nrow(ass_latest) == 0){ ass <- ass[1, ] } else { ass <- ass_latest[1,] }
    
    out$iucn_category <- ass$red_list_category_code
    out$iucn_year <- ass$year_published
    out$iucn_latest <- ass$latest
    out$iucn_url <- ass$url
    
    full <- rredlist::rl_assessment(ass$assessment_id)
    out$iucn_category_name <- full$red_list_category$description$en
    out$iucn_criteria <- full$criteria
    if(!is.null(full$population_trend)) {
      out$iucn_population_trend <- full$population_trend$description$en
      out$iucn_population_trend_code <- full$population_trend$code
    }
    
    locs <- extract_locations(full$locations)
    out$iucn_locations <- locs$names
    out$iucn_location_codes <- locs$codes
    out$iucn_num_locations <- locs$n
    out$iucn_endemic_locations <- locs$endemic_names
    out$iucn_endemic_location_codes <- locs$endemic_codes
    out$iucn_num_endemic_locations <- locs$endemic_n
    out$iucn_is_endemic <- locs$is_endemic
    
    habs <- extract_desc_code(full$habitats)
    out$iucn_habitat_names <- habs$names
    out$iucn_habitat_codes <- habs$codes
    
    threats <- extract_desc_code(full$threats)
    out$iucn_threat_names <- threats$names
    out$iucn_threat_codes <- threats$codes
    
    systems <- extract_desc_code(full$systems)
    out$iucn_system_names <- systems$names
    out$iucn_system_codes <- systems$codes
    
    realms <- extract_desc_code(full$biogeographical_realms)
    out$iucn_realm_names <- realms$names
    out$iucn_realm_codes <- realms$codes
    
    growth <- extract_desc_code(full$growth_forms)
    out$iucn_growth_form_names <- growth$names
    out$iucn_growth_form_codes <- growth$codes
    
    trade <- extract_desc_code(full$use_and_trade)
    out$iucn_use_trade_names <- trade$names
    out$iucn_use_trade_codes <- trade$codes
    
    if(!is.null(full$documentation)){
      out$iucn_doc_range <- full$documentation$range
      out$iucn_doc_habitats <- full$documentation$habitats
      out$iucn_doc_threats <- full$documentation$threats
      out$iucn_doc_rationale <- full$documentation$rationale
    }
    
    Sys.sleep(1) # Sleep to avoid rate limiting
    out
  }, error = function(e) {
    out$query_status <<- "error"
    out$error_message <<- e$message
    message("IUCN query failed for ", sp_name, ": ", e$message)
    out
  })
}

#' Reduce the IUCN Red List table of a run to the category of each species and to counts
#'
#' Reads the species table that Module 11 writes with the assessments of `get_iucn_data()`
#' (`9_Visualization/tables/TABLE_species_all_data.csv`) and returns three tables that hold no
#' assessment content: the species registry with `iucn_found`, `iucn_category` and
#' `iucn_category_name` as its only IUCN columns, and two count tables. Under the IUCN Red List Terms
#' of Use (version 3.1), IUCN places no restriction on the use of the category of a named taxon
#' (section 3), whereas the assessments themselves are not to be redistributed (section 4); counts are
#' published with the citation of the Red List and its version (section 6). The figures of endemism,
#' threats and habitats of Tutorial 3 are drawn from the count tables.
#'
#' The counts cover the ingroup species in the maximum-likelihood tree (`is_ingroup` and `in_ml_tree`
#' true), as the figures of Tutorial 3 do. Countries, threats and habitats are split on `"; "`, the
#' separator `get_iucn_data()` writes, so a species is counted once in every country, threat class or
#' habitat class it has.
#'
#' @param species_table A data frame, or the path to the CSV file, with the columns of the species
#'   registry and of `get_iucn_data()`.
#' @param redlist_version Character. Version of the Red List the assessments were taken from, as
#'   `rredlist::rl_version()` reports it (e.g., `"2026-1"`).
#' @param access_date Character. Date the assessments were retrieved, `YYYY-MM-DD`.
#' @param output_dir Character or `NULL`. Folder where the three tables are written as
#'   `TABLE_species_registry.csv`, `TABLE_iucn_endemic_by_country.csv` and `TABLE_iucn_classes.csv`.
#'   `NULL` (the default) writes nothing.
#' @return A list of three tibbles:
#'   * `registry`: the species table without its `iucn_` columns other than `iucn_found`,
#'     `iucn_category` and `iucn_category_name`;
#'   * `endemic_by_country`: `country`, `iucn_category`, `n` (endemic species of each country in each
#'     category), `redlist_version`, `access_date`;
#'   * `classes`: `field` (`threat`, `habitat` or `endemic_species_in_tree`), `class`, `n`,
#'     `redlist_version`, `access_date`. The row of `endemic_species_in_tree` (class `all`) is the
#'     number of endemic species, each counted once.
#' @references
#' IUCN (2026). *The IUCN Red List of Threatened Species*. Version 2026-1. <https://www.iucnredlist.org>
#' @examples
#' \dontrun{
#' summarise_iucn_registry("9_Visualization/tables/TABLE_species_all_data.csv",
#'                         redlist_version = "2026-1", access_date = "2026-10-15",
#'                         output_dir = "9_Visualization/tables")
#' }
#' @export
summarise_iucn_registry <- function(species_table, redlist_version, access_date, output_dir = NULL) {
  if (missing(redlist_version) || !is.character(redlist_version) || length(redlist_version) != 1L ||
      !nzchar(redlist_version)) {
    stop("`redlist_version` is required: the version of the Red List of the assessments (rredlist::rl_version()).",
         call. = FALSE)
  }
  if (missing(access_date) || length(access_date) != 1L || is.na(as.Date(access_date, optional = TRUE))) {
    stop("`access_date` is required: the date the assessments were retrieved, as YYYY-MM-DD.", call. = FALSE)
  }
  access_date <- format(as.Date(access_date))
  sp <- if (is.character(species_table)) {
    readr::read_csv(species_table, show_col_types = FALSE, guess_max = 10000)
  } else {
    tibble::as_tibble(species_table)
  }
  needed <- c("in_ml_tree", "is_ingroup", "iucn_category", "iucn_is_endemic", "iucn_endemic_locations",
              "iucn_threat_names", "iucn_habitat_names")
  absent <- setdiff(needed, names(sp))
  if (length(absent)) stop("The species table lacks the columns: ", paste(absent, collapse = ", "), call. = FALSE)

  kept <- c("iucn_found", "iucn_category", "iucn_category_name")
  registry <- sp[, !grepl("^iucn_", names(sp)) | names(sp) %in% kept, drop = FALSE]

  in_tree <- sp[sp$in_ml_tree %in% TRUE & sp$is_ingroup %in% TRUE, , drop = FALSE]
  endemic <- in_tree[in_tree$iucn_is_endemic %in% TRUE, , drop = FALSE]
  stamp <- function(d) {
    d$n <- as.integer(d$n)
    d$redlist_version <- redlist_version
    d$access_date <- access_date
    d
  }

  endemic_by_country <- endemic |>
    dplyr::select(country = "iucn_endemic_locations", "iucn_category") |>
    tidyr::separate_rows("country", sep = "; ") |>
    dplyr::filter(!is.na(.data$country), nzchar(.data$country)) |>
    dplyr::count(dplyr::across(dplyr::all_of(c("country", "iucn_category"))), name = "n") |>
    dplyr::arrange(.data$country, .data$iucn_category) |>
    stamp()

  count_field <- function(column, field) {
    in_tree |>
      dplyr::select(class = dplyr::all_of(column)) |>
      dplyr::filter(!is.na(.data$class)) |>
      tidyr::separate_rows("class", sep = "; ") |>
      dplyr::filter(nzchar(.data$class)) |>
      dplyr::count(dplyr::across(dplyr::all_of("class")), name = "n") |>
      dplyr::mutate(field = field, .before = 1)
  }
  classes <- dplyr::bind_rows(
    count_field("iucn_threat_names", "threat"),
    count_field("iucn_habitat_names", "habitat"),
    tibble::tibble(field = "endemic_species_in_tree", class = "all", n = nrow(endemic))
  ) |>
    dplyr::arrange(match(.data$field, c("threat", "habitat", "endemic_species_in_tree")),
                   dplyr::desc(.data$n), .data$class) |>
    stamp()

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(registry, file.path(output_dir, "TABLE_species_registry.csv"), na = "NA")
    readr::write_csv(endemic_by_country, file.path(output_dir, "TABLE_iucn_endemic_by_country.csv"))
    readr::write_csv(classes, file.path(output_dir, "TABLE_iucn_classes.csv"))
  }
  list(registry = registry, endemic_by_country = endemic_by_country, classes = classes)
}

#' Standardized IUCN Category Color Scale for ggplot2
#' 
#' Applies a standardized, publication-grade hex color palette representing official IUCN Red List threat categories
#' (EX, EW, CR, EN, VU, NT, LC, DD, NE) onto `ggplot2` comparative phylogenetic graphics.
#'
#' @param ... Additional arguments passed directly to `ggplot2::scale_fill_manual`.
#' @return A `ggplot2` manual scale fill layer object.
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(df, aes(x = species, fill = iucn_category)) +
#'   geom_bar() +
#'   scale_fill_iucn()
#' }
#' @export
scale_fill_iucn <- function(...) {
  ggplot2::scale_fill_manual(values = phylocactus_colours("iucn"), drop = FALSE,
                             na.value = phylocactus_colours("neutral")[["line"]], ...)
}
