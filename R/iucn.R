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
#' Queries the IUCN Red List API (`rredlist`) to retrieve verified macroevolutionary conservation threat categories (e.g., CR, EN, VU),
#' population trends, habitat descriptions, and geographic endemism metrics for a focal species.
#' Enriching phylogenetic datasets with IUCN conservation statuses enables evolutionary distinctness and extinction vulnerability analyses.
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
  }, error = function(e){
    message("FAILED: ", sp_name)
    out
  })
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
  iucn_colors <- c(
    "EX" = "#000000",
    "EW" = "#5C155A",
    "CR" = "#D81E05",
    "EN" = "#FC7F3F",
    "VU" = "#F9E814",
    "NT" = "#CCE226",
    "LC" = "#60C659",
    "DD" = "#D1D1C6",
    "NE" = "#FFFFFF"
  )
  ggplot2::scale_fill_manual(values = iucn_colors, drop = FALSE, na.value = "grey50", ...)
}
