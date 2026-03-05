#' Data inputs preparation fxns
#' The functions here read in the base data files of the pipeline
#' and prepare them for analysis.


# PREPARE GIS FILES -------------------------------------------------------

prepare_regions <- function(filepath) {
  regions <- sf::st_read(filepath)
  # This assumes the gpkg polygons are already in correct order. 
  # Assign factor levels to the existing row order.
  regions$region <- factor(regions$region, levels = forcats::fct_inorder(regions$region))
  return(regions)
}

# This version doesn't have the uncertainty col
prepare_nests <- function(filepath, regions = regions) {
  nests <- sf::st_read(filepath)
  # Intersect to get region the nest falls within,
  # then condense repeat visits to a nest into one geometry record
  nests <- nests |> 
    sf::st_intersection(regions) |>
    dplyr::select(MMCR_NAME, region) |>
    unique()
  return(nests)
}


# PREPARE SURVEYS ---------------------------------------------------------

# Read survey data, clean up, convert to sf object, and assign pipeline regions
prepare_surveys <- function(filepath, regions = regions) {
  s <- read.csv(filepath)
  
  # Basic dataframe tidying
  s <- janitor::clean_names(s)
  s$survey_date <- lubridate::make_date(year = s$year, month = s$month, day = s$day)
  s$doy <- lubridate::yday(s$survey_date)
  s$site <- stringr::str_trim(s$site)
  s$site <- stringr::str_squish(s$site)
  
  # Create spatial object
  s <- s[!is.na(s$lon), ]
  s <- sf::st_as_sf(s, coords = c("lon", "lat"), 
                    crs = 4326, 
                    remove = FALSE) |>
    sf::st_transform(3005)
  
  # Merge in regions from pipeline
  # Anything already stored in the pre-existing `region` col
  # will be moved to `region_old`
  s$region_old <- s$region
  s <- dplyr::select(s, -region)
  regions <- sf::st_make_valid(regions)
  s <- sf::st_intersection(s, regions)
  return(s)
}
