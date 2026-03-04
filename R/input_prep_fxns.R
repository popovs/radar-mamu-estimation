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