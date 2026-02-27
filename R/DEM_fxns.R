#' DEM fxns
#' Functions to download DEM data and extract
#' elevation information

# Returns the filepath of the VRT
query_cded <- function(regions, output_dir) {
  # Function health checks
  stopifnot("`regions` must be an `sf` object with POLYGON geometry." = all(sf::st_is(regions, "POLYGON")))
  # Create output directory for DEM VRT
  # (actual raster tiles downloaded in secret bcmaps folder deep in drive)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  # Download DEM data
  out <- bcmaps::cded(aoi = regions,
                      dest_vrt = file.path(output_dir, "CDED_VRT.vrt"))
  return(out)
}