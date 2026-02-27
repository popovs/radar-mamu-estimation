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

# Resample a DEM to a given target resolution
resample_rast <- function(rast, res) {
  # CHECK resolution of dem - only 'resample' if it doesn't match
  if (!all(res == terra::res(rast))) {
    message("Resampling DEM to target resolution of ", res, " m...")
    r <- rast
    terra::res(r) <- res
    rast <- terra::resample(rast, r)
  }
  return(rast)
}