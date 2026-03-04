#' DEM fxns
#' Functions to download DEM data and extract
#' elevation information


# Query rnaturalearth for USA land polygons
get_usa_land <- function(extent) {
  usa <- rnaturalearth::ne_states("united states of america")
  usa <- usa[usa$name %in% c("Alaska", "Washington"), ]
  usa <- sf::st_as_sf(usa)
  usa <- sf::st_transform(usa, 3005)
  usa <- sf::st_crop(usa, extent)
  usa <- terra::vect(usa)
  return(usa)
}

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

reclass_rast <- function(rast, cutoffs, min_col, max_col, region_name) { 
  cutoffs <- cutoffs[cutoffs$region == region_name,]
  message("Reclassifying ", cutoffs$region, "...")
  min <- cutoffs[[min_col]]
  max <- cutoffs[[max_col]]
  rast <- terra::ifel(rast > 0, rast, NA)
  rast <- terra::ifel(rast < max, rast, NA)
  rast <- terra::ifel(rast < min, 1, 2)
  return(rast)
}

crop_rast <- function(rast, region_name, regions) {
  message("Cropping ", region_name)
  regions <- regions[regions$region == region_name, ]
  tmp <- terra::crop(rast, regions)
  tmp <- terra::mask(tmp, regions) # This can be done in one step with terra::crop(mask = T), but was resulting in buggy raster values - so doing it in two steps here.
  return(tmp)
}
