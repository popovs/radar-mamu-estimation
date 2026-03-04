#' Raster preparation fxns
#' Functions to: download DEM data and extract elevation information
#'

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

create_canvas <- function(target_rast) {
  # Extract res of target
  # Going to assume all data products are square res... 
  # so just take the first number and assume it's equal to the second
  res <- terra::res(target_rast)[1]
  extent <- terra::ext(target_rast)
  
  # If the resolution is too high, we're going to
  # make it a courser resolution by a factor of 10 - so each
  # pixel will be e.g. 250m x 250m rather than 25m x 25m - otherwise,
  # many of the subsequent raster calculations will fail or 
  # take far far too long. This means that we will be calculating
  # distance from the coast to a maximum of 1/4 km accuracy.
  
  # Create canvas
  if (res >= 100) {
    canvasrow <- nrow(target_rast)
    canvascol <- ncol(target_rast)
  } else {
    canvasrow <- nrow(target_rast) / 10
    canvascol <- ncol(target_rast) / 10
  }
  
  canvas <- terra::rast(extent = extent, # `extent` defined in 03-1 above 
                        nrow = canvasrow, 
                        ncol = canvascol, 
                        nlyr = 1)
  terra::values(canvas) <- 0 # assume everything is the sea (0)
  terra::crs(canvas) <- "epsg:3005"
  return(canvas)
}

block_land <- function(dem, canvas, land) {
  canvas <- terra::mask(canvas, land, inverse = TRUE) # Chunky interior land blocks
  canvas <- terra::ifel(is.na(canvas), 1, canvas) # masked areas == land == 1
  
  land <- terra::ifel(dem > 0, 1, 0) # now extract land areas from DEM
  land <- terra::resample(land, canvas) # resample `land` to match `canvas` extent/resolution
  
  out <- terra::merge(land, canvas, first = TRUE) # now merge the two
  return(out)
}
