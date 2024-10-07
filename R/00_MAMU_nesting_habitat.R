#' 0. MAMU NESTING HABITAT AREA

#' This script contains functions that will download the BC 
#' Digital Elevation Model (DEM) for our study region, then 
#' clip it down to the areas that are highly likely to contain 
#' suitable marbled murrelet nesting habitat area (regions 
#' that are 30 km from shore and under <1500 m elevation). 
#' The exact maximum elevation cutoff varies by our six 
#' conservation regions of interest (North,  Central, and South 
#' Mainland Coast, Haida Gwaii, North + West Vancouver Island, 
#' and East Vancouver Island).


# PREPARE GIS FILES -------------------------------------------------------

prepare_regions <- function(filepath = "GIS/regions.gpkg") {
  regions <- sf::st_read(filepath)
  regions$region <- factor(regions$region, 
                           levels = c("AKB", "HG", "NC", "CC", "SC", "NVI", "MWVI", "SWVI", "EVI"))
  return(regions)
}

prepare_nests <- function(filepath = "GIS/MAMU_nests.gpkg", 
                          uncertainty_col = "LOC_UNCE_1",
                          regions = regions) {
  nests <- sf::st_read(filepath)
  nests <- sf::st_intersection(nests, regions)
  # Add buffers of uncertainty around the nests
  if (any(is.na(nests[[uncertainty_col]]))) stop("Nest locations must have an uncertainty associated with them (", uncertainty_col, ").")
  message("Assuming nest uncertainty is in meters.")
  #stopifnot("Nest location uncertainly calculations currently only support uncertainty expressed in meters (LOC_UNCE_2 column)." = all(tolower(unique(nests$LOC_UNCE_2)) %in% c("m", "meters", "meter")))
  nests <- st_buffer(nests, dist = nests[[uncertainty_col]])
  return(nests)
}



# PREPARE DEM -------------------------------------------------------------


download_dem_tile <- function(tile,
                              output_dir,
                              save_output = TRUE,
                              overwrite = TRUE) {
  # Create output dir
  output_dir <- file.path(output_dir, "DEM_tiles")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  # Download the tile
  MAMU::BC_DEM(letterblock = tile, 
               save_output = save_output,
               overwrite = overwrite,
               output_dir = output_dir)
  out <- list.files(output_dir, full.names = TRUE)
  out <- out[grep(paste0(tile, ".dem$"), out)]
  return(out)
}

merge_dem <- function(vrt_path,
                      output_file,
                      overwrite = FALSE) {
  # Re-project VRT
  vrt <- terra::rast(vrt_path)
  dem_3005 <- terra::project(vrt, "EPSG:3005")
  # Save it
  terra::writeRaster(dem_3005,
                     filename = output_file, 
                     overwrite = overwrite)
  return(output_file)
}

resample_dem <- function(dem_path, res) {
  dem <- terra::rast(dem_path)
  # CHECK resolution of dem - only 'resample' if it doesn't match
  if (!all(res == terra::res(dem))) {
    message("Resampling DEM to target resolution of ", res, " m...")
    r <- dem
    terra::res(r) <- res
    dem <- terra::resample(dem, r)
  }
  return(dem)
}
