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


# REGIONAL ELEVATION CUTOFFS ----------------------------------------------

crop_dem <- function(dem, region_name, regions) {
  message("Cropping and saving ", region_name)
  regions <- regions[regions$region == region_name, ]
  tmp <- terra::crop(dem, regions)
  tmp <- terra::mask(tmp, regions) # This can be done in one step with terra::crop(mask = T), but was resulting in buggy raster values - so doing it in two steps here.
  return(tmp)
}


nest_elev_quantile <- function(nests, elev_data, quantiles = c(0.025, 0.975)) {
  # TODO: if this overwrites the data, it might trigger an endless pipeline reassessment loop
  # TODO: regions within nest data might not necessarily line up with regions_map in pipeline
  nests$elev_m <- elev_data
  # Calculate quantiles
  quants <- aggregate(elev_m ~ region, nests, FUN = quantile, quantiles[1], na.rm = TRUE)
  quants[3] <- aggregate(elev_m ~ region, nests, FUN = quantile, quantiles[2], na.rm = TRUE)[2]
  names(quants) <- c("region", "elev_m_min", "elev_m_max")
  # Round to nearest 10 
  quants[,2:3] <- round(quants[,2:3], -1)
  # Fill in missing values, if they're missing
  # If NVI is NULL, use mean of the other 3 regions of VI
  # If AKB and NC are NULL, use CC cutoffs
  if (!"NVI" %in% quants$region) quants <- rbind(quants, c("NVI", round(mean(quants[grep("VI", quants$region), "elev_m_min"])), round(mean(quants[grep("VI", quants$region), "elev_m_max"]))))
  if (!"AKB" %in% quants$region) quants <- rbind(quants, c("AKB", quants[["elev_m_min"]][quants$region == "CC"], quants[["elev_m_max"]][quants$region == "CC"]))
  if (!"NC" %in% quants$region)  quants <- rbind(quants, c("NC", quants[["elev_m_min"]][quants$region == "CC"], quants[["elev_m_max"]][quants$region == "CC"]))
  # rbind converts numerics to character... convert them back
  quants$elev_m_min <- as.numeric(quants$elev_m_min)
  quants$elev_m_max <- as.numeric(quants$elev_m_max)
  return(quants)
}

reclass_elevation <- function(dem, elev_cutoffs, region, ...) { # ... param to ignore other cols in the dataframe when it gets passed in
  elev_cutoffs <- elev_cutoffs[elev_cutoffs$region == region,]
  message("Reclassifying ", elev_cutoffs$region, "...")
  elev_min <- elev_cutoffs$elev_m_min
  elev_max <- elev_cutoffs$elev_m_max
  dem <- terra::ifel(dem > 0, dem, NA)
  dem <- terra::ifel(dem < elev_max, dem, NA)
  dem <- terra::ifel(dem < elev_min, 1, 2)
  return(dem)
}


# DISTANCE FROM COAST CUTOFF ----------------------------------------------

get_usa_land <- function(extent) {
  usa <- rnaturalearth::ne_states("united states of america")
  usa <- usa[usa$name %in% c("Alaska", "Washington"), ]
  usa <- sf::st_as_sf(usa)
  usa <- sf::st_transform(usa, 3005)
  usa <- sf::st_crop(usa, extent)
  usa <- terra::vect(usa)
  return(usa)
}

get_bc_land <- function(extent) {
  # Rather than block off inland BC with a rough BC-shaped
  # polygon from Natural Earth, we're just going to use crude
  # rectangles. The reason for this is the DEM coastline is far
  # more detailed than the course NE polygon, and we would risk
  # losing coastline data if we masked the DEM with such a 
  # course polygon. So rectangles it is.
  # Specify a few rectangles w coordinates from bottom left corner clockwise to bottom right corner
  x1 <- rbind(c(1280841, extent[3]), c(1280841, extent[4]), c(extent[2], extent[4]), c(extent[2], extent[3]))
  x2 <- rbind(c(1109405, 651763), c(1109405, extent[4]), c(extent[2], extent[4]), c(extent[2], 651763))
  x3 <- rbind(c(979227, 877974), c(979227, extent[4]), c(extent[2], extent[4]), c(extent[2], 877974))
  x4 <- rbind(c(827708, 1154691), c(827708, extent[4]), c(extent[2], extent[4]), c(extent[2], 1154691))
  x5 <- rbind(c(610000, 1300000), c(610000, extent[4]), c(extent[2], extent[4]), c(extent[2], 1300000))
  # lol this is terribly inefficient but it works
  land <- list(x1, x2, x3, x4, x5)
  land <- lapply(land, terra::vect, type = "polygons", crs = "epsg:3005")
  land1 <- terra::union(land[[1]], land[[2]])
  land2 <- terra::union(land1, land[[3]])
  land3 <- terra::union(land2, land[[4]])
  land4 <- terra::union(land3, land[[5]])
  land <- terra::aggregate(land4)
  return(land)
}

merge_land <- function(usa_land, bc_land) {
  land <- terra::aggregate(terra::union(usa_land, bc_land))
  return(land)
}
