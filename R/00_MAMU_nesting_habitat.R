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


nest_quantiles <- function(nests, quant_data, prefix, quantiles = c(0.025, 0.975)) {
  # TODO: if this overwrites the data, it might trigger an endless pipeline reassessment loop
  # TODO: regions within nest data might not necessarily line up with regions_map in pipeline
  nests$quant_col <- quant_data
  # Calculate quantiles
  quants <- aggregate(quant_col ~ region, nests, FUN = quantile, quantiles[1], na.rm = TRUE)
  quants[3] <- aggregate(quant_col ~ region, nests, FUN = quantile, quantiles[2], na.rm = TRUE)[2]
  # Set up output names
  min_col <- paste0(prefix, "_min")
  max_col <- paste0(prefix, "_max")
  names(quants) <- c("region", min_col, max_col)
  # Round to nearest 10 
  quants[,2:3] <- round(quants[,2:3], -1)
  # Fill in missing values, if they're missing
  # If NVI is NULL, use mean of the other 3 regions of VI
  # If AKB and NC are NULL, use CC cutoffs
  if (!"NVI" %in% quants$region) quants <- rbind(quants, c("NVI", round(mean(quants[grep("VI", quants$region), min_col])), round(mean(quants[grep("VI", quants$region), max_col]))))
  if (!"AKB" %in% quants$region) quants <- rbind(quants, c("AKB", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  if (!"NC" %in% quants$region)  quants <- rbind(quants, c("NC", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  # rbind converts numerics to character... convert them back
  quants[[min_col]] <- as.numeric(quants[[min_col]])
  quants[[max_col]] <- as.numeric(quants[[max_col]])
  return(quants)
}


# Helper functions to extract whisker boxplot values
# Anything above that is plotted as outlier dots
upper_whisker <- function(x) {
  x <- x[!is.na(x)]
  min(max(x), as.numeric(quantile(x, 0.75)) + (IQR(x) * 1.5))
}

lower_whisker <- function(x) {
  x <- x[!is.na(x)]
  max(min(x), as.numeric(quantile(x, 0.25)) - (IQR(x) * 1.5))
}

# This is an alternative function to nest_quantiles - instead of taking
# the 95% (or whatever) quantiles of X nest data, this function instead
# finds the min and max of each data group and removes any outliers
# from the interquartile range (i.e., boxplot outliers).
nest_minmax_sans_outliers <- function(nests, quant_data, prefix) {
  # TODO: if this overwrites the data, it might trigger an endless pipeline reassessment loop
  # TODO: regions within nest data might not necessarily line up with regions_map in pipeline
  nests$quant_col <- quant_data
  
  # Remove overall nest quant outliers first
  upper_outliers <- upper_whisker(quant_data)
  nests$quant_outlier_yn <- nests$quant_col >= upper_outliers
  # IQR method
  #upper_outliers <- boxplot.stats(nests$quant_col)$out[boxplot.stats(nests$quant_col)$out > median(nests$quant_col, na.rm = TRUE)]
  # Quantile method
  #upper_outliers <- quantile(nests[["quant_col"]], 0.95, na.rm = TRUE)
  #nests$quant_outlier_yn <- nests$quant_col >= min(upper_outliers)
  
  # Pull minimum value by group after overall quant outliers are cut out
  quants <- aggregate(quant_col ~ region, nests, FUN = min, na.rm = TRUE)
  # Pull maximum value by group after overall quant outliers are cut out
  quants$max <- aggregate(quant_col ~ region, 
                          nests[nests$quant_outlier_yn == FALSE,], 
                          FUN = max, 
                          na.rm = TRUE)[[2]]

  # Set up output names
  min_col <- paste0(prefix, "_min")
  max_col <- paste0(prefix, "_max")
  names(quants) <- c("region", min_col, max_col)
  
  # Round to nearest 10
  quants$cost_min <- floor(quants$cost_min / 10) * 10 # round DOWN to nearest 10
  quants$cost_max <- ceiling(quants$cost_max / 10) * 10 # round UP to nearest 10
  
  # Fill in missing values, if they're missing
  # If NVI is NULL, use mean of the other 3 regions of VI
  # If AKB and NC are NULL, use CC cutoffs
  if (!"NVI" %in% quants$region) quants <- rbind(quants, c("NVI", round(mean(quants[grep("VI", quants$region), min_col])), round(mean(quants[grep("VI", quants$region), max_col]))))
  if (!"AKB" %in% quants$region) quants <- rbind(quants, c("AKB", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  if (!"NC" %in% quants$region)  quants <- rbind(quants, c("NC", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  # rbind converts numerics to character... convert them back
  quants[[min_col]] <- as.numeric(quants[[min_col]])
  quants[[max_col]] <- as.numeric(quants[[max_col]])
  return(quants)
  
}

# Finally, this is what seems to be the best approach - the 8 outliers
# across the entire dataset - nests with costs > 7k - are cut out; 
# then rather than taking the maximum per region, the IQR is taken
# per region (ie, cutting out regional outliers). This approach is 
# defendable because our regional groupings are fairly arbitrary
# and don't necessarily reflect actual regional behavioral patterns
# in MAMU activity. So, there may be a few nests in one region (e.g.,
# MWVI) that could reasonably just as easily fall within another region
# (e.g., NVI). 
nest_iqr_sans_outliers <- function(nests, quant_data, prefix, hg_exception) {
  # TODO: if this overwrites the data, it might trigger an endless pipeline reassessment loop
  # TODO: regions within nest data might not necessarily line up with regions_map in pipeline
  nests$quant_col <- quant_data
  
  # Remove overall nest quant outliers first
  upper_outliers <- upper_whisker(quant_data)
  nests$quant_outlier_yn <- nests$quant_col >= upper_outliers
  
  # Pull minimum value by group after overall quant outliers are cut out
  quants <- aggregate(quant_col ~ region, nests, FUN = min, na.rm = TRUE)
  # Now check if any individual regions have IQR outliers, after overall quant outliers are cut out
  quants$upr_whisker <- aggregate(quant_col ~ region, 
                                  nests[nests$quant_outlier_yn == FALSE,], 
                                  FUN = upper_whisker)[[2]]
  
  # Set up output names
  min_col <- paste0(prefix, "_min")
  max_col <- paste0(prefix, "_max")
  names(quants) <- c("region", min_col, max_col)
  
  # Round to nearest 10
  quants$cost_min <- floor(quants$cost_min / 10) * 10 # round DOWN to nearest 10
  quants$cost_max <- ceiling(quants$cost_max / 10) * 10 # round UP to nearest 10
  
  # HG has so little data and birds tend to fly differently there,
  # as there aren't distinctly well defined watersheds in the same way 
  # as the mainland or VI, which have more rugose coastlines. 
  # So, as an option, simply take the maximum value for HG rather than
  # cutting out any outliers. 
  if (hg_exception) {
    hg_max <- max(nests[["quant_col"]][nests$quant_outlier_yn == FALSE & nests$region == "HG"])
    hg_max <- ceiling(hg_max / 10) * 10
    quants[["cost_max"]][quants$region == "HG"] <- hg_max
  }
  
  # Fill in missing values, if they're missing
  # If NVI is NULL, use mean of the other 3 regions of VI
  # If AKB and NC are NULL, use CC cutoffs
  if (!"NVI" %in% quants$region) quants <- rbind(quants, c("NVI", round(mean(quants[grep("VI", quants$region), min_col])), round(mean(quants[grep("VI", quants$region), max_col]))))
  if (!"AKB" %in% quants$region) quants <- rbind(quants, c("AKB", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  if (!"NC" %in% quants$region)  quants <- rbind(quants, c("NC", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  
  # rbind converts numerics to character... convert them back
  quants[[min_col]] <- as.numeric(quants[[min_col]])
  quants[[max_col]] <- as.numeric(quants[[max_col]])
  
  return(quants)
}


nest_isoforest <- function(nests, quant_data, prefix) {
  nests$quant_col <- quant_data
  # Create isolation forest model and calcuate outlier scores
  data <- sf::st_drop_geometry(nests[,c("region", "quant_col")])
  if_model <- isotree::isolation.forest(data, sample_size = 3, ndim=1, ntrees=10, nthreads=1)
  scores <- isotree::predict.isolation_forest(if_model, data, type="avg_depth")
  # Choose cutoff for what counts as an 'outlier score'. 
  # Choose scores that 99% of the data fall into
  outlier_threshold <- quantile(scores, 0.01)[[1]]
  # Add that info back to `nests`
  nests$scores <- scores
  nests$outlier_yn <- nests$scores <= outlier_threshold # less than or = to threshold
  # Choose max by group after cutting out regional outliers
  nests <- nests[nests$outlier_yn == FALSE, ]
  quants <- aggregate(quant_col ~ region, nests, FUN = "max")
  # Set minimum to be zero with this method
  quants$min <- 0
  
  # Reorder!
  # NOTE this order is switched from other methods above!
  quants <- quants[,c(1,3,2)]
  
  # Set up output names
  min_col <- paste0(prefix, "_min")
  max_col <- paste0(prefix, "_max")
  names(quants) <- c("region", min_col, max_col) 
  
  # Round to nearest 10
  quants[[min_col]] <- floor(quants[[min_col]] / 10) * 10 # round DOWN to nearest 10
  quants[[max_col]] <- ceiling(quants[[max_col]] / 10) * 10 # round UP to nearest 10
  
  # Fill in missing values, if they're missing
  # If NVI is NULL, use mean of the other 3 regions of VI
  # If AKB and NC are NULL, use CC cutoffs
  if (!"NVI" %in% quants$region) quants <- rbind(quants, c("NVI", round(mean(quants[grep("VI", quants$region), min_col])), round(mean(quants[grep("VI", quants$region), max_col]))))
  if (!"AKB" %in% quants$region) quants <- rbind(quants, c("AKB", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  if (!"NC" %in% quants$region)  quants <- rbind(quants, c("NC", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  
  # rbind converts numerics to character... convert them back
  quants[[min_col]] <- as.numeric(quants[[min_col]])
  quants[[max_col]] <- as.numeric(quants[[max_col]])
  
  return(quants)
}

reclass_raster <- function(dem, cutoffs, min_col, max_col, region, ...) { # ... param to ignore other cols in the dataframe when it gets passed in
  cutoffs <- cutoffs[cutoffs$region == region,]
  message("Reclassifying ", cutoffs$region, "...")
  min <- cutoffs[[min_col]]
  max <- cutoffs[[max_col]]
  dem <- terra::ifel(dem > 0, dem, NA)
  dem <- terra::ifel(dem < max, dem, NA)
  dem <- terra::ifel(dem < min, 1, 2)
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
  
coast_distance <- function(elev, sea, dist_km, exclude_islands, ...) {
  if (terra::ext(elev) != terra::ext(sea)) elev <- terra::resample(elev, sea)
  
  # OPTIONAL: exclude HG and VI from coast distance calculation?
  # We have evidence birds fly and nest within all regions of VI.
  if (exclude_islands) {
    # Unpack dots
    dots <- list(...)
    regions <- dots$region
    # Prep `regions` to create two masking blocks for HG and VI
    islands <- dplyr::summarise(regions[regions$region %in% c("HG", "NVI", "MWVI", "EVI", "SWVI", "WNVI"),]) |> nngeo::st_remove_holes()
    # Mask and save islands as raster == 1
    islands <- terra::crop(elev, islands, mask = TRUE)
    islands <- terra::ifel(islands > 0, 1, NA)
  } 
  
  coast_distance <- terra::merge(elev, sea)
  coast_distance <- terra::ifel(coast_distance > 0 , 1, coast_distance) # we don't care about the elevation cutoff minimums here - so reclassify them to only keep the max cutoffs
  coast_distance <- terra::gridDist(coast_distance, target = 0) # calculate distance from the sea! (values of 0 on the raster == sea)
  
  # Cut down to only include 30km distance and clip to land areas
  dist_m <- dist_km * 1000
  c_dist <- coast_distance <= dist_m
  # Merge w sea and/or islands
  if (exclude_islands) {
    c_dist <- terra::merge(islands, c_dist)
  } else {
    c_dist <- terra::merge(sea, c_dist)
  }
  c_dist <- terra::ifel(c_dist == 1, 1, NA) # set any non-valid nesting areas == NA
  
  return(c_dist)
}


# NEST COST DISTANCE ------------------------------------------------------

cost_distance <- function(dem) {
  ## Prepare cost layer ##
  
  # Extract res of target
  # Going to assume all data products are square res... 
  # so just take the first number and assume it's equal to the second
  res <- terra::res(dem)[1]
  
  # First, similar to above, lower the resolution of high res
  # DEM by a factor of 10, or else the calculations will fail.
  if (res < 100) {
    c <- dem
    terra::res(c) <- res * 10
    c <- terra::resample(dem, c)
  } else {
    c <- dem
  }
  
  # Ensure no negative values in the raster, or the costDist
  # function will fail.
  c <- terra::ifel(c < 0, 0, c)
  
  ## Calculate cost ##
  # The cost distance function calculates distance from shore (i.e.,
  # the `gridDist` function we just used above) and multiplies it
  # by the 'cost' layer (elevation). Higher elevations are more
  # costly to fly over. 
  cost <- terra::costDist(c, target = 0, scale = 1000) # divide values by 1000 so output numbers are smaller
  return(cost)
}



# FOREST COVER CUTOFF -----------------------------------------------------

download_forest_tiles <- function(url, output_dir) {
  # Create output dir
  output_dir <- file.path(output_dir)
  dir.create(output_dir, showWarnings = FALSE)
  output_file <- file.path(output_dir, basename(url))
  # These files will take longer than one minute to download, so
  # increase the download timeout to 5 mins
  options(timeout = 300)
  # Download % forest cover data (year 2000) for BC
  utils::download.file(url = url, destfile = output_file)
  return(output_file)
}

