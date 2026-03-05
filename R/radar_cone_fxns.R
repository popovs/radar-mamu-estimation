#' Radar cone fxns
#' 
#' Functions to calculate mean flight headings and 
#' generate cones about the mean.
#' 


# HEADINGS ----------------------------------------------------------------


# CircStats::circ.mean
# Define circular mean function
# where x is your heading angle IN RADIANS
circmean <- function(x) { 
  sinr <- sum(sin(x))
  cosr <- sum(cos(x))
  circmean <- atan2(sinr, cosr)
  circmean
}

# Define circular bootstrapping function
# For a given list `x` of radian measurements,
# estimate the circular mean on a subsample of x.
# Repeat `n` number of times to derive the bootstrapped
# sample; derive quantile breaks `alpha` from the 
# bootstrapped sample and output the result.
# Thanks to Ben Bolker (again): https://stackoverflow.com/a/53916042/1454785
# Quantiles are appropriate for getting CIs for bootstrapped samples: https://stats.stackexchange.com/questions/515057/bootstrapping-mean-difference-standard-error-versus-quantiles 
circboot <- function(x, n, alpha) {
  bootsample <- replicate(n, circmean(sample(x, replace = TRUE)))
  out <- setNames(c(mean(bootsample), quantile(bootsample, 
                                               c(alpha/2, 1-alpha/2))), 
                  c("mean", "lower", "upper"))
  return(out)
}

calc_polar_mean <- function(headings, n_reps, alpha) {
  # For each station, take the POLAR MEAN INCOMING & OUTGOING heading.
  # We need to use POLAR/CIRCULAR stats to deal with 'wraparound' angles -
  # e.g. the mean of 45° and 315° is ~due north, but the cartesian mean 
  # of the two -- 180° -- is directly due south!
  
  # First convert degrees to radians
  headings$rad <- headings$heading * pi / 180
  
  # Calculate bootstrapped means and confidence intervals
  # for each station
  message("Calculating bootstrapped headings and 95 CI for bird headings at each station...")
  h <- aggregate(rad ~ site, headings, FUN = circboot, n = n_reps, alpha = alpha)
  message("Done bootstrapping.")
  
  # Tidy it up...
  h <- cbind(h[1], data.frame(h$rad))
  
  # Calculate cone width
  h$theta <- h$upper - h$lower
  
  # Convert radians back to degrees
  h <- h |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ .x * 180 / pi)) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ ifelse(. < 0, . + 360, .)))
  
  return(h)
}


# GENERATE CONES -----------------------------------------------------------


generate_cones <- function(h, stn, radius, res) {
  # Merge `h` and `stn`
  h <- merge(h, stn, by = "site")
  
  h <- sf::st_as_sf(h) |>
    sf::st_transform(3005)
  
  # Create cones from bootstrapped wedge values in `h`
  cones <- lapply(1:nrow(h), function(x){
    message("Generating cone for ", h$site[[x]], "...")
    radar_cone(pt = h[x,],
               radius = radius,
               theta = h$theta[[x]], # using 95CI heading boundaries
               heading = h$mean[[x]],
               res = res) 
  })
  
  # Check that all of them actually rendered correctly - for some
  # reason the above sometimes results in NULL rasters.... suspected 
  # memory issue.
  isNaN <- lapply(cones, terra::minmax)
  isNaN <- data.table::transpose(as.data.frame(isNaN))
  isNaN <- row.names(isNaN[is.na(isNaN$V1),]) # Extract all the cones with null values
  isNaN <- as.numeric(isNaN)
  
  while (length(isNaN) > 0) {
    # Re-run the radar_cone function for NaN rasters...
    for (i in isNaN) {
      message("Error with ", h$site[[i]], ". Re-generating cone for ", h$site[[i]], "...")
      cones[[i]] <- radar_cone(pt = h[i,],
                               radius = radius,
                               theta = h$theta[[i]],
                               heading = h$mean[[i]],
                               res = 25)
    }
    
    isNaN <- lapply(cones, terra::minmax)
    isNaN <- data.table::transpose(as.data.frame(isNaN))
    isNaN <- row.names(isNaN[is.na(isNaN$V1),]) # Extract all the cones with null values
    isNaN <- as.numeric(isNaN)
  }
  
  # Vectorize cones
  for (i in 1:nrow(h)) {
    message("Vectorizing ", h[i,][["site"]])
    x <- cones[[i]]
    x <- x == 1
    x <- terra::as.polygons(x, crs = "epsg:3005")
    x <- sf::st_as_sf(x)
    x <- x[x$lyr.1 == 1, ] # Keep only area == 1
    # Next draw a convex hull around the cone AND the station coordinate.
    # Some cones are so narrow the station gets dropped
    x <- sf::st_geometry(x) |> 
      sf::st_cast("MULTIPOINT") |>
      sf::st_union(h[i,]) |>
      sf::st_convex_hull() |>
      sf::st_as_sf()
    x$site <- h[i,][["site"]]
    cones[[i]] <- x
    rm(x)
  }
  
  cones <- dplyr::bind_rows(cones)
  
  cones$cone_area_ha <- units::set_units(sf::st_area(cones), "ha")
  
  sf::st_geometry(cones) <- "geometry"
  
  return(cones)
}
