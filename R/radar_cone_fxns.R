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


# Function to create a 'cone' of interpolated values around a given point.

#' Create a conical gradient of interpolated values around a point
#'
#' Using a given coordinate, this function will generate a circle
#' of a given radius around the point and fill the circle with a
#' gradient of values in a conical pattern. By default, the interpolated
#' values range from zero to one, with values of one point in the direction
#' of the specified heading.
#'
#' When using this function with spatial features,
#' ensure your data are projected in an equal-area coordinate system. This
#' function assumes meters as the default unit.
#'
#' @details This function was originally developed with the intention of
#' generating 'cones of probability' for estimating the dispersal
#' of birds flying across a landscape. Radar stations based across British
#' Columbia, Canada, collected data on the mean headings of birds flying
#' across a radar station. These headings indicate the more likely direction
#' birds will be flying to, while headings in the opposite direction
#' indicate areas the birds are flying away from. Thus, the gradient of
#' values from 0 (flying away) to 1 (flying towards) was developed around
#' each given radar station using this function.
#'
#' @param pt A single (X,Y) coordinate of class `sfc_POINT` (created by the `sf` package), where the CRS is a projected coordinate system.
#' @param radius Numeric. Distance in meters from the point `pt` to draw the circle radius.
#' @param heading Numeric. Heading, in compass degrees, which to interpolate values from.
#' @param theta Numeric (0 to 360). Angle, in degrees, indicating the amount of the circle that should be taken up by the `maxvalue`. For example, if you want to create a wedge where 1/3 of the circle == 1, set `theta = 120`. Default `0`.
#' @param res Numeric. Resolution, in meters, of the output raster. Default `500`.
#' @param mask Boolean (TRUE/FALSE). Should the output be masked to a circle shape (rather than square)? Default `TRUE`.
#' @param invert Boolean (TRUE/FALSE). Should the values of the output cone be inverted, such that the values in the direction of the heading are 0 rather than 1? Default `FALSE`.
#' @param maxvalue Numeric. Maximum value of the interpolated values of the cone. The interpolated values run from 0 to `maxvalue`. Default `1`.
#'
#' @return A `SpatRaster` object
#' @export
#'
#' @examples
#' library(magrittr)
#' pt <- sf::st_point(c(-126, 53)) %>%
#'   sf::st_sfc(crs = 4326) %>%
#'   sf::st_transform(crs = 3005)
#'
#' radius <- 30000 # 30 km
#' heading <- 45 # NE direction
#'
#' terra::plot(radar_cone(pt = pt,
#'                        radius = radius,
#'                        heading = heading))
#'
#' # Adjust heading
#' terra::plot(radar_cone(pt = pt,
#'                        radius = radius,
#'                        heading = 120))
#'
#' # Check values ≈ 1
#' terra::plot(radar_cone(pt = pt,
#'                        radius = radius,
#'                        heading = 120) > 0.99)
#'
#' # Increase wedge area where values ≈ 1 to 1/4 of the circle
#' terra::plot(radar_cone(pt = pt,
#'                        radius = radius,
#'                        heading = 120,
#'                        theta = 90) > 0.99)
#'
#' # Increase wedge area where values ≈ 1 to 1/3 of the circle
#' terra::plot(radar_cone(pt = pt,
#'                        radius = radius,
#'                        heading = 120,
#'                        theta = 120) > 0.99)
#'
#' Invert the cone
#' terra::plot(radar_cone(pt = pt,
#'                        radius = radius,
#'                        heading = 120,
#'                        theta = 90,
#'                        invert = TRUE))
#'
#' terra::plot(radar_cone(pt = pt,
#'                        radius = radius,
#'                        heading = 120,
#'                        theta = 120,
#'                        invert = TRUE))
radar_cone <- function(pt, radius, heading,
                       theta = 0, res = 500, mask = TRUE,
                       invert = FALSE, maxvalue = 1) {
  # Data health checks
  # Check that pt is a sf object
  stopifnot("`pt` must be a sf feature." = inherits(pt, c("sf", "sfc", "sfc_POINT")) == TRUE)
  # Check that pt is POINT geometry type
  stopifnot("`pt` must be POINT type geometry." = all(sf::st_geometry_type(pt) == "POINT"))
  # Check that radius > 0
  stopifnot("`radius` must have a positive value." = radius > 0)
  # Check that theta >= 0
  stopifnot("`theta` must have a positive value." = theta >= 0)
  # Check that theta <= 360
  stopifnot("`theta` must be between 0-360°." = theta <= 360)
  # Check that res > 0
  stopifnot("`res` must have a positive value." = res > 0)
  
  # Convert degrees to radians
  # This is sloppy trig, but it works, so......
  # 1) Multiply heading by -1 to ensure radians go clockwise (radians are measured counter-clockwise, by default)
  # 2) Multiply by pi/180 to convert to radians
  # 3) Subtract from 180° to rotate the angle (radians by default are measured from the x-axis,
  #     i.e. "east", rather than "north" on a compass-rose)
  heading <- (180-heading) * (-pi/180)
  p_buff <- sf::st_buffer(pt, radius)
  p_buff <- terra::vect(p_buff)
  
  ## Rasterize it to 500m resolution
  cone <- terra::rast(p_buff, resolution = res)
  xy <- suppressWarnings(terra::crds(cone))
  
  xy[,1] <- xy[,1] - sf::st_coordinates(pt)[1] # x col becomes x distance from origin
  xy[,2] <- xy[,2] - sf::st_coordinates(pt)[2] # y col becomes y distance from origin
  
  # Take the arctan of each point
  # Arctan: given a vector (x, y), returns angle `v`
  # from the x-axis in radians such that `-π ≤ θ ≤ π`
  v <- atan2(xy[,2], xy[,1])
  
  # Adjust values to fall along our heading of interest
  v2 <- sin(v + heading) * -1
  
  # Re-scale values from 0-1 instead of -1-1
  # If theta = 0 (default), the output produces a smooth
  # linear interpolation from 0-1 around the whole circle.
  # If theta > 0, a wedge with an angle of theta degrees
  # will be set to 1.
  theta <- theta * (pi/180) # Multiply by pi/180 to convert to radians
  v2 <- (v2 - min(v2)) / (cos(theta/2) - min(v2))
  v2[v2 > 1] <- 1
  
  # Fill in cone values
  cone[] <- v2
  
  # Invert such that 0 = heading and 1 = opposite heading
  if (invert) cone[] <- (cone - 1) * -1
  
  # Adjust scale that is output
  cone[] <- cone * maxvalue
  
  # Mask to circular output
  if (mask) cone <- terra::mask(cone, p_buff)
  
  return(cone)
}


# For supplied heading, station, radius, and resolution, create cones
generate_cones <- function(h, stn, radius, res) {
  # Merge `h` and `stn`
  h <- merge(h, stn, by = "site")
  
  h <- sf::st_as_sf(h)
  
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
