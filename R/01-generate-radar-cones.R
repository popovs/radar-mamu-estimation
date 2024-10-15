#' 01 GENERATE RADAR CONES

#' This script contains functions that will extract
#' radar station points, extract and calculate mean
#' bird flight headings, and combine that information
#' with watershed data and nest cost distance to
#' create MAMU population catchments. That is, if
#' a MAMU flies into a given watershed mouth, it is
#' assumed it is targeting a nest within the corresponding
#' catchment delineation.


# PREPARE SURVEY DATA -----------------------------------------------------

prepare_surveys <- function(path) {
  s <- MAMU::process_radar_data(path)
  s$site <- s$new_name # `process_radar_data()` stores cleaned up/consolidated site names in the `new_name` col
  # Create spatial object
  s <- s[!is.na(s$lon), ]
  s <- sf::st_as_sf(s, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  return(s)
}

# Extract individual stations from survey data
# bc station lat/lon can change by a few meters each year
# (impossible to set up the radar station in EXACTLY the
# same spot each year), take the mean lat/lon by
# site.
prepare_stn <- function(s, regions) {
  # Grab region from the `regions` sf
  # First drop existing regions col
  s <- dplyr::select(s, -region)
  s <- sf::st_transform(s, 3005)
  s <- suppressWarnings(sf::st_intersection(s, regions))
  stn <- sf::st_drop_geometry(s) # trying to get mean sf coords by group is a pain. 
  stn <- stn |> 
    dplyr::select(site, region, loc, lat, lon) |>
    dplyr::group_by(site, region, loc) |>
    dplyr::summarise(lat = mean(lat),
                     lon = mean(lon)) |>
    sf::st_as_sf(crs = 4326, 
                 coords = c("lon", "lat"), 
                 remove = FALSE)
  return(stn)
}



# GENERATE CONES -----------------------------------------------------------


prepare_headings <- function(path) {
  headings <- readxl::read_excel(path,
                                 na = c("", "NA", "#N/A", "N/A"),
                                 col_types = "text")
  headings <- janitor::clean_names(headings)
  headings <- headings[!is.na(headings$heading),]
  headings$lat <- as.numeric(headings$lat)
  headings$lon <- as.numeric(headings$lon)
  headings$speed <- as.numeric(headings$speed)
  headings$distance <- as.numeric(headings$distance)
  headings$heading <- as.numeric(headings$heading)
  
  # Flag if there's any likely data errors
  stopifnot("You have headings >360 that are likely typos." = !any(headings$heading > 360))
  
  # Clean up station names as needed
  stn_lookup <- MAMU::rs # rename station lookup table
  headings <- merge(headings, stn_lookup, by.x = "name", by.y = "original_name", all.x = TRUE)
  headings$name <- ifelse(is.na(headings$new_name), headings$name, headings$new_name)
  
  # TODO: Cut out stations with <5 headings?
  
  # Split out incoming and outgoing headings
  inc <- headings[grep("Incoming", headings$flightpath_type, ignore.case = T),c("name", "heading")]
  inc$flightpath_type <- "Incoming"
  out <- headings[grep("Outgoing", headings$flightpath_type, ignore.case = T),c("name", "heading")]
  out$flightpath_type <- "Outgoing"
  
  inc <- inc[complete.cases(inc),]
  out <- out[complete.cases(out),]
  
  # Make `out` the opposite heading
  # That is, we assume that 180° (polar opposite) direction
  # of the outgoing headings is equivalent to incoming ones
  out$heading <- out$heading - 180
  out$heading <- ifelse(out$heading < 0, out$heading + 360, out$heading)
  
  # Now merge in and out back together
  # `h` now functionally contains *'incoming' headings only*
  h <- rbind(inc, out)
  
  return(h)
}

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
  h <- aggregate(rad ~ name, headings, FUN = circboot, n = n_reps, alpha = alpha)
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


generate_cones <- function(h, stn, radius, res) {
  # Merge `h` and `stn`
  h <- merge(h, stn, by.x = "name", by.y = "site")
  names(h)[1] <- "site"
  
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
                               radius = 30000,
                               theta = 90,
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



# PLOTTING FXNS -----------------------------------------------------------

plot_headings <- function(site, headings, h) {
  # Subset to needed data
  headings <- headings[headings$name == site, ]
  h <- h[h$name == site, ]
  # Incoming headings plot
  p_inc <- ggplot2::ggplot(data = headings[headings$flightpath_type == "Incoming",]) + 
    ggplot2::geom_histogram(ggplot2::aes(x = heading)) + 
    ggplot2::geom_vline(xintercept = h[["mean"]],
                        color = "red") +
    ggplot2::geom_vline(xintercept = h[["lower"]],
                        color = "grey",
                        linetype = "dashed") +
    ggplot2::geom_vline(xintercept = h[["upper"]],
                        color = "grey",
                        linetype = "dashed") +
    ggplot2::xlim(0, 360) + 
    ggplot2::ggtitle("Incoming") + 
    ggplot2::theme(axis.title = ggplot2::element_blank()) +
    ggplot2::theme_minimal()
  # Outgoing headings plot
  p_out <- ggplot2::ggplot(data = headings[headings$flightpath_type == "Outgoing",]) + 
    ggplot2::geom_histogram(ggplot2::aes(x = heading)) + 
    ggplot2::geom_vline(xintercept = h[["mean"]],
                        color = "red") +
    ggplot2::geom_vline(xintercept = h[["lower"]],
                        color = "grey",
                        linetype = "dashed") +
    ggplot2::geom_vline(xintercept = h[["upper"]],
                        color = "grey",
                        linetype = "dashed") +
    ggplot2::xlim(0, 360) + 
    ggplot2::ggtitle("Outgoing") + 
    ggplot2::theme(axis.title = ggplot2::element_blank()) +
    ggplot2::theme_minimal()
  
  p <- ggpubr::ggarrange(p_inc, p_out, ncol = 1)
  print(p)
}


plot_watersheds <- function(site, watersheds, cones, stn) {
  # Subset to needed data
  watersheds <- watersheds[watersheds$site == site, ]
  cones <- cones[cones$site == site, ]
  stn <- stn[stn$site == site, ]
  # Plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = watersheds,
                     ggplot2::aes(color = keep_yn),
                     show.legend = FALSE) +
    ggplot2::geom_sf(data = cones,
                     fill = NA) +
    ggplot2::geom_sf(data = stn) +
    ggplot2::geom_text(data = watersheds,
                       ggplot2::aes(label = round(prct_coverage, 2),
                                    x = centroid_x,
                                    y = centroid_y),
                       size = 2) +
    ggplot2::ggtitle(site) +
    ggplot2::theme(axis.title = ggplot2::element_blank())
  print(p)
}


plot_cost <- function(site, cost_catchment, watersheds, 
                      cones, stn, nests, cost_cutoffs) {
  # Subset to needed data
  cost_catchment <- cost_catchment[[site]]
  watersheds <- watersheds[watersheds$site == site, ]
  cones <- cones[cones$site == site, ]
  stn <- stn[stn$site == site, ]
  nests <- suppressWarnings(sf::st_intersection(nests, watersheds))
  # Merge watersheds, stn, and cost_cutoffs to get all
  # attributes in one sf object
  watersheds <- merge(watersheds, sf::st_drop_geometry(stn[,c("site", "region")]))
  watersheds <- merge(watersheds, cost_cutoffs)
  # Plot
  p <- ggplot2::ggplot() + 
    tidyterra::geom_spatraster_contour_filled(data = cost_catchment,
                                              show.legend = FALSE) +
    tidyterra::geom_spatraster_contour(data = cost_catchment,
                                       breaks = c(#watersheds[["cost_min"]][watersheds$site == x],
                                         watersheds[["cost_max"]])) +
    tidyterra::scale_fill_whitebox_d() +
    ggplot2::geom_sf(data = cones,
                     fill = NA) +
    ggplot2::geom_sf(data = stn) +
    ggplot2::geom_sf(data = nests, 
                     color = "red", 
                     fill = "red") +
    ggplot2::ggtitle(site, subtitle = paste(watersheds[["region"]], "max cost =", watersheds[["cost_max"]])) +
    ggplot2::theme(axis.title = ggplot2::element_blank())
  print(p)
}




# WATERSHEDS --------------------------------------------------------------


select_watersheds <- function(watersheds, 
                              cones, 
                              min_cone_coverage = 0.02,
                              output_plots = TRUE,
                              output_dir,
                              ...) {
  # Get rid of tiny tiny watersheds + river deltas
  # These mess up the pathfinding algorithm when assuming 
  # MAMU can or cannot cross from one watershed to another
  # EXCEPT for Haida Gwaii. Nest data shows they'll nest in 
  # these tiny scraps of watersheds.
  watersheds$AREA_SQM <- sf::st_area(watersheds)
  watersheds <- watersheds[watersheds$AREA_SQM > units::as_units(150000, "m^2") | watersheds$region == "HG",]
  #watersheds <- watersheds[!is.na(watersheds$WTRSHDCD_2),]
  
  # Get a matrix of which watersheds intersect with cones
  # Each column is a cone, each row a watershed
  wat_cat <- sf::st_intersects(watersheds, cones, sparse = FALSE)
  
  # Loop through each column to extract out the watersheds that 
  # touch each cone
  # Replace the watersheds df with the subsetted list
  watersheds <- lapply(1:ncol(wat_cat), function(x){
    watersheds[wat_cat[,x],]
  })
  names(watersheds) <- cones$site
  
  # Bind it all into one df; keep the names as the id column
  watersheds <- dplyr::bind_rows(watersheds, .id = 'site')
  
  # Remove watersheds with <2% cone coverage -- EXCEPT
  # for Haida Gwaii!! HG contains many, many small 
  # cliffside catchments that empty directly into the
  # sea -- and nesting data shows they will readily
  # nest in these tiny catchments, unlike the mainland. 
  # So with HG just choose everything.
  
  # I.e., remove any watershed slivers
  intersect_area <- sf::st_intersection(cones, watersheds)
  intersect_area$intersect_area <- units::set_units(sf::st_area(intersect_area), "ha")
  intersect_area <- intersect_area |>
    dplyr::filter(site == site.1) |> # we don't care about cone slivers intersecting with other site watersheds
    dplyr::select(site, WSD_ID, intersect_area) |>
    sf::st_drop_geometry() |>
    dplyr::group_by(site, WSD_ID) |>
    dplyr::summarize(intersect_area = sum(intersect_area), .groups = "keep") |>
    dplyr::group_by(site) |>
    dplyr::mutate(total_area = sum(intersect_area),
                  prct_coverage = intersect_area / total_area * 100)
  
  #watersheds <- merge(watersheds, st_drop_geometry(cones), by = "site") # merge cone area in
  watersheds <- merge(watersheds, intersect_area, by = c("site", "WSD_ID"))
  
  #watersheds$prct_cone_overlap <- (watersheds$intersect_area / watersheds$cone_area_ha) * 100
  #watersheds$prct_cone_overlap <- units::drop_units(watersheds$prct_cone_overlap)
  watersheds$prct_coverage <- units::drop_units(watersheds$prct_coverage)
  
  watersheds$centroid_x <- sf::st_coordinates(sf::st_centroid(watersheds))[,1]
  watersheds$centroid_y <- sf::st_coordinates(sf::st_centroid(watersheds))[,2]
  
  watersheds$keep_yn <- (watersheds$prct_coverage > (100 * min_cone_coverage) | watersheds$region == "HG")
  
  # Save plots to inspect each one
  # TODO: ideally split this out and make it its own 
  # target that uses the selected watersheds output. It's time
  # consuming to recreate all the selected watersheds just because
  # of a plotting error.
  if (output_plots == TRUE) {
    # Unpack dots
    plot_data <- list(...)
    stn <- plot_data$stn
    headings <- plot_data$headings
    h <- plot_data$h
    # Create plots
    dir.create(output_dir, showWarnings = F)
    sites <- unique(cones$site)
    for (i in sites) {
      out_path <- file.path(output_dir,
                            fs::path_sanitize(paste0(i, ".png"), replacement = "-"))
      message("Saving plot for ", i, " to '", out_path, "'")
      p1 <- plot_watersheds(site = i, 
                            watersheds = watersheds,
                            cones = cones,
                            stn = stn)
      p2 <- plot_headings(site = i, headings = headings, h = h)
      p_all <- ggpubr::ggarrange(p1, p2, nrow = 1)
      ggplot2::ggsave(filename = out_path, plot = p_all)
    }
  }
  
  # Drop anything that failed keep_yn tests
  watersheds <- watersheds[watersheds$keep_yn,]
  watersheds <- dplyr::select(watersheds, -keep_yn)
  
  return(watersheds)
  
}



# COST DISTANCE WITHIN CATCHMENTS -----------------------------------------

watershed_cost <- function(watersheds, 
                           dem, 
                           cones, 
                           stn, 
                           cost_cutoffs,
                           output_plots = TRUE,
                           output_dir,
                           ...) {
  # Dissolve watersheds together by site
  watersheds <- watersheds |>
    dplyr::select(site) |>
    dplyr::group_by(site) |>
    dplyr::summarize()
  
  # Reproject stn
  stn <- sf::st_transform(stn, 3005)
  
  sites <- unique(watersheds$site)
  # Make 'cost catchments'
  catchments <- lapply(sites, function(x) {
    message("Creating catchment for ", x)
    tmp <- watersheds[watersheds$site == x,]
    tmp <- terra::crop(dem, tmp, mask = TRUE)
    tmp2 <- terra::crop(dem, cones[cones$site == x, ], mask = TRUE) # also extract anything in the cone path - e.g. water - so birds can cross over water areas in the cone's path
    tmp <- terra::merge(tmp, tmp2)
    # Choose the point closest on the raster to the radar
    # station as the origin point
    # Extract origin coordinate
    i <- stn[stn$site == x,]
    j <- sf::st_union(sf::st_as_sf(terra::as.polygons(tmp)))
    origin <- sf::st_nearest_points(i, j)
    origin <- sf::st_coordinates(origin)
    origin <- origin[,1:2] # drop 'L1' column
    i <- sf::st_coordinates(i)
    origin <- rbind(origin, i)
    # drop the duplicated rows - that's our point `stn`, and we only
    # care about origin point in/on the polygon
    if (nrow(unique(origin)) == 1) { # vanilla `ifelse` doesn't play nicely with matrices
      origin <- unique(origin)
    } else if (nrow(unique(origin == 2))) {
      origin <- origin[!(duplicated(origin) | duplicated(origin, fromLast = TRUE)), ]
    } else {
      message("Something weird going on with ", x)
    }
    origin <- matrix(origin, ncol = 2)
    # Assign -1 to origin
    if (is.na(terra::cellFromXY(tmp, origin))) message("Unable to find origin for ", x)
    tmp[terra::cellFromXY(tmp, origin)] <- -1
    tmp <- terra::costDist(tmp, -1, scale = 1000, maxiter = 100)
    tmp <- terra::ifel(tmp > 0, tmp, NA)
    tmp <- terra::crop(tmp, watersheds[watersheds$site == x,], mask = TRUE) # now crop to watersheds shape
    return(tmp)
  })
  
  names(catchments) <- sites
  
  # Save plots to inspect...
  # TODO: ideally split this out and make it its own 
  # target that uses the catchments output. It's time
  # consuming to recreate all the catchments just because
  # of a plotting error.
  if (output_plots == TRUE) {
    # Unpack dots
    plot_data <- list(...)
    headings <- plot_data$headings
    h <- plot_data$h
    nests <- plot_data$nests
    # Create plots
    dir.create(output_dir, showWarnings = F)
    for (i in sites) {
      out_path <- file.path(output_dir,
                            fs::path_sanitize(paste0(i, ".png"), replacement = "-"))
      message("Saving plot for ", i, " to '", out_path, "'")
      p1 <- plot_cost(site = i, 
                      cost_catchment = catchments,
                      watersheds = watersheds,
                      cones = cones,
                      stn = stn, 
                      nests = nests,
                      cost_cutoffs = cost_cutoffs)
      p2 <- plot_headings(site = i, headings = headings, h = h)
      p_all <- ggpubr::ggarrange(p1, p2, nrow = 1)
      ggplot2::ggsave(filename = out_path, plot = p_all)
    }
  }
  
  # Merge watersheds, stn, and cost_cutoffs to get all
  # attributes in one sf object
  watersheds <- merge(watersheds, sf::st_drop_geometry(stn[,c("site", "region")]))
  watersheds <- merge(watersheds, cost_cutoffs)
  # Apply the regional cost cutoff to each individual
  # cost catchment (i.e., draw boundaries on each cost
  # catchment)
  catchments <- lapply(sites, function(x) {
    message("Setting max boundaries of cost catchment for ", x)
    # Extract catchment + apply cutoff
    tmp <- catchments[[x]]
    cutoff <- watersheds[["cost_max"]][watersheds$site == x]
    tmp <- terra::ifel(tmp > cutoff, NA, 1)
    # Vectorize
    tmp <- terra::as.polygons(tmp, crs = "epsg:3005")
    tmp <- sf::st_as_sf(tmp)
  })
  
  # Clean up output
  catchments <- dplyr::bind_rows(catchments)
  catchments$site <- sites
  catchments <- catchments[, "site"]
  
  return(catchments)
  
}

