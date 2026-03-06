#' 01 GENERATE RADAR CATCHMENTS

#' This script contains functions that will extract
#' radar station points, extract and calculate mean
#' bird flight headings, and combine that information
#' with watershed data and nest cost distance to
#' create MAMU population catchments. That is, if
#' a MAMU flies into a given watershed mouth, it is
#' assumed it is targeting a nest within the corresponding
#' catchment delineation.






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
    ggplot2::ggtitle("180° - Outgoing") + 
    ggplot2::theme(axis.title = ggplot2::element_blank()) +
    ggplot2::theme_minimal()
  
  p <- ggpubr::ggarrange(p_inc, p_out, ncol = 1)
  return(p)
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
  return(p)
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
  return(p)
}


plot_catchment <- function(site, cost_catchment, accessible_catchment,
                           watersheds, cones, stn, nests) {
  # Subset to needed data
  cost_catchment <- cost_catchment[cost_catchment$site == site, ]
  accessible_catchment <- accessible_catchment[accessible_catchment$site == site, ]
  watersheds <- watersheds[watersheds$site == site, ]
  cones <- cones[cones$site == site, ]
  stn <- stn[stn$site == site, ]
  nests <- suppressWarnings(sf::st_intersection(nests, watersheds))
  # Plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = watersheds,
                     show.legend = FALSE) +
    ggplot2::geom_sf(data = cost_catchment,
                     fill = "grey",
                     color = NA,
                     alpha = 0.7,
                     show.legend = FALSE) +
    ggplot2::geom_sf(data = accessible_catchment,
                     color = NA, 
                     fill = "#26D1EA",
                     alpha = 0.3,
                     show.legend = FALSE) +
    ggplot2::geom_sf(data = cones,
                     fill = NA) +
    ggplot2::geom_sf(data = stn) +
    ggplot2::geom_sf(data = nests, 
                     color = "red", 
                     fill = "red") +
    ggplot2::ggtitle(site, subtitle = stn[["loc"]][stn$site == site])
  return(p)
}

plot_final_map <- function(site, 
                           catchments,
                           watersheds,
                           suitable_habitat,
                           cones, 
                           stn, 
                           nests,
                           apikey) {
  x <- site
  
  # Pull catchment bbox and add 5km plotting buffer
  bbox <- sf::st_bbox(catchments[catchments$site == x, ])
  bbox[1:2] <- bbox[1:2] - 5000 # add 5km buffer for visualizing
  bbox[3:4] <- bbox[3:4] + 5000
  
  # Prepare data for mapping
  # stn
  p_stn <- sf::st_transform(stn, 3005)
  p_stn <- cbind(p_stn, sf::st_coordinates(p_stn))
  p_stn <- sf::st_crop(p_stn, bbox)
  
  # watersheds
  watersheds <- sf::st_collection_extract(watersheds, "POLYGON")
  watersheds <- sf::st_crop(watersheds, bbox)
  
  # suitable habitat
  if (inherits(suitable_habitat, "sf")) {
    suitable_habitat <- sf::st_collection_extract(suitable_habitat, "POLYGON")
    suitable_habitat <- sf::st_crop(suitable_habitat, bbox)
  } else if (inherits(suitable_habitat, "SpatRaster")) {
    suitable_habitat <- terra::crop(suitable_habitat, bbox)
    suitable_habitat <- suitable_habitat |> 
      terra::as.polygons() |> 
      sf::st_as_sf()
  }
  
  
  # nests
  nests <- sf::st_centroid(nests)
  
  # Pull maptile for the site
  jawg_terrain <- maptiles::create_provider(
    name = "Jawg.Terrain",
    url = "https://tile.jawg.io/jawg-terrain/{z}/{x}/{y}.png?access-token={apikey}",
    citation = "© Jawg Maps"
  )
  
  tile <- maptiles::get_tiles(x = bbox,
                              provider = jawg_terrain,
                              apikey = apikey,
                              zoom = 12,
                              crop = TRUE)
  
  p <- ggplot2::ggplot() +
    tidyterra::geom_spatraster_rgb(data = tile) +
    ggplot2::geom_sf(data = watersheds,
                     color = "red",
                     fill = NA,
                     show.legend = FALSE) +
    ggplot2::geom_sf(data = suitable_habitat,
                     color = "#E92162",#"#B41347",
                     fill = "#FF296E",
                     linewidth = 0.05,
                     alpha = 0.3) + 
    ggplot2::geom_sf(data = catchments[catchments$site == x,],
                     #ggplot2::aes(fill = site),
                     #fill = "#C6E921",
                     #color = "#A1BE19",
                     lwd = 0.1,
                     fill = "#26D1EA",
                     color = "#0996AB",
                     alpha = 0.5,
                     show.legend = FALSE) +
    ggplot2::geom_sf(data = nests,
                     color = "#333333",
                     shape = 18) +
    ggrepel::geom_text_repel(data = p_stn,
                             ggplot2::aes(label = site,
                                          x = X,
                                          y = Y),
                             size = 3,
                             nudge_x = 0,
                             nudge_y = 250,
                             color = "black",
                             bg.color = "white",
                             bg.r = 0.15) +
    ggplot2::geom_sf(data = cones[cones$site == x, ],
                     fill = "orange",
                     color = "#D85426",
                     alpha = 0.15) +
    ggplot2::geom_sf(data = p_stn,
                     color = "#222222") +
    ggplot2::coord_sf(xlim = c(bbox[1], bbox[3]),
                      ylim = c(bbox[2], bbox[4]),
                      expand = FALSE) +
    ggplot2::ggtitle(x) +
    ggplot2::theme(axis.title = ggplot2::element_blank())
  
    return(p)
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
  
  watersheds <- merge(watersheds, intersect_area, by = c("site", "WSD_ID"))
  
  watersheds$prct_coverage <- units::drop_units(watersheds$prct_coverage)
  
  watersheds$centroid_x <- sf::st_coordinates(sf::st_centroid(watersheds))[,1]
  watersheds$centroid_y <- sf::st_coordinates(sf::st_centroid(watersheds))[,2]
  
  #watersheds$keep_yn <- (watersheds$prct_coverage > (100 * min_cone_coverage) | watersheds$region == "HG")
  watersheds$keep_yn <- (watersheds$prct_coverage > (100 * min_cone_coverage))
  
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
                           cost,
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
    # IMPORTANT! Our cut distance cutoffs assume the origin is from
    # a point at sea. For inland stations, we need to add the base
    # cost of *how much it costs to fly further from that station.*
    # Therefore, we need to load up the cost of the station and add 
    # that value to the `tmp` cost raster.
    # if (stn[stn$site == x,]$loc == "Inland") {
    #   stn_cost <- terra::extract(cost, origin)[[1]]
    #   tmp <- tmp + stn_cost
    # }
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
  catchments2 <- lapply(sites, function(x) {
    message("Setting max boundaries of cost catchment for ", x)
    # Extract catchment + apply cutoff
    tmp <- catchments[[x]]
    cutoff <- watersheds[["cost_max"]][watersheds$site == x]
    tmp <- terra::ifel(tmp > cutoff, NA, 1)
    # Vectorize
    tmp <- terra::as.polygons(tmp, crs = "epsg:3005")
    tmp <- sf::st_as_sf(tmp)
    if(nrow(tmp) > 0) tmp$site <- x # after applying inland cutoffs, some sites might have NULL catchments!
    return(tmp)
  })
  
  # Clean up output
  catchments2 <- dplyr::bind_rows(catchments2)
  # catchments2$site <- sites
  catchments2 <- catchments2[, "site"]
  
  return(catchments2)
  
}



# DIRECTIONALITY CUTOFF ---------------------------------------------------

# We've extracted the watershed regions that contain birds,
# but a few of them could be whittled down further - e.g.
# see Brittain or Kwinamass. Birds won't be flying *backwards*
# into our catchment areas. 

# Let's assume birds won't fly backwards. 
# Cut away any catchment area directly behind 
# the flight heading.

directionality_crop <- function(cost_catchments,
                                stn,
                                h, 
                                watersheds, 
                                cones, 
                                res) {
  # Dissolve watersheds together by site
  watersheds <- watersheds |>
    dplyr::select(site) |>
    dplyr::group_by(site) |>
    dplyr::summarize()
  
  # Merge `h` and `stn`
  stn <- merge(stn, h, by = "site")
  stn <- sf::st_transform(stn, 3005)
  
  sites <- unique(cost_catchments$site)
  catchments2 <- lapply(sites, function(x) {
    message("Removing areas behind cone for ", x)
    # Create the cut line, perpendicular to h. The line
    # should be long enough to be sure to cut any weird
    # catchment danglies off. 50km each side should be
    # long enough. 
    # NOTE!! This will cause a BUG if the string isn't
    # wide enough to chop full_cc polygon in half!!
    stn <- stn[stn$site == x, ]
    x0 <- sf::st_coordinates(stn)[1]
    y0 <- sf::st_coordinates(stn)[2]
    alpha <- (180 - stn$mean) * (pi / 180) # convert degrees to radians + rotate 90°
    d <- 50000 # 30km
    x1 <- x0 + (d * cos(alpha))
    y1 <- y0 + (d * sin(alpha))
    x2 <- x0 - (d * cos(alpha))
    y2 <- y0 - (d * sin(alpha))
    string <- data.frame(geom = NA)
    string$geom <- sprintf("LINESTRING(%s %s, %s %s, %s %s)", x1, y1, x0, y0, x2, y2)
    string <- sf::st_as_sf(string, wkt = "geom", crs = 3005)
    string <- sf::st_as_sfc(string)
    # Cut the catchment with the string
    c <- cost_catchments[cost_catchments$site == x, ]
    c <- lwgeom::st_split(c, string)
    c <- sf::st_make_valid(c)
    c <- sf::st_collection_extract(c, "POLYGON")
    # Create a line directly under (one unit of resolution) the cut
    beta <- (90 - stn$mean) * (pi / 180)
    xx <- x0 - (res * cos(beta))
    yy <- y0 - (res * sin(beta))
    x1 <- xx + (5000 * cos(alpha)) # let's make this line much shorter, 10km total
    y1 <- yy + (5000 * sin(alpha))
    x2 <- xx - (5000 * cos(alpha))
    y2 <- yy - (5000 * sin(alpha))
    #xxyy <- st_point(c(xx, yy)) %>% st_sfc(crs = 3005)
    string <- data.frame(geom = NA)
    string$geom <- sprintf("LINESTRING(%s %s, %s %s, %s %s)", x1, y1, xx, yy, x2, y2)
    string <- sf::st_as_sf(string, wkt = "geom", crs = 3005)
    string <- sf::st_as_sfc(string)
    # Delete the section directly under the the cut
    #c <- c[!(st_contains(c, xxyy, sparse = FALSE)),]
    c <- c[!(sf::st_intersects(c, string, sparse = F)),]
    c <- c |> 
      dplyr::group_by(site) |>
      dplyr::summarise() # merge any pieces that may have been cut
    # Clean up the edges a bit
    c <- sf::st_buffer(c, res) |>
      smoothr::smooth() |>
      sf::st_intersection(watersheds[watersheds$site == x, ]) |>
      dplyr::filter(site == site.1) |>
      dplyr::select(-site.1)
    # Multipart polygon to singlepart
    c <- sf::st_cast(c, "POLYGON", warn = FALSE)
    # Now select the piece closest to the station
    #c <- c[st_nearest_feature(stn, c),]
    c <- c[sf::st_intersects(c, cones[cones$site == x,], sparse = FALSE), ]
  })
  
  catchments2 <- dplyr::bind_rows(catchments2)
  
  # Also, let's drop any catchments with tiny areas that might
  # confuse the wat_maz selection algorithm (e.g., see Port Chanal)
  catchments2$area_m2 <- units::drop_units(sf::st_area(catchments2))
  catchments2 <- catchments2[catchments2$area_m2 > 2000, ]
  
  return(catchments2)
  
}


# ACCESSIBLE AREAS --------------------------------------------------------



# Intersect the cost catchments with the 'MAMU accessible zone'
# area. If the area is inaccessible to MAMU, it should be cut out,
# even if it's within the cost catchment!

access_catchments <- function(cost_catchments, maz, stn, cones, 
                              raster_stats = TRUE, # extract summary statistics from other eg cost layers? T/F
                              output_plots = TRUE,
                              output_dir,
                              ...) {
  # Prep `maz`
  # Load up as a `stars` object - faster to vectorize
  maz <- stars::st_as_stars(maz)
  # Polygonize
  maz <- sf::st_as_sf(maz,
                      as_points = FALSE,
                      merge = TRUE,
                      na.rm = TRUE)
  maz <- sf::st_union(maz) # merge into one polygon
  # Simplify this giant maz polygon
  maz <- smoothr::smooth(maz, method = "chaikin")
  
  # Intersect with catchments
  cc_maz <- sf::st_intersection(cost_catchments, maz)
  
  # Prep `stn`
  stn <- sf::st_transform(stn, 3005)
  
  # Select pieces within a 10km radius of the origin
  sites <- unique(cost_catchments$site)
  cc_maz <- lapply(sites, function(x) {
    message("Selecting pieces accessible from the origin for ", x)
    origin <- stn[stn$site == x, ]
    tmp <- cc_maz[cc_maz$site == x, ]
    tmp <- sf::st_make_valid(tmp) |>
      sf::st_collection_extract("POLYGON") |>
      sf::st_cast("POLYGON", warn = FALSE)
    if (origin$region == "HG") {
      tmp <- tmp[sf::st_intersects(tmp, sf::st_buffer(origin, 10000), sparse = FALSE),]
    } else {
      #tmp <- tmp[sf::st_nearest_feature(origin, tmp),]
      x_cone <- cones[cones$site == x, ]
      tmp <- tmp[sf::st_intersects(x_cone, tmp)[[1]], ]
    }
    return(tmp)
  })
  
  #cc_maz <- sf::st_collection_extract(cc_maz, "POLYGON")
  cc_maz <- dplyr::bind_rows(cc_maz)
  
  # Drop empty geometries
  # E.g., Upper Campbell catchment doesn't
  # intersect with any `maz` area...
  cc_maz <- cc_maz[!is.na(cc_maz$site),]
  
  # Final catchment cleanup
  cc_maz <- cc_maz |>
    dplyr::select(site) |>
    dplyr::group_by(site) |>
    dplyr::summarise()
  
  # Recalculate catchment area
  cc_maz$area_ha <- units::set_units(sf::st_area(cc_maz), "ha")
  
  # Merge with `stn` to get region and loc
  cc_maz <- merge(cc_maz, sf::st_drop_geometry(stn[,c("site", "region", "loc")]), by = "site")
  
  # Optional: extract cost info?
  if (raster_stats == TRUE) {
    # TODO: this fails if one of the dots args is not supplied.
    # Unpack dots
    dots <- list(...)
    cost <- dots$cost
    # Extract raster values
    cc_maz$mean_cost <- exactextractr::exact_extract(cost, cc_maz, "mean")
    rm(dots)
  }
  
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
    watersheds <- plot_data$watersheds
    cones <- plot_data$cones
    nests <- plot_data$nests
    # Create plots
    dir.create(output_dir, showWarnings = F)
    for (i in sites) {
      out_path <- file.path(output_dir,
                            fs::path_sanitize(paste0(i, ".png"), replacement = "-"))
      message("Saving plot for ", i, " to '", out_path, "'")
      p1 <- plot_catchment(site = i, 
                           cost_catchment = cost_catchments,
                           accessible_catchment = cc_maz,
                           watersheds = watersheds,
                           cones = cones,
                           stn = stn, 
                           nests = nests)
      p2 <- plot_headings(site = i, headings = headings, h = h)
      p_all <- ggpubr::ggarrange(p1, p2, nrow = 1)
      ggplot2::ggsave(filename = out_path, plot = p_all)
    }
  }
  
  cc_maz <- sf::st_collection_extract(cc_maz, "POLYGON")
  return(cc_maz)
}


