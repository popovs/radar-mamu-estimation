#' Catchment fxns
#' 
#' This script contains functions that combine information 
#' about flight headings (radar_cone_fxns.R) 
#' with watershed data and nest cost distance to
#' create MAMU population catchments. That is, if
#' a MAMU flies into a given watershed mouth, it is
#' assumed it is targeting a nest within the corresponding
#' catchment delineation.



# SELECT TARGETED WATERSHEDS ----------------------------------------------

# Select watersheds that incoming MAMU may be targeting by
# selecting all watersheds that overlap with at least 2% of
# the generated radar cone.
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

# Calculate the cost of flying through the selected watersheds 
# from the radar entry point (mouth of the watershed)
# This fxn is designed for use for one site at a time
st_cost_catchment <- function(site,
                           watersheds,
                           dem, # note this is land + sea merged
                           cones,
                           stn,
                           nest_likelihood,
                           cost_function = "e" # cost fxn, one of movecost::mc_cost_functions
                           ) {
  message("Calculating cost catchment for ", site)
  
  # Dissolve watersheds together by site
  watersheds <- watersheds |>
    dplyr::select(site) |>
    dplyr::group_by(site) |>
    dplyr::summarize()
  
  # Filter all inputs down to the appropriate site
  w <- watersheds[watersheds$site == site, ]
  cone <- cones[cones$site == site, ]
  stn <- stn[stn$site == site, ]
  nl <- terra::crop(nest_likelihood, w, mask = TRUE) # crop nest likelihood raster to our watershed
  
  # Prepare watershed DEM
  tmp <- terra::crop(dem, w, mask = TRUE)
  tmp2 <- terra::crop(dem, cone, mask = TRUE) # also extract anything in the cone path - e.g. water - so birds can cross over water areas in the cone's path
  tmp <- terra::merge(tmp, tmp2)
  
  # Calculate resistance surface for the selected watersheds
  surf <- movecost::mc_surface(dtm = tmp, funct = cost_function)
  
  # Calculate the origin pixel of the surface raster
  # (In many cases, the station coordinate does not actually
  # exactly overlap the raster)
  p <- terra::as.points(tmp, values = FALSE) # extract raster centroids
  p <- sf::st_as_sf(p)
  origin <- p[sf::st_nearest_feature(stn, p), ] # pull closest raster centroid to our station
  
  # Calculate accumulated cost surface from origin
  acc <- movecost::mc_accum(surf, origin = origin)
  
  # Now let's make 100 'pseudo nests' within the watershed, all
  # falling within the most likely habitat areas, and calculate
  # path lengths of traveling from the origin to each pseudo-nest
  # Sample 10 'nests' from nlf
  pseudo_nests <- terra::spatSample(nl,
                                    size = 100, 
                                    method = "weights", # incorporate the raster probability into the sampling
                                    na.rm = TRUE,
                                    as.points = TRUE) |>
    sf::st_as_sf()
  
  # Calculate least-cost paths between origin and the pseudonests
  paths <- movecost::mc_paths(surf, 
                              origin = origin,
                              destin = pseudo_nests)
  
  # We assume throughout the next that MAMU are willing to fly at 
  # most 30km in a straight line to their nest; be definition, we
  # assume that same cutoff for meandering (i.e. longer) pathways.
  # If a MAMU takes a very meandering, long path from this origin
  # to a pseudonest, the assumption is that they in fact will prefer
  # to take a more efficient route to said nest via a different
  # watershed mouth. Therefore, we apply the same 30km distance 
  # cutoff to the least cost paths routes. 
  # So, extract out the raster cost of the longest route, up to 30km
  # in length.
  # The exception is Haida Gwaii: it's a very narrow island and
  # all known nests are within 5km of the coast. Let's assume a max
  # length route of 10km to be conservative.
  max_path <- ifelse(stn$region == "HG", 15000, 30000)
  max_cost <- paths$paths[paths$paths$length <= units::as_units(max_path, "m"), ] |> 
    dplyr::pull(cost) |> 
    max()
  
  # Draw cost catchment boundary
  cat <- movecost::mc_boundary(surface = surf, 
                               origin = origin,
                               limit = max_cost)
  
  # Prepare our function output
  out <- cat$boundaries
  out$site <- site
  out$region <- stn$region

  # Crop out any sea-areas (both for aesthetics and
  # also so sea areas don't get counted within the 
  # catchment areas)
  out <- sf::st_intersection(out, w)
  
  # Recalculate area; rearrange cols
  out <- out |>
    dplyr::rename(max_cost = limit) |>
    dplyr::mutate(area = sf::st_area(geometry),
                  area_ha = units::set_units(area, "ha")) |>
    dplyr::select(site, region, max_cost, area_ha, perimeter)
  
  return(out)
    
}


