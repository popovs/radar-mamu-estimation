#' 2. CALCULATE MAMU DENSITY AND POPULATION
#' 
#' Following the creation of the catchments and standardization
#' of the radar survey data, we can now calculate the density
#' of birds within each catchment. From there, we will group 
#' each catchment by conservation region (based on the 
#' conservation region that the survey station was based in)
#' and calculate the mean density of birds per conservation
#' region. We will then multiply this mean density of birds
#' per region times the amount of available habitat area to 
#' come up with a population estimate across all of BC.


# MAX MAMU SURVEY ---------------------------------------------------------

# TO USE IN LIEU OF GLMM STANDARDIZED SURVEYS:
max_mamu <- function(s, stn, CI_level = 95) {
  s <- s |> dplyr::group_by(site) |> dplyr::mutate(total_effort = dplyr::n())
  s <- aggregate(mamuinpd ~ site + year + total_effort, s, FUN = "max")
  names(s)[4] <- "max_mamu_count"
  # Select most recent max count of each
  # s <- s |>
  #   dplyr::arrange(site, year) |>
  #   dplyr::group_by(site) |>
  #   dplyr::slice(dplyr::n()) |>
  #   dplyr::select(site, region, loc, year, total_effort, mamu_count) # rearrange cols
  # Take the mean of the maximum mamu count across all years
  s <- s |> 
    dplyr::group_by(site) |> 
    dplyr::summarise(year_min = min(year),
                     year_max = max(year),
                     n_surveys = max(total_effort),
                     mean_max_mamu = round(mean(max_mamu_count)),
                     sd = round(sd(max_mamu_count))) # while the overall distribution is non-normal, the max count for each station follows a normal dist.
  # Calculate 95% CI
  if (CI_level > 1) CI_level <- CI_level / 100
  z_score <- qnorm(1 - ((1 - CI_level) / 2)) # get the z-score for the CI
  s$CI <- z_score * s$sd / sqrt(s$n_surveys)
  # Merge in region info
  stn <- sf::st_drop_geometry(stn)
  stn <- stn[,c("site", "region", "loc")]
  s <- merge(s, stn, by = "site")
  # Rearrange cols
  s <- s[,c("site", "region", "loc", "year_min", "year_max", 
            "n_surveys", "mean_max_mamu", "sd", "CI")]
  s$CI_level <- CI_level
  # Return
  return(s)
}


# PREPARE MAMU HABITAT ----------------------------------------------------


prepare_mamu_habitat_gpkg <- function(path, maz, regions) {
  # Read files
  gpkgs <- list.files(path, full.names = TRUE)
  gpkgs <- gpkgs[grep(".gpkg$", gpkgs)]
  gpkgs <- lapply(gpkgs, sf::st_read)
  sh <- dplyr::bind_rows(gpkgs)
  sh <- janitor::clean_names(sh)
  #sh <- sh[,c("suit_hab_cl", "geom")]
  
  # Mask the suitable habitat with the MAMU-accessible zone (`maz`)
  # There's so much habitat deep in the interior that might technically be 
  # appropriate trees for nesting, but are >100km from the coast.
  # It's highly unlikely that birds are nesting there.
  # The suitable habitat layer did not incorporate distance from
  # the coast. 
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
  # Intersect suitable habitat with MAMU-accessible zone
  sh <- sf::st_intersection(sh, maz)
  
  # Calc area
  sh$sh_area_ha <- units::set_units(sf::st_area(sh), "ha")
  # Merge with regions
  sh <- sf::st_intersection(sh, regions)
  return(sh)
}

prepare_mamu_habitat_tiff <- function(path, maz, band = NA) {
  # Read files
  tiff <- terra::rast(path)
  # Extract the band of interest, if applicable
  if (!is.na(band)) {
    tiff <- tiff[band]
  }
  # Mask the suitable habitat with the MAMU-accessible zone (`maz`)
  # There's so much habitat deep in the interior that might technically be 
  # appropriate trees for nesting, but are >100km from the coast.
  # It's highly unlikely that birds are nesting there.
  # The suitable habitat layer did not incorporate distance from
  # the coast. 
  maz <- terra::resample(maz, tiff) # resample `maz` to be same extent as habitat
  tiff <- (maz == 1) * tiff # only choose habitat areas where maz == TRUE
  
  return(tiff)
}


# Nest probability decay function
# Certain habitats may be more or less suitable for nesting.
# However, while habitat may be suitable for MAMU nesting in terms
# of tree species composition, the suitable habitat layers do not
# explicitly take into account the fact that ~99% of nests occur
# within 30km of the coastline, and, crucially, that the further
# from the coast you are, the less likely the nests are likely to 
# occur. The nest data follow a gamma distribution of likelihood
# vs distance from shore. So, apply a gamma distribution decay
# curve to the suitable habitat layer such that distances <30km
# from shore are more likely, while distances >30km are less so.
nest_gamma_decay <- function(nests, # `nests` target 
                             coast, # `sea` target
                             habitat # `suitable_habitat` target
) {
  # Max nest dist 
  max_nest_dist <- ceiling(max(nests$dist_km, na.rm = TRUE))
  # 1) fit gamma function to the nest data
  # First fit the gamma distribution to the nest data
  fit <- fitdistrplus::fitdist(nests[["dist_km"]], 
                               distr = "gamma", 
                               method = "mle")
  
  # 2) derive a raster of distance from coast
  # Fill in the NA values of raster with `0`;
  # Replace any areas already == 0 with `2`
  # (i.e. all sea == 2, while all land == 0)
  coast <- terra::ifel(is.na(coast), 0, 2)
  coast_dist <- terra::gridDist(coast, target = 2) # generate raster with all distances from cells == 2
  coast_dist <- terra::ifel(coast_dist == 0, NA, coast_dist) # turn any sea areas to NA
  coast_dist <- coast_dist / 1000 # convert to m
  
  # 3) replace distance from coast with gamma probabilities
  # Derive a table of every distance from 0-62 km with the 
  # gamma density function value at each distance, at 10m intervals
  distances <- terra::values(coast_dist)
  # Predict values at each distance following the fit gamma distr 
  g <- dgamma(distances, 
              shape = fit$estimate[[1]], 
              scale =  1 / fit$estimate[[2]])
  # Normalize to 0-1
  p <- g / max(g, na.rm = TRUE)
  
  # Now put the nest probability values into a raster
  p <- terra::rast(vals = p, 
                   crs = terra::crs(coast_dist), 
                   terra::ext(coast_dist), 
                   res = terra::res(coast_dist))
  
  # 4) rasterize the suitable habitat (if not a raster already)
  if (!inherits(habitat, "SpatRaster")) {
    temp <- terra::rast(terra::vect(sf::st_geometry(habitat)),
                        res = terra::res(coast_dist))
    habitat <- terra::rasterize(terra::vect(habitat), temp)
  }
  
  habitat <- terra::resample(habitat, p) # resample `habitat` to be same extent as `p`
  
  # 5) ensure raster runs from 0-1, if not already
  h_vals <- terra::minmax(habitat)
  if (any(h_vals > 1)) {
    habitat <- habitat / max(terra::values(habitat), na.rm = TRUE)
  }
  
  # 6) multiply dist from coast nest probability raster * habitat raster
  nest_prob <- habitat * p
  
  return(nest_prob)
}

# use_probability: if the raster is a layer of habitat probabilities,
# should that be incorporated into the habitat estimation? TRUE or FALSE.
# If TRUE, it multiplies the number of non-NA pixels * the resolution * 
# the probability that that pixel is habitat to derive an estimate of habitat area.
# If FALSE, it multiplies the number of non-NA pixels * the resolution
# to derive an estimate of habitat area, under the assumption that any
# non-NA pixel is habitat (i.e. a binary yes/no habitat raster)
# NOTE habitat must be a singleband raster.
habitat_in_catchments <- function(catchments, habitat, use_probability = FALSE) {
  # Run one set of functions if `habitat` is supplied
  # as a `sf` vs raster.
  if (inherits(habitat, "sf")) {
    # Drop any attributes from `habitat`. Some catchments
    # cross over regional boundaries, so it can cause issues
    # with later aggregating density by region. So, instead,
    # we will rely on the region col that is in `catchments`.
    habitat <- dplyr::select(habitat, suit_hab_cl)
    
    # Intersect
    c_hab <- sf::st_intersection(habitat, catchments)
    
    # Specify what the area col actually references
    names(c_hab)[grep("^area_ha$", names(c_hab))] <- "cat_area_ha"
    
    # Re-calc hab area (some habitat polygons may have been clipped
    # by the catchment)
    c_hab$sh_area_ha <- units::set_units(sf::st_area(c_hab), "ha")
    
  } else if (inherits(habitat, "SpatRaster")) {
    
    # Specify what the area col actually references
    names(catchments)[grep("^area_ha$", names(catchments))] <- "cat_area_ha"
    
    # Now extract the amount of habitat area within the catchments,
    # using the exact extract method.
    res <- unique(terra::res(habitat)) # assuming square cells here
    
    if (use_probability) {
      habitat_probs <- exactextractr::exact_extract(habitat, catchments, "sum")
      habitat_m2 <- habitat_probs * res^2 # each raster cell is res m x res m (e.g., 25m x 25m), aka res meters squared
    } else {
      habitat_count <- exactextractr::exact_extract(habitat, catchments, "count")
      habitat_m2 <- habitat_count * res^2 # each raster cell is res m x res m (e.g., 25m x 25m), aka res meters squared
    }
    
    habitat_m2 <- units::set_units(habitat_m2, "m2") # assuming resolution is in m2
    catchments$sh_area_ha <- units::set_units(habitat_m2, "ha") # set to hectares and add as a column to the data
    
    c_hab <- catchments
  }
  
  return(c_hab)
}


# use_probability: if the raster is a layer of habitat probabilities,
# should that be incorporated into the habitat estimation? TRUE or FALSE.
# If TRUE, it multiplies the number of non-NA pixels * the resolution * 
# the probability that that pixel is habitat to derive an estimate of habitat area.
# If FALSE, it multiplies the number of non-NA pixels * the resolution
# to derive an estimate of habitat area, under the assumption that any
# non-NA pixel is habitat (i.e. a binary yes/no habitat raster)
# NOTE habitat must be a singleband raster.
total_habitat_area <- function(habitat, use_probability = FALSE) {
  if (inherits(habitat, "sf")) {
    out <- sum(habitat$sh_area_ha)
    
  } else if (inherits(habitat, "SpatRaster")) {
    res <- unique(terra::res(habitat)) # assuming square cells here
    
    if (use_probability) {
      habitat_probs <- sum(terra::values(habitat)[,1], na.rm = TRUE)
      habitat_m2 <- habitat_probs * res^2 # each raster cell is res m x res m (e.g., 25m x 25m), aka res meters squared
    } else {
      habitat_count <- length(terra::cells(habitat))
      habitat_m2 <- habitat_count * res^2 # each raster cell is res m x res m (e.g., 25m x 25m), aka res meters squared
    }
    
    habitat_m2 <- units::set_units(habitat_m2, "m2")
    out <- units::set_units(habitat_m2, "ha") # set to hectares and add as a column to the data
    
  }
  return(out)
}


regional_habitat_area <- function(habitat, regions = NA, use_probability = TRUE) {
  if (inherits(habitat, "sf")) {
    out <- aggregate(sh_area_ha ~ region, habitat, FUN = "sum")
    
  } else if (inherits(habitat, "SpatRaster")) {
    # Now extract the amount of habitat area within the catchments,
    # using the exact extract method.
    res <- unique(terra::res(habitat)) # assuming square cells here
    
    if (use_probability) {
      habitat_probs <- exactextractr::exact_extract(habitat, regions, "sum")
      habitat_m2 <- habitat_probs * res^2 # each raster cell is res m x res m (e.g., 25m x 25m), aka res meters squared
    } else {
      habitat_count <- exactextractr::exact_extract(habitat, regions, "count")
      habitat_m2 <- habitat_count * res^2 # each raster cell is res m x res m (e.g., 25m x 25m), aka res meters squared
    }
    
    habitat_m2 <- units::set_units(habitat_m2, "m2")
    regions$sh_area_ha <- units::set_units(habitat_m2, "ha") # set to hectares and add as a column to the data
    
    out <- sf::st_drop_geometry(regions[,c("region", "sh_area_ha")])
  }
  
  return(out)
}



catchment_density <- function(catchment_habitat, mamu_station_count, 
                              mamu_count_col, CI_col) {
  summary <- catchment_habitat |> 
    sf::st_drop_geometry() |> 
    dplyr::group_by(site, region) |> 
    dplyr::summarize(sh_area_ha = sum(sh_area_ha),
                     cat_area_ha = sum(cat_area_ha)) |> 
    dplyr::mutate(catchment_habitat_density = sh_area_ha / cat_area_ha)
  summary <- merge(summary, mamu_station_count, by = c("site", "region"), all.x = TRUE) # some MAMU stations have radar counts, but no headings :(
  summary$mamu_sh_density <- summary[[mamu_count_col]] / summary$sh_area_ha
  summary$density_CI <- summary[[CI_col]] / summary$sh_area_ha
  return(summary)
}

extrapolate_density <- function(mamu_density, regional_sh_area, min_ss) {
  #reg_dens <- aggregate(mamu_sh_density ~ region, mamu_density, FUN = "mean")
  # Remove density estimates with fewer than minimum sample size total
  mamu_density <- mamu_density[mamu_density$n_surveys >= min_ss,]
  reg_dens <- mamu_density |> 
    dplyr::mutate(density_min = mamu_sh_density - density_CI,
                  density_max = mamu_sh_density + density_CI) |>
    dplyr::summarise(.by = region,
                     mamu_sh_density = round(mean(mamu_sh_density), 3),
                     density_min = mean(density_min), # ideally, cutting out minimum sample size will prevent NA's sneaking in here
                     density_max = mean(density_max),
                     n_catchments = dplyr::n(),
                     CI_level = mean(CI_level))
  reg_dens <- merge(reg_dens, regional_sh_area, by = "region")
  # Mean MAMU count
  reg_dens$mamu_count <- reg_dens$mamu_sh_density * reg_dens$sh_area_ha
  reg_dens$mamu_count <- round(reg_dens$mamu_count)
  # Min MAMU count
  reg_dens$min_count <- reg_dens$density_min * reg_dens$sh_area_ha
  reg_dens$min_count <- round(reg_dens$min_count)
  # Max MAMU count
  reg_dens$max_count <- reg_dens$density_max * reg_dens$sh_area_ha
  reg_dens$max_count <- round(reg_dens$max_count)
  # Reorder
  reg_dens <- reg_dens[order(reg_dens$region),]
  # Round other cols
  reg_dens$sh_area_ha <- round(reg_dens$sh_area_ha, 0)
  # Select cols
  reg_dens <- reg_dens[,c("region", "n_catchments", "mamu_count", "min_count", "max_count", "mamu_sh_density", "sh_area_ha", "CI_level")]
  return(reg_dens)
}
  
  
  