#' 3. CALCULATE MAMU DENSITY AND POPULATION
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


prepare_mamu_habitat_gpkg <- function(path, regions) {
  # Read files
  gpkgs <- list.files(path, full.names = TRUE)
  gpkgs <- gpkgs[grep(".gpkg$", gpkgs)]
  gpkgs <- lapply(gpkgs, sf::st_read)
  sh <- dplyr::bind_rows(gpkgs)
  sh <- janitor::clean_names(sh)
  #sh <- sh[,c("suit_hab_cl", "geom")]
  # Calc area
  sh$sh_area_ha <- units::set_units(sf::st_area(sh), "ha")
  # Merge with regions
  sh <- sf::st_intersection(sh, regions)
  return(sh)
}

prepare_mamu_habitat_tiff <- function(path, regions) {
  # Read files
  tiff <- terra::rast(path)
  return(tiff)
}

habitat_in_catchments <- function(catchments, habitat) {
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
    
    habitat_count <- exactextractr::exact_extract(habitat, catchments, "count")
    habitat_m2 <- habitat_count * res^2 # each raster cell is res m x res m (e.g., 25m x 25m), aka res meters squared
    habitat_m2 <- units::set_units(habitat_m2, "m2")
    catchments$sh_area_ha <- units::set_units(habitat_m2, "ha") # set to hectares and add as a column to the data
    
    c_hab <- catchments
  }
  
  return(c_hab)
}


total_habitat_area <- function(habitat) {
  if (inherits(habitat, "sf")) {
    out <- sum(suitable_habitat$sh_area_ha)
    
  } else if (inherits(habitat, "SpatRaster")) {
    res <- unique(terra::res(habitat)) # assuming square cells here
    habitat_count <- length(terra::cells(habitat))
    habitat_m2 <- habitat_count * res^2 # each raster cell is res m x res m (e.g., 25m x 25m), aka res meters squared
    habitat_m2 <- units::set_units(habitat_m2, "m2")
    out <- units::set_units(habitat_m2, "ha") # set to hectares and add as a column to the data
    
  }
  return(out)
}


regional_habitat_area <- function(habitat, regions = NA) {
  if (inherits(habitat, "sf")) {
    out <- aggregate(sh_area_ha ~ region, suitable_habitat, FUN = "sum")
    
  } else if (inherits(habitat, "SpatRaster")) {
    # Now extract the amount of habitat area within the catchments,
    # using the exact extract method.
    res <- unique(terra::res(habitat)) # assuming square cells here
    
    habitat_count <- exactextractr::exact_extract(habitat, regions, "count")
    habitat_m2 <- habitat_count * res^2 # each raster cell is res m x res m (e.g., 25m x 25m), aka res meters squared
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
                     mamu_sh_density = mean(mamu_sh_density),
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
  # Select cols
  reg_dens <- reg_dens[,c("region", "n_catchments", "mamu_count", "min_count", "max_count", "mamu_sh_density", "sh_area_ha", "CI_level")]
  return(reg_dens)
}
  
  
  