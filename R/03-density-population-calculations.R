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
max_mamu <- function(s, stn) {
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
                     mean_max_mamu = mean(max_mamu_count))
  # Merge in region info
  stn <- sf::st_drop_geometry(stn)
  stn <- stn[,c("site", "region", "loc")]
  s <- merge(s, stn, by = "site")
  # Rearrange cols
  s <- s[,c("site", "region", "loc", "year_min", "year_max", "n_surveys", "mean_max_mamu")]
  # Return
  return(s)
}


# PREPARE MAMU HABITAT ----------------------------------------------------


prepare_mamu_habitat <- function(path, regions) {
  # Read files
  gpkgs <- list.files(path, full.names = TRUE)
  gpkgs <- gpkgs[grep(".gpkg$", gpkgs)]
  gpkgs <- lapply(gpkgs, sf::st_read)
  sh <- dplyr::bind_rows(gpkgs)
  sh <- janitor::clean_names(sh)
  sh <- sh[,c("suit_hab_cl", "geom")]
  # Calc area
  sh$sh_area_ha <- units::set_units(sf::st_area(sh), "ha")
  # Merge with regions
  sh <- sf::st_intersection(sh, regions)
  return(sh)
}

habitat_in_catchments <- function(catchments, habitat) {
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
  
  return(c_hab)
}

catchment_density <- function(catchment_habitat, mamu_station_count, mamu_count_col) {
  summary <- catchment_habitat |> 
    sf::st_drop_geometry() |> 
    dplyr::group_by(site, region) |> 
    dplyr::summarize(sh_area_ha = sum(sh_area_ha),
                     cat_area_ha = sum(cat_area_ha)) |> 
    dplyr::mutate(catchment_habitat_density = sh_area_ha / cat_area_ha)
  summary <- merge(summary, mamu_station_count, by = c("site", "region"), all.x = TRUE) # some MAMU stations have radar counts, but no headings :(
  summary$mamu_sh_density <- summary[[mamu_count_col]] / summary$sh_area_ha
  summary$mamu_cat_density <- summary[[mamu_count_col]] / summary$cat_area_ha
  return(summary)
}

extrapolate_density <- function(mamu_density, regional_sh_area) {
  reg_dens <- aggregate(mamu_sh_density ~ region, mamu_density, FUN = "mean")
  reg_dens <- merge(reg_dens, regional_sh_area, by = "region")
  reg_dens$mamu_count <- reg_dens$mamu_sh_density * reg_dens$sh_area_ha
  reg_dens <- reg_dens[order(reg_dens$region),]
  return(reg_dens)
}
  
  
  