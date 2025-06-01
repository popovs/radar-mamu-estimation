#' 2. CALCULATE MAMU DENSITY AND POPULATION
#' 
#' Following the creation of the inxn and standardization
#' of the radar survey data, we can now calculate the density
#' of birds within each catchment. From there, we will group 
#' each catchment by conservation region (based on the 
#' conservation region that the survey station was based in)
#' and calculate the mean density of birds per conservation
#' region. We will then multiply this mean density of birds
#' per region times the amount of available habitat area to 
#' come up with a population estimate across all of BC.


# MAX MAMU SURVEY ---------------------------------------------------------

# A small internal function to pass to the bootstrap function
bootmean_internal <- function(x, indices, na.rm = TRUE) {
  d <- x[indices] # allows boot to select sample 
  return(mean(d, na.rm = na.rm))
  }

# TO USE IN LIEU OF GLMM STANDARDIZED SURVEYS:
# Bootstrapped annual mean MAMU count + bootstrapped CIs
# Using this generic fxn which calculates the bootstrapped
# mean maximum by group
bootmean <- function(dat, 
                     group_by = "site", 
                     dat_col,
                     CI_level = 0.95) {
  #dat <- units::drop_units(dat)
  grps <- unique(dat[[group_by]])
  
  boots <- lapply(grps, function(x) {
    bootdat <- dat[which(dat[[group_by]] == x), ]
    bootdat <- bootdat[[dat_col]]
    out <- boot::boot(data = bootdat,
                      statistic = bootmean_internal,
                      R = 1000)
    return(out)
  })
  
  boots_ci <- lapply(boots, boot::boot.ci, type = "perc", conf = CI_level)
  
  # Extract means from boots as a vector
  out <- data.frame(grp = grps,
                    bootmean = unlist(lapply(boots, function(b) b$t0)))
  
  # Extract bootstrapped min/max 95% CIs
  bmin <- lapply(1:length(boots_ci), function(b) boots_ci[[b]]$percent[4])
  bmax <- lapply(1:length(boots_ci), function(b) boots_ci[[b]]$percent[5])
  
  bmin[lengths(bmin) == 0] <- NA # replace NULLs with NA
  bmax[lengths(bmax) == 0] <- NA # replace NULLs with NA
  
  # Now we can `unlist(bmin)` or `unlist(bmax)` without losing the NULLs
  out$boot_min <- unlist(bmin)
  out$boot_max <- unlist(bmax)
  
  # Rename `grp` col
  names(out)[1] <- group_by
  
  # Merge in sample size per group
  out <- dat |> 
    sf::st_drop_geometry() |>
    dplyr::group_by_at(group_by) |>
    dplyr::summarise(N = dplyr::n()) |>
    merge(out)
  
  return(out)
}


# SUPERCEDED. Using bootstrapped mean + CIs instead.
# max_mamu <- function(s, stn, CI_level = 95) {
#   s <- s |> dplyr::group_by(site) |> dplyr::mutate(total_effort = dplyr::n())
#   s <- aggregate(mamuinpd ~ site + year + total_effort, s, FUN = "max")
#   names(s)[4] <- "max_mamu_count"
#   # Select most recent max count of each
#   # s <- s |>
#   #   dplyr::arrange(site, year) |>
#   #   dplyr::group_by(site) |>
#   #   dplyr::slice(dplyr::n()) |>
#   #   dplyr::select(site, region, loc, year, total_effort, mamu_count) # rearrange cols
#   # Take the mean of the maximum mamu count across all years
#   s <- s |>
#     dplyr::group_by(site) |>
#     dplyr::summarise(year_min = min(year),
#                      year_max = max(year),
#                      n_surveys = max(total_effort),
#                      mean_max_mamu = round(mean(max_mamu_count)),
#                      sd = round(sd(max_mamu_count))) # while the overall distribution is non-normal, the max count for each station follows a normal dist.
#   # Calculate 95% CI
#   if (CI_level > 1) CI_level <- CI_level / 100
#   z_score <- qnorm(1 - ((1 - CI_level) / 2)) # get the z-score for the CI
#   s$CI <- z_score * s$sd / sqrt(s$n_surveys)
#   # Merge in region info
#   stn <- sf::st_drop_geometry(stn)
#   stn <- stn[,c("site", "region", "loc")]
#   s <- merge(s, stn, by = "site")
#   # Rearrange cols
#   s <- s[,c("site", "region", "loc", "year_min", "year_max",
#             "n_surveys", "mean_max_mamu", "sd", "CI")]
#   s$CI_level <- CI_level
#   # Return
#   return(s)
# }


# PREPARE MAMU HABITAT ----------------------------------------------------


prepare_mamu_habitat_gpkg <- function(path, maz, regions) {
  # Read files
  gpkgs <- list.files(path, full.names = TRUE)
  gpkgs <- gpkgs[grep(".gpkg$", gpkgs)]
  gpkgs <- lapply(gpkgs, sf::st_read)
  sh <- dplyr::bind_rows(gpkgs)
  sh <- janitor::clean_names(sh)
  #sh <- sh[,c("suit_hab_cl", "geom")]
  
  # 2025-05 UPDATE: Instead of intersecting with MAZ, 
  # we will be applying a gamma decay function of 
  # distance from shore. This gamma decay curve was
  # fit to nest distances from shore data.
  
  # # Mask the suitable habitat with the MAMU-accessible zone (`maz`)
  # # There's so much habitat deep in the interior that might technically be 
  # # appropriate trees for nesting, but are >100km from the coast.
  # # It's highly unlikely that birds are nesting there.
  # # The suitable habitat layer did not incorporate distance from
  # # the coast. 
  # # Prep `maz`
  # # Load up as a `stars` object - faster to vectorize
  # maz <- stars::st_as_stars(maz)
  # # Polygonize
  # maz <- sf::st_as_sf(maz,
  #                     as_points = FALSE,
  #                     merge = TRUE,
  #                     na.rm = TRUE)
  # maz <- sf::st_union(maz) # merge into one polygon
  # # Simplify this giant maz polygon
  # maz <- smoothr::smooth(maz, method = "chaikin")
  # # Intersect suitable habitat with MAMU-accessible zone
  # sh <- sf::st_intersection(sh, maz)
  
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
  
  # 2025-05 UPDATE: Instead of intersecting with MAZ, 
  # we will be applying a gamma decay function of 
  # distance from shore. This gamma decay curve was
  # fit to nest distances from shore data.
  
  # # Mask the suitable habitat with the MAMU-accessible zone (`maz`)
  # # There's so much habitat deep in the interior that might technically be 
  # # appropriate trees for nesting, but are >100km from the coast.
  # # It's highly unlikely that birds are nesting there.
  # # The suitable habitat layer did not incorporate distance from
  # # the coast. 
  # maz <- terra::resample(maz, tiff) # resample `maz` to be same extent as habitat
  # tiff <- (maz == 1) * tiff # only choose habitat areas where maz == TRUE
  
  return(tiff)
}

# use_probability: if the raster is a layer of habitat probabilities,
# should that be incorporated into the habitat estimation? TRUE or FALSE.
# If TRUE, it multiplies the number of non-NA pixels * the resolution * 
# the probability that that pixel is habitat to derive an estimate of habitat area.
# If FALSE, it multiplies the number of non-NA pixels * the resolution
# to derive an estimate of habitat area, under the assumption that any
# non-NA pixel is habitat (i.e. a binary yes/no habitat raster)
# NOTE habitat must be a singleband raster.
# `inxn` = "intersection" - ie a sf to intersect with
st_habitat_in_sf <- function(sf, habitat, use_probability = FALSE) {
  inxn <- sf
  sf::st_geometry(inxn) <- "geometry" # rename geometry col to be consistent
  # Run one set of functions if `habitat` is supplied
  # as a `sf` vs raster.
  if (inherits(habitat, "sf")) {
    # Drop any attributes from `habitat`
    habitat <- sf::st_geometry(habitat)
    
    # Intersect
    out <- sf::st_intersection(inxn, habitat)
    
    # Clean up output
    out <- sf::st_collection_extract(out, "POLYGON")
    out <- sf::st_make_valid(out)
    
    # Dissolve by non-geom attribute cols
    # (i.e. one polygon per original input feature)
    out <- out |>
      dplyr::group_by(dplyr::across(c(-geometry))) |>
      dplyr::summarise()
    
    # Calculate suitable habitat area per polygon
    out$sh_area_ha <- sf::st_area(out)
    out$sh_area_ha <- units::set_units(out$sh_area_ha, "ha")
    
  } else if (inherits(habitat, "SpatRaster")) {
    
    # Now extract the amount of habitat area within the inxn,
    # using the exact extract method.
    res <- unique(terra::res(habitat)) # assuming square cells here
    
    if (use_probability) {
      habitat_probs <- exactextractr::exact_extract(habitat, inxn, "sum")
      habitat_m2 <- habitat_probs * res^2 # each raster cell is res m x res m (e.g., 25m x 25m), aka res meters squared
    } else {
      habitat_count <- exactextractr::exact_extract(habitat, inxn, "count")
      habitat_m2 <- habitat_count * res^2 # each raster cell is res m x res m (e.g., 25m x 25m), aka res meters squared
    }
    
    habitat_m2 <- units::set_units(habitat_m2, "m2") # assuming resolution is in m2
    inxn$sh_area_ha <- units::set_units(habitat_m2, "ha") # set to hectares and add as a column to the data
    
    out <- inxn
  }
  
  return(out)
}



# CALCULATE DENSITY -------------------------------------------------------


# MAMU density within catchments
catchment_density <- function(mm = mean_max_mamu_cc, # mean_max_mamu_cc is the default
                              catchment_habitat = cc_habitat,
                              area_col = "sh_area_ha",
                              dat_cols = c("bootmean", "boot_min", "boot_max")) {
  c_hab <- sf::st_drop_geometry(catchment_habitat)
  mm <- merge(mm, c_hab, by = "site")
  
  # Calculate density
  mm$density <- mm[[dat_cols[1]]] / mm[[area_col]]
  mm$density_lwr <- mm[[dat_cols[2]]] / mm[[area_col]]
  mm$density_upr <- mm[[dat_cols[3]]] / mm[[area_col]]
  
  # Double check density_lwr < density_upr
  lwr_lt_upr <- all(mm$density_lwr < mm$density_upr, na.rm = TRUE)
  if(!lwr_lt_upr) warning("Did you mix up your upper and lower bounds cols? `density_lwr` is greater than `density_upr`.")
  
  # Return
  return(mm)
}

# From individual catchment densities, calculate a regional estimate
# with 95% bootstrapped CIs.
regional_density <- function(catchment_density, 
                           group_by = "region", 
                           dat_col = "density", 
                           CI_level = 0.95,
                           add_AKB = TRUE) {
  
  out <- bootmean(dat = catchment_density, 
           group_by = group_by,
           dat_col = dat_col,
           CI_level = CI_level) |>
    dplyr::rename(density = bootmean,
                  density_lwr = boot_min,
                  density_upr = boot_max) |>
    # Units are sometimes annoying, but set them here
    # to maintain consistency btwn `cc_density` & `reg_density`
    dplyr::mutate(dplyr::across(c(density, density_lwr, density_upr),
                                ~units::set_units(., "1/ha")))
  
  # Add in Alaska Border Region?
  if (add_AKB) {
    out <- rbind(out, 
                 data.frame(region = "AKB",
                            N = 0, 
                            density = out[["density"]][out$region == "NC"],
                            density_lwr = out[["density_lwr"]][out$region == "NC"],
                            density_upr = out[["density_upr"]][out$region == "NC"]))
  }
  
  return(out)
}



# DENSITY RASTER MANIPULATION ---------------------------------------------



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
                             habitat = NULL # `suitable_habitat` target
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
  
  if (is.null(habitat)) {
    
    return(p)
    
  } else {
    
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
  
  
}

# Given a table of densities and a shapefile of areas, 
# create a raster with values of `density` within the shapes of `sf`.
# The `merge_by` arg provides the column with which to merge the density
# table values to the sf. 
rasterize_density <- function(density, 
                              sf, # TODO: this will fail if sf is a raster 
                              merge_by = "region",
                              res
                              ) {
  # Health check
  if (any(plyr::count(density$region)$freq > 1)) stop(paste0("Please provide one unique density value per ", merge_by, "."))
  
  # Convert MAMU density per hectare to MAMU density per
  # unit cell in the raster pipeline. E.g., this density
  # in `density` is {x1 mamu}/1 ha; but all the rasters in the 
  # pipeline are NOT at 1/ha resolution. So we'd be underestimating
  # the amount of MAMU per raster cell. We need to get
  # {x2 mamu}/{res ha}. 
  cell_res <- res^2
  cell_res <- units::set_units(cell_res, "m2")
  cell_res <- units::set_units(cell_res, "ha") # so we need density per 6.25ha
  
  # The density currently is {x1 mamu}/1 ha.
  # We need it to be {x2 mamu}/6.25 ha.
  density$density <- density$density * cell_res
  density$density_lwr <- density$density_lwr * cell_res
  density$density_upr <- density$density_upr * cell_res
  
  # Merge to regions
  sf <- merge(sf, density, "region")
  
  # DROP UNITS! Otherwise terra interprets the values as 
  # CATEGORICAL INTEGERS (factors)
  sf <- units::drop_units(sf)
  sf <- terra::vect(sf)
  
  # Convert to density raster by polygon
  rdens <- terra::rast(terra::ext(sf),
                       res = res,
                       crs = terra::crs(sf))
  
  rdens <- terra::rasterize(sf, 
                            rdens, 
                            field = "density")
  
  rdens_lwr <- terra::rasterize(sf, 
                            rdens, 
                            field = "density_lwr")
  
  rdens_upr <- terra::rasterize(sf, 
                                rdens, 
                                field = "density_upr")
  
  # Stack em up
  rdens <- c(rdens, rdens_lwr, rdens_upr)
  
  # NOTE!! The extent of this raster matches `sf`.
  # It might not match the extents of other rasters in the dataset.
  return(rdens)
}


# Given a target mean density value and a nest_likelihood
# raster, calculate the value that needs to be added to
# the (density * nest_likelihood) raster to shift the 
# mean density value to reach the target density value.
# This needs to be done because otherwise we'd be underestimating
# regional density - we need to ensure it follows a gamma
# distribution, but still sums up to a mean regional MAMU
# density that lines up with the reg_density table!
adj_gamma_mean <- function(nest_likelihood, 
                           target_density) {
  # Extract the values of the nest_likelihood raster
  nls <- terra::values(nest_likelihood)
  nls <- nls[!is.na(nls)]
  # Calculate the adjustment value
  adj <- ((target_density * length(nls)) / sum(nls)) - target_density
  return(adj)
}


# Given:
# - a nest_likelihood raster
# - a rdens raster
# - a target_density
# Produce:
# A correctly adjusted gamma density raster
# This output will be a raster of MAMU
# habitat that has a mean density value 
# as calculated in the `reg_density`. 
# However, instead of being a flat raster
# assuming equal MAMU density across the 
# whole region, it is instead adjusted 
# such that the density is higher closer to
# shore and less dense further from shore.
# I.e... "gammify" (gammafy?) it!
gammify <- function(nl_rast, 
                    d_rast, 
                    target_density) {
  # Health checks
  # TODO: might be more efficient to do it in one go?
  # `d_rast` should be singleband
  #stopifnot("`d_rast` needs to be a singleband raster input." = length(names(d_rast)) == 1)
  
  # Make sure the extents match
  if (!(terra::ext(nl_rast) == terra::ext(d_rast))) {
    d_rast <- terra::resample(d_rast, nl_rast)
  }
  
  # Cut down nl_rast to exclude non-habitat
  # Otherwise we will underestimate our gamma
  # adjustment
  nl_rast <- (d_rast > 0) * nl_rast
  
  # While d_nl0 follows the probability distribution of 
  # the nest gamma decay raster now, it unfortunately now
  # has a lower mean density value than it originally 
  # started with. So we need to adjust it to match our
  # target mean raster value.
  adj <- adj_gamma_mean(nl_rast, target_density = target_density)
  
  # Apply the adjustment, then multiple density * likelihood
  # See the 'regional density with nest distance.R' script
  # for the algebra deriving this
  d_nl <- (d_rast + adj) * nl_rast
  
  # Check that the mean value approximates the target_density
  if (!all.equal(mean(terra::values(d_nl), na.rm = TRUE), target_density)) {
    warning("The gamma adjustment failed to produce a raster with a mean density value matching the target mean density value.")
  }
  
  return(d_nl)
  
}



# TODO: change this so it outputs a density map for the whole study
# area. Then this is the data product that gets summed up for total 
# MAMU. 
# Inputs: density table, nest_likelihood_full, and habitat.
extrapolate_density <- function(cc_density, regional_sh_area, min_ss) {
  #reg_dens <- aggregate(mamu_sh_density ~ region, cc_density, FUN = "mean")
  # Remove density estimates with fewer than minimum sample size total
  cc_density <- cc_density[cc_density$n_surveys >= min_ss,]
  reg_dens <- cc_density |> 
    dplyr::mutate(density_min = mamu_sh_density - density_CI,
                  density_max = mamu_sh_density + density_CI) |>
    dplyr::summarise(.by = region,
                     mamu_sh_density = round(mean(mamu_sh_density), 3),
                     density_min = mean(density_min), # ideally, cutting out minimum sample size will prevent NA's sneaking in here
                     density_max = mean(density_max),
                     n_inxn = dplyr::n(),
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
  reg_dens <- reg_dens[,c("region", "n_inxn", "mamu_count", "min_count", "max_count", "mamu_sh_density", "sh_area_ha", "CI_level")]
  return(reg_dens)
}
  
  
  