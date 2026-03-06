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
# NOTE 2025-11 This is an UNWEIGHTED MEAN. The density from a catchment 
# 1 ha in size is weighted the same as one from a 50 ha catchment.
# It is the equivalent error to taking the mean of several means.
# For this reason we've switched to weighted mean density for the regional
# density estimates. This is now obsolete.
# regional_density <- function(catchment_density, 
#                            group_by = "region", 
#                            dat_col = "density", 
#                            CI_level = 0.95,
#                            add_AKB = TRUE) {
#   
#   out <- bootmean(dat = catchment_density, 
#            group_by = group_by,
#            dat_col = dat_col,
#            CI_level = CI_level) |>
#     dplyr::rename(density = bootmean,
#                   density_lwr = boot_min,
#                   density_upr = boot_max) |>
#     # Units are sometimes annoying, but set them here
#     # to maintain consistency btwn `cc_density` & `reg_density`
#     dplyr::mutate(dplyr::across(c(density, density_lwr, density_upr),
#                                 ~units::set_units(., "1/ha")))
#   
#   # Add in Alaska Border Region?
#   if (add_AKB) {
#     out <- rbind(out, 
#                  data.frame(region = "AKB",
#                             N = 0, 
#                             density = out[["density"]][out$region == "NC"],
#                             density_lwr = out[["density_lwr"]][out$region == "NC"],
#                             density_upr = out[["density_upr"]][out$region == "NC"]))
#   }
#   
#   # For convenience... rearrange table such
#   # that records listed from North to South
#   if ("region" %in% names(out)) {
#     out$region <- factor(out$region,
#                          levels = c("AKB", "HG", "NC", "CC", "SC", "WNVI", "NVI", "MWVI", "SWVI", "EVI"))
#     out <- out[order(out$region), ]
#   }
#   
#   return(out)
# }



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
  max_nest_dist <- ceiling(max(nests$nest_dist_km, na.rm = TRUE))
  # 1) fit gamma function to the nest data
  # First fit the gamma distribution to the nest data
  fit <- fitdistrplus::fitdist(nests[["nest_dist_km"]], 
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
                              res) {
  # Health check
  if (any(plyr::count(density$region)$freq > 1)) stop(paste0("Please provide one unique density value per ", merge_by, "."))
  
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
  # Also equivalent to:
  # adj <- ((target_density) * exactextractr::ext_extract(nls_raster, region, "mean")) - target_density
  # But it's inconvenient to require a `regions` input as well for this fxn.
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

# Loop through each region and 'gammify', then merge into one big raster.
# Inputs: density table, nest_likelihood_full, and habitat.
extrapolate_density <- function(regions, reg_density, nest_likelihood) {
  # First loop through each region and calculated the 
  # adjustment factor that needs to be added to the 
  # density field overall so that after multiplying by
  # the nest_likelihood layer the mean density is still
  # the same value as in the `reg_density` table
  adj_tbl <- lapply(1:nrow(regions), function(x) {
    reg_x <- regions[x,] # Pull out the region in question
    reg_x <- merge(reg_x, reg_density) # Merge in the densities for that region
    nl_x <- terra::crop(nest_likelihood, reg_x, mask = TRUE) # Intersect nest_likelihood with region
    # Drop geometry of region to stop annoying geom from popping 
    # up in all subsequent operations
    reg_x <- sf::st_drop_geometry(reg_x)
    # Pull out each density value for that region
    dvals <- reg_x[,c("density", "density_lwr", "density_upr")]
    dvals <- as.numeric(dvals)
    # Now, for each density val in reg_x, run `adj_gamma_mean`
    adj_vals <- lapply(dvals, adj_gamma_mean, nest_likelihood = nl_x)
    adj_vals <- unlist(adj_vals)
    # Put into df to return from fxn
    out <- data.frame(region = reg_x[["region"]],
                      adj = adj_vals[1],
                      adj_lwr = adj_vals[2],
                      adj_upr = adj_vals[3])
    return(out)
  })
  
  adj_tbl <- dplyr::bind_rows(adj_tbl)
  
  # Drop AKB row. The density adjustment assumes the rasters contain
  # the FULL habitat of the region when doing the adjustment. In
  # the case of AKB, simply use North Coast values.
  adj_tbl <- adj_tbl[which(adj_tbl$region != "AKB"), ]
  
  # Add in AKB using NC vals.
  adj_tbl <- rbind(adj_tbl, 
                   data.frame(region = "AKB",
                              adj = adj_tbl[["adj"]][adj_tbl$region == "NC"],
                              adj_lwr = adj_tbl[["adj_lwr"]][adj_tbl$region == "NC"],
                              adj_upr = adj_tbl[["adj_upr"]][adj_tbl$region == "NC"]))
  
  adj_tbl
  
  # Great. Now for each region we have the density value and
  # desired adjustment value. Rasterize it with the adjustments.
  
  # DESIRED OUTPUT: raster map with the actual densities of birds
  # (N MAMU/1 ha) per cell
  
  # First update our density table with the adjustments.
  adj_reg_density <- merge(reg_density, adj_tbl, by = "region")
  
  adj_reg_density$density <- adj_reg_density$density + units::set_units(adj_reg_density$adj, "1/ha")
  adj_reg_density$density_lwr <- adj_reg_density$density_lwr + units::set_units(adj_reg_density$adj_lwr, "1/ha")
  adj_reg_density$density_upr <- adj_reg_density$density_upr + units::set_units(adj_reg_density$adj_upr, "1/ha")
  
  # Now rasterize adj_reg_density
  rdens_adj <- rasterize_density(density = adj_reg_density,
                                 sf = regions,
                                 merge_by = "region", 
                                 res = terra::res(nest_likelihood)[1])
  
  # Next actually produce the final density maps. 
  # Make sure extents match
  rdens_adj <- terra::resample(rdens_adj, nest_likelihood)
  # Apply the nest likelihood gamma decay to the adjusted density layer
  final_rdens <- rdens_adj * nest_likelihood
  # Return
  return(final_rdens)
}



# Now calculate the population from the MAMU density map!
calculate_population <- function(density_map, sf = NULL, merge_df = NULL) {
  # Convert MAMU density per hectare to MAMU density per
  # cell in the `density_map` raster. E.g., the density
  # value is {Y mamu}/1 ha; but all the rasters in the
  # pipeline are NOT at 1/ha resolution. So we'd be underestimating
  # the amount of MAMU per raster cell. We need to get
  # {Z mamu}/{res ha}.
  res <- terra::res(density_map)[1] # grab the resolution in m of the map cells
  cell_res <- res^2
  cell_res <- units::set_units(cell_res, "m2")
  cell_res <- units::set_units(cell_res, "ha") # so we need density per 6.25ha
  cell_res <- units::drop_units(cell_res) # now drop the units so it plays nice w terra
  
  # The density currently is {Y mamu}/1 ha.
  # We need it to be {Z mamu}/6.25 ha.
  density_map <- density_map * cell_res 
  
  # Now the map shows N MAMU per pixel, rather than density 
  # of MAMU per 1 ha. We can simply add the pixels up to get
  # our population of MAMU.
  
  # If a region (or other) sf was supplied, calculate the MAMU
  # count per polygon. Otherwise, simply count up the total
  # MAMU population across the whole raster.
  if (!is.null(sf)) {
    
    out <- cbind(sf::st_drop_geometry(sf),
                 exactextractr::exact_extract(density_map, sf, "sum"))
    
    n <- ncol(out)
    names(out)[(n-2):n] <- c("mamu", "mamu_lwr", "mamu_upr")
    
    out[,(n-2):n] <- lapply(out[,(n-2):n], round)
    
    # Merge in other df cols if they are supplied
    if (!is.null(merge_df)) {
      out <- merge(out, merge_df)
      out <- dplyr::select(out, region, dplyr::one_of("N"), dplyr::everything())
    }
    
  } else {
    # Else calculate total raster population
    out <- colSums(terra::values(density_map), na.rm = TRUE)
    names(out) <- c("mamu", "mamu_lwr", "mamu_upr")
    out <- round(out)
  }
  
  # For convenience... rearrange table such
  # that records listed from North to South
  if ("region" %in% names(out)) {
    out$region <- factor(out$region,
                         levels = c("AKB", "HG", "NC", "CC", "SC", "WNVI", "NVI", "MWVI", "SWVI", "EVI"))
    out <- out[order(out$region), ]
  }
  
  return(out)
  
}



# Outdated
# extrapolate_density <- function(cc_density, regional_sh_area, min_ss) {
#   #reg_dens <- aggregate(mamu_sh_density ~ region, cc_density, FUN = "mean")
#   # Remove density estimates with fewer than minimum sample size total
#   cc_density <- cc_density[cc_density$n_surveys >= min_ss,]
#   reg_dens <- cc_density |> 
#     dplyr::mutate(density_min = mamu_sh_density - density_CI,
#                   density_max = mamu_sh_density + density_CI) |>
#     dplyr::summarise(.by = region,
#                      mamu_sh_density = round(mean(mamu_sh_density), 3),
#                      density_min = mean(density_min), # ideally, cutting out minimum sample size will prevent NA's sneaking in here
#                      density_max = mean(density_max),
#                      n_inxn = dplyr::n(),
#                      CI_level = mean(CI_level))
#   reg_dens <- merge(reg_dens, regional_sh_area, by = "region")
#   # Mean MAMU count
#   reg_dens$mamu_count <- reg_dens$mamu_sh_density * reg_dens$sh_area_ha
#   reg_dens$mamu_count <- round(reg_dens$mamu_count)
#   # Min MAMU count
#   reg_dens$min_count <- reg_dens$density_min * reg_dens$sh_area_ha
#   reg_dens$min_count <- round(reg_dens$min_count)
#   # Max MAMU count
#   reg_dens$max_count <- reg_dens$density_max * reg_dens$sh_area_ha
#   reg_dens$max_count <- round(reg_dens$max_count)
#   # Reorder
#   reg_dens <- reg_dens[order(reg_dens$region),]
#   # Round other cols
#   reg_dens$sh_area_ha <- round(reg_dens$sh_area_ha, 0)
#   # Select cols
#   reg_dens <- reg_dens[,c("region", "n_inxn", "mamu_count", "min_count", "max_count", "mamu_sh_density", "sh_area_ha", "CI_level")]
#   return(reg_dens)
# }
#   
  
  