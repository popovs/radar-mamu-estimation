#' Data inputs preparation fxns
#' 
#' The functions here read in the base data files of the pipeline
#' and prepare them for analysis.


# PREPARE GIS FILES -------------------------------------------------------

prepare_regions <- function(filepath) {
  regions <- sf::st_read(filepath)
  # This assumes the gpkg polygons are already in correct order. 
  # Assign factor levels to the existing row order.
  regions$region <- factor(regions$region, levels = forcats::fct_inorder(regions$region))
  return(regions)
}

# This version doesn't have the uncertainty col
prepare_nests <- function(filepath, regions = regions) {
  nests <- sf::st_read(filepath)
  # Intersect to get region the nest falls within,
  # then condense repeat visits to a nest into one geometry record
  nests <- nests |> 
    sf::st_intersection(regions) |>
    dplyr::select(MMCR_NAME, region) |>
    unique()
  return(nests)
}


# PREPARE SURVEYS ---------------------------------------------------------

# Read survey data, clean up, convert to sf object, and assign pipeline regions
#' The function first performs a few housekeeping tasks:
#' * Make R-friendly column names
#' * Create the `survey_date` column from the year, month, and day columns
#' * Re-calculate day-of-year (`doy`)
#'
#' It then filters down the data with the following criteria:
#' * Include only "Complete" surveys
#' * Include only pre-dawn surveys (i.e., after 2 am and before 8 am)
#' * Exclude August surveys
#' * Exclude Alaska border region surveys
#' * Exclude 'Interior' (i.e. only keep coastal) stations
prepare_surveys <- function(filepath, 
                            regions = regions,
                            N_years_min = 2) {
  s <- read.csv(filepath)
  
  # Basic dataframe tidying
  s <- janitor::clean_names(s)
  s$survey_date <- lubridate::make_date(year = s$year, month = s$month, day = s$day)
  s$doy <- lubridate::yday(s$survey_date)
  s$site <- stringr::str_trim(s$site)
  s$site <- stringr::str_squish(s$site)
  
  # Calculate N unique years within the data
  s <- s |> 
    dplyr::group_by(site) |> 
    dplyr::mutate(N_years = dplyr::n_distinct(year))
  
  # Filter down to appropriate data
  # Use only complete surveys
  s <- s[which(s$status == "Complete"), ]
  # Use only pre-dawn surveys
  s <- s[which(s$start_hr < 8), ]
  # Exclude August surveys
  s <- s[which(s$month != 8), ]
  # Exclude Alaska border region
  s <- s[which(s$region != "AS"), ]
  # Exclude Inland sites
  s <- s[which(s$loc == "Coastal"), ]
  # Minimum sample size (N years of data) cutoff
  s <- s[which(s$N_years >= N_years_min), ]
  
  # Create spatial object
  s <- s[!is.na(s$lon), ]
  s <- sf::st_as_sf(s, coords = c("lon", "lat"), 
                    crs = 4326, 
                    remove = FALSE) |>
    sf::st_transform(3005)
  
  # Merge in regions from pipeline
  # Anything already stored in the pre-existing `region` col
  # will be moved to `region_old`
  s$region_old <- s$region
  s <- dplyr::select(s, -region)
  regions <- sf::st_make_valid(regions)
  s <- sf::st_intersection(s, regions)
  return(s)
}



# PREPARE HEADINGS --------------------------------------------------------


prepare_headings <- function(filepath) {
  headings <- read.csv(filepath)
  headings <- janitor::clean_names(headings)
  headings <- headings[!is.na(headings$heading),]
  headings$lat <- as.numeric(headings$lat)
  headings$lon <- as.numeric(headings$lon)
  headings$speed <- as.numeric(headings$speed)
  headings$distance <- as.numeric(headings$distance)
  headings$heading <- as.numeric(headings$heading)
  
  # Flag if there's any likely data errors
  stopifnot("You have headings >360 that are likely typos." = !any(headings$heading > 360))
  
  # First clean up station lat/longs. While the radar 
  # station is on land, the birds are flying in off
  # to the side. So, for each observation at each station,
  # take the mean "bird entry point" - that is, the mean
  # location *relative to the radar unit* that the 
  # birds are flying into. All the headings are then
  # estimated from *that point* on the radar screen.
  headings$initial_direction <- tolower(headings$initial_direction)
  compass_directions <- data.frame(initial_direction = c("n", "nne", "ne", "ene",
                                                         "e", "ese", "se", "sse",
                                                         "s", "ssw", "sw", "wsw",
                                                         "w", "wnw", "nw", "nnw"),
                                   initial_degrees = seq(0, 360, 22.5)[1:16])
  headings <- merge(headings, compass_directions, all.x = TRUE)
  # Now, using the bearing from the station and the distance from the 
  # station, find the true "mean bird entry point" on to the radar
  # field, relative to the radar station.
  relative_pts <- as.data.frame(geosphere::destPoint(p = headings[,c("lon", "lat")], 
                                                     b = headings$initial_degrees, 
                                                     d = headings$distance))
  names(relative_pts) <- c("rel_lon", "rel_lat")
  headings <- cbind(headings, relative_pts)
  
  # TODO: Cut out stations with <5 headings?
  
  # Split out incoming and outgoing headings
  inc <- headings[grep("Incoming", headings$flightpath_type, ignore.case = T),
                  c("site", "lat", "lon", "rel_lat", "rel_lon", "initial_direction", "initial_degrees", "distance", "heading")]
  inc$flightpath_type <- "Incoming"
  out <- headings[grep("Outgoing", headings$flightpath_type, ignore.case = T),
                  c("site", "lat", "lon", "rel_lat", "rel_lon", "initial_direction", "initial_degrees", "distance", "heading")]
  out$flightpath_type <- "Outgoing"
  
  inc <- inc[!is.na(inc$heading),]
  out <- out[!is.na(out$heading),]
  
  # Make `out` the opposite heading
  # That is, we assume that 180° (polar opposite) direction
  # of the outgoing headings is equivalent to incoming ones
  out$heading <- out$heading - 180
  out$heading <- ifelse(out$heading < 0, out$heading + 360, out$heading)
  
  # Now merge in and out back together
  # `h` now functionally contains *'incoming' headings only*
  h <- rbind(inc, out)
  h <- h[order(h$site),]
  
  return(h)
}


# Extract individual stations from survey data
# bc station lat/lon can change by a few meters each year
# (impossible to set up the radar station in EXACTLY the
# same spot each year), take the mean lat/lon by
# site. Additionally, for radar stations with headings,
# use the *relative* latitude and longitude. Radar stations
# are typically set up on land, while birds are flying
# over water several hundred m away. So, extract the
# mean lat/lon of where birds are actually flying over
# when entering a watershed mouth to act as the 'true'
# radar station coordinate.
prepare_stn <- function(s, headings, regions) {
  # Grab region from the `regions` sf
  # First drop existing regions col and
  # replace with region the point falls into in the
  # project shapefile
  s <- dplyr::select(s, -region)
  s <- sf::st_transform(s, 3005)
  s <- suppressWarnings(sf::st_intersection(s, regions))
  
  # Get the mean coordinate
  stn <- sf::st_drop_geometry(s) # trying to get mean sf coords by group is a pain. 
  stn <- stn |> 
    dplyr::select(site, region, loc, lat, lon) |>
    dplyr::group_by(site, region, loc) |>
    dplyr::summarise(lat = mean(lat),
                     lon = mean(lon))
  
  # Next do the same thing with radar stations
  h_stn <- headings |>
    dplyr::select(site, rel_lat, rel_lon) |>
    dplyr::group_by(site) |>
    dplyr::summarise(rel_lat = mean(rel_lat, na.rm = TRUE),
                     rel_lon = mean(rel_lon, na.rm = TRUE))
  
  # Merge with stn, then grab the rel_lat/lon if it's available
  stn <- merge(stn, h_stn, by = "site", all.x = TRUE)
  stn$lat <- ifelse(is.na(stn$rel_lat), stn$lat, stn$rel_lat)
  stn$lon <- ifelse(is.na(stn$rel_lon), stn$lon, stn$rel_lon)
  stn$rel_coord <- !is.na(stn$rel_lat) # keep track of which stations used relative coords
  
  # Drop the two `rel` cols and turn into a spatial object
  stn <- stn |>
    dplyr::select(-rel_lat, -rel_lon) |>
    sf::st_as_sf(crs = 4326, 
                 coords = c("lon", "lat"), 
                 remove = FALSE) |>
    sf::st_transform(3005)
  
  return(stn)
}



# PREPARE WATERSHEDS ------------------------------------------------------

prepare_watersheds <- function(filepath, regions) {
  ws <- sf::st_read(filepath)
  ws <- sf::st_transform(ws, 3005)
  ws <- sf::st_intersection(ws, regions)
  return(ws)
}




# PREPARE SUITABLE HABITAT ------------------------------------------------


prepare_mamu_habitat_gpkg <- function(filepath, maz, regions) {
  # Read files
  sh <- sf::st_read(filepath)
  sh <- janitor::clean_names(sh)
  
  # Calc area
  sh$sh_area_ha <- units::set_units(sf::st_area(sh), "ha")
  # Merge with regions
  sh <- sf::st_intersection(sh, regions)
  return(sh)
}

prepare_mamu_habitat_tiff <- function(filepath, maz, band = NA) {
  # Read files
  tiff <- terra::rast(path)
  # Extract the band of interest, if applicable
  if (!is.na(band)) {
    tiff <- tiff[band]
  }
  
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



# Combine one or the other - detect if the suitable habitat file
# is a gpkg or raster and go from there. 
prepare_habitat <- function(filepath, regions) {
  gpkg_yn <- grepl(".gpkg", filepath)
  
  if (gpkg_yn) {
    out <- prepare_mamu_habitat_gpkg(filepath = filepath,
                                     regions = regions)
  } else {
    out <- prepare_mamu_habitat_tiff(filepath = filepath,
                                     band = "X1")
  }
  
  return(out)
  
}
