# 1. GENERATE RADAR CATCHMENTS

# In this script, radar survey 'catchments' containing
# marbled murrelet nesting habitat will be generated, using
# a few basic assumptions about bird habitat + BC Digital
# Elevation Model (DEM) data. 

# First, we read in bird headings data. We're using a simplified
# assumption here and taking the mean of all the bird flight
# headings here. While the heading on the radar screen might
# not necessarily translate 1-to-1 to the real world flight
# direction, the error in this will be captured by the cone
# we generate around the mean heading.

# DISCLAIMER 1: This script will take well over an hour to run.
# If you wish to simply test the script, modify the raster
# resolutions to a courser scale (e.g., 500m) to run through
# the script quickly.

# DISCLAIMER 2: If you are running this script on MacOS, there
# is a known bug with the default MacOS graphics device that 
# makes plotting really reaaallly slow. Install `ragg` and 
# change the graphics device to "AGG" in Options > General >
# Graphics. See here:
# https://gis.stackexchange.com/questions/424897/why-is-sf-struggling-to-run-plots-quickly-on-uk-datasets?rq=1

# TODO: set it so if overwrite == FALSE, still creates the dude in memory

# 01 SETUP ----------------------------------------------------------------

#library(devtools)
#devtools::install_github("popovs/MAMU")
library(MAMU) # MUST be >= 0.2.1
library(dplyr)
library(sf)
library(ggplot2)

dir.create("temp", showWarnings = FALSE)

res <- 250 # resolution in meters of your cones
overwrite <- FALSE # wether you wish to overwrite previous run files

# 01-1 Create `stn` object ----

# Read in survey data
s <- MAMU::process_radar_data("data/ECCC_FLNR_MAMU-RadarData-20240307.xlsx")

s$site <- s$new_name

# Create survey spatial object
s$lat <- s$lat
s$lon <- s$lon
s <- s[!is.na(s$lon), ]
s <- st_as_sf(s, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Create object of mean coordinate per survey station
# Also include count column to get n surveys at each
# unique site
stn <- s %>%
  st_transform(crs = 4326) %>%
  select(site, loc) %>%
  group_by(site) %>%
  aggregate(.,
            by = list(.$site),
            function(x) x = x[1]) %>%
  st_centroid() %>%
  select(-Group.1) %>%
  mutate("lat" = st_coordinates(.)[,2],
         "lon" = st_coordinates(.)[,1]) %>%
  select(site, loc, lat, lon)

# 01-2 Process headings ----

# Read in headings & clean up
headings <- readxl::read_excel("data/headings.xlsx",
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
stn_lookup <- MAMU::rs
headings <- merge(headings, stn_lookup, by.x = "name", by.y = "original_name", all.x = TRUE)
headings$name <- ifelse(is.na(headings$new_name), headings$name, headings$new_name)

# TODO: Cut out stations with <5 headings?

# Split out incoming and outgoing headings
inc <- headings[grep("Incoming", headings$flightpath_type, ignore.case = T),c("name", "heading")]
out <- headings[grep("Outgoing", headings$flightpath_type, ignore.case = T),c("name", "heading")]

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
rm(inc, out)

# 01-3 Calculate polar mean ----

# For each station, take the POLAR MEAN INCOMING & OUTGOING heading.
# We need to use POLAR/CIRCULAR stats to deal with 'wraparound' angles -
# e.g. 45° and 315° is ~due north, but the cartesian mean of the two -
# 180° - is directly due south!

# First convert degrees to radians
h$rad <- h$heading * pi / 180

# Define circular mean function
# where x is your heading angle IN RADIANS
circmean <- function(x) { 
  sinr <- sum(sin(x))
  cosr <- sum(cos(x))
  circmean <- atan2(sinr, cosr)
  circmean
}

# Calculate bootstrapped means and confidence intervals
# for each station

# For a given list `x` of radian measurements,
# estimate the circular mean on a subsample of x.
# Repeat `n` number of times to derive the bootstrapped
# sample; derive quantile breaks `alpha` from the 
# bootstrapped sample and output the result.
# Thanks to Ben Bolker (again): https://stackoverflow.com/a/53916042/1454785
circboot <- function(x, n = 1000, alpha = 0.05) {
  bootsample <- replicate(n, circmean(sample(x, replace = TRUE)))
  out <- setNames(c(mean(bootsample), quantile(bootsample, 
                                               c(alpha/2, 1-alpha/2))), 
                  c("mean", "lower", "upper"))
  return(out)
}

message("Calculating bootstrapped headings and 95 CI for bird headings at each station...")
h <- aggregate(rad ~ name, h, FUN = circboot)
message("Done bootstrapping.")

# Tidy it up...
h <- cbind(h[1], data.frame(h$rad))

# Calculate cone width
h$theta <- h$upper - h$lower

# Convert radians back to degrees
h <- h %>% 
  mutate(across(where(is.numeric), ~ .x * 180 / pi)) %>%
  mutate(across(where(is.numeric), ~ ifelse(. < 0, . + 360, .)))

# Merge with stations
h <- merge(h, stn, by.x = "name", by.y = "site")
names(h)[1] <- "site"

h <- st_as_sf(h) %>%
  st_transform(3005)

# Merge in region 
regions <- st_read("GIS/regions.gpkg")
h <- st_intersection(h, regions)
h <- h[order(h$site),]

rm(stn_lookup, s)

# 02 GENERATE CONES -------------------------------------------------------

# Create cones out of them with bootstrapped wedges
# Note your points MUST be in a projected coordinate system in meters!!
cones <- lapply(1:nrow(h), function(x){
  message("Generating cone for ", h$site[[x]], "...")
  radar_cone(pt = h[x,],
             radius = 30000,
             theta = h$theta[[x]], # using 95CI heading boundaries
             heading = h$mean[[x]],
             res = res) 
})

# Check that all of them actually rendered correctly - for some
# reason the above just results in NULL rasters.... suspected 
# memory issue.
isNaN <- lapply(cones, terra::minmax)
isNaN <- data.table::transpose(as.data.frame(isNaN))
isNaN <- row.names(isNaN[is.na(isNaN$V1),]) # Extract all the cones with null values
isNaN <- as.numeric(isNaN)

while (length(isNaN) > 0) {
  # Re-run the radar_cone function for NaN rasters...
  for (i in isNaN) {
    message("RE-Generating cone for ", h$site[[i]], "...")
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


# 02-1 Vectorize ----

# Next, turn all these catchment raster cones into vectors.
# It is far faster to vectorize first -> then mask with the 
# MAMU habitat area, rather than the other way around.

#catchments <- list()
for (i in 1:nrow(h)) {
  message("Vectorizing ", h[i,][["site"]])
  x <- cones[[i]]
  x <- x == 1
  x <- terra::as.polygons(x, crs = "epsg:3005")
  x <- sf::st_as_sf(x)
  x <- x[x$lyr.1 == 1, ] # Keep only area == 1
  # Next draw a convex hull around the cone AND the station coordinate.
  # Some cones are so narrow the station gets dropped
  x <- st_geometry(x) %>% 
    st_cast("MULTIPOINT") %>% 
    st_union(h[i,]) %>%
    st_convex_hull() %>%
    st_as_sf()
  x$site <- h[i,][["site"]]
  cones[[i]] <- x
  rm(x)
}

cones <- dplyr::bind_rows(cones)

cones$cone_area_ha <- units::set_units(sf::st_area(cones), "ha") 
mean(cones$cone_area_ha)
#hist(cones$cone_area_ha, breaks = 100)
plot(cones[1]) # Visually inspect

gc()


# 03 SELECT WATERSHED CATCHMENTS ------------------------------------------

# 03-1 Select watersheds that fall within cones ----

# Evidence indicates that MAMU generally fly within a given
# watershed when flying to a nesting site. 
# TODO: CITE!!!!!\
watersheds <- st_read("GIS/Watersheds/code_2_watersheds.gpkg")

# Get rid of tiny tiny watersheds + river deltas
# These mess up the pathfinding algorithm when assuming 
# MAMU can or cannot cross from one watershed to another
watersheds$AREA_SQM <- st_area(watersheds)
watersheds <- watersheds[watersheds$AREA_SQM > units::as_units(150000, "m^2"),]
#watersheds <- watersheds[!is.na(watersheds$WTRSHDCD_2),]

# Get a matrix of which watersheds intersect with cones
# Each column is a cone, each row a watershed
wat_cat <- st_intersects(watersheds, cones, sparse = FALSE)

# Loop through each column to extract out the watersheds that 
# touch each cone
# Replace the watersheds df with the subsetted list
watersheds <- lapply(1:ncol(wat_cat), function(x){
  watersheds[wat_cat[,x],]
})
names(watersheds) <- cones$site

# Bind it all into one df; keep the names as the id column
watersheds <- bind_rows(watersheds, .id = 'site')

ggplot() +
  geom_sf(data = watersheds,
          aes(color = site),
          show.legend = FALSE) +
  geom_sf(data = cones,
          fill = NA) +
  geom_sf(data = h)

dev.off()

# 03-2 Remove watersheds with <2% cone coverage ----

# Remove watersheds with <2% cone coverage -- EXCEPT
# for Haida Gwaii!! HG contains many, many small 
# cliffside catchments that empty directly into the
# sea -- and nesting data shows they will readily
# nest in these tiny catchments, unlike the mainland. 
# So HG minimum will be <1% cone coverage.

# I.e., remove any watershed slivers
intersect_area <- st_intersection(cones, watersheds) %>%
  mutate(intersect_area = units::set_units(st_area(.), "ha")) %>%
  filter(site == site.1) %>% # we don't care about cone slivers intersecting with other site watersheds
  dplyr::select(site, WSD_ID, intersect_area) %>%
  st_drop_geometry() %>%
  group_by(site, WSD_ID) %>%
  summarize(intersect_area = sum(intersect_area), .groups = "keep") %>%
  group_by(site) %>%
  mutate(total_area = sum(intersect_area),
         prct_coverage = intersect_area / total_area * 100)

#watersheds <- merge(watersheds, st_drop_geometry(cones), by = "site") # merge cone area in
watersheds <- merge(watersheds, intersect_area, by = c("site", "WSD_ID"))

#watersheds$prct_cone_overlap <- (watersheds$intersect_area / watersheds$cone_area_ha) * 100
#watersheds$prct_cone_overlap <- units::drop_units(watersheds$prct_cone_overlap)
watersheds$prct_coverage <- units::drop_units(watersheds$prct_coverage)

watersheds$centroid_x <- st_coordinates(st_centroid(watersheds))[,1]
watersheds$centroid_y <- st_coordinates(st_centroid(watersheds))[,2]

#watersheds$keep_yn <- watersheds$prct_cone_overlap > 3.5
watersheds$keep_yn <- (watersheds$prct_coverage > 2 | (watersheds$region == "HG" & watersheds$prct_coverage > 1))

# 03-2A Inspect each one... ----
dir.create("temp/cone_inspection", showWarnings = F)
lapply(h$site, function(x){
  message(x)
  p1 <- ggplot() +
    geom_sf(data = watersheds[watersheds$site == x, ],
            aes(color = keep_yn),
            show.legend = FALSE) +
    geom_sf(data = cones[cones$site == x, ],
            fill = NA) +
    geom_sf(data = h[h$site == x, ]) +
    geom_text(data = watersheds[watersheds$site == x, ],
              aes(label = round(prct_coverage, 2),
                  x = centroid_x,
                  y = centroid_y),
              size = 1) +
    ggtitle(x) +
    theme(axis.title = element_blank())
  p2 <- ggplot(data = headings[headings$name == x & headings$flightpath_type == "Incoming",]) + 
    geom_histogram(aes(x = heading)) + 
    geom_vline(xintercept = h[["mean"]][h$site == x],
               color = "red") +
    geom_vline(xintercept = h[["lower"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    geom_vline(xintercept = h[["upper"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    xlim(0, 360) + 
    ggtitle("Incoming") + 
    theme(axis.title = element_blank()) +
    theme_minimal()
  p3 <- ggplot(data = headings[headings$name == x & headings$flightpath_type == "Outgoing",]) + 
    geom_histogram(aes(x = heading)) + 
    geom_vline(xintercept = h[["mean"]][h$site == x],
               color = "red") +
    geom_vline(xintercept = h[["lower"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    geom_vline(xintercept = h[["upper"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    xlim(0, 360) + 
    ggtitle("Outgoing") + 
    theme(axis.title = element_blank()) +
    theme_minimal()
  p_all <- ggpubr::ggarrange(p1, ggpubr::ggarrange(p2, p3, ncol = 1), nrow = 1)
  ggsave(paste0("temp/cone_inspection/", x, ".png"),
         p_all)
})

# Drop anything w less than 2% of the cone overlapping it
watersheds <- watersheds[watersheds$keep_yn,]

rm(wat_cat)

# 03-3 Merge watersheds into catchments ----

# Save all watersheds for plotting
all_watersheds <- watersheds

# Dissolve watersheds together by site
watersheds <- watersheds %>%
  select(site) %>%
  group_by(site) %>%
  dplyr::summarize()

# 04 COST DISTANCE WITHIN CATCHMENTS --------------------------------------

# From the mouth of the selected watersheds, run a cost
# distance simulation as if a MAMU were flying up the 
# watershed. Then apply the 95% maximum cost cutoff as
# derived from the nest data for each region.


# 04-1 Read DEM ----

message("Reading in existing DEM file...")
dem_3005 <- terra::rast("GIS/DEM/BC_DEM_EPSG3005.tiff")
# Resample DEM to match the resolution specified above
if (all(terra::res(dem_3005) != res)) {
  message("Resampling DEM to target resolution (", res, "m)...")
  r <- dem_3005
  terra::res(r) <- res
  dem_3005 <- terra::resample(dem_3005, r)
  rm(r)
}

# 04-2 Intersect DEM with watersheds and run cost distance ----

catchments <- lapply(h$site, function(x) {
  message("Creating catchment for ", x)
  tmp <- watersheds[watersheds$site == x,]
  # TODO: find better way to get birds to cross water. 
  # perhaps if cone crosses over NA area, assume it's == 0
  
  # tmp <- st_coordinates(tmp)[,1:2] %>%  # extract concave hull of watershed, so we can include ocean areas in the raster
  #   st_multipoint() %>% 
  #   st_sfc(crs = 3005) %>% 
  #   st_concave_hull(ratio = 1) %>%
  #   st_as_sf()
  tmp <- terra::crop(dem_3005, tmp, mask = TRUE)
  tmp2 <- terra::crop(dem_3005, cones[cones$site == x, ], mask = TRUE) # also extract anything in the cone path - e.g. water - so birds can cross over water areas in the cone's path
  tmp <- terra::merge(tmp, tmp2)
  # Choose the point closest on the raster to the radar
  # station as the origin point
  # Extract origin coordinate
  i <- h[h$site == x,]
  j <- st_union(st_as_sf(terra::as.polygons(tmp)))
  origin <- st_nearest_points(i, j)
  origin <- st_coordinates(origin)
  origin <- origin[,1:2] # drop 'L1' column
  i <- st_coordinates(i)
  origin <- rbind(origin, i)
  # drop the duplicated rows - that's our point `h`, and we only
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
})

names(catchments) <- h$site

#rm(dem_3005)

# 04-3 Extract regional nest costs ----

nests <- st_read("GIS/MAMU_nests.gpkg")
nests <- st_buffer(nests, dist = nests$LOC_UNCE_1)

regions <- st_read("GIS/regions.gpkg")
nests <- st_intersection(nests, regions)

cost <- terra::rast("GIS/Cost_cutoffs/cost_layer.tiff")
# Resample to match the resolution specified above
if (all(terra::res(cost) != res)) {
  r <- cost
  terra::res(r) <- res
  cost <- terra::resample(cost, r)
  rm(r)
}

nests$cost <- exactextractr::exact_extract(cost, nests, 'mean')

cost_cutoffs <- setNames(data.frame(regions$region,
                                    unlist(lapply(regions$region, function(x){quantile(nests[["cost"]][nests$region == x], 0.025, na.rm = TRUE)[[1]]})),
                                    unlist(lapply(regions$region, function(x){quantile(nests[["cost"]][nests$region == x], 0.975, na.rm = TRUE)[[1]]}))
                                    ),
                         c("region", "cost_min", "cost_max"))
# If northern mainland coast and alaska border are NA, use CMC elevation cutoff
# Minimum
if (is.na(cost_cutoffs[["cost_min"]][cost_cutoffs$region == "NC"])) cost_cutoffs[["cost_min"]][cost_cutoffs$region == "NC"] <- cost_cutoffs[["cost_min"]][cost_cutoffs$region == "CC"]
if (is.na(cost_cutoffs[["cost_min"]][cost_cutoffs$region == "AKB"])) cost_cutoffs[["cost_min"]][cost_cutoffs$region == "AKB"] <- cost_cutoffs[["cost_min"]][cost_cutoffs$region == "CC"]
if (is.na(cost_cutoffs[["cost_min"]][cost_cutoffs$region == "NVI"])) cost_cutoffs[["cost_min"]][cost_cutoffs$region == "NVI"] <- mean(cost_cutoffs[["cost_min"]][cost_cutoffs$region %in% c("EVI", "SWVI", "MWVI")], na.rm = TRUE)
# Maximum
if (is.na(cost_cutoffs[["cost_max"]][cost_cutoffs$region == "NC"])) cost_cutoffs[["cost_max"]][cost_cutoffs$region == "NC"] <- cost_cutoffs[["cost_max"]][cost_cutoffs$region == "CC"]
if (is.na(cost_cutoffs[["cost_max"]][cost_cutoffs$region == "AKB"])) cost_cutoffs[["cost_max"]][cost_cutoffs$region == "AKB"] <- cost_cutoffs[["cost_max"]][cost_cutoffs$region == "CC"]
if (is.na(cost_cutoffs[["cost_max"]][cost_cutoffs$region == "NVI"])) cost_cutoffs[["cost_max"]][cost_cutoffs$region == "NVI"] <- mean(cost_cutoffs[["cost_max"]][cost_cutoffs$region %in% c("EVI", "SWVI", "MWVI")], na.rm = TRUE)

# Round to nearest 10 meter
cost_cutoffs$cost_min <- round(cost_cutoffs$cost_min, -1)
cost_cutoffs$cost_max <- round(cost_cutoffs$cost_max, -1)

rm(cost, regions)

# 04-3A Inspect... ----

# Grab region attributes
watersheds <- merge(watersheds, st_drop_geometry(h[,c("site", "region")]), by = "site")

# Grab regional cost cutoffs
watersheds <- merge(watersheds, cost_cutoffs)

library(tidyterra)

dir.create("temp/cost_inspection", showWarnings = F)
lapply(names(catchments), function(x){
  message(x)
  p1 <- ggplot() + 
    geom_spatraster_contour_filled(data = catchments[[x]],
                                   show.legend = FALSE) +
    geom_spatraster_contour(data = catchments[[x]],
                            breaks = c(#watersheds[["cost_min"]][watersheds$site == x],
                                       watersheds[["cost_max"]][watersheds$site == x])) +
    scale_fill_whitebox_d() +
    geom_sf(data = cones[cones$site == x, ],
            fill = NA) +
    geom_sf(data = h[h$site == x, ]) +
    geom_sf(data = st_intersection(nests, watersheds[watersheds$site == x,]),
            color = "red", 
            fill = "red") +
    ggtitle(x, subtitle = paste(watersheds[["region"]][watersheds$site == x], "max cost =", watersheds[["cost_max"]][watersheds$site == x])) +
    theme(axis.title = element_blank())
  p2 <- ggplot(data = headings[headings$name == x & headings$flightpath_type == "Incoming",]) + 
    geom_histogram(aes(x = heading)) + 
    geom_vline(xintercept = h[["mean"]][h$site == x],
               color = "red") +
    geom_vline(xintercept = h[["lower"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    geom_vline(xintercept = h[["upper"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    xlim(0, 360) + 
    ggtitle("Incoming") + 
    theme(axis.title = element_blank()) +
    theme_minimal()
  p3 <- ggplot(data = headings[headings$name == x & headings$flightpath_type == "Outgoing",]) + 
    geom_histogram(aes(x = heading)) + 
    geom_vline(xintercept = h[["mean"]][h$site == x],
               color = "red") +
    geom_vline(xintercept = h[["lower"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    geom_vline(xintercept = h[["upper"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    xlim(0, 360) + 
    ggtitle("Outgoing") + 
    theme(axis.title = element_blank()) +
    theme_minimal()
  p_all <- ggpubr::ggarrange(p1, ggpubr::ggarrange(p2, p3, ncol = 1), nrow = 1)
  ggsave(paste0("temp/cost_inspection/", x, ".png"),
         p_all)
})


# 04-4 Apply cost cutoff ----

catchments <- lapply(h$site, function(x) {
  message("Extracting max-cost catchment for ", x)
  # Extract catchment + apply cutoff
  tmp <- catchments[[x]]
  cutoff <- watersheds[["cost_max"]][watersheds$site == x]
  tmp <- terra::ifel(tmp > cutoff, NA, 1)
  # Vectorize
  tmp <- terra::as.polygons(tmp, crs = "epsg:3005")
  tmp <- sf::st_as_sf(tmp)
})

catchments <- dplyr::bind_rows(catchments)
catchments$site <- h$site
catchments <- catchments[, "site"]


# 05 DIRECTIONALITY CUTOFF ------------------------------------------------

# We've extracted the watershed regions that contain birds,
# but a few of them could be whittled down further - e.g.
# see Brittain or Kwinamass. Birds won't be flying *backwards*
# into our catchment areas. 

# 05-1 Cut out pieces behind heading ----
# Let's assume birds won't fly backwards. 
# Cut away any catchment area directly behind 
# the flight heading. 

catchments2 <- lapply(h$site, function(x) {
  message("Removing areas behind cone for ", x)
  # Create the cut line, perpendicular to h. The line
  # should be long enough to be sure to cut any weird
  # catchment danglies off. 10km each side should be
  # long enough. 
  stn <- h[h$site == x, ]
  x0 <- st_coordinates(stn)[1]
  y0 <- st_coordinates(stn)[2]
  alpha <- (180 - stn$mean) * (pi / 180) # convert degrees to radians + rotate 90°
  d <- 30000 # 30km
  x1 <- x0 + (d * cos(alpha))
  y1 <- y0 + (d * sin(alpha))
  x2 <- x0 - (d * cos(alpha))
  y2 <- y0 - (d * sin(alpha))
  string <- data.frame(geom = NA)
  string$geom <- sprintf("LINESTRING(%s %s, %s %s, %s %s)", x1, y1, x0, y0, x2, y2)
  string <- st_as_sf(string, wkt = "geom", crs = 3005)
  string <- st_as_sfc(string)
  # Cut the catchment with the string
  c <- catchments[catchments$site == x, ]
  c <- lwgeom::st_split(c, string)
  c <- st_make_valid(c)
  c <- st_collection_extract(c, "POLYGON")
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
  string <- st_as_sf(string, wkt = "geom", crs = 3005)
  string <- st_as_sfc(string)
  # Delete the section directly under the the cut
  #c <- c[!(st_contains(c, xxyy, sparse = FALSE)),]
  c <- c[!(st_intersects(c, string, sparse = F)),]
  c <- c %>% group_by(site) %>% dplyr::summarise() # merge any pieces that may have been cut
  # Clean up the edges a bit
  c <- st_buffer(c, res) %>%
    smoothr::smooth() %>%
    st_intersection(watersheds[watersheds$site == x, ]) %>%
    filter(site == site.1) %>%
    select(site, region)
  # Multipart polygon to singlepart
  c <- st_cast(c, "POLYGON", warn = FALSE)
  # Now select the piece closest to the station
  #c <- c[st_nearest_feature(stn, c),]
  c <- c[st_intersects(c, cones[cones$site == x,], sparse = FALSE), ]
})

catchments2 <- dplyr::bind_rows(catchments2)

# 05-1A Inspect... ----
dir.create("temp/catchment_inspection", showWarnings = F)
lapply(catchments$site, function(x){
  message(x)
  p1 <- ggplot() +
    geom_sf(data = all_watersheds[all_watersheds$site == x, ],
            aes(color = x),
            show.legend = FALSE) +
    geom_sf(data = catchments2[catchments2$site == x, ],
            fill = "grey",
            color = NA,
            alpha = 0.7,
            show.legend = FALSE) +
    geom_sf(data = cones[cones$site == x, ],
            fill = NA) +
    geom_sf(data = h[h$site == x, ]) +
    ggtitle(x, subtitle = h[["loc"]][h$site == x])
  p2 <- ggplot(data = headings[headings$name == x & headings$flightpath_type == "Incoming",]) + 
    geom_histogram(aes(x = heading)) + 
    geom_vline(xintercept = h[["mean"]][h$site == x],
               color = "red") +
    geom_vline(xintercept = h[["lower"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    geom_vline(xintercept = h[["upper"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    xlim(0, 360) + 
    ggtitle("Incoming") + 
    theme(axis.title = element_blank()) +
    theme_minimal()
  p3 <- ggplot(data = headings[headings$name == x & headings$flightpath_type == "Outgoing",]) + 
    geom_histogram(aes(x = heading)) + 
    geom_vline(xintercept = h[["mean"]][h$site == x],
               color = "red") +
    geom_vline(xintercept = h[["lower"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    geom_vline(xintercept = h[["upper"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    xlim(0, 360) + 
    ggtitle("Outgoing") + 
    theme(axis.title = element_blank()) +
    theme_minimal()
  p_all <- ggpubr::ggarrange(p1, ggpubr::ggarrange(p2, p3, ncol = 1), nrow = 1)
  ggsave(paste0("temp/catchment_inspection/", x, ".png"),
         p_all)
})

catchments <- catchments2
rm(catchments2)

# 06 ACCESSIBLE AREAS -----------------------------------------------------

# 06-1 Intersect with accessible areas ----

# Load up accessible nesting habitat tiff
# This time using stars - far more efficient to vectorize
# a stars object
if (any(grepl("\bmaz\b", ls()))) {
  maz <- stars::st_as_stars(maz)
} else {
  #maz <- stars::read_stars("GIS/MAMU_accessible_zone.tiff")
  maz <- stars::read_stars("GIS/Elevation_cutoffs/elevation_cutoffs.tiff")
}

# Polygonize
# This will take ~5 minutes
maz <- st_as_sf(maz,
                as_points = FALSE,
                merge = TRUE,
                na.rm = TRUE)
maz <- st_union(maz) # merge into one polygon

# Simplify this giant maz polygon
maz <- smoothr::smooth(maz, method = "chaikin")

# This will take just under ~1 minute
wat_maz <- st_intersection(catchments, maz)
rm(maz)

# 06-2 Select piece closest to origin ----

wat_maz <- lapply(h$site, function(x) {
  message("Selecting pieces accessible from the origin for ", x)
  origin <- h[h$site == x, ]
  tmp <- wat_maz[wat_maz$site == x, ]
  tmp <- st_make_valid(tmp) %>%
    st_collection_extract("POLYGON") %>%
    st_cast("POLYGON", warn = FALSE)
  tmp <- tmp[st_nearest_feature(origin, tmp),]
})

wat_maz <- dplyr::bind_rows(wat_maz)

# 06-1A Inspect each one... ----
# This will write over previous catchment inspection
dir.create("temp/catchment_inspection", showWarnings = F)
lapply(h$site, function(x){
  message(x)
  p1 <- ggplot() +
    geom_sf(data = all_watersheds[all_watersheds$site == x, ],
            aes(color = x),
            show.legend = FALSE) +
    geom_sf(data = catchments[catchments$site == x, ],
            color = NA,
            fill = "grey",
            alpha = 0.5,
            show.legend = FALSE) +
    geom_sf(data = wat_maz[wat_maz$site == x, ],
            color = NA, 
            fill = "#26D1EA",
            alpha = 0.3,
            show.legend = FALSE) +
    geom_sf(data = cones[cones$site == x, ],
            fill = NA) +
    geom_sf(data = h[h$site == x, ]) +
    ggtitle(x)
  # p1 <- ggplot() +
  #   geom_spatraster_contour_filled(data = catchments[[x]],
  #                                  show.legend = FALSE) +
  #   geom_sf(data = wat_maz[wat_maz$site == x,],
  #           color = "#26D1EA",
  #           fill = "white",
  #           alpha = 0.3) +
  #   geom_spatraster_contour(data = catchments[[x]],
  #                           breaks = c(#watersheds[["cost_min"]][watersheds$site == x],
  #                             watersheds[["cost_max"]][watersheds$site == x]),
  #                           linewidth = 1) +
  #   scale_fill_whitebox_d() +
  #   geom_sf(data = cones[cones$site == x, ],
  #           fill = NA) +
  #   geom_sf(data = h[h$site == x, ]) +
  #   geom_sf(data = st_intersection(nests, watersheds[watersheds$site == x,]),
  #           color = "red",
  #           fill = "red") +
  #   ggtitle(x, subtitle = paste(watersheds[["region"]][watersheds$site == x], "max cost =", watersheds[["cost_max"]][watersheds$site == x])) +
  #   theme(axis.title = element_blank())
  p2 <- ggplot(data = headings[headings$name == x & headings$flightpath_type == "Incoming",]) +
    geom_histogram(aes(x = heading)) +
    geom_vline(xintercept = h[["mean"]][h$site == x],
               color = "red") +
    geom_vline(xintercept = h[["lower"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    geom_vline(xintercept = h[["upper"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    xlim(0, 360) +
    ggtitle("Incoming") +
    theme(axis.title = element_blank()) +
    theme_minimal()
  p3 <- ggplot(data = headings[headings$name == x & headings$flightpath_type == "Outgoing",]) +
    geom_histogram(aes(x = heading)) +
    geom_vline(xintercept = h[["mean"]][h$site == x],
               color = "red") +
    geom_vline(xintercept = h[["lower"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    geom_vline(xintercept = h[["upper"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    xlim(0, 360) +
    ggtitle("Outgoing") +
    theme(axis.title = element_blank()) +
    theme_minimal()
  p_all <- ggpubr::ggarrange(p1, ggpubr::ggarrange(p2, p3, ncol = 1), nrow = 1)
  ggsave(paste0("temp/catchment_inspection/", x, ".png"),
         p_all)
})

catchments <- wat_maz
rm(wat_maz)

# 06 MINIMUM CUTOFF ZONES -------------------------------------------------

# Finally, we can cut out areas that MAMU can techincally
# *access*, e.g. <10m elevation, but are highly *unlikely*
# to nest in (i.e., the bottom 2.5% values for nest elev-
# ation and cost). 

# Load up accessible nesting habitat tiff
# This time using stars - far more efficient to vectorize
# a stars object
if (any(grepl("mnh", ls()))) {
  mnh <- stars::st_as_stars(mnh)
} else {
  mnh <- stars::read_stars("GIS/MAMU_nesting_habitat.tiff")
}

# Polygonize
# This will take ~5 minutes
mnh <- st_as_sf(mnh,
                as_points = FALSE,
                merge = TRUE,
                na.rm = TRUE)
mnh <- st_union(mnh) # merge into one polygon

# Simplify this giant mnh polygon
mnh <- smoothr::smooth(mnh, method = "chaikin")

# This will take just under ~1 minute
catchments2 <- st_intersection(catchments, mnh)

# Inspect final catchments...
dir.create("temp/final_catchment_inspection", showWarnings = F)
lapply(h$site, function(x){
  message(x)
  p1 <- ggplot() +
    geom_sf(data = all_watersheds[all_watersheds$site == x, ],
            aes(color = x),
            show.legend = FALSE) +
    geom_sf(data = catchments[catchments$site == x, ],
            color = NA, 
            fill = "#26D1EA",
            alpha = 0.3,
            show.legend = FALSE) +
    geom_sf(data = catchments2[catchments2$site == x, ],
            color = NA,
            fill = "#343434",
            alpha = 0.3,
            show.legend = FALSE) +
    geom_sf(data = cones[cones$site == x, ],
            fill = NA) +
    geom_sf(data = h[h$site == x, ]) +
    ggtitle(x)
  p2 <- ggplot(data = headings[headings$name == x & headings$flightpath_type == "Incoming",]) + 
    geom_histogram(aes(x = heading)) + 
    geom_vline(xintercept = h[["mean"]][h$site == x],
               color = "red") +
    geom_vline(xintercept = h[["lower"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    geom_vline(xintercept = h[["upper"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    xlim(0, 360) + 
    ggtitle("Incoming") + 
    theme(axis.title = element_blank()) +
    theme_minimal()
  p3 <- ggplot(data = headings[headings$name == x & headings$flightpath_type == "Outgoing",]) + 
    geom_histogram(aes(x = heading)) + 
    geom_vline(xintercept = h[["mean"]][h$site == x],
               color = "red") +
    geom_vline(xintercept = h[["lower"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    geom_vline(xintercept = h[["upper"]][h$site == x],
               color = "grey",
               linetype = "dashed") +
    xlim(0, 360) + 
    ggtitle("Outgoing") + 
    theme(axis.title = element_blank()) +
    theme_minimal()
  p_all <- ggpubr::ggarrange(p1, ggpubr::ggarrange(p2, p3, ncol = 1), nrow = 1)
  ggsave(paste0("temp/final_catchment_inspection/", x, ".png"),
         p_all)
})


# 06 EXTRACT CATCHMENT VALUES ---------------------------------------------


# Recalculate catchment area
catchments$area_ha <- units::set_units(sf::st_area(catchments), "ha")


# Extract mean cost and mean forest cover per catchment
cost <- terra::rast("GIS/Cost_cutoffs/cost_layer.tiff")
forest <- terra::rast("GIS/Forest_cover/forest_cover.tiff")

catchments$mean_cost <- exactextractr::exact_extract(cost, catchments, "mean")
catchments$mean_forest <- exactextractr::exact_extract(forest, catchments, "mean")

rm(cost, forest)

# Save em
if (!(any(grepl("radar_derived_catchments.gpkg", list.files("GIS"))))|overwrite == TRUE) st_write(catchments, "GIS/radar_derived_catchments.gpkg", append = FALSE)
st_write(cones, "GIS/flight_headings.gpkg")

# 03 CLEAN UP -------------------------------------------------------------

#if (!(any(grepl("MAMU_nesting_habitat.gpkg", list.files("GIS"))))|overwrite == TRUE) st_write(mnh, "GIS/MAMU_nesting_habitat.gpkg", append = FALSE)

rm(list = ls())
dev.off()
