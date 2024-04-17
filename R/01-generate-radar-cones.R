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
# TODO: remove ggplots

# 01 SETUP ----------------------------------------------------------------

#library(devtools)
#devtools::install_github("popovs/MAMU")
library(MAMU) # MUST be >= 0.2.1
library(dplyr)
library(sf)

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

#library(Hmisc)
# h <- h %>%
#   select(name, rad) %>%
#   group_by(name) %>%
#   summarise(data = list(smean.cl.boot(pick(everything()), conf.int = .95, B = 1000, na.rm = TRUE))) %>%
#   tidyr::unnest_wider(data) %>%
#   as.data.frame()

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

rm(stn_lookup)

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
                             heading = h$inc_out2[[i]],
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
  #if (st_geometry_type(x) == "MULTIPOLYGON") x <- st_convex_hull(x) # merge broken pixels together into single polygon if multipart
  # TODO:: st_convex_hull should include the original station coordinate
  cones[[i]] <- x
  # if (i == 1) {
  #   catchments <- x
  # } else {
  #   catchments <- rbind(catchments, x)
  # }
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
# TODO: if this works, move to GIS folder!
watersheds <- st_read("../Watersheds/conservation_region_watersheds.gpkg")

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

library(ggplot2)
ggplot() +
  geom_sf(data = watersheds,
          aes(color = site),
          show.legend = FALSE) +
  geom_sf(data = cones,
          fill = NA) +
  geom_sf(data = h)

# 03-2 Remove watersheds with <5% cone coverage ----

# I.e., remove any watershed slivers
intersect_area <- st_intersection(cones, watersheds) %>%
  mutate(intersect_area = units::set_units(st_area(.), "ha")) %>%
  filter(site == site.1) %>% # we don't care about cone slivers intersecting with other site watersheds
  dplyr::select(site, WSD_ID, intersect_area) %>%
  st_drop_geometry() %>%
  group_by(site, WSD_ID) %>%
  summarize(intersect_area = sum(intersect_area), .groups = "keep")

watersheds <- merge(watersheds, st_drop_geometry(cones), by = "site") # merge cone area in
watersheds <- merge(watersheds, intersect_area, by = c("site", "WSD_ID"))

watersheds$prct_cone_overlap <- (watersheds$intersect_area / watersheds$cone_area_ha) * 100
watersheds$prct_cone_overlap <- units::drop_units(watersheds$prct_cone_overlap)

# Drop anything w less than 5% of the cone overlapping it
watersheds <- watersheds[watersheds$prct_cone_overlap > 5,]

# 03-2A Inspect each one... ----
dir.create("temp/cone_inspection", showWarnings = F)
lapply(h$site, function(x){
  message(x)
  p1 <- ggplot() +
    geom_sf(data = watersheds[watersheds$site == x, ],
            aes(color = x),
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
  ggsave(paste0("temp/cone_inspection/", x, ".png"),
         p_all)
})


# Dissolve watersheds together by site
# watersheds <- watersheds %>%
#   select(site) %>%
#   group_by(site) %>%
#   dplyr::summarize()

rm(wat_cat)

# 03-3 Intersect with accessible areas ----

# Load up accessible nesting habitat tiff
# This time using stars - far more efficient to vectorize
# a stars object
if (any(grepl("maz", ls()))) {
  maz <- stars::st_as_stars(maz)
} else {
  maz <- stars::read_stars("GIS/MAMU_accessible_zone.tiff")
}

# Polygonize
# This will take ~5 minutes
maz <- st_as_sf(maz,
                as_points = FALSE,
                merge = TRUE,
                na.rm = TRUE)
maz <- st_union(maz) # merge into one polygon

# Simplify this giant maz polygon
maz <- smooth(maz, method = "chaikin")


# This will take just under ~1 minute
catchments <- st_intersection(watersheds, maz)
rm(maz)

# 03-3A Inspect each one... ----
dir.create("temp/catchment_inspection", showWarnings = F)
lapply(h$site, function(x){
  message(x)
  p1 <- ggplot() +
    geom_sf(data = watersheds[watersheds$site == x, ],
            aes(color = x),
            show.legend = FALSE) +
    geom_sf(data = catchments[catchments$site == x, ],
            color = NA,
            aes(fill = WSD_ID),
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
  ggsave(paste0("temp/catchment_inspection/", x, ".png"),
         p_all)
})


# 03-4 Select accessible areas ----

# Now we've cut out the inaccessible zones, we can see some
# of our catchments are cut up into pieces. (e.g., see 
# Aaltanhash). We need to choose the pieces that are actually
# reachable from the radar station. This will entail two 
# assumptions:
#   1) Any catchment pieces that are adjacent to/touch the 
#      origin piece are accessible to a flying MAMU.
#   2) The catchment piece closest to the radar station is thh
#      'origin' catchment piece

# 03-4A Merge adjacent catchment pieces ----
# Dissolve catchments together by site, then split multipart
# polygons into singlepart 
# Need to cast to multipart first bc of bug: https://github.com/r-spatial/sf/issues/763
catchments <- catchments %>%
  select(site) %>%
  group_by(site) %>%
  dplyr::summarize() %>%
  st_cast("MULTIPOLYGON") %>%
  st_cast("POLYGON", warn = FALSE)

catchments$id <- rownames(catchments)

# 03-4B Select pieces closest to origin ----

# If a flying bird can enter the origin piece, logic dictates
# that it can then subsequently fly into any neighboring areas
# that are also accessible. Therefore, select any pieces that
# touch the origin piece.

# For each radar station site and each corresponding catchment,
# choose the origin piece (catchment piece that is closest to
# the radar station)
catchments2 <- list()
catchments2 <- lapply(h$site, function(x) {
  message("Finding origin site for ", x, "...")
  hx <- h[h$site == x,]
  cx <- catchments[catchments$site == x,]
  cx$origin <- FALSE
  cx[["origin"]][st_nearest_feature(hx, cx)] <- TRUE
  cx
})
catchments2 <- dplyr::bind_rows(catchments2)

# Inspect again...
lapply(h$site, function(x){
  message(x)
  p1 <- ggplot() +
    geom_sf(data = watersheds[watersheds$site == x, ],
            aes(color = x),
            show.legend = FALSE) +
    geom_sf(data = catchments2[catchments2$site == x, ],
            color = NA,
            aes(fill = origin),
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
  ggsave(paste0("temp/catchment_inspection/", x, ".png"),
         p_all)
})

# Replace catchments w catchments2
catchments2 <- catchments2[catchments2$origin == TRUE,]
catchments <- catchments2
rm(catchments2)



# 03-3 xxxCUT SECTION? Select correct watershed ----

# The catchments are now essentially chopped up cones divided
# up by either the non-nesting zones from `mnh` or into 
# multiple watersheds. Now we need to choose the correct
# piece of the catchment that relates to our radar station and
# discard the rest of the miscellaneous values.

# Next, choose whichever watershed ID falls in line/intersects 
# with the heading
# https://stackoverflow.com/questions/69820690/finding-coordinates-from-heading-and-distance-in-r

# h2 <- st_transform(h, crs = 4326)
# h2$lon.1 <- NA
# h2$lat.1 <- NA
# 
# for (i in 1:nrow(h2)) {
#   ic <- st_coordinates(h2[i,])
#   b <- sf::st_drop_geometry(h2[i, "inc_out2"])[[1]]
#   nc <- geosphere::destPoint(ic, b = b, d = 30000)
#   h2[i, "lon.1"] <- nc[1]
#   h2[i, "lat.1"] <- nc[2]
# }
# rm(ic, b, nc, i)
# 
# h2 <- st_drop_geometry(h2)
# h2$geom <- sprintf("LINESTRING(%s %s, %s %s)",
#                    h2$lon, h2$lat,
#                    h2$lon.1, h2$lat.1)
# h2 <- st_as_sf(h2, wkt = "geom")
# st_crs(h2) <- 4326
# h2 <- st_transform(h2, 3005)
# 
# 
# # Check it out...
# x <- "Southgate"
# x <- "Aaltanhash"
# 
# ggplot() +
#   geom_sf(data = catchments[catchments$site == x,],
#           aes(color = WSD_ID)) +
#   geom_sf(data = h[h$site == x,]) +
#   geom_sf(data = h2[h2$site == x,])
# 
# ggplot() +
#   geom_sf(data = catchments[catchments$site == x,][st_intersects(catchments[catchments$site == x,], h2[h2$site == x,], sparse = F),],
#           aes(color = WSD_ID)) +
#   geom_sf(data = h[h$site == x,]) +
#   geom_sf(data = h2[h2$site == x,])

# `h2` is now our dataframe of headings (lines that are 30km 
# long oriented in the direction of each station's mean
# heading). 
# The intersection between each line and the catchment pieces
# it touches are the 'true' catchment. All other pieces can be
# dropped.
# Now cycle through each catchment and extract out the pieces
# that intersect with it's corresponding `h2` line.


# catchments$area_ha <- units::set_units(sf::st_area(catchments), "ha")
# catchments <- catchments[order(catchments$site, catchments$area_ha),]
# catchments$keep_yn <- NA
# 
# for (i in 1:nrow(h2)) {
#   site <- h2[["site"]][i]
#   message("Cleaning up ", site, "...")
#   catchments[["keep_yn"]][catchments$site == site] <- st_intersects(catchments[catchments$site == site,], h2[h2$site == site,], sparse = F)
# }

# 
# keep_cat_yn <- lapply(h2$site, function(x) {
#   st_intersects(catchments[catchments$site == x,], h2[h2$site == x,], sparse = F)
# })
# catchments$keep_yn <- unlist(keep_cat_yn)

# Check it out...

# x <- "Aaltanhash"
# x <- "Southgate"
# x <- "Watta"
# 
# ggplot() +
#   geom_sf(data = catchments[catchments$site == x,][st_intersects(catchments[catchments$site == x,], h2[h2$site == x,], sparse = F),],
#           aes(color = WSD_ID)) +
#   geom_sf(data = h[h$site == x,]) +
#   geom_sf(data = h2[h2$site == x,]) +
#   ggtitle(x)
# 
# ggplot() +
#   geom_sf(data = catchments[catchments$site == x & catchments$keep_yn == TRUE,]) +
#   geom_sf(data = h[h$site == x,]) +
#   ggtitle(x)


# Remove non-intersecting catchment pieces
# catchments <- catchments[catchments$keep_yn == TRUE, ]



# There should now be 95 features, matching the original 95 headings...
nrow(h) == nrow(catchments)



# 04 DIRECTIONALITY CUTOFF ------------------------------------------------

# We've extracted the watershed regions that contain birds,
# but a few of them could be whittled down further - e.g.
# see Brittain or Kwinamass. Birds won't be flying *backwards*
# into our catchment areas. 

# 04-1 Calculate cost cones ----
# Cost cones
cost_cones <- lapply(1:nrow(h), function(x){
  message("Calculating cost cone for ", h[x,]$site)
  tmp <- radar_cone(pt = h[x,],
                    radius = 30000,
                    theta = h$theta[[x]], # using 95CI heading boundaries
                    heading = h$mean[[x]],
                    invert = TRUE, # we want our bird flight direction to be zero - zero difficult flying through it
                    res = res) 
  # Now we're going to assume the exact *inverse* of the heading
  # cone is where the bird will NOT go. We'll set a 1/3 cone with
  # CONFIDENCE as no-go
  tmp2 <- radar_cone(pt = h[x,],
                     radius = 30000,
                     theta = 120,
                     heading = h$mean[[x]] - 180,
                     res = res) 
  tmp2 <- terra::ifel(tmp2 < 1, tmp, NA) # If tmp2 < 1, REPLACE WITH VALUES FROM TMP
  # Now do our cost difficulty!
  # It will be NO DIFFICULTY AT ALL to fly within the cone area (tmp2 == 0).
  # It will be INCREASINGLY DIFFICULT to fly away from cone area (tmp > 0).
  # It will be IMPOSSIBLE to fly directly behind the cone area (1/3 wedge where tmp is NA).
  # Our (x,y) origin point will be -1 (where we are flying FROM).
  xy <- sf::st_coordinates(h[x,])
  xy <- terra::cellFromXY(tmp2, xy) # extract the cell number that contains our origin point
  tmp2[xy] <- -1 # set our origin cell value to -1
  cost <- terra::costDist(tmp2, target = -1)
  # Next we need to choose our cutoff point. We're going to assume our bird
  # can access 50% of the remaining area.
  cutoff <- quantile(terra::values(cost), 0.5, na.rm = TRUE)
  # Now vectorize
  cost <- terra::ifel(cost < cutoff, 1, NA)
  cost <- terra::as.polygons(cost, crs = "epsg:3005")
  cost <- sf::st_as_sf(cost)
  cost <- smoothr::smooth(cost, method = "ksmooth")
  cost$site <- h[x,]$site
  cost <- cost[,c("site", "geometry")]
  return(cost)
})

cost_cones <- bind_rows(cost_cones)
plot(cost_cones)

# 04-2 Intersect with catchments ----

# Subset each catchment and intersect with it's 
# corresponding cost cone.
catchments2 <- lapply(catchments$site, function(x) {
  message("Intersecting ", x)
  i <- st_make_valid(catchments[catchments$site == x, ])
  j <- st_make_valid(cost_cones[cost_cones$site == x, ])
  st_intersection(i, j)
})

catchments2 <- bind_rows(catchments2)

# Inspect again...
lapply(h$site, function(x){
  message(x)
  p1 <- ggplot() +
    geom_sf(data = watersheds[watersheds$site == x, ],
            aes(color = x),
            show.legend = FALSE) +
    geom_sf(data = catchments[catchments$site == x, ],
            color = NA,
            fill = "#a9a9a9",
            alpha = 0.3,
            show.legend = FALSE) +
    geom_sf(data = catchments2[catchments2$site == x, ],
            color = NA,
            aes(fill = origin),
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
  ggsave(paste0("temp/catchment_inspection/", x, ".png"),
         p_all)
})


# Replace catchments w catchments2
catchments2 <- catchments2[catchments2$origin == TRUE,]
catchments <- catchments2
rm(catchments2, cost_cones)



# 05 xxxCUT SECTION? DISTANCE FROM ORIGIN -------------------------------------------------

# Now we have our feasible MAMU travel catchments. If a MAMU
# had an infinite distance it could travel, it could reach 
# anywhere within the network. However, we will cut off flights
# with a 30km distance from the origin (for inland areas, this
# will likely hit the edge of the allowable MAMU zone quicker).

# 05-1 Rasterize catchments ----

# raster_cats <- lapply(1:nrow(catchments), function(x) {
#   message("Rasterizing ", catchments[["site"]][x])
#   tmp <- catchments[x,]
#   template <- terra::rast(terra::vect(tmp), res = ifelse(res > 25, res/10, res))
#   out <- terra::rasterize(terra::vect(tmp), y = template, touches = TRUE) # give it a slight buffer so we don't lose any area unintentionally
#   return(out)
#   })
# 
# # 05-2 Calculate distance from origin ----
# 
# catchments2 <- lapply(1:nrow(h), function(x){
#   message("Calculating distance from ", h[["site"]][x], " radar station...")
#   # Extract origin coordinate
#   i <- h[x,]
#   j <- catchments[x,]
#   origin <- st_nearest_points(i, j)
#   origin <- st_coordinates(origin)
#   origin <- origin[,1:2] # drop 'L1' column
#   i <- st_coordinates(i)
#   origin <- rbind(origin, i)
#   # drop the duplicated rows - that's our point `h`, and we only 
#   # care about origin point in/on the polygon
#   if (nrow(unique(origin)) == 1) { # vanilla `ifelse` doesn't play nicely with matrices
#     origin <- unique(origin)
#   } else {
#     origin <- origin[!(duplicated(origin) | duplicated(origin, fromLast = TRUE)), ]
#   }
#   origin <- matrix(origin, ncol = 2)
#   # Assign origin cell in raster a value of '0'
#   tmp <- raster_cats[[x]]
#   tmp[terra::cellFromXY(tmp, origin)] <- 0 # Set the origin value == 0
#   # Compute grid distance from origin
#   tmp <- terra::gridDist(tmp, 0)
#   # Extract anything < 30km
#   tmp <- terra::ifel(tmp <= 30000, 1, NA)
#   # Vectorize
#   tmp <- terra::as.polygons(tmp, crs = "epsg:3005")
#   tmp <- sf::st_as_sf(tmp)
#   tmp$site <- h[["site"]][x]
#   tmp <- tmp[,"site"] # drop 'layer' column
#   tmp <- smoothr::smooth(tmp, method = "ksmooth")
#   tmp <- st_make_valid(tmp)
#   return(tmp)
# })
# 
# catchments2 <- dplyr::bind_rows(catchments2)
# 
# 
# # Save them
# catchments <- catchments2
# rm(catchments2)



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
mnh <- smooth(mnh, method = "chaikin")


# This will take just under ~1 minute
catchments <- st_intersection(catchments, mnh)

# Inspect final catchments...
dir.create("temp/final_catchment_inspection", showWarnings = F)
lapply(h$site, function(x){
  message(x)
  p1 <- ggplot() +
    geom_sf(data = watersheds[watersheds$site == x, ],
            aes(color = x),
            show.legend = FALSE) +
    geom_sf(data = catchments[catchments$site == x, ],
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
