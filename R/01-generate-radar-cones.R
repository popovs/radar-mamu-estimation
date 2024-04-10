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

# Create cones out of them with 1/4 circle wedges
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
  x$site <- h[i,][["site"]]
  x <- x[x$lyr.1 == 1, ]
  x <- x[,"site"] # drop 'lyr.1' column
  if (st_geometry_type(x) == "MULTIPOLYGON") x <- st_convex_hull(x) # merge broken pixels together into single polygon if multipart
  # TODO:: st_convex_hull should include the original station coordinate
  cones[[i]] <- x
  # if (i == 1) {
  #   catchments <- x
  # } else {
  #   catchments <- rbind(catchments, x)
  # }
}

cones <- dplyr::bind_rows(cones)

# Smooth out jagged edges
library(smoothr)
cones <- smooth(cones, method = "chaikin")

cones$area_ha <- units::set_units(sf::st_area(cones), "ha") 
mean(cones$area_ha)
#hist(cones$area_ha, breaks = 100)
plot(cones[1]) # Visually inspect

gc()

# 03 CLIP watersheds? TO CATCHMENTS ---------------------------------------------

# 03-1 Intersect cones with watersheds ----

# Evidence indicates that MAMU generally fly within a given
# watershed when flying to a nesting site. 
# TODO: CITE!!!!!\
# TODO: if this works, move to GIS folder!
watersheds <- st_read("../Watersheds/conservation_region_watersheds.gpkg")

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

# Inspect each one...
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
watersheds <- watersheds %>% 
  select(site) %>%
  group_by(site) %>%
  dplyr::summarize()

ggplot() +
  geom_sf(data = watersheds,
          aes(color = site),
          show.legend = FALSE) +
  geom_sf(data = cones,
          fill = NA) +
  geom_sf(data = h)



# 03-2 Mask to accessible areas ----

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


# 03-3 Select correct watershed ----

# The catchments are now essentially chopped up cones divided
# up by either the non-nesting zones from `mnh` or into 
# multiple watersheds. Now we need to choose the correct
# piece of the catchment that relates to our radar station and
# discard the rest of the miscellaneous values.

# Next, choose whichever watershed ID falls in line/intersects 
# with the heading
# https://stackoverflow.com/questions/69820690/finding-coordinates-from-heading-and-distance-in-r

h2 <- st_transform(h, crs = 4326)
h2$lon.1 <- NA
h2$lat.1 <- NA

for (i in 1:nrow(h2)) {
  ic <- st_coordinates(h2[i,])
  b <- sf::st_drop_geometry(h2[i, "inc_out2"])[[1]]
  nc <- geosphere::destPoint(ic, b = b, d = 30000)
  h2[i, "lon.1"] <- nc[1]
  h2[i, "lat.1"] <- nc[2]
}
rm(ic, b, nc, i)

h2 <- st_drop_geometry(h2)
h2$geom <- sprintf("LINESTRING(%s %s, %s %s)",
                   h2$lon, h2$lat,
                   h2$lon.1, h2$lat.1)
h2 <- st_as_sf(h2, wkt = "geom")
st_crs(h2) <- 4326
h2 <- st_transform(h2, 3005)


# Check it out...
x <- "Southgate"
x <- "Aaltanhash"

ggplot() +
  geom_sf(data = catchments[catchments$site == x,],
          aes(color = WSD_ID)) +
  geom_sf(data = h[h$site == x,]) +
  geom_sf(data = h2[h2$site == x,])

ggplot() +
  geom_sf(data = catchments[catchments$site == x,][st_intersects(catchments[catchments$site == x,], h2[h2$site == x,], sparse = F),],
          aes(color = WSD_ID)) +
  geom_sf(data = h[h$site == x,]) +
  geom_sf(data = h2[h2$site == x,])

# `h2` is now our dataframe of headings (lines that are 30km 
# long oriented in the direction of each station's mean
# heading). 
# The intersection between each line and the catchment pieces
# it touches are the 'true' catchment. All other pieces can be
# dropped.
# Now cycle through each catchment and extract out the pieces
# that intersect with it's corresponding `h2` line.


catchments$area_ha <- units::set_units(sf::st_area(catchments), "ha")
catchments <- catchments[order(catchments$site, catchments$area_ha),]
catchments$keep_yn <- NA

for (i in 1:nrow(h2)) {
  site <- h2[["site"]][i]
  message("Cleaning up ", site, "...")
  catchments[["keep_yn"]][catchments$site == site] <- st_intersects(catchments[catchments$site == site,], h2[h2$site == site,], sparse = F)
}

# 
# keep_cat_yn <- lapply(h2$site, function(x) {
#   st_intersects(catchments[catchments$site == x,], h2[h2$site == x,], sparse = F)
# })
# catchments$keep_yn <- unlist(keep_cat_yn)

# Check it out...

x <- "Aaltanhash"
x <- "Southgate"
x <- "Watta"

ggplot() +
  geom_sf(data = catchments[catchments$site == x,][st_intersects(catchments[catchments$site == x,], h2[h2$site == x,], sparse = F),],
          aes(color = WSD_ID)) +
  geom_sf(data = h[h$site == x,]) +
  geom_sf(data = h2[h2$site == x,]) +
  ggtitle(x)

ggplot() +
  geom_sf(data = catchments[catchments$site == x & catchments$keep_yn == TRUE,]) +
  geom_sf(data = h[h$site == x,]) +
  ggtitle(x)


# Remove non-intersecting catchment pieces
catchments <- catchments[catchments$keep_yn == TRUE, ]

# 03-4 Clean up catchments ----

# There are still some catchments with multiple watershed IDs
# per radar station. The final task here will be to choose the
# watershed ID that is geographically closest to the radar
# station and toss any extraneous ones.
catchments$keep_yn <- NA
for (i in 1:nrow(h2)) {
  site <- h2[["site"]][i]
  message("Selecting catchment pieces closest to ", site, "...")
  tmp <- catchments[catchments$site == site, ]
  keep_index <- st_nearest_feature(h[h$site == site,], tmp)
  catchments[["keep_yn"]][catchments$site == site][keep_index] <- TRUE
}


# Check it out...

x <- "Aaltanhash"
x <- "Southgate"
x <- "Watta"

ggplot() +
  geom_sf(data = catchments[catchments$site == x & catchments$keep_yn == TRUE,]) +
  geom_sf(data = h[h$site == x,]) +
  ggtitle(x)

# Before tossing the 'non-closest' catchment pieces, let's check
# that it's not physically touching the closest piece. If it's
# touching the closest piece, we can safely assume a MAMU is 
# capable of flying through to the next watershed piece (e.g.,
# Watta!). Otherwise, if the other watershed pieces are *not* 
# touching our closest 'index' piece, toss them.


# Now, you can finally chuck the off ones.
catchments <- catchments[which(catchments$keep_yn == TRUE), ]

# There should now be 95 features, matching the original 95 headings...
nrow(h) == nrow(catchments)

# 04 EXTRACT CATCHMENT VALUES ---------------------------------------------


# Recalculate catchment area
catchments$area_ha <- units::set_units(sf::st_area(catchments), "ha")


# Extract mean cost and mean forest cover per catchment
catchments$mean_cost <- exactextractr::exact_extract(cost, catchments, "mean")
catchments$mean_forest <- exactextractr::exact_extract(forest, catchments, "mean")

# Save em
if (!(any(grepl("radar_derived_catchments.gpkg", list.files("GIS"))))|overwrite == TRUE) st_write(catchments, "GIS/radar_derived_catchments.gpkg", append = FALSE)


# 03 CLEAN UP -------------------------------------------------------------

if (!(any(grepl("MAMU_nesting_habitat.gpkg", list.files("GIS"))))|overwrite == TRUE) st_write(mnh, "GIS/MAMU_nesting_habitat.gpkg", append = FALSE)

rm(h, i, isNaN)
