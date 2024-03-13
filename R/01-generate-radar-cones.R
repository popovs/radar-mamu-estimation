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

# 01 SETUP ----------------------------------------------------------------

#library(devtools)
#devtools::install_github("popovs/MAMU")
library(MAMU) # MUST be >= 0.2.1
library(dplyr)
library(sf)

dir.create("temp", showWarnings = FALSE)

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

# Clean up station names as needed
stn_lookup <- MAMU::rs
headings <- merge(headings, stn_lookup, by.x = "name", by.y = "original_name", all.x = TRUE)
headings$name <- ifelse(is.na(headings$new_name), headings$name, headings$new_name)

# For each station, take the mean INCOMING & OUTGOING heading.
inc <- headings[grep("Incoming", headings$flightpath_type, ignore.case = T),c("name", "heading")]
out <- headings[grep("Outgoing", headings$flightpath_type, ignore.case = T),c("name", "heading")]

inc <- inc[complete.cases(inc),]
out <- out[complete.cases(out),]

inc <- aggregate(heading ~ name, inc, mean)
out <- aggregate(heading ~ name, out, mean)

# Make `out` the opposite heading
out2 <- out
out2$heading <- out2$heading - 180
out2$heading <- ifelse(out2$heading < 0, out2$heading + 360, out2$heading)

h <- merge(inc, out, by = "name")
h <- merge(h, out2, by = "name")
names(h) <- c("name", "inc", "out", "out2")

h$inc_out2 <- rowMeans(h[,c("inc", "out2")])

h <- merge(h, stn, by.x = "name", by.y = "site")
names(h)[1] <- "site"

h <- st_as_sf(h) %>%
  st_transform(3005)

rm(headings, inc, out, out2, stn_lookup)

# 02 GENERATE CONES -------------------------------------------------------

# Create cones out of them with 1/4 circle wedges
# Note your points MUST be in a projected coordinate system in meters!!
cones <- lapply(1:nrow(h), function(x){
  message("Generating cone for ", h$site[[x]], "...")
  radar_cone(pt = h[x,],
             radius = 30000,
             theta = 90,
             heading = h$inc_out2[[x]],
             res = 25)
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

catchments <- list()
for (i in 1:nrow(h)) {
  message("Vectorizing ", h[i,][["site"]])
  x <- cones[[i]]
  x <- x == 1
  x <- terra::as.polygons(x, crs = "epsg:3005")
  x <- sf::st_as_sf(x)
  x$site <- h[i,][["site"]]
  if (i == 1) {
    catchments <- x
  } else {
    catchments <- rbind(catchments, x)
  }
}

catchments <- catchments[catchments$lyr.1 == 1,]
catchments <- catchments[,c("site", "geometry")]

catchments$area_ha <- units::set_units(sf::st_area(catchments), "ha") # Check they all have the same full size area
plot(catchments[1]) # Visually inspect

rm(cones)

# 02-2 Mask to nesting areas ----

# Load up viable nesting habitat tiff
# This time using stars - far more efficient to vectorize
# a stars object
mnh <- stars::read_stars("GIS/MAMU_nesting_habitat.tiff")

# Polygonize
# This will take ~5 minutes
mnh <- st_as_sf(mnh, 
                as_points = FALSE, 
                merge = TRUE, 
                na.rm = TRUE)
mnh <- st_union(mnh) # merge into one polygon

# Simplify this giant mnh polygon
library(smoothr)
mnh <- smooth(mnh, method = "chaikin")

st_write(mnh, "GIS/MAMU_nesting_habitat.gpkg")

# This will take just under ~1 minute
catchments <- st_intersection(catchments, mnh)

# Recalculate catchment area
catchments$area_ha <- units::set_units(sf::st_area(catchments), "ha")

plot(catchments[1])

# 02-3 Save em ----
st_write(catchments, "GIS/radar_derived_catchments.gpkg")


# 03 CLEAN UP -------------------------------------------------------------

rm(h, i, isNaN)
