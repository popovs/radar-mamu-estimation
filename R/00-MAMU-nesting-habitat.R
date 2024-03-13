# 0. MAMU NESTING HABITAT AREA

# This script will download the BC Digital Elevation Model
# (DEM) for our study region, then clip it down to the areas
# that are highly likely to contain suitable marbled murrelet
# nesting habitat area (regions that are 30 km from shore and
# under <1500 m elevation). The exact maximum elevation cutoff
# varies by our six conservation regions of interest (North, 
# Central, and South Mainland Coast, Haida Gwaii, North + West 
# Vancouver Island, and East Vancouver Island).

# TODO: update with actual time 
# DISCLAIMER 1: THIS WHOLE SCRIPT WILL TAKE ABOUT AN HOUR TO RUN.
# There are likely faster or more efficient ways of doing this,
# but the primary goal of this script is reproducibility and
# ease of understanding.

# Disclaimer 2: this was run with 16 GB ram and will likely 
# fail with less - that said, many of these steps can be
# easily replicated in QGIS to the same effect. The benefit of
# this code is the exact reproducibility.

# Create temporary directory to store large scratch files
# This makes it a bit easier on the RStudio memory
dir.create("temp")

# 01 DOWNLOAD DEM ---------------------------------------------------------

# Prior to running any sort of analysis, we first need to
# download the BC Digital Elevation Model (DEM) dataset.
# This script downloads the DEM data and then clips it to 
# our study region. 

# The results of this script are stored in the `GIS` folder. 
# Because the DEM rasters are so large, they are not commited 
# to this repository. Therefore, if you are cloning this repository 
# to replicate this analysis on your own machine, you will need 
# to run this script once.

# We will make use of the BC_DEM function in the MAMU library to
# download this DEM.
# https://popovs.github.io/MAMU/articles/BC_DEM.html

#library(devtools)
#devtools::install_github("popovs/MAMU")
library(MAMU)

tiles_to_download <- c("103k", "103j", "103f", "103g", "103c", "103b", "102o", # haida gwaii
                       "103o", "103p", "103j", "103i", "103g", "103h", "103a", "93e", "93d", "93c", "93l", "93m", "102p", "92m", "92n", "104a", "104b", # north/central coast
                       "102i", "92l", "92k", "92e", "92f", "92c", "92b", # vancouver island
                       "92j", "92g", "92h" # lower mainland
                       )
tiles_to_download <- tolower(tiles_to_download)
tiles_to_download <- unique(tiles_to_download)

# Create a new directory, `GIS/DEM/DEM_tiles`, where the dem data will be stored
dir.create("GIS/DEM/DEM_tiles", recursive = TRUE)

# Download the tiles one at a time and save to the DEM folder
# (R will crash if you try to download your DEM's all in one go
# from within the BC_DEM function - so do it in a loop).
# This will take approximately 15 minutes to complete.
sapply(tiles_to_download, 
       BC_DEM,
       save_output = TRUE,
       overwrite = FALSE,
       output_dir = "GIS/DEM/DEM_tiles/")
beepr::beep()

rm(tiles_to_download)

# Create virtual raster layer of all these tiles
# This will throw a warning for all the xml files - safe to ignore.
vrt <- make_vrt(path = "GIS/DEM/DEM_tiles/", filename = "GIS/DEM/BC_DEM_VRT.vrt")

# Reproject the DEM VRT to BC Albers (EPSG 3005)
# This will take approximately 30 minutes
# `stars` package may be faster, though I have had issues getting it to work
dem_3005 <- terra::project(vrt, "EPSG:3005")
beepr::beep()

# Save your hard-earned DEM
# This will take approximately 1-2 minutes
terra::writeRaster(dem_3005, "GIS/DEM/BC_DEM_EPSG3005.tiff")
beepr::beep()



# 02 REGIONAL ELEVATION CUTOFFS -------------------------------------------

# The DEM now needs to be divided into conservation regions 
# for the elevational cutoff step.

# Read in the conservation regions shapefile
library(sf)
cons_reg <- st_read("GIS/cons_reg.shp")
cons_reg <- st_transform(cons_reg, 3005) # Set BC Albers projection

# 02-1 Intersect DEM with each conservation region ----

# The `stars` package renders images much quicker within R, but
# I have found `terra::writeRaster` to be much faster than
# `stars::write_stars` for actually saving the rasters to disk.
# This will take approximately 10 minutes.
dir.create("GIS/DEM/Regional_DEM")
library(terra)
for (i in 1:nrow(cons_reg)) {
  message("Cropping and saving ", cons_reg$MMCR_NA[i])
  tmp <- terra::crop(dem_3005, cons_reg[i,])
  tmp <- terra::mask(tmp, cons_reg[i,]) # This can be done in one step with terra::crop(mask = T), but was resulting in buggy raster values - so doing it in two steps here.
  filename <- paste0(file.path("GIS/DEM/Regional_DEM", gsub('\\b(\\pL)\\pL{3,}|.','\\U\\1', cons_reg$MMCR_NA[i], perl = TRUE)), "_3005.tiff")
  terra::writeRaster(tmp, filename)
  gc() # Clean up RAM & misc session garbo
}
beepr::beep()
rm(i, tmp, filename)

# To quickly visualize
#plot(stars::read_stars("GIS/DEM/Regional_DEM/NMC_3005.tiff")) # etc.

# 02-2 Apply elevation cutoffs ----

# The exact elevational cutoffs vary by region and are primarily
# dictated by the growing conditions of the trees that MAMU
# prefer to nest in. Large coniferous trees that can support
# MAMU nests can grow in higher elevations at lower latitudes.
# Conversely, the further north you go, the lower the maximum
# MAMU nesting elevation.

# | REGION                | ELEVATION |
# +-----------------------+-----------+
# | Central & Northern MC |  <  800m  |
# | Haida Gwaii           |  <  700m  |
# | East VI               |  < 1100m  |
# | West & North VI       |  < 1200m  |
# | South MC              |  < 1500m  |

# Even though Alaska is not being used in this present analysis,
# I've assigned it a cutoff of 800m in line with Northern MC so 
# the file is handy in the future if needed.

elevation_cutoffs <- setNames(data.frame(cons_reg$MMCR_NA,
                                         gsub('\\b(\\pL)\\pL{3,}|.','\\U\\1', cons_reg$MMCR_NA, perl = TRUE),
                                         c(800, 700, 800, 1500, 1200, 1100, 800)),
                              c("region", "abbreviation", "elevation_m"))


dir.create("GIS/Elevation_cutoffs")

# This will take approximately 5 minutes
for (i in 1:nrow(elevation_cutoffs)) {
  message("Reclassifying ", elevation_cutoffs$region[i])
  tmp <- terra::rast(paste0("GIS/DEM/Regional_DEM/", elevation_cutoffs$abbreviation[i], "_3005.tiff"))
  tmp <- terra::ifel(tmp < elevation_cutoffs$elevation_m[i] & tmp > 0, 1, NA)
  filename <- file.path("GIS/Elevation_cutoffs", paste0(elevation_cutoffs$abbreviation[i], "_", elevation_cutoffs$elevation_m[i], "m.tiff"))
  terra::writeRaster(tmp, filename)
}
beepr::beep()
rm(i, tmp, filename, elevation_cutoffs)


# 02-3 Merge elevation cutoff rasters ----

# Now we can merge all these rasters into one single raster
# for ease of use.

# Read them into environment
f <- list.files("GIS/Elevation_cutoffs", full.names = T)
ec <- lapply(f, terra::rast) # 'ec' for elevational cutoffs
names(ec) <- basename(f)

# Mosaic all the ec rasters together
# This will take approximately 1-2 minutes
ec <- make_vrt("GIS/Elevation_cutoffs/",
               filename = "GIS/Elevation_cutoffs/elevation_cutoffs.vrt")
ec_collection <- terra::sprc(ec)
elev <- terra::mosaic(ec_collection,
                      filename = "GIS/Elevation_cutoffs/elevation_cutoffs.tiff")
terra::plot(elev)

rm(ec, ec_collection, f)


# 03 DISTANCE FROM COAST CUTOFF -------------------------------------------

# Next we will create a separate raster of all points wit
# 30 km distance from the coast, as BC nest survey data
# indicates that 99% of MAMU nests are within 30 km of the
# coastline. We will then intersect the two rasters to get
# our suitable MAMU habitat raster.

# 03-1 Create blank canvas of study area ----

# It's faster to calculate distance from coast using a blank
# canvas raster rather than reclassifying the thousands of 
# cell values in the DEM. Additionally, we're going to make
# it a courser resolution by a factor of 10 here - so each
# pixel will be 250m x 250m rather than 25m x 25m - otherwise,
# many of the subsequent raster calculations will fail or 
# take far far too long. This means that we will be calculating
# distance from the coast to a 1/4 km accuracy.
canvas <- terra::rast(extent = terra::ext(elev), 
                      nrow = nrow(elev)/10, 
                      ncol = ncol(elev)/10, 
                      nlyr = 1)
terra::values(canvas) <- 1
terra::crs(canvas) <- "epsg:3005"

# Clip canvas to cons_reg extent to cut down on memory
canvas <- terra::crop(canvas, cons_reg, mask = TRUE)

# Read in coastline (high water mark)
coast <- sf::st_read("GIS/CHS_HWM_S_line.shp")
sf::st_crs(coast) # check projection

# Extract coastline points
xy <- sf::st_coordinates(sf::st_cast(coast, "POINT"))

rm(coast)

# Change all coastline values of the raster to '2'
# Extract out all the cells in `canvas` that contain our coastline
# (x,y) coordinates, then change the value of those cells to '2'
# This will take ~3 minutes
rcoast <- terra::cellFromXY(canvas, xy)
rcoast <- unique(rcoast)
rcoast <- rcoast[!is.na(rcoast)]
canvas[rcoast] <- 2

rm(rcoast, xy)

# If you visualize it, it looks like not all the coastline is '2',
# but that's more of a rendering issue - if you zoom in the coastline
# is all there.

# 03-2 Pull out mainland from canvas ----

# Next, we're going to split apart the islands from the mainland. 
# The mainland is the only region we're applying the 30km buffer to, 
# whereas we know with confidence MAMU are flying into the most interior 
# ~35-40km bits of Vancouver Island (and Haida Gwaii is <30 km across
# anyway). We are going to include the most interior bits of the islands
# in our analysis and *not* apply the 30 km buffer rule to them.

# Weird terra bug where you sometimes have to separately mask the object
# or else it messes with the raster values?
mainland <- terra::crop(canvas, 
                        cons_reg[!(cons_reg$MMCR_NA %in% c("West and North Vancouver Island", "East Vancouver Island", "Haida Gwaii")),])
mainland <- terra::mask(mainland,
                        cons_reg[!(cons_reg$MMCR_NA %in% c("West and North Vancouver Island", "East Vancouver Island", "Haida Gwaii")),])

rm(canvas)
gc()

# For the islands themselves, we're just going to take the dem, pull 
# out the island, and select any cells with an elevation greater than 0.
# (See section 03-4 below.)

# 03-3 Calculate mainland coast distance ----

# Now compute the distance of every cell != 2 to every cell == 2
# Note that this is a very simple distance calculation - we are
# assuming the birds can fly a straight-line path in any direction
# from the coast. (To assume barriers, you can use the 'costDist'
# function.)
distance <- terra::gridDist(mainland, target = 2)

# Extract out only distances <= 30km
distance <- terra::ifel(distance <= 30000, 1, NA) # projection is in meters

# Re-aggregate to the same resolution as original DEM so we can
# later extract out any land area plus + merge with the islands
# This will take 1-2 minutes
distance <- terra::resample(distance, dem_3005)

# Crop the distance to land areas only
# (This step is purely cosmetic and can be skipped)
# Takes approximately 5 minutes
distance <- (dem_3005 > 0) * distance
distance <- terra::ifel(distance == 1, 1, NA)
beepr::beep()

# 03-4 Pull out islands from DEM ----

# This will take approximately 10 minutes
islands <- terra::crop(dem_3005,
                       cons_reg[(cons_reg$MMCR_NA %in% c("West and North Vancouver Island", "East Vancouver Island", "Haida Gwaii")),])
islands <- terra::mask(dem_3005,
                       cons_reg[(cons_reg$MMCR_NA %in% c("West and North Vancouver Island", "East Vancouver Island", "Haida Gwaii")),])
islands <- terra::ifel(islands > 0, 1, NA)
beepr::beep()

# Resample islands to match extent of the DEM, so we 
# can merge with `distance`
# This will take approximately 2 minutes
islands <- terra::resample(islands, dem_3005)
beepr::beep()

# 03-5 Merge mainland and islands ----

# First saving `distance` and `islands` to temp folder
# in case anything goes awry
terra::writeRaster(distance, "temp/distance.tiff")
terra::writeRaster(islands, "temp/islands.tiff")
beepr::beep()

gc()

# This will take approximately 10 minutes
dir.create("GIS/Distance_cutoffs")
distance <- sum(distance, islands, na.rm = TRUE)
distance <- terra::resample(distance, elev) # set to same spatial extent as `elev`
beepr::beep()

terra::writeRaster(distance, "GIS/Distance_cutoffs/distance_cutoffs.tiff")
beepr::beep()

rm(islands)


# 04 MERGE ELEVATION AND COAST DISTANCE CUTOFFS ---------------------------

# Finally, merge the coast distance + elevation cutoff rasters
# together to create a "MAMU containment zone" area that has a
# high probability of containing high quality MAMU nesting 
# habitat. This raster will be used as our maximum MAMU area
# within BC that we will extrapolate our population estimates
# to.

# `mnh` for "MAMU nesting habitat"
# This will take approximately 3-4 minutes
mnh <- elev + distance

# 04-1 Adjust raster resolution ----

# Per the BC DEM website, this raster product is at a 25m 
# resolution, but gridded to a 0.75 arc-second scale. 
# https://www2.gov.bc.ca/gov/content/data/geographic-data-services/topographic-data/elevation/digital-elevation-model
# This means all our raster data is at a somewhat odd 
# 17.37227 x 17.37227 resolution:
res(dem_3005)
res(elev)
res(mnh)

# To make calculations easier down the line, we will resample
# the final `mnh` resolution to 25m x 25m grid cells.
# This will take 1-2 minutes.
r <- mnh
terra::res(r) <- 25
mnh <- terra::resample(mnh, r)
beepr::beep()

res(mnh)
rm(r)

# 04-2 Calculate habitat area within each conservation region ----

# This serves to both update our conservation region shapefile
# with the correct habitat areas, but also to double check that 
# our numbers roughly line up with what we expect - the Haida 
# Gwaii habitat area should be just under 10k km squared.
cons_reg$mnh_cell_count <- exactextractr::exact_extract(mnh, cons_reg, 'count') # count the number of cells within each cons region
cons_reg$mnh_m <- cons_reg$mnh_cell_count * 625 # each raster cell is 25m x 25m, aka 625 meters squared
cons_reg$mnh_km <- cons_reg$mnh_m / 1000000 # go from m2 to km2
cons_reg$mnh_ha <- cons_reg$mnh_km / 100 # divide by 100 to go from km2 to hectares

# Checks out nicely!
cons_reg[["mnh_km"]][cons_reg$MMCR_NA == "Haida Gwaii"] # HG should be just under 10k
sum(cons_reg[["mnh_km"]][cons_reg$MMCR_NA %in% c("West and North Vancouver Island", "East Vancouver Island")]) # Vancouver Island (incl. Nootka + Quadra + Gulf islands) should be ~32k

# 04-3 Save MAMU nesting habitat ----

terra::writeRaster(mnh, "GIS/MAMU_nesting_habitat.tiff")
st_write(cons_reg, "GIS/cons_reg.shp", append = FALSE) # overwrite the conservation region shapefile with a shp containing the habitat area numbers


# 05 CLEAN UP -------------------------------------------------------------

# Clean up any temp files
unlink("temp", recursive = TRUE)
