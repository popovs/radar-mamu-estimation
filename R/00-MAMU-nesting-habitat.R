# 0. MAMU NESTING HABITAT AREA

# This script will download the BC Digital Elevation Model
# (DEM) for our study region, then clip it down to the areas
# that are highly likely to contain suitable marbled murrelet
# nesting habitat area (regions that are 30 km from shore and
# under <1500 m elevation). The exact maximum elevation cutoff
# varies by our six conservation regions of interest (North, 
# Central, and South Mainland Coast, Haida Gwaii, North + West 
# Vancouver Island, and East Vancouver Island).

# DISCLAIMER 1: THIS WHOLE SCRIPT WILL TAKE ABOUT TWO HOURS TO RUN.
# There are likely faster or more efficient ways of doing this,
# but the primary goal of this script is reproducibility and
# ease of understanding. If you wish to simply test the script, 
# modify the raster resolutions to a courser scale (e.g., 500m) 
# to run through the script quickly.

# DISCLAIMER 2: this was run with 16 GB ram and will likely 
# fail with less - that said, many of these steps can be
# easily replicated in QGIS to the same effect. The benefit of
# this code is the exact reproducibility.

# DISCLAIMER 3: If you are running this script on MacOS, there
# is a known bug with the default MacOS graphics device that 
# makes plotting really reaaallly slow. Install `ragg` and 
# change the graphics device to "AGG" in Options > General >
# Graphics. See here:
# https://gis.stackexchange.com/questions/424897/why-is-sf-struggling-to-run-plots-quickly-on-uk-datasets?rq=1


# 00 SETUP ----------------------------------------------------------------

# Create temporary directory to store large scratch files
# This makes it a bit easier on the RStudio memory
dir.create("temp", showWarnings = FALSE)

# Create `res` variable - the resolution in meters of all 
# GIS calculations. The max resolution for the dataset is
# 25m x 25m, but the runtime for the script will be 3+ hours.
res <- 500 # res MUST be >25, as the DEM goes does to 25m accuracy.
stopifnot("`res` must be >=25, as the DEM goes down to 25m accuracy." = res > 25)
overwrite <- TRUE

#library(devtools)
#devtools::install_github("popovs/MAMU")
library(MAMU)

# Read in the conservation regions shapefile
library(sf)
cons_reg <- st_read("GIS/cons_reg.shp")
cons_reg <- st_transform(cons_reg, 3005) # Set BC Albers projection

# Read in nest data
nests <- st_read("GIS/MAMU_nests.gpkg")

# TODO: replace all the "if length files in DEM folder == 0
# then run this code" w "if overwrite == TRUE" piece above -
# this will make it targets-friendly

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

# Create a new directory, `GIS/DEM/DEM_tiles`, where the dem data will be stored
dir.create("GIS/DEM/DEM_tiles", recursive = TRUE, showWarnings = FALSE)

# And finally, ONLY RUN THIS CODE IF THE DEM DOESN'T EXIST!
if (!("BC_DEM_EPSG3005.tiff" %in% list.files("GIS/DEM"))) {
  
  tiles_to_download <- c("103k", "103j", "103f", "103g", "103c", "103b", "102o", # haida gwaii
                         "103o", "103p", "103j", "103i", "103g", "103h", "103a", "93e", "93d", "93c", "93l", "93m", "102p", "92m", "92n", "104a", "104b", # north/central coast
                         "102i", "92l", "92k", "92e", "92f", "92c", "92b", # vancouver island
                         "92j", "92g", "92h" # lower mainland
                         )
  tiles_to_download <- tolower(tiles_to_download)
  tiles_to_download <- unique(tiles_to_download)
  
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

} else {
  dem_3005 <- terra::rast("GIS/DEM/BC_DEM_EPSG3005.tiff")
  # Resample DEM to match the resolution specified above
  r <- dem_3005
  terra::res(r) <- res
  dem_3005 <- terra::resample(dem_3005, r)
  rm(r)
}



# 02 REGIONAL ELEVATION CUTOFFS -------------------------------------------

# The DEM now needs to be divided into conservation regions 
# for the elevational cutoff step.

# 02-1 Intersect DEM with each conservation region ----

# The `stars` package renders images much quicker within R, but
# I have found `terra::writeRaster` to be much faster than
# `stars::write_stars` for actually saving the rasters to disk.
# This will take approximately 10 minutes.
dir.create("GIS/DEM/Regional_DEM", showWarnings = FALSE)
library(terra)

# ONLY RUN THIS IF FILES DO NOT EXIST YET
if (length(list.files("GIS/DEM/Regional_DEM")) == 0) {
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
}

# To quickly visualize
#plot(stars::read_stars("GIS/DEM/Regional_DEM/NMC_3005.tiff")) # etc.

# 02-2 Apply elevation cutoffs ----

# The exact elevational cutoffs vary by region and are primarily
# dictated by the growing conditions of the trees that MAMU
# prefer to nest in. Large coniferous trees that can support
# MAMU nests can grow in higher elevations at lower latitudes.
# Conversely, the further north you go, the lower the maximum
# MAMU nesting elevation. The nest data is the main source of 
# cutoffs. Adding more nest data will trigger a re-run of this
# analysis.

# Extract nest elevations from DEM
nests$elev_m <- terra::extract(dem_3005, nests)[[2]]

ggplot2::ggplot(data = nests, 
                ggplot2::aes(x = cons_reg, 
                             y = elev_m)) + 
  ggplot2::geom_boxplot() + 
  ggplot2::geom_jitter() +
  ggplot2::theme_minimal()

# | REGION                | ELEVATION* |
# +-----------------------+------------+
# | Central & Northern MC |  <  730m   |
# | Haida Gwaii           |  <  610m   |
# | East VI               |  < 1050m   |
# | West & North VI       |  < 1050m   |
# | South MC              |  < 1220m   |
# * Derived from 95 percentile nest elevations by region

# Even though Alaska is not being used in this present analysis,
# I've assigned it a cutoff of 800m in line with Northern MC so 
# the file is handy in the future if needed.

elevation_cutoffs <- setNames(data.frame(cons_reg$MMCR_NA,
                                         gsub('\\b(\\pL)\\pL{3,}|.','\\U\\1', cons_reg$MMCR_NA, perl = TRUE),
                                         #c(700, 600, 700, 1200, 1000, 1000, 800)
                                         unlist(lapply(cons_reg$MMCR_NA, function(x){quantile(nests[["elev_m"]][nests$cons_reg == x], 0.95, na.rm = TRUE)[[1]]}))
                                         ),
                              c("region", "abbreviation", "elevation_m"))
# If northern mainland coast and alaska border are NA, use CMC elevation cutoff
if (is.na(elevation_cutoffs[["elevation_m"]][elevation_cutoffs$abbreviation == "NMC"])) elevation_cutoffs[["elevation_m"]][elevation_cutoffs$abbreviation == "NMC"] <- elevation_cutoffs[["elevation_m"]][elevation_cutoffs$abbreviation == "CMC"]
if (is.na(elevation_cutoffs[["elevation_m"]][elevation_cutoffs$abbreviation == "AB"])) elevation_cutoffs[["elevation_m"]][elevation_cutoffs$abbreviation == "AB"] <- elevation_cutoffs[["elevation_m"]][elevation_cutoffs$abbreviation == "CMC"]

# Round to nearest 10 meter
elevation_cutoffs$elevation_m <- round(elevation_cutoffs$elevation_m, -1)

if (overwrite == TRUE) unlink("GIS/Elevation_cutoffs", recursive = TRUE)
dir.create("GIS/Elevation_cutoffs", showWarnings = F)

# This will take approximately 5 minutes
# RUN IF FILES DON'T EXIST OR OVERWRITE == TRUE AT TOP OF SCRIPT
if (length(list.files("GIS/Elevation_cutoffs/")) == 0|overwrite == TRUE){
  for (i in 1:nrow(elevation_cutoffs)) {
    message("Reclassifying ", elevation_cutoffs$region[i])
    tmp <- terra::rast(paste0("GIS/DEM/Regional_DEM/", elevation_cutoffs$abbreviation[i], "_3005.tiff"))
    # TODO: change elevation cutoffs to the `res` resolution
    tmp <- terra::ifel(tmp < elevation_cutoffs$elevation_m[i] & tmp > 0, 1, NA)
    filename <- file.path("GIS/Elevation_cutoffs", paste0(elevation_cutoffs$abbreviation[i], "_", elevation_cutoffs$elevation_m[i], "m.tiff"))
    terra::writeRaster(tmp, filename, overwrite = overwrite)
  }
  beepr::beep()
  rm(i, tmp, filename, elevation_cutoffs)
} 


# 02-3 Merge elevation cutoff rasters ----

# Now we can merge all these rasters into one single raster
# for ease of use.

# Read them into environment
f <- list.files("GIS/Elevation_cutoffs", full.names = T)
ec <- lapply(f, terra::rast) # 'ec' for elevational cutoffs
names(ec) <- basename(f)

# Mosaic all the ec rasters together
# This will take approximately 1-2 minutes
# TODO: wrap in 'if overwrite = TRUE' situation
ec <- make_vrt("GIS/Elevation_cutoffs/",
               filename = "GIS/Elevation_cutoffs/elevation_cutoffs.vrt",
               overwrite = overwrite)
ec_collection <- terra::sprc(ec)
elev <- terra::mosaic(ec_collection,
                      filename = "GIS/Elevation_cutoffs/elevation_cutoffs.tiff",
                      overwrite = overwrite)
list.files("GIS/Elevation_cutoffs/")
terra::plot(elev)

# TODO: ELSE, if overwrite = F, read in elev and set it to the
# resolution set above

rm(ec, ec_collection, f)


# 03 DISTANCE FROM COAST CUTOFF -------------------------------------------

# Next we will create a separate raster of all points wit
# 30 km distance from the coast, as BC nest survey data
# indicates that 99% of MAMU nests are within 30 km of the
# coastline. We are going to assume a 30 km distance from 
# shore flying around mountain barriers  but allowing for 
#flight over water.

# 03-1 Get USA land areas ----

# The DEM data doesn't include the USA or inland BC.
# These land areas need to be blocked off so the distance
# algorithm can differentiate between land areas and
# water areas. 

# Pull Alaska + Washington state land area from Natural Earth
library(rnaturalearth)
usa <- ne_states("united states of america")
usa <- usa[usa$name %in% c("Alaska", "Washington"), ]
usa <- st_as_sf(usa)
usa <- st_transform(usa, crs = 3005)

extent <- terra::ext(dem_3005)
usa <- st_crop(usa, extent)

# 03-2 Block of inland BC areas ----

# Rather than block off inland BC with a rough BC-shaped
# polygon from Natural Earth, we're just going to use crude
# rectangles. The reason for this is the DEM coastline is far
# more detailed than the course NE polygon, and we would risk
# losing coastline data if we masked the DEM with such a 
# course polygon. So rectangles it is.
# Specify a few rectangles w coordinates from bottom left corner clockwise to bottom right corner
x1 <- rbind(c(1280841, extent[3]), c(1280841, extent[4]), c(extent[2], extent[4]), c(extent[2], extent[3]))
x2 <- rbind(c(1109405, 651763), c(1109405, extent[4]), c(extent[2], extent[4]), c(extent[2], 651763))
x3 <- rbind(c(979227, 877974), c(979227, extent[4]), c(extent[2], extent[4]), c(extent[2], 877974))
x4 <- rbind(c(827708, 1154691), c(827708, extent[4]), c(extent[2], extent[4]), c(extent[2], 1154691))
x5 <- rbind(c(610000, 1300000), c(610000, extent[4]), c(extent[2], extent[4]), c(extent[2], 1300000))
# lol this is terribly inefficient but it works
land <- list(x1, x2, x3, x4, x5)
land <- lapply(land, terra::vect, type = "polygons", crs = "epsg:3005")
land1 <- terra::union(land[[1]], land[[2]])
land2 <- terra::union(land1, land[[3]])
land3 <- terra::union(land2, land[[4]])
land4 <- terra::union(land3, land[[5]])
land <- terra::aggregate(land4)
land <- terra::aggregate(terra::union(terra::vect(usa), land))
terra::plot(land)
rm(land1, land2, land3, land4, x1, x2, x3, x4, x5)


# 03-3 Create canvas of target area ----

# Now we need to combine our land data with the `elev` data to
# create the raster with which we are going to do distance
# calculations with. The `elev` data will cut off barriers
# that the algorithm will have to travel around when 
# calculating distance, while the `canvas` area will supply
# land and sea data.
#   - Mountain barriers -> must travel around
#   - Land areas -> traveling through
#   - Sea areas -> traveling from
# We will need to employ some if/else logic to correctly combine
# these three data layers.

# Additionally, if the resolution is too high, we're going to
# make it a courser resolution by a factor of 10 - so each
# pixel will be e.g. 250m x 250m rather than 25m x 25m - otherwise,
# many of the subsequent raster calculations will fail or 
# take far far too long. This means that we will be calculating
# distance from the coast to a maximum of 1/4 km accuracy.

if (res >= 100) {
  canvasrow <- nrow(dem_3005)
  canvascol <- ncol(dem_3005)
} else {
  canvasrow <- nrow(dem_3005) / 10
  canvascol <- ncol(dem_3005) / 10
}

canvas <- terra::rast(extent = extent, # `extent` defined in 03-1 above 
                      nrow = canvasrow, 
                      ncol = canvascol, 
                      nlyr = 1)
terra::values(canvas) <- 0 # assume everything is the sea (0)
terra::crs(canvas) <- "epsg:3005"

canvas <- terra::mask(canvas, land, inverse = TRUE)
canvas <- terra::ifel(is.na(canvas), 1, canvas) # masked areas == land == 1
terra::plot(canvas)

land <- terra::ifel(dem_3005 > 0, 1, 0) # now extract land areas from DEM
land <- terra::resample(land, canvas) # resample `land` to match `canvas` extent/resolution
terra::plot(land)

land <- terra::merge(land, canvas, first = TRUE) # now merge the two
sea <- terra::ifel(land == 0, 0, NA) # now we can extract the sea

# 03-4 Calculate coast distance ----
# NOTE: IF WE INCLUDE ALASKA BORDER REGION LATER, we will need to 
# include the Alaska DEM in this to correctly measure distance from
# coast for the Alaska Border region. 
coast_distance <- terra::merge(terra::resample(elev, sea), sea)
coast_distance <- terra::gridDist(coast_distance, target = 0) # calculate distance from the sea!
terra::plot(coast_distance)
terra::plot(coast_distance <= 30000)

# Cut down to only include 30km distance and clip to land areas
c30 <- coast_distance <= 30000
c30 <- terra::merge(sea, c30)
c30 <- terra::ifel(c30 == 1, 1, NA) # set any non-valid nesting areas == NA

# Set it to match the resolution specified above
r <- c30
terra::res(r) <- res
c30 <- terra::resample(c30, r)

# Save it
if (length(list.files("GIS/Distance_cutoffs/")) == 0|overwrite == TRUE) terra::writeRaster(c30, "GIS/Distance_cutoffs/distance_cutoffs.tiff", overwrite = overwrite)

# Clean up
rm(canvas, coast_distance, extent, land, sea, r, usa, canvascol, canvasrow)
gc()

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
cons_reg$mnh_count <- exactextractr::exact_extract(mnh, cons_reg, 'count') # count the number of cells within each cons region
cons_reg$mnh_m <- cons_reg$mnh_count * 625 # each raster cell is 25m x 25m, aka 625 meters squared
cons_reg$mnh_km <- cons_reg$mnh_m / 1000000 # go from m2 to km2
cons_reg$mnh_ha <- cons_reg$mnh_km * 100 # multiply by 100 to go from km2 to hectares

# Checks out nicely!
cons_reg[["mnh_km"]][cons_reg$MMCR_NA == "Haida Gwaii"] # HG should be just under 10k
sum(cons_reg[["mnh_km"]][cons_reg$MMCR_NA %in% c("West and North Vancouver Island", "East Vancouver Island")]) # Vancouver Island (incl. Nootka + Quadra + Gulf islands) should be ~32k

# 04-3 Save MAMU nesting habitat ----

terra::writeRaster(mnh, "GIS/MAMU_nesting_habitat.tiff")
st_write(cons_reg, "GIS/cons_reg.shp", append = FALSE) # overwrite the conservation region shapefile with a shp containing the habitat area numbers


# 05 CLEAN UP -------------------------------------------------------------

# Clean up any temp files
unlink("temp", recursive = TRUE)
