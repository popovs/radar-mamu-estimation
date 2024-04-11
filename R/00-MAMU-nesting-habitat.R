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
res <- 250 # res MUST be >25, as the DEM goes does to 25m accuracy.
stopifnot("`res` must be >=25, as the DEM goes down to 25m accuracy." = res > 25)
overwrite <- FALSE

#library(devtools)
#devtools::install_github("popovs/MAMU")
library(MAMU)

# Read in the conservation regions shapefile
library(sf)
cons_reg <- st_read("GIS/cons_reg.shp")
cons_reg <- st_transform(cons_reg, 3005) # Set BC Albers projection

# 00-1 Prep nest data ----
# Read in nest data
nests <- st_read("GIS/MAMU_nests.gpkg")

# Add buffers of uncertainty around the nests
stopifnot("Nest locations must have an uncertainty associated with them (LOC_UNCE_1 column)." = !(any(is.na(nests$LOC_UNCE_1))))
stopifnot("Nest location uncertainly calculations currently only support uncertainty expressed in meters (LOC_UNCE_2 column)." = all(tolower(unique(nests$LOC_UNCE_2)) %in% c("m", "meters", "meter")))

nests <- st_buffer(nests, dist = nests$LOC_UNCE_1)

# TODO: replace all the "if length files in DEM folder == 0
# then run this code" w "if overwrite == TRUE" piece above -
# this will make it targets-friendly
# TODO: if files exist but overwrite == FALSE, just run the
# code in-memory for all chunks

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

# 02-2 Extract elevation cutoffs ----

# The exact elevational cutoffs vary by region and are primarily
# dictated by the growing conditions of the trees that MAMU
# prefer to nest in. Large coniferous trees that can support
# MAMU nests can grow in higher elevations at lower latitudes.
# Conversely, the further north you go, the lower the maximum
# MAMU nesting elevation. The nest data is the main source of 
# cutoffs. Adding more nest data will trigger a re-run of this
# analysis.

# Extract nest elevations from DEM
# Using exactextract to extract the *mean* elevation within the
# nest uncertainty radius
nests$elev_m <- exactextractr::exact_extract(dem_3005, nests, 'mean')

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

# We're now going to set elevation cutoffs that contain 95%
# of the nests, from both directions (min and max elevations).
# So, the cutoffs will be @ the 2.5%ile on each end.
elevation_cutoffs <- setNames(data.frame(cons_reg$MMCR_NA,
                                         gsub('\\b(\\pL)\\pL{3,}|.','\\U\\1', cons_reg$MMCR_NA, perl = TRUE),
                                         unlist(lapply(cons_reg$MMCR_NA, function(x){quantile(nests[["elev_m"]][nests$cons_reg == x], 0.025, na.rm = TRUE)[[1]]})),
                                         unlist(lapply(cons_reg$MMCR_NA, function(x){quantile(nests[["elev_m"]][nests$cons_reg == x], 0.975, na.rm = TRUE)[[1]]}))
                                         ),
                              c("region", "abbreviation", "elev_min_m", "elev_max_m"))

# If northern mainland coast and alaska border are NA, use CMC elevation cutoff
# Minimum
if (is.na(elevation_cutoffs[["elev_min_m"]][elevation_cutoffs$abbreviation == "NMC"])) elevation_cutoffs[["elev_min_m"]][elevation_cutoffs$abbreviation == "NMC"] <- elevation_cutoffs[["elev_min_m"]][elevation_cutoffs$abbreviation == "CMC"]
if (is.na(elevation_cutoffs[["elev_min_m"]][elevation_cutoffs$abbreviation == "AB"])) elevation_cutoffs[["elev_min_m"]][elevation_cutoffs$abbreviation == "AB"] <- elevation_cutoffs[["elev_min_m"]][elevation_cutoffs$abbreviation == "CMC"]
# Maximum
if (is.na(elevation_cutoffs[["elev_max_m"]][elevation_cutoffs$abbreviation == "NMC"])) elevation_cutoffs[["elev_max_m"]][elevation_cutoffs$abbreviation == "NMC"] <- elevation_cutoffs[["elev_max_m"]][elevation_cutoffs$abbreviation == "CMC"]
if (is.na(elevation_cutoffs[["elev_max_m"]][elevation_cutoffs$abbreviation == "AB"])) elevation_cutoffs[["elev_max_m"]][elevation_cutoffs$abbreviation == "AB"] <- elevation_cutoffs[["elev_max_m"]][elevation_cutoffs$abbreviation == "CMC"]

# Round to nearest 10 meter
elevation_cutoffs$elev_min_m <- round(elevation_cutoffs$elev_min_m, -1)
elevation_cutoffs$elev_max_m <- round(elevation_cutoffs$elev_max_m, -1)

# 02-3 Apply elevation cutoffs ----
if (overwrite == TRUE) unlink("GIS/Elevation_cutoffs", recursive = TRUE)
dir.create("GIS/Elevation_cutoffs", showWarnings = F)

# This will take approximately 5 minutes
# RUN IF FILES DON'T EXIST OR OVERWRITE == TRUE AT TOP OF SCRIPT
if (length(list.files("GIS/Elevation_cutoffs/")) == 0|overwrite == TRUE){
  for (i in 1:nrow(elevation_cutoffs)) {
    message("Reclassifying ", elevation_cutoffs$region[i])
    tmp <- terra::rast(paste0("GIS/DEM/Regional_DEM/", elevation_cutoffs$abbreviation[i], "_3005.tiff"))
    if (all(terra::res(tmp) != res)) {
      r <- tmp 
      terra::res(r) <- res
      tmp <- terra::resample(tmp, r)
      rm(r)
    }
    tmp <- terra::ifel(tmp > 0, tmp, NA)
    tmp <- terra::ifel(tmp < elevation_cutoffs$elev_max_m[i], tmp, NA)
    tmp <- terra::ifel(tmp < elevation_cutoffs$elev_min_m[i], 1, 2)
    #tmp <- terra::ifel(tmp < elevation_cutoffs$elev_max_m[i] & tmp > elevation_cutoffs$elev_min_m[i], 1, NA)
    filename <- file.path("GIS/Elevation_cutoffs", paste0(elevation_cutoffs$abbreviation[i], "_", elevation_cutoffs$elev_max_m[i], "m.tiff"))
    terra::writeRaster(tmp, filename, overwrite = overwrite)
  }
  beepr::beep()
  rm(i, tmp, filename)
} else {
  # TODO: make it run vs read in stuff by choice
  # Else just run this in memory
  ec <- list()
  for (i in 1:nrow(elevation_cutoffs)) {
    message("Reclassifying ", elevation_cutoffs$region[i])
    tmp <- terra::rast(paste0("GIS/DEM/Regional_DEM/", elevation_cutoffs$abbreviation[i], "_3005.tiff"))
    if (all(terra::res(tmp) != res)) {
      r <- tmp
      terra::res(r) <- res
      tmp <- terra::resample(tmp, r)
      rm(r)
    }
    tmp <- terra::ifel(tmp > 0, tmp, NA)
    tmp <- terra::ifel(tmp < elevation_cutoffs$elev_max_m[i], tmp, NA)
    tmp <- terra::ifel(tmp < elevation_cutoffs$elev_min_m[i], 1, 2)
    #tmp <- terra::ifel(tmp < elevation_cutoffs$elev_max_m[i] & tmp > elevation_cutoffs$elev_min_m[i], 1, NA)
    ec[[i]] <- tmp
  }
  beepr::beep()
  rm(i, tmp)
}


# 02-4 Merge elevation cutoff rasters ----

# Now we can merge all these rasters into one single raster
# for ease of use.

if (!("elevation_cutoffs.vrt" %in% list.files("GIS/Elevation_cutoffs"))|overwrite == TRUE) {
  # Read them into environment
  f <- list.files("GIS/Elevation_cutoffs", full.names = T)
  f <- f[!grepl("elevations_cutoffs", f)]
  ec <- lapply(f, terra::rast) # 'ec' for elevational cutoffs
  names(ec) <- basename(f)
  
  # Mosaic all the ec rasters together
  # This will take approximately 1-2 minutes
  ec <- make_vrt("GIS/Elevation_cutoffs/",
                 filename = "GIS/Elevation_cutoffs/elevation_cutoffs.vrt",
                 overwrite = overwrite)
  ec_collection <- terra::sprc(ec)
  elev <- terra::mosaic(ec_collection,
                        filename = "GIS/Elevation_cutoffs/elevation_cutoffs.tiff",
                        overwrite = overwrite)
  #list.files("GIS/Elevation_cutoffs/")
  #terra::plot(elev)
  rm(ec, ec_collection, f)
} else {
  message("Reading in existing elevation cutoffs file...")
  elev <- terra::rast("GIS/Elevation_cutoffs/elevation_cutoffs.tiff")
  # Resample to match the resolution specified above
  if (all(terra::res(elev) != res)) {
    message("Resampling elevation cutoffs to target resolution (", res, "m)...")
    r <- elev
    terra::res(r) <- res
    elev <- terra::resample(elev, r)
    rm(r)
  }
}
# TODO: add another 'else if' if you want to re-run w new res and not read in old file

gc()


# 03 DISTANCE FROM COAST CUTOFF -------------------------------------------

# Next we will create a separate raster of all points wit
# 30 km distance from the coast, as BC nest survey data
# indicates that 99% of MAMU nests are within 30 km of the
# coastline. We are going to assume a 30 km distance from 
# shore flying around mountain barriers  but allowing for 
#flight over water.

# Only run if files are empty or overwrite == TRUE
if (length(list.files("GIS/Distance_cutoffs/")) == 0|overwrite == TRUE){

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
  coast_distance <- terra::ifel(coast_distance > 0 , 1, coast_distance) # we don't care about the elevation cutoff minimums here - so reclassify them to only keep the max cutoffs
  coast_distance <- terra::gridDist(coast_distance, target = 0) # calculate distance from the sea!
  terra::plot(coast_distance)
  terra::plot(coast_distance <= 30000)
  
  # Cut down to only include 30km distance and clip to land areas
  c30 <- coast_distance <= 30000
  c30 <- terra::merge(sea, c30)
  c30 <- terra::ifel(c30 == 1, 1, NA) # set any non-valid nesting areas == NA
  
  # Set it to match the resolution specified above
  if (any(terra::res(c30) != res)) {
    r <- c30
    terra::res(r) <- res
    c30 <- terra::resample(c30, r)
    rm(r)
  }
  
  # Save it
  if (length(list.files("GIS/Distance_cutoffs/")) == 0|overwrite == TRUE) terra::writeRaster(c30, "GIS/Distance_cutoffs/distance_cutoffs.tiff", overwrite = overwrite)
  
  # Clean up
  rm(canvas, coast_distance, extent, land, usa, canvascol, canvasrow)
  gc()
  
} else {
  message("Reading in existing coast distance file...")
  c30 <- terra::rast("GIS/Distance_cutoffs/distance_cutoffs.tiff")
  # Resample to match the resolution specified above
  if (all(terra::res(c30) != res)) {
    message("Resampling coast distance to target resolution (", res, "m)...")
    r <- c30
    terra::res(r) <- res
    c30 <- terra::resample(c30, r)
    rm(r)
  }
}

gc()


# 04 NEST COST DISTANCE ---------------------------------------------------

# Next, we calculate the flight cost values of each nest. 
# Evidence shows that MAMU take the least-cost flightpaths to
# their nests; that is, they tend to fly along valley contours
# rather than straight across ridges (even if the ridges are
# below their nest cutoff elevations). Here we will calculate a
# raster of the flight cost (elevation * distance from coast)
# to generate a landscape where birds are more or less likely to
# nest. This will be used to: 1) delineate the nesting catchment 
# boundaries later and 2) eliminate "high cost" nesting areas 
# that may still be included within the elevation cutoff and 
# coast distance rasters.

if (overwrite == TRUE) unlink("GIS/Cost_cutoffs", recursive = TRUE)
dir.create("GIS/Cost_cutoffs", showWarnings = F)

# Only run this if file does not exist or overwrite == TRUE

if (!("cost_layer.tiff" %in% list.files("GIS/Cost_cutoffs"))|overwrite == TRUE) {
  # 04-1 Prepare cost layer ----
  # First, similar to above, lower the resolution of high res
  # DEM by a factor of 10, or else the calculations will fail.
  if (res < 100) {
    c <- dem_3005
    terra::res(c) <- res * 10
    c <- terra::resample(dem_3005, c)
  } else {
    c <- dem_3005
  }
  
  # Ensure no negative values in the raster, or the costDist
  # function will fail.
  c <- terra::ifel(c < 0, 0, c)
  
  # 04-2 Calculate cost ----
  # The cost distance function calculates distance from shore (i.e.,
  # the `gridDist` function we just used above) and multiplies it
  # by the 'cost' layer (elevation). Higher elevations are more
  # costly to fly over. 
  cost <- terra::costDist(c, target = 0, scale = 1000) # divide values by 1000 so output numbers are smaller
  rm(c)
  
  # Save it - we will still use this layer later to extract the 
  # mean nesting cost per catchment
  if (!("cost_layer.tiff" %in% list.files("GIS/Cost_cutoffs"))|overwrite == TRUE) terra::writeRaster(cost, "GIS/Cost_cutoffs/cost_layer.tiff")

} else {
  cost <- terra::rast("GIS/Cost_cutoffs/cost_layer.tiff")
  # Resample to match the resolution specified above
  if (all(terra::res(cost) != res)) {
    r <- cost
    terra::res(r) <- res
    cost <- terra::resample(cost, r)
    rm(r)
  }
}

terra::plot(cost)

# 04-3 Extract nest costs ----
# Using exactextract to extract the *mean* cost within the
# nest uncertainty radius
nests$cost <- exactextractr::exact_extract(cost, nests, 'mean')

ggplot2::ggplot(data = nests, 
                ggplot2::aes(x = cons_reg, 
                             y = cost)) + 
  ggplot2::geom_boxplot() + 
  ggplot2::geom_jitter() +
  ggplot2::theme_minimal()


cost_cutoffs <- setNames(data.frame(cons_reg$MMCR_NA,
                                    gsub('\\b(\\pL)\\pL{3,}|.','\\U\\1', cons_reg$MMCR_NA, perl = TRUE),
                                    unlist(lapply(cons_reg$MMCR_NA, function(x){quantile(nests[["cost"]][nests$cons_reg == x], 0.025, na.rm = TRUE)[[1]]})),
                                    unlist(lapply(cons_reg$MMCR_NA, function(x){quantile(nests[["cost"]][nests$cons_reg == x], 0.975, na.rm = TRUE)[[1]]}))
                                    ),
                         c("region", "abbreviation", "cost_min", "cost_max"))
# If northern mainland coast and alaska border are NA, use CMC elevation cutoff
# Minimum
if (is.na(cost_cutoffs[["cost_min"]][cost_cutoffs$abbreviation == "NMC"])) cost_cutoffs[["cost_min"]][cost_cutoffs$abbreviation == "NMC"] <- cost_cutoffs[["cost_min"]][cost_cutoffs$abbreviation == "CMC"]
if (is.na(cost_cutoffs[["cost_min"]][cost_cutoffs$abbreviation == "AB"])) cost_cutoffs[["cost_min"]][cost_cutoffs$abbreviation == "AB"] <- cost_cutoffs[["cost_min"]][cost_cutoffs$abbreviation == "CMC"]
# Maximum
if (is.na(cost_cutoffs[["cost_max"]][cost_cutoffs$abbreviation == "NMC"])) cost_cutoffs[["cost_max"]][cost_cutoffs$abbreviation == "NMC"] <- cost_cutoffs[["cost_max"]][cost_cutoffs$abbreviation == "CMC"]
if (is.na(cost_cutoffs[["cost_max"]][cost_cutoffs$abbreviation == "AB"])) cost_cutoffs[["cost_max"]][cost_cutoffs$abbreviation == "AB"] <- cost_cutoffs[["cost_max"]][cost_cutoffs$abbreviation == "CMC"]

# Round to nearest 10 meter
cost_cutoffs$cost_min <- round(cost_cutoffs$cost_min, -1)
cost_cutoffs$cost_max <- round(cost_cutoffs$cost_max, -1)

# 04-4 Apply cost cutoffs ----

# This will take approximately 5 minutes
# RUN IF FILES DON'T EXIST OR OVERWRITE == TRUE AT TOP OF SCRIPT
if (length(list.files("GIS/Cost_cutoffs/")) == 0|overwrite == TRUE){
  for (i in 1:nrow(cost_cutoffs)) {
    message("Reclassifying ", cost_cutoffs$region[i])
    # First split apart/mask the main `cost` raster into regions
    tmp <- terra::mask(cost, cons_reg[cons_reg$MMCR_NA == cost_cutoffs$region[i],], touches = TRUE)
    tmp <- terra::ifel(tmp > 0, tmp, NA)
    tmp <- terra::ifel(tmp < cost_cutoffs$cost_max[i], tmp, NA)
    tmp <- terra::ifel(tmp < cost_cutoffs$cost_min[i], 1, 2)
    #tmp <- terra::ifel(tmp < cost_cutoffs$cost[i] & tmp > 0, 1, NA)
    filename <- file.path("GIS/Cost_cutoffs", paste0(cost_cutoffs$abbreviation[i], "_", cost_cutoffs$cost_max[i], ".tiff"))
    terra::writeRaster(tmp, filename, overwrite = overwrite)
  }
  beepr::beep()
  rm(i, tmp, filename)
} #else {
#   # TODO: make this run if you don't want to overwrite but want to produce new vals 
#   # jsut like in 'apply elevation cutoffs'
#   cc <- list()
#   for (i in 1:nrow(cost_cutoffs)) {
#     message("Reclassifying ", cost_cutoffs$region[i])
#     # First split apart/mask the main `cost` raster into regions
#     tmp <- terra::mask(cost, cons_reg[cons_reg$MMCR_NA == cost_cutoffs$region[i],], touches = TRUE)
#     tmp <- terra::ifel(tmp > 0, tmp, NA)
#     tmp <- terra::ifel(tmp < cost_cutoffs$cost_max[i], tmp, NA)
#     tmp <- terra::ifel(tmp < cost_cutoffs$cost_min[i], 1, 2)
#     #tmp <- terra::ifel(tmp < cost_cutoffs$cost[i] & tmp > 0, 1, NA)
#     cc[[i]] <- tmp
#   }
#   beepr::beep()
#   rm(i, tmp)
# }

# To quickly visualize
list.files("GIS/Cost_cutoffs")
#terra::plot(terra::rast("GIS/Cost_cutoffs/CMC_2670.tiff")) # etc.
gc()

# 04-5 Merge cost cutoff rasters ----

# Now we can merge all these rasters into one single raster
# for ease of use.

if (!("cost_cutoffs.tiff" %in% list.files("GIS/Cost_cutoffs"))|overwrite == TRUE) {
  # Only do this step if its the first time running or overwrite == TRUE
  # Read them into environment
  f <- list.files("GIS/Cost_cutoffs", full.names = T)
  f <- f[!grepl("cost_cutoffs.vrt|cost_cutoffs.tiff|cost_layer.tiff", f)]
  cc <- lapply(f, terra::rast) # 'cc' for cost cutoffs
  names(cc) <- basename(f)
  
  # Mosaic all the ec rasters together
  # This will take approximately 1-2 minutes
  cc_collection <- terra::sprc(cc)
  cc <- terra::mosaic(cc_collection,
                      filename = "GIS/Cost_cutoffs/cost_cutoffs.tiff",
                      overwrite = overwrite)
  rm(cc_collection, f)
} else {
  message("Reading in existing cost cutoffs files...")
  cc <- terra::rast("GIS/Cost_cutoffs/cost_cutoffs.tiff")
  # Resample to match the resolution specified above
  if (all(terra::res(cc) != res)) {
    message("Resampling cost cutoffs to target resolution (", res, ")...")
    r <- cc
    terra::res(r) <- res
    cc <- terra::resample(cc, r)
    rm(r)
  }
}
# TODO: add another 'else if' if you want to re-run w new res and not read in old file

gc()

# 05 FOREST COVER ---------------------------------------------------------


# MAMU are almost certainly not nesting in urban or other completely
# treeless areas. Cut those out. 

# Download the University of Maryland global forest cover dataset
# https://glad.earthengine.app/view/global-forest-change#bl=off;old=0;dl=off;lon=-486.9726343690056;lat=51.450799512228464;zoom=6;

dir.create("GIS/Tree_cover/", showWarnings = FALSE)

# 05-1 Download forest cover data ---- 
# This will take ~1 hour to run
# This can also be done in QGIS/ArcGIS - mosaic the 3 map tiles
# into one raster -> project/warp to EPSG3005 -> crop/mask to 
# study area
if (length(list.files("GIS/Forest_cover/")) == 0) {
  urls <- c("https://storage.googleapis.com/earthenginepartners-hansen/GFC-2023-v1.11/Hansen_GFC-2023-v1.11_treecover2000_60N_140W.tif",
            "https://storage.googleapis.com/earthenginepartners-hansen/GFC-2023-v1.11/Hansen_GFC-2023-v1.11_treecover2000_60N_130W.tif",
            "https://storage.googleapis.com/earthenginepartners-hansen/GFC-2023-v1.11/Hansen_GFC-2023-v1.11_treecover2000_50N_130W.tif"
  )
  
  # These files will take longer than one minute to download, so
  # increase the download timeout to 5 mins
  options(timeout = 300)
  
  # Download % forest cover data (year 2000) for BC
  lapply(urls, function(x) {
    download.file(url = x, 
                  destfile = file.path("GIS/Forest_cover", basename(x))
                  )
    })
  
  options(timeout = 60) # reset default
  rm(urls)

  # Read files into memory
  forest <- lapply(list.files("GIS/Forest_cover", full.names = TRUE), terra::rast)
  
  # Mosaic
  forest <- terra::sprc(forest)
  forest <- terra::mosaic(forest)

  # Crop to conservation regions
  forest <- terra::project(forest, "EPSG:3005") # This will take the longest
  forest <- terra::mask(forest, cons_reg)
  forest <- terra::crop(forest, cons_reg)
  gc()
  
  # Save
  terra::writeRaster(forest, "GIS/Forest_cover/forest_cover.tiff")
  beepr::beep()
  
} else {
  message("Reading in existing forest cover file...")
  forest <- terra::rast("GIS/Forest_cover/forest_cover.tiff")
  if (all(terra::res(forest) != res)) {
    message("Resampling forest cover to target resolution (", res, "m)...")
    r <- forest
    terra::res(r) <- res
    forest <- terra::resample(forest, r)
    rm(r)
  }
}

# 05-2 Extract nest forest cover ----
nests$forest <- exactextractr::exact_extract(forest, nests, 'mean')

ggplot2::ggplot(data = nests, 
                ggplot2::aes(x = cons_reg, 
                             y = forest)) + 
  ggplot2::geom_boxplot() + 
  ggplot2::geom_jitter() +
  ggplot2::theme_minimal()

hist(nests$forest, breaks = 100)

# 05-3 Apply forest cutoff ----
# In the case of forests, we know with 100% confidence that
# MAMU will nest in areas with 100% tree cover. As such, we
# don't need to define a 'maximum tree cover' cutoff for 
# the forested areas layer. Instead, we need to choose a 
# minimum cutoff. In the same vein, setting a variable minimum
# forest cover cutoff by region might be cutting out too much
# habitat area. This forest data doesn't contain information
# on species composition, whereas elevation cutoffs are done
# with the max elevation of prefered nest tree species by 
# region in mind.
# Instead, we'll just cut out the bottom 2.5% tree cover (6
# nests out of 242 in the dataset).
min_forest <- quantile(nests$forest, 0.025)

# `fc` for forest cutoff
fc <- terra::ifel(forest < min_forest, NA, 1)

if (!any(grepl("forest_cutoff.tiff", list.files("GIS/Forest_cover/")))|overwrite == TRUE){
  terra::writeRaster(fc, "GIS/Forest_cover/forest_cutoff.tiff", overwrite = overwrite)
}

gc()


# 06 MERGE CUTOFFS --------------------------------------------------------

# Finally, merge the cutoff rasters together to create a 
# "MAMU containment zone" area that has a high probability of
# containing high quality MAMU nesting habitat. This raster 
# will be used as our maximum MAMU area within BC that we will 
# extrapolate our population estimates to.

# 06-1 Resample raster resolutions ----

# The raster datasets are all various resolutions as they come
# from different base datasets (e.g. the DEM versus the forest
# cover data) or the resolution was increased for calculation
# to even be feasible (e.g. cost distance). 

# Per the BC DEM website, this raster product is at a 25m 
# resolution, but gridded to a 0.75 arc-second scale. 
# https://www2.gov.bc.ca/gov/content/data/geographic-data-services/topographic-data/elevation/digital-elevation-model
# This means all our raster data is at a somewhat odd 
# 17.37227 x 17.37227 resolution:
#res(dem_3005)
#res(elev)
# As such, the maximum resolution we can safely go for is 25m.

# Resample raster resolutions to the `res` value specified above,
# to a maximum of 25m. In theory the minimum should be set in the
# script above but just in case it's not here's another failsafe.
res <- ifelse(res < 25, 25, res)

# Resample `elev`
if (all(terra::res(elev) != res)) {
  r <- elev
  terra::res(r) <- res
  elev <- terra::resample(elev, r)
  beepr::beep()
  rm(r)
}

# Resample `c30`
if (all(terra::res(c30) != res)) {
  r <- c30
  terra::res(r) <- res
  c30 <- terra::resample(c30, r)
  beepr::beep()
  rm(r)
}

# Resample `cc`
if (all(terra::res(cc) != res)) {
  r <- cc
  terra::res(r) <- res
  cc <- terra::resample(cc, r)
  beepr::beep()
  rm(r)
}

# Resample `fc`
if (all(terra::res(fc) != res)) {
  r <- fc
  terra::res(r) <- res
  cc <- terra::resample(fc, r)
  beepr::beep()
  rm(r)
}

# Confirm that the raster extents align
terra::ext(c30) == terra::ext(cc)
terra::ext(c30) == terra::ext(elev)
terra::ext(cc) == terra::ext(elev)
terra::ext(fc) == terra::ext(elev)

# Resample everything to match c30 and cc
elev <- terra::resample(elev, cc)
fc <- terra::resample(fc, cc)

terra::ext(c30) == terra::ext(cc)
terra::ext(c30) == terra::ext(elev)
terra::ext(cc) == terra::ext(elev)
terra::ext(fc) == terra::ext(elev)


# 06-2 Merge nest layers ----

# Now we can merge our rasters into a 'suitable nesting
# habitat' final layer.

# `mnh` for "MAMU nesting habitat"
# This will take approximately 3-4 minutes
mnh <- sum(cc, c30, elev, na.rm = TRUE)
terra::plot(mnh)

# We're going to save two layers - first, 'mamu accessible',
# i.e. they can reach it at all to nest (it's below elevation,
# within 30km, and within cost). Second, we're going to save
# 'mamu likely', which will include accessible regions but
# then also cut out places they are unlikely to nest (it's
# below minimum required elevation, cost, or forest cover).

# First layer to save - MAMU accessible areas
# Anything accesible was stored as >= 1 in each raster layer.
# Because 3 raster layers went into mnh, we can choose any
# cells >= 3 for the "mamu accessible zone".
maz <- terra::ifel(mnh >= 3, 1, NA)
terra::plot(maz)
if (!any(grepl("MAMU_accessible_zone.tiff", list.files("GIS")))|overwrite == TRUE) terra::writeRaster(maz, "GIS/MAMU_accessible_zone.tiff", overwrite = overwrite)

# Second layer to save - MAMU nesting areas
# Here we add in the forest cutoff layer, as we want to 
# exclude any non-forested areas from the nesting habitat.
mnh <- sum(mnh, fc, na.rm = TRUE)
terra::plot(mnh)
max_mnh <- terra::minmax(mnh, compute = FALSE)[2] # extract max value
mnh <- terra::ifel(mnh == max_mnh, 1, NA)
terra::plot(mnh)
if (!(any(grepl("MAMU_nesting_habitat.tiff", list.files("GIS"))))|overwrite == TRUE) terra::writeRaster(mnh, "GIS/MAMU_nesting_habitat.tiff", overwrite = overwrite)


# 06-3 Calculate habitat area within each conservation region ----

# This serves to both update our conservation region shapefile
# with the correct habitat areas, but also to double check that 
# our numbers roughly line up with what we expect - the Haida 
# Gwaii habitat area should be just under 10k km squared.
cons_reg$mnh_count <- exactextractr::exact_extract(mnh, cons_reg, 'count') # count the number of cells within each cons region
cons_reg$mnh_m <- cons_reg$mnh_count * res^2 # each raster cell is res m x res m (e.g., 25m x 25m), aka 625 meters squared
cons_reg$mnh_km <- cons_reg$mnh_m / 1000000 # go from m2 to km2
cons_reg$mnh_ha <- cons_reg$mnh_km * 100 # multiply by 100 to go from km2 to hectares

# Checks out nicely!
cons_reg[["mnh_km"]][cons_reg$MMCR_NA == "Haida Gwaii"] # HG should be just under 10k
sum(cons_reg[["mnh_km"]][cons_reg$MMCR_NA %in% c("West and North Vancouver Island", "East Vancouver Island")]) # Vancouver Island (incl. Nootka + Quadra + Gulf islands) should be under ~32k

# 06-4 Save MAMU nesting habitat ----

if (!(any(grepl("MAMU_nesting_habitat.tiff", list.files("GIS"))))|overwrite == TRUE) {
  st_write(cons_reg, "GIS/cons_reg.shp", append = FALSE) # overwrite the conservation region shapefile with a shp containing the habitat area numbers
}


# 07 CLEAN UP -------------------------------------------------------------

rm(list = ls())

# Clean up any temp files
unlink("temp", recursive = TRUE)
