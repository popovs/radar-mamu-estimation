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
# THIS WHOLE SCRIPT WILL TAKE ABOUT AN HOUR TO RUN.


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
rm(i, tmp, filename, dem_3005)

# To quickly visualize
#plot(stars::read_stars("GIS/DEM/Regional_DEM/NMC_3005.tiff")) # etc.

# 02-2 Apply elevation cutoffs ----

elevation_cutoffs <- setNames(data.frame(cons_reg$MMCR_NA,
                                         gsub('\\b(\\pL)\\pL{3,}|.','\\U\\1', cons_reg$MMCR_NA, perl = TRUE),
                                         c(800, 700, 800, 1500, 1200, 1100, 800)),
                              c("region", "abbreviation", "elevation_m"))


dir.create("GIS/Elevation_cutoffs")

# This will take approximately 5 minutes
for (i in 1:nrow(elevation_cutoffs)) {
  message("Reclassifying ", elevation_cutoffs$region[i])
  tmp <- terra::rast(paste0("GIS/DEM/Regional_DEM/", elevation_cutoffs$abbreviation[i], "_3005.tiff"))
  tmp <- terra::ifel(tmp < elevation_cutoffs$elevation_m[i], 1, NA)
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
m <- terra::mosaic(ec_collection, # probably need to sum and ignore NA
                   filename = "GIS/Elevation_cutoffs/elevation_cutoffs.tiff")
