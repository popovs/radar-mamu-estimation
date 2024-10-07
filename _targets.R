# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Needed for `tar_map()`

# DISCLAIMER 1: THIS WHOLE SCRIPT WILL TAKE SEVERAL HOURS TO RUN.
# There are likely faster or more efficient ways of doing this,
# but the primary goal of this script is reproducibility and
# ease of understanding. If you wish to simply test the script, 
# modify the raster resolutions to a courser scale (e.g., 500m) 
# to run through the script quickly.

# DISCLAIMER 2: this was run with 16 GB ram and will likely 
# fail with less - that said, many of these steps can be
# easily replicated in QGIS to the same effect. The benefit of
# this code is the exact reproducibility.

# Set resolution (in meterse) for all raster analysis
# smaller number = MUCH SLOWER
res <- 250 # res MUST be >25, as the DEM goes does to 25m accuracy.

# Set target options:
tar_option_set(
  packages = c("MAMU", # remotes::install_github("popovs/MAMU")
               "sf",
               "terra"),
  format = tar_format( # Default qs is superceded by qs2. Install qs2 and specify read & write fxns
    read = function(path) { qs2::qs_read(path) }, 
    write = function(object, path) { qs2::qs_save(object = object, file = path) }
    ), 
  memory = "transient", # unload memory for each target line after it's completed
  garbage_collection = TRUE, # run gc() prior to each target
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Objects to reference in pipeline
tiles_to_download <- data.frame(tiles = c("103k", "103j", "103f", "103g", "103c", "103b", "102o", # haida gwaii
                                          "103o", "103p", "103j", "103i", "103g", "103h", "103a", "93e", "93d", "93c", "93l", "93m", "102p", "92m", "92n", "104a", "104b", # north/central coast
                                          "102i", "92l", "92k", "92e", "92f", "92c", "92b", # vancouver island
                                          "92j", "92g", "92h" # lower mainland
                                          ))
DEM_dir <- "GIS/DEM"

# Run tar_make() to execute the pipeline
list(
  tar_target(regions, prepare_regions()),
  tar_target(nests, prepare_nests(regions = regions)),
  # PREPARE DEM 
  # All the DEM data products will be saved as raster files on the
  # disk rather than simply as _targets objects, so that QGIS can
  # access and plot them. The targets pipeline will track changes
  # to the files. 
  # Download DEM tiles
  tar_map(
    values = tiles_to_download, # params need to be passed as a df/tibble, defined OUTSIDE the pipeline
    # Download DEM tiles
    tar_target(DEM_tiles,
               download_dem_tile(tile = tiles, # `tiles` in this case refers to the `tiles` column in `tiles_to_download` df
                                 save_output = TRUE,
                                 overwrite = TRUE,
                                 output_dir = DEM_dir),
               format = "file")
    ),
  # Track files within the DEM_dir
  tar_target(DEM_tile_files, list.files(file.path(DEM_dir, "DEM_tiles"), full.names = TRUE)),
  # Make VRT
  tar_target(DEM_VRT, 
             terra::vrt(x = DEM_tile_files,
                        filename = file.path(DEM_dir, "BC_DEM_VRT.vrt"),
                        overwrite = TRUE,
                        return_filename = TRUE), # for format "file" targets, the output MUST be a filepath
             format = "file"),
  # Merge DEM tiles and reproject to 3005 (~30 mins)
  tar_target(BC_DEM_3005, 
             merge_dem(vrt_path = DEM_VRT,
                       output_file = file.path(DEM_dir, "BC_DEM_EPSG3005.tiff"),
                       overwrite = TRUE),
             format = "file"),
  # Resample to pipeline resolution (but no need to save it as its own tiff file)
  tar_target(DEM, resample_dem(dem_path = BC_DEM_3005, res = res))
)
