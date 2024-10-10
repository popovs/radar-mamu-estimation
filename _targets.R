# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Needed for `tar_map()`
library(geotargets) # Needed to save `terra` SpatRaster targets

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
               "terra",
               "rnaturalearth"),
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
# BC DEM Ortho tiles
tiles_to_download <- data.frame(tiles = c("103k", "103j", "103f", "103g", "103c", "103b", "102o", # haida gwaii
                                          "103o", "103p", "103j", "103i", "103g", "103h", "103a", "93e", "93d", "93c", "93l", "93m", "102p", "92m", "92n", "104a", "104b", # north/central coast
                                          "102i", "92l", "92k", "92e", "92f", "92c", "92b", # vancouver island
                                          "92j", "92g", "92h" # lower mainland
                                          ))
# File directory to store DEM data
DEM_dir <- "GIS/DEM"
# Regions to iterate GIS operations over
regions_map <- sf::st_drop_geometry(sf::st_read("GIS/regions.gpkg"))

# Run tar_make() to execute the pipeline
list(
  # TODO: track regions and nests files themselves to track changes & execute pipeline if necessary
  tar_target(regions, prepare_regions()),
  tar_target(nests, prepare_nests(regions = regions)),
  #### PREPARE DEM ####
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
  tar_terra_rast(DEM, resample_dem(dem_path = BC_DEM_3005, res = res)),
  #### REGIONAL ELEVATION CUTOFFS ####
  # Extract nest elevations
  # The exact elevational cutoffs vary by region and are primarily
  # dictated by the growing conditions of the trees that MAMU
  # prefer to nest in. Large coniferous trees that can support
  # MAMU nests can grow in higher elevations at lower latitudes.
  # Conversely, the further north you go, the lower the maximum
  # MAMU nesting elevation. The nest data is the main source of 
  # cutoffs. Adding more nest data will trigger a re-run of this
  # analysis.
  # Using exactextract to extract the mean elevation within the 
  # nest uncertainty radius (nests were buffered in `prepare_nests()`)
  tar_target(nest_elev_m, exactextractr::exact_extract(DEM, nests, 'mean')),
  # Calculate and store regional elevation cutoffs in a table
  tar_group_by(elevation_cutoffs,
               nest_quantiles(nests, 
                              quant_data = nest_elev_m,
                              prefix = "elev_m"),
               region), # group by region col so targets later knows to run `reclass_elevation()` by row-level grouping
  #### NEST COST DISTANCE ####
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
  tar_terra_rast(cost, cost_distance(DEM)),
  # Extract nest cost
  # Same process as extracting elevation values + applying cutoffs
  tar_target(nest_cost, exactextractr::exact_extract(cost, nests, 'mean')),
  # Calculate and store regional cost cutoffs in a table
  tar_group_by(cost_cutoffs,
               nest_quantiles(nests, 
                              quant_data = nest_cost,
                              prefix = "cost"),
               region),
  #### ITERATE OVER REGIONAL CUTOFFS ####
  # For each region, iterate the following:
  # 1. Chop up the DEM into regions
  # 2. Apply elevation cutoffs to each regional DEM
  # 3. Apply cost distance cutoffs to each regional DEM
  mapped <- tar_map(
    values = regions_map, # params need to be passed as a df/tibble, defined OUTSIDE the pipeline
    # Intersect DEM with each conservation region
    tar_terra_rast(regional_DEM, 
                   crop_dem(
                     dem = DEM,
                     region_name = region, # `region` in this case refers to the `region` column in `regions_map` df
                     regions = regions # `regions` is the regions sf object (very first target)
                   )),
    # Apply elevation cutoffs for each region
    tar_terra_rast(regional_elev_cutoffs,
                   reclass_raster(
                     dem = regional_DEM,
                     region = region,
                     cutoffs = elevation_cutoffs,
                     min_col = "elev_m_min",
                     max_col = "elev_m_max"
                     )),
    # Apply cost cutoffs for each region 
    tar_terra_rast(regional_cost_cutoffs,
                   reclass_raster(
                     dem = regional_DEM,
                     region = region,
                     cutoffs = cost_cutoffs,
                     min_col = "cost_min",
                     max_col = "cost_max"
                   ))
    ),
  # Combine into single raster
  # `elev` for 'elevation cutoffs'
  tar_combine(elev_sprc, mapped[[2]], command = terra::sprc(list(!!!.x)) |> terra::wrap()), # need to wrap it to prevent "invalid pointer" error: https://stackoverflow.com/questions/74855695/load-raster-data-with-terra-in-targets-pipeline
  tar_terra_rast(elev, terra::unwrap(elev_sprc) |> 
                   terra::mosaic(fun = "max") |> # choose 'max' - by default assume value is 2, or accessible, in cases where mosaicing rasters results in some 1 or 2 overlap
                   terra::resample(DEM)), # resample to match DEM resolution
  # `cc` for 'cost cutoffs'
  tar_combine(cc_sprc, mapped[[1]], command = terra::sprc(list(!!!.x)) |> terra::wrap()), # need to wrap it to prevent "invalid pointer" error: https://stackoverflow.com/questions/74855695/load-raster-data-with-terra-in-targets-pipeline
  tar_terra_rast(cc, terra::unwrap(cc_sprc) |>
                   terra::mosaic(fun = "max") |> # choose 'max' - by default assume value is 2, or accessible, in cases where mosaicing rasters results in some 1 or 2 overlap
                   terra::resample(DEM)),
  #### DISTANCE FROM COAST CUTOFF ####
  # Next we will create a separate raster of all points with
  # 30 km distance from the coast, as BC nest survey data
  # indicates that 99% of MAMU nests are within 30 km of the
  # coastline. We are going to assume a 30 km distance from 
  # shore flying around mountain barriers  but allowing for 
  # flight over water.
  # The DEM data doesn't include the USA or inland BC.
  # These land areas need to be blocked off so the distance
  # algorithm can differentiate between land areas and
  # water areas.
  # Get USA land areas
  tar_terra_vect(usa_land, get_usa_land(extent = terra::ext(elev))),
  # Block off inland BC areas
  tar_terra_vect(bc_land, get_bc_land(extent = terra::ext(elev))),
  tar_terra_vect(land_mask, merge_land(usa_land, bc_land)),
  # Create canvas of target area
  tar_terra_rast(canvas, create_canvas(target_rast = elev)),
  # Block off land and extract `land` and `sea` areas
  tar_terra_rast(land, block_land(dem = DEM, 
                                  canvas = canvas, 
                                  land = land_mask)),
  tar_terra_rast(sea, terra::ifel(land == 0, 0, NA)),
  # Calculate coast distance
  # Now we need to combine our land data with the `elev` data to
  # create the raster with which we are going to do distance
  # calculations with. The `elev` data will cut off barriers
  # that the algorithm will have to travel around when 
  # calculating distance, while the `canvas` area will supply
  # land and sea data.
  #   - Mountain barriers -> traveling around
  #   - Land areas -> traveling through
  #   - Sea areas -> traveling from
  # We will need to employ some if/else logic to correctly combine
  # these three data layers.
  # NOTE: IF WE INCLUDE ALASKA BORDER REGION LATER, we will need to 
  # include the Alaska DEM in this to correctly measure distance from
  # coast for the Alaska Border region. 
  tar_terra_rast(coast_dist, coast_distance(elev = elev, sea = sea, dist_km = 30))
)
