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

# Set target options:
tar_option_set(
  packages = c("MAMU", # remotes::install_github("popovs/MAMU")
               "janitor",
               "units",
               "sf",
               "terra",
               "stars",
               "rnaturalearth",
               "dplyr",
               "readxl",
               "smoothr"),
  format = tar_format( # Default qs is superceded by qs2. Install qs2 and specify read & write fxns
    read = function(path) { qs2::qs_read(path) }, 
    write = function(object, path) { qs2::qs_save(object = object, file = path) }
    ), 
  memory = "transient", # unload memory for each target line after it's completed
  garbage_collection = TRUE, # run gc() prior to each target
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

#### STATIC PIPELINE OBJECTS ####

# Set resolution (in meterse) for all raster analysis
# Per the BC DEM website, this raster product is at a 25m 
# resolution, but gridded to a 0.75 arc-second scale. 
# https://www2.gov.bc.ca/gov/content/data/geographic-data-services/topographic-data/elevation/digital-elevation-model
# This means all our raster data is at a somewhat odd 
# 17.37227 x 17.37227 resolution:
#res(BC_DEM_3005)
# As such, the highest resolution we can safely go for is 25m.
# Note that a smaller number = MUCH SLOWER SCRIPT
res <- 250 # res MUST be >25, as the DEM goes does to 25m accuracy.

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

# Download the University of Maryland global forest cover dataset
# https://glad.earthengine.app/view/global-forest-change#bl=off;old=0;dl=off;lon=-486.9726343690056;lat=51.450799512228464;zoom=6;
# Forest cover tiles
forest_urls <- data.frame(names = c("60N_140W", "60N_130W", "50N_130W"),
                          url = c("https://storage.googleapis.com/earthenginepartners-hansen/GFC-2023-v1.11/Hansen_GFC-2023-v1.11_treecover2000_60N_140W.tif",
                                  "https://storage.googleapis.com/earthenginepartners-hansen/GFC-2023-v1.11/Hansen_GFC-2023-v1.11_treecover2000_60N_130W.tif",
                                  "https://storage.googleapis.com/earthenginepartners-hansen/GFC-2023-v1.11/Hansen_GFC-2023-v1.11_treecover2000_50N_130W.tif"))
# File directory to store forest data
forest_dir <- "GIS/Forest_cover"

#### PIPELINE ####
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
               # nest_quantiles(nests, 
               #                quant_data = nest_cost,
               #                prefix = "cost"),
               # A few outliers really skew the quantiles. Instead,
               # grab the simple min/max of each region once outliers
               # are removed. Outliers are define w a simple boxplot
               # whisker - if they are outside the 1.5*IQR range, it's 
               # an outlier.
               nest_minmax_sans_outliers(nests, 
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
  tar_combine(cc_sprc, mapped[[3]], command = terra::sprc(list(!!!.x)) |> terra::wrap()), # need to wrap it to prevent "invalid pointer" error: https://stackoverflow.com/questions/74855695/load-raster-data-with-terra-in-targets-pipeline
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
  tar_terra_rast(coast_dist, coast_distance(elev = elev, sea = sea, dist_km = 30)),
  #### FOREST COVER CUTOFF ####
  # MAMU are almost certainly not nesting in urban or other completely
  # treeless areas. Cut those out.
  # Download the University of Maryland global forest cover dataset
  # https://glad.earthengine.app/view/global-forest-change#bl=off;old=0;dl=off;lon=-486.9726343690056;lat=51.450799512228464;zoom=6;
  tar_map(
    values = forest_urls, # params need to be passed as a df/tibble, defined OUTSIDE the pipeline
    names = names,
    # Download DEM tiles
    tar_target(forest_tiles,
               download_forest_tiles(url = url,
                                     output_dir = forest_dir),
               format = "file")
  ),
  # Track files within the forest_dir
  tar_target(forest_tile_files, list.files(file.path(forest_dir), full.names = TRUE)),
  # Make VRT
  tar_target(forest_VRT, 
             terra::vrt(x = forest_tile_files,
                        filename = file.path(forest_dir, "forest.vrt"),
                        overwrite = TRUE,
                        return_filename = TRUE), # for format "file" targets, the output MUST be a filepath
             format = "file"),
  # Merge forest tiles and reproject to 3005 (~30 mins)
  # TODO: something wrong with output.
  tar_target(forest_cover, 
             merge_dem(vrt_path = forest_VRT,
                       output_file = file.path(forest_dir, "forest_cover.tiff"),
                       overwrite = TRUE),
             format = "file"),
  # Resample to pipeline resolution (but no need to save it as its own tiff file)
  tar_terra_rast(forest, resample_dem(dem_path = forest_cover, res = res)),
  # Extract nest forest cover
  tar_target(nest_forest, exactextractr::exact_extract(forest, nests, 'mean')),
  # Apply forest cutoff
  # In the case of forests, we know with 100% confidence that
  # MAMU will nest in areas with 100% tree cover. As such, we
  # don't need to define a 'maximum tree cover' cutoff for 
  # the forested areas layer. Instead, we need to choose a 
  # minimum cutoff. In the same vein, setting a variable minimum
  # forest cover cutoff by region might be cutting out too much
  # habitat area. 
  # Instead, we'll just cut out the bottom 2.5% tree cover (6
  # nests out of 242 in the dataset).
  # `fc` for forest cutoff
  tar_terra_rast(fc, terra::ifel(forest < quantile(nest_forest, 0.025, na.rm = TRUE), NA, 1)),
  #### MERGE CUTOFFS ####
  # Finally, merge the cutoff rasters together to create a 
  # "MAMU containment zone" area that has a high probability of
  # containing high quality MAMU nesting habitat. This raster 
  # will be used as our maximum MAMU area within BC that we will 
  # extrapolate our population estimates to.
  # All the raster prep functions above set the resolution to match
  # `res` and the extent to match `DEM`. Therefore we can safely
  # layer all our raster layers.
  # Create MAMU Accessible Zone (MAZ)
  # First layer - 'MAMU accessible', i.e. it's below elevation cutoff,
  # within 30km of ocean, and within cost distance. Does not 
  # necessarily imply it's all suitable nesting habitat; rather, 
  # this area is assumed to encompass ~95% of all suitable
  # nesting habitat within each region.
  tar_terra_rast(maz, ((elev > 0 ) * (cc > 0) * (coast_dist > 0))), # Anything accessible was stored as either '1' or '2' in each raster layer.
  #### GENERATE RADAR CONES ####
  # Here, radar survey 'catchments' containing marbled murrelet 
  # nesting habitat will be generated, using a few basic 
  # assumptions about bird habitat + BC Digital Elevation Model 
  # (DEM) data. 
  # First, we read in bird headings data. We're using a simplified
  # assumption here and taking the mean of all the bird flight
  # headings here. While the heading on the radar screen might
  # not necessarily translate 1-to-1 to the real world flight
  # direction, the error in this will be captured by the cone
  # we generate around the mean heading.
  # Read in MAMU radar survey data
  tar_target(s_file, "data/ECCC_FLNR_MAMU-RadarData-20240307.xlsx", format = "file"),
  tar_target(s, prepare_surveys(s_file)),
  # Extract individual station coords
  tar_target(stn, prepare_stn(s, regions)),
  # Read in flight headings data
  tar_target(h_file, "data/headings.xlsx", format = "file"),
  tar_target(h_0, prepare_headings(h_file)),
  # Calculate polar mean flight headings
  tar_target(h, calc_polar_mean(headings = h_0, n_reps = 1000, alpha = 0.05)),
  # Calculate cones
  tar_target(cones, generate_cones(h = h, stn = stn, radius = 30000, res = res)),
  tar_target(cones_gpkg, 
             save_sf(sf = cones, output_path = "GIS/flight_headings.gpkg"), 
             format = "file"),
  #### SELECT TARGETED WATERSHED CATCHMENTS ####
  # Now, based on the MAMU flight headings at each radar station,
  # select the watershed catchments the birds are targeting (i.e.,
  # the watersheds the radar cones overlap with).
  tar_target(watersheds_file, "GIS/Watersheds/code_2_watersheds.gpkg", format = "file"),
  tar_target(watersheds_raw, sf::st_read(watersheds_file)),
  tar_target(watersheds, select_watersheds(watersheds = watersheds_raw, 
                                           cones = cones, 
                                           min_cone_coverage = 0.02,
                                           output_plots = TRUE,
                                           output_dir = "temp/cone_inspection",
                                           stn = stn,
                                           headings = h_0,
                                           h = h)),
  # Calculate how much it costs to fly within the selected watersheds
  tar_target(full_cc, watershed_cost(watersheds = watersheds,
                                             dem = DEM,
                                             cones = cones,
                                             stn = stn, 
                                             cost_cutoffs = cost_cutoffs,
                                             output_dir = "temp/cost_inspection",
                                             headings = h_0,
                                             h = h,
                                             nests = nests)),
  # Cut out pieces behind heading 
  # Birds aren't flying backwards from radar station
  tar_target(cropped_cc, directionality_crop(cost_catchments = full_cc,
                                             stn = stn, 
                                             h = h,
                                             watersheds = watersheds,
                                             cones = cones, 
                                             res = res)),
  # Intersect with MAMU-accessible areas
  tar_target(final_cc, access_catchments(cost_catchments = cropped_cc, 
                                         maz = maz, 
                                         stn = stn,
                                         raster_stats = TRUE,
                                         forest = forest,
                                         cost = cost,
                                         output_dir = "temp/final_catchment_inspection",
                                         headings = h_0,
                                         h = h,
                                         watersheds = watersheds,
                                         cones = cones,
                                         nests = nests)),
  tar_target(cc_gpkg, 
             save_sf(sf = final_cc, output_path = "GIS/radar_derived_catchments.gpkg"), 
             format = "file")
) 
