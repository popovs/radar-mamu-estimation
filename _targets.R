# HEAD --------------------------------------------------------------------

# This is the 'mothership' script of this R repository. From here, you
# track all your data inputs, manage the data workflow (the 'pipeline'),
# and create all your outputs. This pipeline is built using the {targets}
# package: https://books.ropensci.org/targets/

# Note this repository DOES NOT contain all the raw data files
# needed to successfully reproduce the data pipeline on another machine.
# You will need to contact Sarah Popov (sarah.popov@gov.bc.ca) to receive
# all the raw data (large GIS files, nest locations, etc.) to run this
# successfully on your machine.

# BASIC USAGE
# Assuming you have set up all the raw data files and directories, ...
# 1) Run `renv::restore()` to install all necessary R packages the pipeline
#     depends on.
# 2) Load the {targets} library. Run `tar_make()` to run the pipeline.
#     NOTE the full pipeline takes 5-6 hours to run if running from scratch.
#     You can comment out sections of the pipeline that you are not
#     interested in recreating if you wish to skip the creation of them,
#     assuming nothing downstream of the pipeline depends upon it.
#     TIP: you can run `tar_visnetwork()` to see how various targets depend
#     upon each other.
# 3) Once you have created all your targets, you are ready to play with
#     the data outputs. In a separate R script, you can run 
#     `tar_load(<target name>)` to quickly load up the target in your R
#     session and manipulate it from there.

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

# SETUP -------------------------------------------------------------------

library(targets)
library(tarchetypes) # Needed for `tar_map()`
library(geotargets) # Needed to save `terra` SpatRaster targets

# Set target options:
tar_option_set(
  packages = c("sf", "terra"), # packages needed for static pipeline obj
  format = tar_format( # Default qs is superceded by qs2. Install qs2 and specify read & write fxns
    read = function(path) { qs2::qs_read(path) }, 
    write = function(object, path) { qs2::qs_save(object = object, file = path) }
  ), 
  memory = "transient", # unload memory for each target line after it's completed
  garbage_collection = TRUE, # run gc() prior to each target
)

# Let {targets} find custom fxns located within the 'R' folder
# when running `tar_make()`:
tar_source()


#### STATIC PIPELINE OBJECTS ####

# Set resolution (in meters) for all raster analysis
# Per the BC DEM website, this raster product is at a 25m 
# resolution, but gridded to a 0.75 arc-second scale. 
# https://www2.gov.bc.ca/gov/content/data/geographic-data-services/topographic-data/elevation/digital-elevation-model
# This means all our raster data is at a somewhat odd 
# 17.37227 x 17.37227 resolution:
#res(BC_DEM_3005)
# As such, the highest resolution we can safely go for is 25m.
# Note that a smaller number = MUCH SLOWER SCRIPT
res <- 250 # Set res to ~250 to run things much faster

# Regions to iterate GIS operations over
#regions_map <- sf::st_drop_geometry(sf::st_read("GIS/regions.gpkg"))

# Bounding box of our entire study area
study_bounds <- sf::st_bbox(c(xmin = 164728, ymin = 333912, xmax = 1387220, ymax = 1778675),
                            crs = 3005)

# Get land area for AK and WA
# This will be used to clean up raster boundaries of the AK raster +
# block off land areas of WA
# usa_land <- rnaturalearth::ne_states("united states of america") |>
#   dplyr::filter(name %in% c("Alaska", "Washington")) |>
#   sf::st_as_sf() |>
#   sf::st_transform(3005) |>
#   sf::st_crop(study_bounds) |>
#   terra::vect() # transform to terra SpatVect obj to play better w later terra raster objects


# API tokens
# Necessary for plotting fxns
#source("temp/apikey.R")

# SAVE PLOTS? Yes or no
save_plots <- FALSE


# PIPELINE ----------------------------------------------------------------

list(
  #### PREPARE DATA FILES ####
  ##### Regions #####
  tar_target(regions_path, "GIS/cons_reg.gpkg", format = "file"), # Track regions gpkg file
  tar_target(regions, prepare_regions(filepath = regions_path)),
  tar_target(regions_region, regions$region), # regions names within regions gpkg
  ##### Nests #####
  tar_target(nests_path, "GIS/MAMU_Nests_BC_CDC.gpkg", format = "file"), # Track nests gpkg file
  # nests_geom target will contain just the lat/longs of nests + region they're in
  # Later, we will add extracted raster values to this dataset.
  tar_target(nests_geom, prepare_nests(filepath = nests_path,
                                       regions = regions)),

  ##### Radar surveys #####
  # TODO: maybe convert this to csv?
  tar_target(s_path, "data/ECCC_FLNR_MAMU-RadarData-20240307.xlsx", format = "file"), # Track surveys excel file
  
  # Alt arrangement: track all paths first in one section, 
  # Then prep all inputs in another section
  # tar_target(regions, prepare_regions(filepath = regions_path)),
  # tar_target(s, prepare_surveys(filepath = s_path,
  #                               regions = regions)),
  # tar_target(nests, prepare_nests(filepath = nests_path,
  #                                 regions = regions)),
  
  ##### Land area masks #####
  # This study assumes MAMU can access land from any sea area on 
  # the map. Therefore, need to effectively split land from sea.
  # Create vector masks to accurately mask land vs sea areas of the DEM.
  # While the BC DEM is quite accurate, the AK DEM in particular has some
  # anomalously high sea areas (e.g., >50 m "above sea level" for some
  # ocean patches). These need to be masked manually.
  ###### AK mask #####
  # First let's query USA areas.
  tar_terra_vect(usa_vect, get_usa_land(extent = study_bounds)), # `study_bounds` is defined in the static pipeline objects
  # Buffer the `regions` polygon by 10km, so there are no gaps when
  # it is merged with the `usa_vect` polygon
  tar_terra_vect(reg_vect, regions |>
                   dplyr::filter(region != "HG") |>
                   sf::st_buffer(10000) |>
                   terra::vect()),
  # Create a polygon to mask the AK DEM to land areas only, and crop
  # out ocean areas.
  # The reason we don't bother doing this for the BC DEM is because
  # it is quite accurate for sea vs land elevation.
  tar_terra_vect(ak_mask, terra::union(usa_vect, reg_vect)),
  ###### BC interior mask ######
  # Next create polygon to mask interior BC land areas
  tar_target(bc_interior_vect_path, "GIS/DEM/land_mask.gpkg", format = "file"),
  tar_terra_vect(bc_interior_vect, sf::st_read(bc_interior_vect_path) |>
                   sf::st_crop(study_bounds) |>
                   terra::vect()),
  tar_terra_vect(land_mask, terra::union(usa_vect, bc_interior_vect)),
  ###### Study region mask ######
  # Finally, create a mask for the entire study region - this is
  # to crop the DEM to our `regions` polygon + AK coast and is 
  # purely for aesthetic purposes, to trim off the jagged DEM 
  # edges. 
  tar_terra_vect(study_region_mask, terra::union(terra::buffer(usa_vect, 10000),
                                                 terra::vect(regions))),
  
  #### PREPARE DEM ####
  ##### Query BC CDED tiles #####
  # Note this target simply points to the DEM VRT filepath - it is not 
  # a raster in and of itself.
  tar_target(CDED_VRT_path, 
             query_cded(regions = regions,
                        output_dir = "GIS/DEM"),
             format = "file"),
  ##### Prepare BC coast DEM #####
  # Resample to target resolution and project to EPSG 3005
  # This takes about ~30 mins at 250m resolution
  tar_terra_rast(BC_DEM, terra::rast(CDED_VRT_path) |> # read raster
                   terra::project("epsg:3005") |> # reproject to m-based projection
                   resample_rast(res = res)), # resample to target resolution
  ##### Prepare Alaska coast DEM #####
  # Resample to target resolution and project to EPSG 3005,
  # then mask/select to land areas only
  tar_target(AK_DEM_path, "GIS/DEM/Alaska_DEM.tiff", format = "file"),
  tar_terra_rast(AK_DEM, terra::rast(AK_DEM_path) |> # read raster
                   terra::project("epsg:3005") |> # reproject to m-based projection
                   resample_rast(res = res) |> # resample to target resolution
                   terra::mask(ak_mask) |> # mask to land areas only
                   {\(.) terra::ifel(. > 0, ., NA)}()), # select only pixels whose elevation > 0. See here for `\(.)` notation syntax: https://stackoverflow.com/a/76422551/1454785
  ##### Full DEM #####
  tar_terra_rast(DEM, terra::merge(BC_DEM, AK_DEM, first = TRUE) |> # merge BC and AK DEMs. In cases where they overlap, take the BC_DEM value.
                   terra::mask(study_region_mask) |> # mask to `regions` polygon + AK coast. Purely for aesthetic purposes, to cut off jagged DEM edges.
                   {\(.) terra::ifel(. > 0, ., NA)}()), # select only pixels whose elevation > 0. See here for `\(.)` notation syntax: https://stackoverflow.com/a/76422551/1454785
  
  #### LAND + SEA AREAS ####
  # Using a combination of the land area masks and DEMs, create
  # rasters of land area and sea area.
  # Create raster canvas of target area
  tar_terra_rast(canvas, create_canvas(target_rast = DEM)),
  ##### Land #####
  tar_terra_rast(land, block_land(dem = DEM,
                                  canvas = canvas,
                                  land = land_mask)),
  ##### Sea #####
  tar_terra_rast(sea, terra::ifel(land == 0, 0, NA)),
  
  #### DISTANCE TO SEA ####
  ##### Sea distance #####
  # Distance from sea, in km
  tar_terra_rast(sea_dist, distance(sea) / 1000),
  ##### Cost distance #####
  # Elevation * Distance from coast
  tar_terra_rast(cost, terra::costDist(terra::merge(DEM, sea), 
                                       target = 0) |> 
                   {\(.) terra::ifel(. > 0, ., NA)}()),
  
  #### EXTRACT NEST VALUES ####
  ##### Elevation, distance from sea, cost #####
  # Elevation
  tar_target(nest_elev_m, terra::extract(DEM, nests_geom, ID = FALSE)[[1]]),
  # Distance from sea
  tar_target(nest_dist_km, terra::extract(sea_dist, nests_geom, ID = FALSE)[[1]]),
  # Cost
  tar_target(nest_cost, terra::extract(cost, nests_geom, ID = FALSE)[[1]]),
  ##### Merge nest data #####
  # Merge into a single dataset
  tar_target(nests, cbind(nests_geom, nest_elev_m, nest_dist_km, nest_cost)),

  #### REGIONAL CUTOFFS ####
  ##### Elevation cutoffs #####
  # Extract nest elevations
  # The exact elevational cutoffs vary by region and are primarily
  # dictated by the growing conditions of the trees that MAMU
  # prefer to nest in. Large coniferous trees that can support
  # MAMU nests can grow in higher elevations at lower latitudes.
  # Conversely, the further north you go, the lower the maximum
  # MAMU nesting elevation. The nest data is the main source of
  # cutoffs. Adding more nest data will trigger a re-run of this
  # analysis.
  # Using the nest elevations extracted in the previous section,
  # calculate and store regional elevation cutoffs in a table
  tar_group_by(elevation_cutoffs,
               # OPTION A: 95% QUANTILES
               # nest_quantiles(nests,
               #                quant_data = nest_elev_m,
               #                prefix = "elev_m"),
               # OPTION B: ISOLATION TREE OUTLIER DETECTION
               # nest_isoforest(nests,
               #                quant_data = nest_elev_m,
               #                prefix = "elev_m"),
               # OPTION C: NO OUTLIERS REMOVED
               aggregate(nest_elev_m ~ region, nests, FUN = "max") |>
                 dplyr::mutate(elev_m_max = ceiling((nest_elev_m) / 100) * 100) |>
                 tidyr::complete(region) |>
                 # If no nests present in certain regions, fill in with data from other regions
                 # If no nests in NC, use the same value as CC
                 # If no nests in AKB, use the same value as CC
                 # If no nests in NVI, use mean value of all other VI regions
                 dplyr::mutate(elev_m_max = dplyr::replace_when(elev_m_max, 
                                                                region == "NC" & is.na(nest_elev_m) ~ max(elev_m_max[region == "CC"]),
                                                                region == "AKB" & is.na(nest_elev_m) ~ max(elev_m_max[region == "CC"]),
                                                                region == "NVI" & is.na(nest_elev_m) ~ mean(elev_m_max[grepl("VI", region)]))) |>
                 dplyr::mutate(elev_m_min = 0) |>
                 dplyr::select(region, elev_m_min, elev_m_max),
               region), # group by region col so targets later knows to run `reclass_elevation()` by row-level grouping
  
  ##### Cost cutoffs #####
  # Next, we calculate the flight cost values of each nest.
  # Evidence shows that MAMU take the least-cost flight paths to
  # their nests; that is, they tend to fly along valley contours
  # rather than straight across ridges (even if the ridges are
  # below their nest cutoff elevations). Here we will calculate a
  # raster of the flight cost (elevation * distance from coast)
  # to generate a landscape where birds are more or less likely to
  # nest. This will be used to: 1) delineate the nesting catchment
  # boundaries later and 2) eliminate "high cost" nesting areas
  # that may still be included within the elevation cutoff and
  # coast distance rasters.
  # Using the nest elevations extracted in the previous section,
  # calculate and store regional cost cutoffs in a table
  tar_group_by(cost_cutoffs,
               # OPTION A: 95% QUANTILES
               # nest_quantiles(nests,
               #                quant_data = nest_cost,
               #                prefix = "cost"),
               # OPTION B: MAX VALUE AFTER GLOBAL OUTLIERS REMOVED
               # A few outliers really skew the quantiles. Instead,
               # grab the simple min/max of each region once outliers
               # are removed. Outliers are define w a simple boxplot
               # whisker - if they are outside the 1.5*IQR range, it's
               # an outlier.
               # nest_minmax_sans_outliers(nests,
               #                           quant_data = nest_cost,
               #                           prefix = "cost"),
               # OPTION C: REGIONAL IQR AFTER GLOBAL OUTLIERS REMOVED
               # nest_iqr_sans_outliers(nests = nests,
               #                        quant_data = nest_cost,
               #                        prefix = "cost",
               #                        hg_exception = TRUE),
               # OPTION D: ISOLATION TREE OUTLIER DETECTION
               # nest_isoforest(nests = nests,
               #                quant_data = nest_cost,
               #                prefix = "cost"),
               # OPTION X: NO OUTLIERS REMOVED
               aggregate(nest_cost ~ region, nests, FUN = "max") |>
                 dplyr::mutate(cost_max = ceiling((nest_cost) / 100) * 100) |>
                 tidyr::complete(region) |>
                 # If no nests present in certain regions, fill in with data from other regions
                 # If no nests in NC, use the same value as CC
                 # If no nests in AKB, use the same value as CC
                 # If no nests in NVI, use mean value of all other VI regions
                 dplyr::mutate(cost_max = dplyr::replace_when(cost_max, 
                                                              region == "NC" & is.na(nest_cost) ~ max(cost_max[region == "CC"]),
                                                              region == "AKB" & is.na(nest_cost) ~ max(cost_max[region == "CC"]),
                                                              region == "NVI" & is.na(nest_cost) ~ mean(cost_max[grepl("VI", region)]))) |>
                 #dplyr::mutate(cost_max = cost_max / 1000 / 1000) |> # Re-express cost in km^2 rather than m^2
                 dplyr::mutate(cost_min = 0) |>
                 dplyr::select(region, cost_min, cost_max),
               region),
  
  #### ITERATE OVER REGIONAL CUTOFFS ####
  # For each region, iterate the following:
  # 1. Chop up the raster into regions
  # 2. Apply elevation cutoffs to each regional raster
  # 3. Merge into single raster

  ##### Regional elevation cutoffs #####
  # Chop up DEM
  tar_terra_rast(regional_DEM, 
                 crop_rast(rast = DEM,
                           region_name = regions_region,
                           regions = regions),
                 pattern = map(regions_region)),
  # Apply elevation cutoffs to each regional DEM
  tar_terra_rast(regional_elev_cutoffs,
                 reclass_rast(rast = regional_DEM,
                              region_name = regions_region,
                              cutoffs = elevation_cutoffs,
                              min_col = "elev_m_min",
                              max_col = "elev_m_max"),
                 pattern = map(regional_DEM, regions_region)),
  # Merge
  # `elev` for 'elevation cutoffs'
  tar_terra_rast(elev, regional_elev_cutoffs |>
                   terra::sprc() |>
                   terra::mosaic() |>
                   terra::resample(DEM)), # resample to same size/resolution as base DEM
  
  ##### Regional cost cutoffs #####
  # Chop up cost raster
  tar_terra_rast(regional_cost,
                 crop_rast(rast = cost,
                           region_name = regions_region,
                           regions = regions),
                 pattern = map(regions_region)),
  # Apply cost cutoffs to each regional cost raster
  tar_terra_rast(regional_cost_cutoffs,
                 reclass_rast(rast = regional_cost,
                              region_name = regions_region,
                              cutoffs = cost_cutoffs,
                              min_col = "cost_min",
                              max_col = "cost_max"),
                 pattern = map(regional_cost, regions_region)),
  # Merge
  # `cc` for 'cost cutoffs'
  tar_terra_rast(cc, regional_cost_cutoffs |>
                   terra::sprc() |>
                   terra::mosaic() |>
                   terra::resample(DEM)), # resample to same size/resolution as base DEM
  
  ##### Coast distance cutoff #####
  # Next we will create a separate raster of all points with
  # 30 km distance from the coast, as BC nest survey data
  # indicates that 99% of MAMU nests are within 30 km of the
  # coastline. We are going to assume a straight-line 30 km
  # distance from the shore, ignoring any barriers.
  # This won't have any regional variation aside from
  # Vancouver Island - we are going to assume that all of the
  # island is regularly accessed for nests by MAMU.
  # `c30km` for 'cutoff - 30km' 
  tar_terra_rast(c30km, sea_dist |>
                   terra::crop(regions[grepl("VI", regions$region), ],
                               mask = TRUE) |>
                   {\(.) terra::ifel(. > 0, 1, NA)}() |>
                   terra::merge(sea_dist, first = TRUE) |>
                   {\(.) terra::ifel(. > 0, sea_dist, NA)}() |>
                   {\(.) terra::ifel(. <= 30, 1, NA)}())
  
  #### MERGE CUTOFFS ####
  # Finally, merge the cutoff rasters together to create a
  # "MAMU containment zone" area that meets minimum threshold criteria
  # for spatially delineating MAMU nesting habitat survey catchments.
  # All the raster prep functions above set the resolution to match
  # `res` and the extent to match `DEM`. Therefore we can safely
  # layer all our raster layers.
  # Create MAMU Accessible Zone (MAZ)
  # 'MAMU accessible', i.e. it's below elevation cutoff,
  # within 30km of ocean, and within cost distance. Does not
  # necessarily imply it's all suitable nesting habitat; rather,
  # this area is assumed to encompass most of the suitable
  # nesting habitat within each region. Finally, we assume that
  # any sea areas are accessible as well.
  # NOTE: cost cutoffs are so conservative they don't make any difference to MAZ
  # i.e. [(elev > 0) * (cc > 0) * (coast_dist > 0)] == [(elev > 0) * (coast_dist > 0)]
  # tar_terra_rast(maz, terra::cover(((elev > 0 ) * (coast_dist > 0)), # Anything accessible was stored as either '1' or '2' in each raster layer.
  #                                  sea + 1)), # and then we add the sea (+1, because it's all stored as 0 or NA)
  
  

  
  
# QUARANTINE --------------------------------------------------------------

# #### GENERATE RADAR CONES ####
  # # Here, radar survey 'catchments' containing marbled murrelet 
  # # nesting habitat will be generated, using a few basic 
  # # assumptions about bird habitat + BC Digital Elevation Model 
  # # (DEM) data. 
  # # First, we read in bird headings data. We're using a simplified
  # # assumption here and taking the mean of all the bird flight
  # # headings here. While the heading on the radar screen might
  # # not necessarily translate 1-to-1 to the real world flight
  # # direction, the error in this will be captured by the cone
  # # we generate around the mean heading.
  # # Read in flight headings data
  # tar_target(h_file, "data/headings.xlsx", format = "file"),
  # tar_target(h_0, prepare_headings(h_file)),
  # # Extract individual station coords
  # tar_target(stn, prepare_stn(s, h_0, regions)),
  # tar_target(stn_gpkg, 
  #            save_sf(sf = stn, output_path = "temp/stn.gpkg"), 
  #            format = "file"),
  # # Calculate polar mean flight headings
  # tar_target(h, calc_polar_mean(headings = h_0, n_reps = 1000, alpha = 0.05)),
  # # Calculate cones
  # tar_target(cones, generate_cones(h = h, 
  #                                  stn = stn, 
  #                                  radius = 30000, # length cones to select appropriate watersheds (in meters)
  #                                  res = res)),
  # tar_target(cones_gpkg, 
  #            save_sf(sf = cones, output_path = "GIS/flight_headings.gpkg"), 
  #            format = "file"),
  # #### SELECT TARGETED WATERSHED CATCHMENTS ####
  # # Now, based on the MAMU flight headings at each radar station,
  # # select the watershed catchments the birds are targeting (i.e.,
  # # the watersheds the radar cones overlap with).
  # #tar_target(watersheds_file, "GIS/Watersheds/code_2_watersheds.gpkg", format = "file"),
  # #tar_target(watersheds_file, "GIS/Watersheds/WSA_WS_SVW_polygon.shp", format = "file"),
  # tar_target(watersheds_file, "GIS/Watersheds/code_2_HG_merge.gpkg", format = "file"),
  # tar_target(watersheds_raw, prepare_watersheds(watersheds_file, regions)),
  # tar_target(watersheds, select_watersheds(watersheds = watersheds_raw, 
  #                                          cones = cones, 
  #                                          min_cone_coverage = 0.01,
  #                                          output_plots = save_plots, # defined at the top of script
  #                                          output_dir = "temp/cone_inspection",
  #                                          stn = stn,
  #                                          headings = h_0,
  #                                          h = h)),
  # # Calculate how much it costs to fly within the selected watersheds
  # tar_target(full_cc, watershed_cost(watersheds = watersheds,
  #                                    dem = DEM,
  #                                    cones = cones,
  #                                    stn = stn, 
  #                                    cost = cost,
  #                                    output_plots = save_plots, # defined at the top of script
  #                                    cost_cutoffs = cost_cutoffs,
  #                                    output_dir = "temp/cost_inspection",
  #                                    headings = h_0,
  #                                    h = h,
  #                                    nests = nests)),
  # # Cut out pieces behind heading 
  # # Birds aren't flying backwards from radar station
  # tar_target(cropped_cc, directionality_crop(cost_catchments = full_cc,
  #                                            stn = stn, 
  #                                            h = h,
  #                                            watersheds = watersheds,
  #                                            cones = cones, 
  #                                            res = res)),
  # # Intersect with MAMU-accessible areas
  # tar_target(final_cc, access_catchments(cost_catchments = cropped_cc, 
  #                                        maz = maz, 
  #                                        stn = stn,
  #                                        raster_stats = TRUE,
  #                                        cost = cost,
  #                                        output_plots = save_plots, # defined at top of script
  #                                        output_dir = "temp/final_catchment_inspection",
  #                                        headings = h_0,
  #                                        h = h,
  #                                        watersheds = watersheds,
  #                                        cones = cones,
  #                                        nests = nests)),
  # tar_target(cc_gpkg, 
  #            save_sf(sf = final_cc, output_path = "GIS/radar_derived_catchments.gpkg"), 
  #            format = "file"),
  # #### CATCHMENTS x HABITAT ####
  # # Intersect with 2024 MAMU suitable habitat layer
  # tar_target(suitable_habitat, prepare_mamu_habitat_gpkg(path = "GIS/2024_suitable_habitat/",
  #                                                   maz = maz,
  #                                                   regions = regions)),
  # # Intersect with 2025 MAMU suitable habitat layer
  # # tar_terra_rast(suitable_habitat, prepare_mamu_habitat_tiff(path = "GIS/2025_suitable_habitat/mamu_predict_2025_feb_03.tif",
  # #                                                            maz = maz,
  # #                                                            band = "X1")),
  # # Intersect suitable habitat in each catchment
  # tar_target(cc_habitat, st_habitat_in_sf(sf = final_cc,
  #                                         habitat = suitable_habitat,
  #                                         use_probability = TRUE)),
  # # Intersect suitable habitat in each region
  # tar_target(reg_habitat, st_habitat_in_sf(sf = regions,
  #                                          habitat = suitable_habitat,
  #                                          use_probability = TRUE)),
  # # Final visualization
  # # tar_render(final_visualization,
  # #            path = "Rmd/catchments_visualization.Rmd",
  # #            output_file = "catchments_visualization.pdf",
  # #            #error = "null",
  # #            quiet = TRUE,
  # #            params = list(catchments = final_cc,
  # #                          watersheds = watersheds_raw,
  # #                          suitable_habitat = suitable_habitat,
  # #                          cones = cones,
  # #                          stn = stn,
  # #                          nests = nests,
  # #                          apikey = jawg_token,
  # #                          headings = h_0,
  # #                          h = h)
  # #            ),
  # #### ANNUAL MAX MAMU COUNTS ####
  # # Select maximum MAMU count per station per year
  # tar_target(annual_max_mamu, s |> 
  #              sf::st_drop_geometry() |>
  #              dplyr::group_by(site, region, year) |> 
  #              dplyr::summarise(max_mamu = max(mamuinpd, na.rm = TRUE),
  #                               N = dplyr::n())),
  # # Calculate bootstrapped mean annual maximum per catchment
  # # (across all years)
  # tar_target(mean_max_mamu_cc, bootmean(annual_max_mamu, 
  #                                       group_by = "site",
  #                                       dat_col = "max_mamu", 
  #                                       CI_level = 0.95) |>
  #              dplyr::mutate(dplyr::across(c(bootmean, boot_min, boot_max), 
  #                                          round))),
  # # Regional mean annual maximum per catchment
  # # (across all years) for paper table purposes
  # tar_target(mean_max_mamu_reg, bootmean(annual_max_mamu, 
  #                                       group_by = "region",
  #                                       dat_col = "max_mamu", 
  #                                       CI_level = 0.95) |>
  #              dplyr::mutate(dplyr::across(c(bootmean, boot_min, boot_max), 
  #                                          round))),
  # #### DENSITY CALCS ####
  # # NOTE suitable habitat =/= even MAMU density across whole layer!
  # # While interior habitat might be equally 'suitable' to coastal habitat, the 
  # # habitat closer to sea will have more MAMU population than more interior sites.
  # # Total habitat (ha) across whole suitable habitat layer
  # tar_target(total_suit_hab_area_ha, sum(reg_habitat$sh_area_ha)), 
  # # Habitat (ha) summarized by region
  # tar_target(regional_suit_hab_area_ha, sf::st_drop_geometry(reg_habitat)), 
  # # Calculate the mean density of birds within each catchment
  # # Using the bootstrapped mean + upper + lower CIs
  # tar_target(cc_density, catchment_density(mm = mean_max_mamu_cc, 
  #                                          catchment_habitat = cc_habitat)),
  # # Calculate the weighted mean density of birds within each region
  # # Catchments that take up more area of the region have higher weight
  # tar_target(reg_density, cc_density |> 
  #              dplyr::filter(region == "NC") |>
  #              dplyr::mutate(region = "AKB") |> # Add in dummy rows for AKB, using NC density
  #              dplyr::bind_rows(cc_density) |>
  #              dplyr::group_by(region) |> 
  #              dplyr::summarise(N = dplyr::n(),
  #                               bootmean = sum(bootmean), 
  #                               boot_min = sum(boot_min, na.rm = TRUE), 
  #                               boot_max = sum(boot_max, na.rm = TRUE), 
  #                               area_ha = sum(area_ha), 
  #                               sh_area_ha = sum(sh_area_ha)) |>
  #              dplyr::mutate(density = bootmean / sh_area_ha, 
  #                            density_lwr = boot_min / sh_area_ha, 
  #                            density_upr = boot_max / sh_area_ha)),
  # # Fit a nest gamma decay function + rasterize it
  # # Certain habitats may meet the criteria for "suitable habitat".
  # # However, while habitat may be suitable for MAMU nesting in terms
  # # of tree species composition, the suitable habitat layers do not
  # # explicitly take into account the fact that ~99% of nests occur
  # # within 30km of the coastline, and, crucially, that the further
  # # from the coast you are, the less likely the nests are likely to
  # # occur. The nest data follow a gamma distribution of likelihood
  # # vs distance from shore. So, fit a gamma distribution to the nest
  # # data, and then map that gamma distribution decay curve to a 
  # # raster. Cells <30km from shore will have a higher probability,
  # # closer to 1, while distances >30km will decay down to 0 probability.
  # # `nest_likelihood` has cut out all non-habitat pieces
  # tar_terra_rast(nest_likelihood, nest_gamma_decay(nests = nests,
  #                                                  coast = sea,
  #                                                  habitat = suitable_habitat)),
  # # `nest_likelihood_full` is primarily for visualization purposes
  # tar_terra_rast(nest_likelihood_full, nest_gamma_decay(nests = nests,
  #                                                       coast = sea)),
  # # Rasterize the regional mean density
  # # Apply the nest probability decay function to regional densities
  # # Now, apply that gamma decay function to the density layer such 
  # # that the mean density per region will remain the same as calculated
  # # in `reg_density`, but the *spatial pattern* of the density follows
  # # the nest gamma distribution.
  # tar_terra_rast(density_map, extrapolate_density(regions, 
  #                                                 reg_density, 
  #                                                 nest_likelihood)),
  # #### POPULATION CALCS ####
  # # Calculate MAMU population!
  # tar_target(regional_population, calculate_population(density_map, 
  #                                                      sf = regions, 
  #                                                      merge_df = reg_density) |>
  #              dplyr::arrange(region)),
  # tar_target(total_population, calculate_population(density_map)),
  # tar_target(total_density, colSums(cc_density[,c("bootmean", "boot_min", "boot_max")], na.rm = TRUE) / sum(cc_density$sh_area_ha)),
  # tar_target(bc_density, total_population / total_suit_hab_area_ha)
) 
