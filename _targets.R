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
               "isotree",
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

# Set resolution (in meters) for all raster analysis
# Per the BC DEM website, this raster product is at a 25m 
# resolution, but gridded to a 0.75 arc-second scale. 
# https://www2.gov.bc.ca/gov/content/data/geographic-data-services/topographic-data/elevation/digital-elevation-model
# This means all our raster data is at a somewhat odd 
# ~17.37227 x 17.37227 resolution (depending on latitude):
#res(BC_DEM_3005)
# As such, the highest resolution we can safely go for is 25m.
# Note that a smaller number = MUCH SLOWER SCRIPT
res <- 250 # res MUST be >25, as the DEM goes does to 25m accuracy.

# File directory to store DEM data
DEM_dir <- "GIS/DEM"

# Regions to iterate GIS operations over
#regions_map <- sf::st_drop_geometry(sf::st_read("GIS/regions.gpkg")) # if using the other region demarcations
regions_map <- sf::st_drop_geometry(sf::st_read("GIS/cons_reg.gpkg"))

# API tokens
source("temp/apikey.R")

# SAVE PLOTS? Yes or no
save_plots <- TRUE

#### PIPELINE ####
# Run tar_make() to execute the pipeline
list(
  # We track the files themselves to automatically re-prepare the
  # `regions` and `nests` targets if changes in the file are detected.
  tar_target(regions_file, "GIS/cons_reg.gpkg", format = "file"),
  tar_target(s_file, "data/ECCC_FLNR_MAMU-RadarData-20240307.xlsx", format = "file"),
  tar_target(nests_file, "GIS/MAMU_nests.gpkg", format = "file"),
  
  # Read in MAMU radar survey data
  tar_target(regions, prepare_regions(filepath = regions_file)),
  tar_target(s, prepare_surveys(filepath = s_file,
                                regions = regions)),
  tar_target(nests, prepare_nests(filepath = nests_file,
                                  regions = regions)),
  
  # tar_target(regions, prepare_regions()),
  # tar_target(nests, prepare_nests(regions = regions)),
  #### PREPARE DEM ####
  # All the DEM data products will be saved as raster files on the
  # disk rather than simply as _targets objects, so that QGIS can
  # access and plot them. The targets pipeline will track changes
  # to the files. 
  # Download DEM tiles that overlap `regions` & create VRT of them all
  tar_target(DEM_VRT, 
             download_dem_tiles(regions = regions,
                                output_dir = file.path(DEM_dir, "BC_DEM_VRT.vrt")),
             format = "file"),
  # Merge DEM tiles and reproject to 3005 (~30 mins)
  tar_target(BC_DEM_3005, 
             merge_vrt(vrt_path = DEM_VRT,
                       output_file = file.path(DEM_dir, "BC_DEM_EPSG3005.tiff"),
                       overwrite = TRUE),
             format = "file"),
  # Prepare Alaska coast DEM
  tar_target(AK_DEM_path, "GIS/DEM/Alaska_coast_DEM.tif", format = "file"),
  tar_terra_rast(AK_DEM, resample_dem(dem_path = AK_DEM_path,
                                  res = res)),
  # Resample to pipeline resolution (but no need to save it as its own tiff file)
  tar_terra_rast(DEM_target_res, resample_dem(dem_path = BC_DEM_3005, res = res)),
  # Merge BC DEM and AK DEM
  tar_terra_rast(DEM, merge_dem(DEM_target_res, AK_DEM)),
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
                 tibble::add_row(region = "NC", elev_m_max = 800) |> # manually specify NC (== CC)
                 #tibble::add_row(region = "NVI", elev_m_max = 1200) |> # manually specify NVI (mean of all other VI areas)
                 dplyr::mutate(elev_m_min = 0) |>
                 dplyr::select(region, elev_m_min, elev_m_max), 
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
  tar_terra_rast(cost, cost_distance(DEM, sea)),
  # Extract nest cost
  # Same process as extracting elevation values + applying cutoffs
  tar_target(nest_cost, exactextractr::exact_extract(cost, nests, 'mean')),
  # Calculate and store regional cost cutoffs in a table
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
                 tibble::add_row(region = "NC", cost_max = 2900) |> # manually specify NC (== CC)
                 #tibble::add_row(region = "NVI", cost_max = 6900) |> # manually specify NVI (mean of all other VI areas)
                 dplyr::mutate(cost_min = 0) |>
                 dplyr::select(region, cost_min, cost_max),
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
  # coastline. We are going to assume a straight-line 30 km 
  # distance from the shore, ignoring any barriers.
  # The DEM data doesn't include inland BC.
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
  #   - Mountain barriers -> NOT traveling around # 2025-03-01: to travel around mountains, use coast_distance_barriers fxn
  #   - Land areas -> traveling through
  #   - Sea areas -> traveling from
  # We will need to employ some if/else logic to correctly combine
  # these three data layers.
  # NOTE: IF WE INCLUDE ALASKA BORDER REGION LATER, we will need to 
  # include the Alaska DEM in this to correctly measure distance from
  # coast for the Alaska Border region. 
  tar_terra_rast(coast_dist, coast_distance(#elev = elev, 
                                            sea = sea, 
                                            dist_km = 30,
                                            exclude_islands = TRUE, # exclude HG and VI from coast distance calcs - we know birds that access inland regions on HG/VI
                                            regions = regions)),
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
  # nesting habitat within each region. Then we assume that
  # anything over/crossing the sea is accessible as well.
  # NOTE: cost cutoffs are so conservative they don't make any difference to MAZ
  # i.e. [(elev > 0) * (cc > 0) * (coast_dist > 0)] == [(elev > 0) * (coast_dist > 0)]
  tar_terra_rast(maz, terra::cover(((elev > 0 ) * (coast_dist > 0)), # Anything accessible was stored as either '1' or '2' in each raster layer.
                                   sea + 1)), # and then we add the sea (+1, because it's all stored as 0 or NA)
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
  # Read in flight headings data
  tar_target(h_file, "data/headings.xlsx", format = "file"),
  tar_target(h_0, prepare_headings(h_file)),
  # Extract individual station coords
  tar_target(stn, prepare_stn(s, h_0, regions)),
  tar_target(stn_gpkg, 
             save_sf(sf = stn, output_path = "temp/stn.gpkg"), 
             format = "file"),
  # Calculate polar mean flight headings
  tar_target(h, calc_polar_mean(headings = h_0, n_reps = 1000, alpha = 0.05)),
  # Calculate cones
  tar_target(cones, generate_cones(h = h, 
                                   stn = stn, 
                                   radius = 30000, # length cones to select appropriate watersheds (in meters)
                                   res = res)),
  tar_target(cones_gpkg, 
             save_sf(sf = cones, output_path = "GIS/flight_headings.gpkg"), 
             format = "file"),
  #### SELECT TARGETED WATERSHED CATCHMENTS ####
  # Now, based on the MAMU flight headings at each radar station,
  # select the watershed catchments the birds are targeting (i.e.,
  # the watersheds the radar cones overlap with).
  #tar_target(watersheds_file, "GIS/Watersheds/code_2_watersheds.gpkg", format = "file"),
  #tar_target(watersheds_file, "GIS/Watersheds/WSA_WS_SVW_polygon.shp", format = "file"),
  tar_target(watersheds_file, "GIS/Watersheds/code_2_HG_merge.gpkg", format = "file"),
  tar_target(watersheds_raw, prepare_watersheds(watersheds_file, regions)),
  tar_target(watersheds, select_watersheds(watersheds = watersheds_raw, 
                                           cones = cones, 
                                           min_cone_coverage = 0.01,
                                           output_plots = save_plots, # defined at the top of script
                                           output_dir = "temp/cone_inspection",
                                           stn = stn,
                                           headings = h_0,
                                           h = h)),
  # Calculate how much it costs to fly within the selected watersheds
  tar_target(full_cc, watershed_cost(watersheds = watersheds,
                                     dem = DEM,
                                     cones = cones,
                                     stn = stn, 
                                     cost = cost,
                                     output_plots = save_plots, # defined at the top of script
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
                                         cost = cost,
                                         output_plots = save_plots, # defined at top of script
                                         output_dir = "temp/final_catchment_inspection",
                                         headings = h_0,
                                         h = h,
                                         watersheds = watersheds,
                                         cones = cones,
                                         nests = nests)),
  tar_target(cc_gpkg, 
             save_sf(sf = final_cc, output_path = "GIS/radar_derived_catchments.gpkg"), 
             format = "file"),
  #### CATCHMENTS x HABITAT ####
  # Intersect with 2024 MAMU suitable habitat layer
  tar_target(suitable_habitat, prepare_mamu_habitat_gpkg(path = "GIS/2024_suitable_habitat/",
                                                    maz = maz,
                                                    regions = regions)),
  # Intersect with 2025 MAMU suitable habitat layer
  # tar_terra_rast(suitable_habitat, prepare_mamu_habitat_tiff(path = "GIS/2025_suitable_habitat/mamu_predict_2025_feb_03.tif",
  #                                                            maz = maz,
  #                                                            band = "X1")),
  # Intersect suitable habitat in each catchment
  tar_target(cc_habitat, st_habitat_in_sf(sf = final_cc,
                                          habitat = suitable_habitat,
                                          use_probability = TRUE)),
  # Intersect suitable habitat in each region
  tar_target(reg_habitat, st_habitat_in_sf(sf = regions,
                                           habitat = suitable_habitat,
                                           use_probability = TRUE)),
  # Final visualization
  # tar_render(final_visualization,
  #            path = "Rmd/catchments_visualization.Rmd",
  #            output_file = "catchments_visualization.pdf",
  #            #error = "null",
  #            quiet = TRUE,
  #            params = list(catchments = final_cc,
  #                          watersheds = watersheds_raw,
  #                          suitable_habitat = suitable_habitat,
  #                          cones = cones,
  #                          stn = stn,
  #                          nests = nests,
  #                          apikey = jawg_token,
  #                          headings = h_0,
  #                          h = h)
  #            ),
  #### ANNUAL MAX MAMU COUNTS ####
  # Select maximum MAMU count per station per year
  tar_target(annual_max_mamu, s |> 
               sf::st_drop_geometry() |>
               dplyr::group_by(site, region, year) |> 
               dplyr::summarise(max_mamu = max(mamuinpd, na.rm = TRUE))),
  # Calculate bootstrapped mean annual maximum per catchment
  # (across all years)
  tar_target(mean_max_mamu_cc, bootmean(annual_max_mamu, 
                                        group_by = "site",
                                        dat_col = "max_mamu", 
                                        CI_level = 0.95) |>
               dplyr::mutate(dplyr::across(c(bootmean, boot_min, boot_max), 
                                           round))),
  # Regional mean annual maximum per catchment
  # (across all years) for paper table purposes
  tar_target(mean_max_mamu_reg, bootmean(annual_max_mamu, 
                                        group_by = "region",
                                        dat_col = "max_mamu", 
                                        CI_level = 0.95) |>
               dplyr::mutate(dplyr::across(c(bootmean, boot_min, boot_max), 
                                           round))),
  #### DENSITY CALCS ####
  # NOTE suitable habitat =/= even MAMU density across whole layer!
  # While interior habitat might be equally 'suitable' to coastal habitat, the 
  # habitat closer to sea will have more MAMU population than more interior sites.
  # Total habitat (ha) across whole suitable habitat layer
  tar_target(total_suit_hab_area_ha, sum(reg_habitat$sh_area_ha)), 
  # Habitat (ha) summarized by region
  tar_target(regional_suit_hab_area_ha, sf::st_drop_geometry(reg_habitat)), 
  # Calculate the mean density of birds within each catchment
  # Using the bootstrapped mean + upper + lower CIs
  tar_target(cc_density, catchment_density(mm = mean_max_mamu_cc, 
                                           catchment_habitat = cc_habitat)),
  # Calculate the bootstrapped mean density of birds within each region
  tar_target(reg_density, regional_density(cc_density, 
                                           group_by = "region", 
                                           dat_col = "density", 
                                           CI_level = 0.95,
                                           add_AKB = TRUE)),
  # Fit a nest gamma decay function + rasterize it
  # Certain habitats may meet the criteria for "suitable habitat".
  # However, while habitat may be suitable for MAMU nesting in terms
  # of tree species composition, the suitable habitat layers do not
  # explicitly take into account the fact that ~99% of nests occur
  # within 30km of the coastline, and, crucially, that the further
  # from the coast you are, the less likely the nests are likely to
  # occur. The nest data follow a gamma distribution of likelihood
  # vs distance from shore. So, fit a gamma distribution to the nest
  # data, and then map that gamma distribution decay curve to a 
  # raster. Cells <30km from shore will have a higher probability,
  # closer to 1, while distances >30km will decay down to 0 probability.
  # `nest_likelihood` has cut out all non-habitat pieces
  tar_terra_rast(nest_likelihood, nest_gamma_decay(nests = nests,
                                                   coast = sea,
                                                   habitat = suitable_habitat)),
  # `nest_likelihood_full` is primarily for visualization purposes
  tar_terra_rast(nest_likelihood_full, nest_gamma_decay(nests = nests,
                                                        coast = sea)),
  # Rasterize the regional mean density
  # Apply the nest probability decay function to regional densities
  # Now, apply that gamma decay function to the density layer such 
  # that the mean density per region will remain the same as calculated
  # in `reg_density`, but the *spatial pattern* of the density follows
  # the nest gamma distribution.
  tar_terra_rast(density_map, extrapolate_density(regions, 
                                                  reg_density, 
                                                  nest_likelihood)),
  #### POPULATION CALCS ####
  # Calculate MAMU population!
  tar_target(regional_population, calculate_population(density_map, 
                                                       sf = regions, 
                                                       merge_df = reg_density) |>
               dplyr::arrange(region)),
  tar_target(total_population, calculate_population(density_map))
) 
