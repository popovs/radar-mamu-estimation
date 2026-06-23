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
source("temp/apikey.R")

# SAVE PLOTS? Yes or no
save_plots <- FALSE


# PIPELINE ----------------------------------------------------------------

list(
  #### PREPARE INPUTS ####
  ##### Track files #####
  tar_target(regions_path, "GIS/cons_reg.gpkg", format = "file"), # Track regions gpkg file
  tar_target(nests_path, "GIS/MAMU_Nests_BC_CDC.gpkg", format = "file"), # Track nests gpkg file
  tar_target(s_path, "data/ECCC_FLNR_MAMU-RadarData-20240307.csv", format = "file"), # Track surveys csv file
  tar_target(h_path, "data/headings.csv", format = "file"), # Track flight headings csv file
  tar_target(ws_path, "GIS/Watersheds/code_2_HG_merge.gpkg", format = "file"), # Track watersheds gpkg file
  tar_target(sh_path, "GIS/Suitable Habitat/2024_provincial_suitable_habitat.gpkg", format = "file"), # Track suitable habitat file
  ##### Read & prepare files #####
  ###### Regions ######
  tar_target(regions, prepare_regions(filepath = regions_path)),
  tar_target(regions_region, regions$region), # track regions names within regions gpkg
  ###### Nests ######
  # The 'nests' target will be prepared down the line, after
  # all relevant raster data has been extracted at the nests.
  # 'nests_0' target contains just the coords of nests + region they're in
  # Later, we will add extracted raster values to this dataset.
  tar_target(nests_0, prepare_nests(filepath = nests_path,
                                       regions = regions)),
  ###### Radar Surveys ######
  tar_target(s, prepare_surveys(filepath = s_path,
                                regions = regions,
                                N_years_min = 3)), # minimum sample size (N unique years) cutoff
  ###### Flight Headings ######
  tar_target(h_0, prepare_headings(h_path) |>
               # Now filter it to only stations/sites that meet our
               # filtering criteria to applied to our main `s` dataset
               dplyr::filter(site %in% s$site)),
  ###### Watersheds ######
  tar_target(watersheds_0, prepare_watersheds(ws_path, regions)),
  ###### Suitable Habitat ######
  tar_target(suitable_habitat, prepare_habitat(sh_path, regions)),
  
  
  
  #### PREPARE DEM ####
  ##### Land area masks #####
  # These masks will be used to clean up the DEM data and split apart
  # land areas from sea areas.
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
  tar_terra_rast(land, mask_land(dem = DEM,
                                 canvas = canvas,
                                 land = land_mask)),
  ##### Sea #####
  tar_terra_rast(sea, terra::ifel(land == 0, 0, NA)),
  
  #### DISTANCE TO SEA ####
  ##### Sea distance #####
  # Distance from sea, in km
  tar_terra_rast(sea_dist, terra::distance(sea) / 1000),
  
  #### EXTRACT NEST VALUES ####
  ##### Elevation, distance from sea, cost #####
  # Elevation
  tar_target(nest_elev_m, terra::extract(DEM, nests_0, ID = FALSE)[[1]]),
  # Distance from sea

  ##### Merge nest data #####
  # Merge into a single dataset
  tar_target(nests, cbind(nests_0, nest_elev_m, nest_dist_km)), # nest_cost)),
  
  ##### Nest probability #####
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
  
  # TODO: look into if habitat is being parsed the best way. 
  # When you simply mask the full raster, the results look a bit
  # different.
  
  # `nest_likelihood` has cut out all non-habitat pieces
  tar_terra_rast(nest_likelihood, nest_gamma_decay(nests = nests,
                                                   sea_dist = sea_dist,
                                                   habitat = suitable_habitat)),
  # `nest_likelihood_full` is primarily for visualization purposes
  tar_terra_rast(nest_likelihood_full, nest_gamma_decay(nests = nests,
                                                        sea_dist = sea_dist)),

  #### GENERATE RADAR CONES ####
  # Here, radar survey 'catchments' containing marbled murrelet
  # nesting habitat will be generated, using a few basic
  # assumptions about bird habitat + BC Digital Elevation Model
  # (DEM) data.
  ##### stn #####
  # Extract individual station coords
  # While radar stations themselves are typically on land,
  # the birds typically are flying over the water in the inlet.
  # Calculate the 'mean flight entry point' where birds enter
  # the radar screen relative to the radar unit itself, and
  # use that as our 'stn' coordinates. I.e. we care where 
  # the *birds* are relative to land entry, don't care about
  # where the radar station is set up.
  tar_target(stn, prepare_stn(s, h_0, regions)),
  tar_target(stn_gpkg,
             save_sf(sf = stn, output_path = "temp/QGIS temp/stn.gpkg"),
             format = "file"),
  
  ##### Polar mean flight headings #####
  # Calculate polar mean flight headings
  # Most birds will be flying in the same general direction
  # as they head inland, with some variation. Calculate the 
  # polar mean flight heading with 95% bootstrapped confidence
  # intervals about the mean.
  # TODO: clean up comments/fxns
  # OLD METHOD:
  # tar_target(h, calc_polar_mean(headings = h_0, n_reps = 1000, alpha = 0.05)),
  # NEW METHOD (using CircStats pkg):
  tar_target(h, h_0 |>
               dplyr::mutate(rad = heading * pi / 180) |>
               dplyr::group_by(site) |>
               dplyr::summarise(n = CircStats::circ.disp(rad)[[1]],
                                r = CircStats::circ.disp(rad)[[2]],
                                rbar = CircStats::circ.disp(rad)[[3]],
                                var = CircStats::circ.disp(rad)[[4]],
                                mean = CircStats::circ.mean(rad),
                                lower = mean - var,
                                upper = mean + var) |>
               # Filter to only include samples with minumum sample size of 30
               dplyr::filter(n >= 30) |>
               dplyr::mutate(mean = mean * 180 / pi,
                             lower = lower * 180 / pi,
                             upper = upper * 180 / pi,
                             theta = upper - lower) |>
               dplyr::mutate(mean = dplyr::if_else(mean < 0, mean + 360, mean),
                             lower = dplyr::if_else(lower < 0, lower + 360, lower),
                             upper = dplyr::if_else(upper < 0, upper + 360, upper),
                             theta = dplyr::if_else(theta > 180, 360 - theta, theta))),
  
  ##### Calculate cones #####
  tar_target(cones, generate_cones(h = h,
                                   stn = stn,
                                   radius = 30000, # length cones to select appropriate watersheds (in meters)
                                   res = res)),
  tar_target(cones_gpkg,
             save_sf(sf = cones, output_path = "temp/QGIS temp/flight_headings.gpkg"),
             format = "file"),
  #### CREATE CATCHMENTS ####
  ##### Select targeted watersheds #####
  # Now, based on the MAMU flight headings at each radar station,
  # select the watershed catchments the birds are targeting (i.e.,
  # the watersheds the radar cones overlap with).
  tar_target(watersheds, select_watersheds(watersheds = watersheds_0,
                                           cones = cones,
                                           min_cone_coverage = 0.01,
                                           output_plots = save_plots, # defined at the top of script
                                           output_dir = "temp/cone_inspection",
                                           stn = stn,
                                           headings = h_0,
                                           h = h)),
  
  ##### Calculate cost within watersheds #####
  tar_target(cost_catchments, lapply(unique(h$site),
                                     st_cost_catchment, # our function that defines cost catchments for one site
                                     watersheds = watersheds, 
                                     dem = terra::merge(DEM, sea), 
                                     cones = cones, 
                                     stn = stn, 
                                     nest_likelihood = nest_likelihood) |>
               dplyr::bind_rows()),
  
  # Save gpkg
  tar_target(cc_gpkg,
             save_sf(sf = cost_catchments, output_path = "temp/QGIS temp/radar_derived_catchments.gpkg"),
             format = "file"),

  #### CATCHMENTS X HABITAT ####
  # Intersect suitable habitat in each catchment
  # TODO: dplyr warning:
  # `summarise()` has regrouped the output.       [4m 39.6s, 3+, 63-]
  # ℹ Summaries were computed grouped by site, area_ha, region, loc, and mean_cost.
  # ℹ Output is grouped by site, area_ha, region, and loc.
  # ℹ Use `summarise(.groups = "drop_last")` to silence this message.
  # ℹ Use `summarise(.by = c(site, area_ha, region, loc, mean_cost))` for
  # per-operation grouping (`?dplyr::dplyr_by`) instead.
  tar_target(cc_habitat, st_habitat_in_sf(sf = cost_catchments,
                                          habitat = suitable_habitat,
                                          use_probability = TRUE)),
  # Intersect suitable habitat in each region
  # TODO: dplyr warning:
  # ℹ Summaries were computed grouped by MMCR_NAME and region.
  # ℹ Output is grouped by MMCR_NAME.
  # ℹ Use `summarise(.groups = "drop_last")` to silence this message.
  # ℹ Use `summarise(.by = c(MMCR_NAME, region))` for per-operation grouping
  # (`?dplyr::dplyr_by`) instead.
  tar_target(reg_habitat, st_habitat_in_sf(sf = regions,
                                           habitat = suitable_habitat,
                                           use_probability = TRUE)),
  # Final visualization
  # TODO: rm viz scripts
  # tar_render(final_visualization,
  #            path = "Rmd/catchments_visualization.Rmd",
  #            output_file = "catchments_visualization.pdf",
  #            #error = "null",
  #            quiet = TRUE,
  #            params = list(catchments = cost_catchments,
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
               dplyr::summarise(max_mamu = max(mamu_in_pd, na.rm = TRUE),
                                N = dplyr::n())),
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
  # Calculate the weighted mean density of birds within each region
  # Catchments that take up more area of the region have higher weight
  tar_target(reg_density, cc_density |>
               dplyr::filter(region == "NC") |>
               dplyr::mutate(region = "AKB") |> # Add in dummy rows for AKB, using NC density
               dplyr::bind_rows(cc_density) |>
               dplyr::group_by(region) |>
               dplyr::summarise(N = dplyr::n(),
                                bootmean = sum(bootmean),
                                boot_min = sum(boot_min, na.rm = TRUE),
                                boot_max = sum(boot_max, na.rm = TRUE),
                                area_ha = sum(area_ha),
                                sh_area_ha = sum(sh_area_ha)) |>
               dplyr::mutate(density = bootmean / sh_area_ha,
                             density_lwr = boot_min / sh_area_ha,
                             density_upr = boot_max / sh_area_ha)),

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
  tar_target(total_population, calculate_population(density_map)),
  tar_target(total_density, colSums(cc_density[,c("bootmean", "boot_min", "boot_max")], na.rm = TRUE) / sum(cc_density$sh_area_ha)),
  tar_target(bc_density, total_population / total_suit_hab_area_ha)
  
) 
