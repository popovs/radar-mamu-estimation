# A novel application of radar survey data to estimate population size of the Marbled Murrelet, a cryptic seabird in British Columbia, Canada

This repository contains all the code needed to estimate the population size of Marbled Murrelets (MAMU) in B.C. from radar monitoring data. This project comes bundled with a `renv` file in order to reproduce the package environment that was in place while running these analyses, so the code should work out-of-the-box. Some data comes pre-bundled within the data repository, while other data will automatically be downloaded onto your machine via the script.

## Quick start

This project relies on the [targets](https://books.ropensci.org/targets/) package to create a reproducible and manageable data pipeline. All objects within the pipeline are created step-by-step with various functions -- just like in a standard R script -- but are managed and tracked by the targets package to keep things organized. 

### BASIC USAGE

Assuming you have cloned this repository and the directory structure is identical on your machine...
1) You will first need to add two GIS files that are too large to track via Git to the following subdirectories:
    - `GIS/Suitable Habitat/2024_provincial_suitable_habitat.gpkg` -> download link will be updated here shortly.
    - TODO: ensure suitable habitat polygons are downloaded via pipeline. 
2) Run `renv::restore()` to install all necessary R packages the pipeline depends on. See the [renv documentation](https://rstudio.github.io/renv/reference/restore.html) for troubleshooting details.
3) Load the `{targets}` library. Open `_targets.R` and run `tar_make()` to run the pipeline. NOTE the full pipeline can take 1-2 hours to run if running from scratch. You can comment out sections of the pipeline that you are not interested in recreating if you wish to skip the creation of them, assuming nothing downstream of the pipeline depends upon it (e.g. Sup Mat generation at the end of the pipeline). TIP: you can run `tar_visnetwork()` to see how various targets depend upon each other.
4) Once you have created all your targets, you are ready to play with the data outputs. In a separate R script, you can run `tar_load(<target_name>)` to quickly load up the target in your R session and manipulate it from there.

 DISCLAIMER 1: THIS WHOLE SCRIPT WILL TAKE ~1 HOUR TO RUN.
 There are likely faster or more efficient ways of doing this,
 but the primary goal of this script is reproducibility and
 ease of understanding. If you wish to simply test the script, 
 modify the raster resolutions to a courser scale (e.g., 500m) 
 to run through the script quickly.

 DISCLAIMER 2: this was run with 16 GB ram and will likely 
 fail with less - that said, many of these steps can be
 easily replicated in QGIS to the same effect. The benefit of
 this code is the exact reproducibility.

## Project directories

### data

This folder contains the data used in this analysis, including the raw radar count survey data and flight headings data. It also includes supplemental data related to the project but not used explicitly within the script.

### GIS

This folder contains the GIS data necessary for the project, including regional boundaries, DEM data of the Alaska panhandle, suitable habitat polygons, and watershed polygons. **Note:** the suitable habitat polygons will need to be downloaded from the BC Data Catalogue (TODO: link will be added to this Github soon).

### R

This folder contains all R functions that the targets pipeline calls. 

### docs

This folder contains R Markdown scripts that generate the Supplemental Material documents to accompany the publication.

## System Requirements

The analysis pipeline requires up to **10 GB** of free space on your machine. This is due to the large file sizes of the DEM rasters and suitable habitat and watershed polygons. To fully clean this project off your machine, in addition to deleting your local copy of this repositry, be sure to *uninstall the bcdata package*, in order to remove the large BC DEM raster.
