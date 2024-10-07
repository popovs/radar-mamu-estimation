# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

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

# Set resolution (in meterse)
res <- 250 # res MUST be >25, as the DEM goes does to 25m accuracy.

# Set target options:
tar_option_set(
  packages = c("MAMU", # remotes::install_github("popovs/MAMU")
               ),
  format = "qs", 
  memory = "transient", # unload memory for each target line after it's completed
  garbage_collection = TRUE, # run gc() prior to each target
)

tar_source()

list(
  
)
