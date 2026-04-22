#' Utility fxns
#' 
#' Miscellaneous helper & target management functions.

# SAVE OUTPUTS ------------------------------------------------------------

# ... and still have it be tracked by `targets`

save_sf <- function(sf, output_path) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  sf::st_write(sf, output_path, append = FALSE)
  return(output_path)
}