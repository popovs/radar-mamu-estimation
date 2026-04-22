#' 0. MAMU NESTING HABITAT AREA

#' This script contains functions that will download the BC 
#' Digital Elevation Model (DEM) for our study region, then 
#' clip it down to the areas that are highly likely to contain 
#' suitable marbled murrelet nesting habitat area (regions 
#' that are 30 km from shore and under <1500 m elevation). 
#' The exact maximum elevation cutoff varies by our six 
#' conservation regions of interest (North,  Central, and South 
#' Mainland Coast, Haida Gwaii, North + West Vancouver Island, 
#' and East Vancouver Island).




# REGIONAL ELEVATION CUTOFFS ----------------------------------------------




nest_quantiles <- function(nests, quant_data, prefix, quantiles = c(0.025, 0.975)) {
  # TODO: if this overwrites the data, it might trigger an endless pipeline reassessment loop
  # TODO: regions within nest data might not necessarily line up with regions_map in pipeline
  nests$quant_col <- quant_data
  # Calculate quantiles
  quants <- aggregate(quant_col ~ region, nests, FUN = quantile, quantiles[1], na.rm = TRUE)
  quants[3] <- aggregate(quant_col ~ region, nests, FUN = quantile, quantiles[2], na.rm = TRUE)[2]
  # Set up output names
  min_col <- paste0(prefix, "_min")
  max_col <- paste0(prefix, "_max")
  names(quants) <- c("region", min_col, max_col)
  # Round to nearest 10 
  quants[,2:3] <- round(quants[,2:3], -1)
  # Fill in missing values, if they're missing
  # If NVI is NULL, use mean of the other 3 regions of VI
  # If AKB and NC are NULL, use CC cutoffs
  if (!"NVI" %in% quants$region) quants <- rbind(quants, c("NVI", round(mean(quants[grep("VI", quants$region), min_col])), round(mean(quants[grep("VI", quants$region), max_col]))))
  if (!"AKB" %in% quants$region) quants <- rbind(quants, c("AKB", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  if (!"NC" %in% quants$region)  quants <- rbind(quants, c("NC", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  # rbind converts numerics to character... convert them back
  quants[[min_col]] <- as.numeric(quants[[min_col]])
  quants[[max_col]] <- as.numeric(quants[[max_col]])
  return(quants)
}


# Helper functions to extract whisker boxplot values
# Anything above that is plotted as outlier dots
upper_whisker <- function(x) {
  x <- x[!is.na(x)]
  min(max(x), as.numeric(quantile(x, 0.75)) + (IQR(x) * 1.5))
}

lower_whisker <- function(x) {
  x <- x[!is.na(x)]
  max(min(x), as.numeric(quantile(x, 0.25)) - (IQR(x) * 1.5))
}

# This is an alternative function to nest_quantiles - instead of taking
# the 95% (or whatever) quantiles of X nest data, this function instead
# finds the min and max of each data group and removes any outliers
# from the interquartile range (i.e., boxplot outliers).
nest_minmax_sans_outliers <- function(nests, quant_data, prefix) {
  # TODO: if this overwrites the data, it might trigger an endless pipeline reassessment loop
  # TODO: regions within nest data might not necessarily line up with regions_map in pipeline
  nests$quant_col <- quant_data
  
  # Remove overall nest quant outliers first
  upper_outliers <- upper_whisker(quant_data)
  nests$quant_outlier_yn <- nests$quant_col >= upper_outliers
  # IQR method
  #upper_outliers <- boxplot.stats(nests$quant_col)$out[boxplot.stats(nests$quant_col)$out > median(nests$quant_col, na.rm = TRUE)]
  # Quantile method
  #upper_outliers <- quantile(nests[["quant_col"]], 0.95, na.rm = TRUE)
  #nests$quant_outlier_yn <- nests$quant_col >= min(upper_outliers)
  
  # Pull minimum value by group after overall quant outliers are cut out
  quants <- aggregate(quant_col ~ region, nests, FUN = min, na.rm = TRUE)
  # Pull maximum value by group after overall quant outliers are cut out
  quants$max <- aggregate(quant_col ~ region, 
                          nests[nests$quant_outlier_yn == FALSE,], 
                          FUN = max, 
                          na.rm = TRUE)[[2]]

  # Set up output names
  min_col <- paste0(prefix, "_min")
  max_col <- paste0(prefix, "_max")
  names(quants) <- c("region", min_col, max_col)
  
  # Round to nearest 10
  quants$cost_min <- floor(quants$cost_min / 10) * 10 # round DOWN to nearest 10
  quants$cost_max <- ceiling(quants$cost_max / 10) * 10 # round UP to nearest 10
  
  # Fill in missing values, if they're missing
  # If NVI is NULL, use mean of the other 3 regions of VI
  # If AKB and NC are NULL, use CC cutoffs
  if (!"NVI" %in% quants$region) quants <- rbind(quants, c("NVI", round(mean(quants[grep("VI", quants$region), min_col])), round(mean(quants[grep("VI", quants$region), max_col]))))
  if (!"AKB" %in% quants$region) quants <- rbind(quants, c("AKB", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  if (!"NC" %in% quants$region)  quants <- rbind(quants, c("NC", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  # rbind converts numerics to character... convert them back
  quants[[min_col]] <- as.numeric(quants[[min_col]])
  quants[[max_col]] <- as.numeric(quants[[max_col]])
  return(quants)
  
}

# Finally, this is what seems to be the best approach - the 8 outliers
# across the entire dataset - nests with costs > 7k - are cut out; 
# then rather than taking the maximum per region, the IQR is taken
# per region (ie, cutting out regional outliers). This approach is 
# defendable because our regional groupings are fairly arbitrary
# and don't necessarily reflect actual regional behavioral patterns
# in MAMU activity. So, there may be a few nests in one region (e.g.,
# MWVI) that could reasonably just as easily fall within another region
# (e.g., NVI). 
nest_iqr_sans_outliers <- function(nests, quant_data, prefix, hg_exception) {
  # TODO: if this overwrites the data, it might trigger an endless pipeline reassessment loop
  # TODO: regions within nest data might not necessarily line up with regions_map in pipeline
  nests$quant_col <- quant_data
  
  # Remove overall nest quant outliers first
  upper_outliers <- upper_whisker(quant_data)
  nests$quant_outlier_yn <- nests$quant_col >= upper_outliers
  
  # Pull minimum value by group after overall quant outliers are cut out
  quants <- aggregate(quant_col ~ region, nests, FUN = min, na.rm = TRUE)
  # Now check if any individual regions have IQR outliers, after overall quant outliers are cut out
  quants$upr_whisker <- aggregate(quant_col ~ region, 
                                  nests[nests$quant_outlier_yn == FALSE,], 
                                  FUN = upper_whisker)[[2]]
  
  # Set up output names
  min_col <- paste0(prefix, "_min")
  max_col <- paste0(prefix, "_max")
  names(quants) <- c("region", min_col, max_col)
  
  # Round to nearest 10
  quants$cost_min <- floor(quants$cost_min / 10) * 10 # round DOWN to nearest 10
  quants$cost_max <- ceiling(quants$cost_max / 10) * 10 # round UP to nearest 10
  
  # HG has so little data and birds tend to fly differently there,
  # as there aren't distinctly well defined watersheds in the same way 
  # as the mainland or VI, which have more rugose coastlines. 
  # So, as an option, simply take the maximum value for HG rather than
  # cutting out any outliers. 
  if (hg_exception) {
    hg_max <- max(nests[["quant_col"]][nests$quant_outlier_yn == FALSE & nests$region == "HG"])
    hg_max <- ceiling(hg_max / 10) * 10
    quants[["cost_max"]][quants$region == "HG"] <- hg_max
  }
  
  # Fill in missing values, if they're missing
  # If NVI is NULL, use mean of the other 3 regions of VI
  # If AKB and NC are NULL, use CC cutoffs
  if (!"NVI" %in% quants$region) quants <- rbind(quants, c("NVI", round(mean(quants[grep("VI", quants$region), min_col])), round(mean(quants[grep("VI", quants$region), max_col]))))
  if (!"AKB" %in% quants$region) quants <- rbind(quants, c("AKB", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  if (!"NC" %in% quants$region)  quants <- rbind(quants, c("NC", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  
  # rbind converts numerics to character... convert them back
  quants[[min_col]] <- as.numeric(quants[[min_col]])
  quants[[max_col]] <- as.numeric(quants[[max_col]])
  
  return(quants)
}


nest_isoforest <- function(nests, quant_data, prefix) {
  nests$quant_col <- quant_data
  # Create isolation forest model and calcuate outlier scores
  data <- sf::st_drop_geometry(nests[,c("region", "quant_col")])
  if_model <- isotree::isolation.forest(data, sample_size = 3, ndim=1, ntrees=10, nthreads=1)
  scores <- isotree::predict.isolation_forest(if_model, data, type="avg_depth")
  # Choose cutoff for what counts as an 'outlier score'. 
  # Choose scores that 99% of the data fall into
  outlier_threshold <- quantile(scores, 0.01)[[1]]
  # Add that info back to `nests`
  nests$scores <- scores
  nests$outlier_yn <- nests$scores <= outlier_threshold # less than or = to threshold
  # Choose max by group after cutting out regional outliers
  nests <- nests[nests$outlier_yn == FALSE, ]
  quants <- aggregate(quant_col ~ region, nests, FUN = "max")
  # Set minimum to be zero with this method
  quants$min <- 0
  
  # Reorder!
  # NOTE this order is switched from other methods above!
  quants <- quants[,c(1,3,2)]
  
  # Set up output names
  min_col <- paste0(prefix, "_min")
  max_col <- paste0(prefix, "_max")
  names(quants) <- c("region", min_col, max_col) 
  
  # Round to nearest 10
  quants[[min_col]] <- floor(quants[[min_col]] / 10) * 10 # round DOWN to nearest 10
  quants[[max_col]] <- ceiling(quants[[max_col]] / 10) * 10 # round UP to nearest 10
  
  # Fill in missing values, if they're missing
  # If NVI is NULL, use mean of the other 3 regions of VI
  # If AKB and NC are NULL, use CC cutoffs
  if (!"NVI" %in% quants$region) quants <- rbind(quants, c("NVI", round(mean(quants[grep("VI", quants$region), min_col])), round(mean(quants[grep("VI", quants$region), max_col]))))
  if (!"AKB" %in% quants$region) quants <- rbind(quants, c("AKB", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  if (!"NC" %in% quants$region)  quants <- rbind(quants, c("NC", quants[[min_col]][quants$region == "CC"], quants[[max_col]][quants$region == "CC"]))
  
  # rbind converts numerics to character... convert them back
  quants[[min_col]] <- as.numeric(quants[[min_col]])
  quants[[max_col]] <- as.numeric(quants[[max_col]])
  
  return(quants)
}






