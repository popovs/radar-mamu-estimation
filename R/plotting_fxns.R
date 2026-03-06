#' Plotting fxns
#' 
#' Self explanatory, I hope.


plot_headings <- function(site, headings, h) {
  # Subset to needed data
  headings <- headings[headings$name == site, ]
  h <- h[h$name == site, ]
  # Incoming headings plot
  p_inc <- ggplot2::ggplot(data = headings[headings$flightpath_type == "Incoming",]) + 
    ggplot2::geom_histogram(ggplot2::aes(x = heading)) + 
    ggplot2::geom_vline(xintercept = h[["mean"]],
                        color = "red") +
    ggplot2::geom_vline(xintercept = h[["lower"]],
                        color = "grey",
                        linetype = "dashed") +
    ggplot2::geom_vline(xintercept = h[["upper"]],
                        color = "grey",
                        linetype = "dashed") +
    ggplot2::xlim(0, 360) + 
    ggplot2::ggtitle("Incoming") + 
    ggplot2::theme(axis.title = ggplot2::element_blank()) +
    ggplot2::theme_minimal()
  # Outgoing headings plot
  p_out <- ggplot2::ggplot(data = headings[headings$flightpath_type == "Outgoing",]) + 
    ggplot2::geom_histogram(ggplot2::aes(x = heading)) + 
    ggplot2::geom_vline(xintercept = h[["mean"]],
                        color = "red") +
    ggplot2::geom_vline(xintercept = h[["lower"]],
                        color = "grey",
                        linetype = "dashed") +
    ggplot2::geom_vline(xintercept = h[["upper"]],
                        color = "grey",
                        linetype = "dashed") +
    ggplot2::xlim(0, 360) + 
    ggplot2::ggtitle("180° - Outgoing") + 
    ggplot2::theme(axis.title = ggplot2::element_blank()) +
    ggplot2::theme_minimal()
  
  p <- ggpubr::ggarrange(p_inc, p_out, ncol = 1)
  return(p)
}


plot_watersheds <- function(site, watersheds, cones, stn) {
  # Subset to needed data
  watersheds <- watersheds[watersheds$site == site, ]
  cones <- cones[cones$site == site, ]
  stn <- stn[stn$site == site, ]
  # Plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = watersheds,
                     ggplot2::aes(color = keep_yn),
                     show.legend = FALSE) +
    ggplot2::geom_sf(data = cones,
                     fill = NA) +
    ggplot2::geom_sf(data = stn) +
    ggplot2::geom_text(data = watersheds,
                       ggplot2::aes(label = round(prct_coverage, 2),
                                    x = centroid_x,
                                    y = centroid_y),
                       size = 2) +
    ggplot2::ggtitle(site) +
    ggplot2::theme(axis.title = ggplot2::element_blank())
  return(p)
}


plot_cost <- function(site, cost_catchment, watersheds, 
                      cones, stn, nests, cost_cutoffs) {
  # Subset to needed data
  cost_catchment <- cost_catchment[[site]]
  watersheds <- watersheds[watersheds$site == site, ]
  cones <- cones[cones$site == site, ]
  stn <- stn[stn$site == site, ]
  nests <- suppressWarnings(sf::st_intersection(nests, watersheds))
  # Merge watersheds, stn, and cost_cutoffs to get all
  # attributes in one sf object
  watersheds <- merge(watersheds, sf::st_drop_geometry(stn[,c("site", "region")]))
  watersheds <- merge(watersheds, cost_cutoffs)
  # Plot
  p <- ggplot2::ggplot() + 
    tidyterra::geom_spatraster_contour_filled(data = cost_catchment,
                                              show.legend = FALSE) +
    tidyterra::geom_spatraster_contour(data = cost_catchment,
                                       breaks = c(#watersheds[["cost_min"]][watersheds$site == x],
                                         watersheds[["cost_max"]])) +
    tidyterra::scale_fill_whitebox_d() +
    ggplot2::geom_sf(data = cones,
                     fill = NA) +
    ggplot2::geom_sf(data = stn) +
    ggplot2::geom_sf(data = nests, 
                     color = "red", 
                     fill = "red") +
    ggplot2::ggtitle(site, subtitle = paste(watersheds[["region"]], "max cost =", watersheds[["cost_max"]])) +
    ggplot2::theme(axis.title = ggplot2::element_blank())
  return(p)
}


plot_catchment <- function(site, cost_catchment, accessible_catchment,
                           watersheds, cones, stn, nests) {
  # Subset to needed data
  cost_catchment <- cost_catchment[cost_catchment$site == site, ]
  accessible_catchment <- accessible_catchment[accessible_catchment$site == site, ]
  watersheds <- watersheds[watersheds$site == site, ]
  cones <- cones[cones$site == site, ]
  stn <- stn[stn$site == site, ]
  nests <- suppressWarnings(sf::st_intersection(nests, watersheds))
  # Plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = watersheds,
                     show.legend = FALSE) +
    ggplot2::geom_sf(data = cost_catchment,
                     fill = "grey",
                     color = NA,
                     alpha = 0.7,
                     show.legend = FALSE) +
    ggplot2::geom_sf(data = accessible_catchment,
                     color = NA, 
                     fill = "#26D1EA",
                     alpha = 0.3,
                     show.legend = FALSE) +
    ggplot2::geom_sf(data = cones,
                     fill = NA) +
    ggplot2::geom_sf(data = stn) +
    ggplot2::geom_sf(data = nests, 
                     color = "red", 
                     fill = "red") +
    ggplot2::ggtitle(site, subtitle = stn[["loc"]][stn$site == site])
  return(p)
}

plot_final_map <- function(site, 
                           catchments,
                           watersheds,
                           suitable_habitat,
                           cones, 
                           stn, 
                           nests,
                           apikey) {
  x <- site
  
  # Pull catchment bbox and add 5km plotting buffer
  bbox <- sf::st_bbox(catchments[catchments$site == x, ])
  bbox[1:2] <- bbox[1:2] - 5000 # add 5km buffer for visualizing
  bbox[3:4] <- bbox[3:4] + 5000
  
  # Prepare data for mapping
  # stn
  p_stn <- sf::st_transform(stn, 3005)
  p_stn <- cbind(p_stn, sf::st_coordinates(p_stn))
  p_stn <- sf::st_crop(p_stn, bbox)
  
  # watersheds
  watersheds <- sf::st_collection_extract(watersheds, "POLYGON")
  watersheds <- sf::st_crop(watersheds, bbox)
  
  # suitable habitat
  if (inherits(suitable_habitat, "sf")) {
    suitable_habitat <- sf::st_collection_extract(suitable_habitat, "POLYGON")
    suitable_habitat <- sf::st_crop(suitable_habitat, bbox)
  } else if (inherits(suitable_habitat, "SpatRaster")) {
    suitable_habitat <- terra::crop(suitable_habitat, bbox)
    suitable_habitat <- suitable_habitat |> 
      terra::as.polygons() |> 
      sf::st_as_sf()
  }
  
  
  # nests
  nests <- sf::st_centroid(nests)
  
  # Pull maptile for the site
  jawg_terrain <- maptiles::create_provider(
    name = "Jawg.Terrain",
    url = "https://tile.jawg.io/jawg-terrain/{z}/{x}/{y}.png?access-token={apikey}",
    citation = "© Jawg Maps"
  )
  
  tile <- maptiles::get_tiles(x = bbox,
                              provider = jawg_terrain,
                              apikey = apikey,
                              zoom = 12,
                              crop = TRUE)
  
  p <- ggplot2::ggplot() +
    tidyterra::geom_spatraster_rgb(data = tile) +
    ggplot2::geom_sf(data = watersheds,
                     color = "red",
                     fill = NA,
                     show.legend = FALSE) +
    ggplot2::geom_sf(data = suitable_habitat,
                     color = "#E92162",#"#B41347",
                     fill = "#FF296E",
                     linewidth = 0.05,
                     alpha = 0.3) + 
    ggplot2::geom_sf(data = catchments[catchments$site == x,],
                     #ggplot2::aes(fill = site),
                     #fill = "#C6E921",
                     #color = "#A1BE19",
                     lwd = 0.1,
                     fill = "#26D1EA",
                     color = "#0996AB",
                     alpha = 0.5,
                     show.legend = FALSE) +
    ggplot2::geom_sf(data = nests,
                     color = "#333333",
                     shape = 18) +
    ggrepel::geom_text_repel(data = p_stn,
                             ggplot2::aes(label = site,
                                          x = X,
                                          y = Y),
                             size = 3,
                             nudge_x = 0,
                             nudge_y = 250,
                             color = "black",
                             bg.color = "white",
                             bg.r = 0.15) +
    ggplot2::geom_sf(data = cones[cones$site == x, ],
                     fill = "orange",
                     color = "#D85426",
                     alpha = 0.15) +
    ggplot2::geom_sf(data = p_stn,
                     color = "#222222") +
    ggplot2::coord_sf(xlim = c(bbox[1], bbox[3]),
                      ylim = c(bbox[2], bbox[4]),
                      expand = FALSE) +
    ggplot2::ggtitle(x) +
    ggplot2::theme(axis.title = ggplot2::element_blank())
  
  return(p)
}