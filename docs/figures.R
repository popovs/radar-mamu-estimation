
# This script creates all figures within the paper. Note that
# the targets pipeline must be run completely from start to finish
# by running `tar_make()` first before any of the targets can be 
# loaded in this script.

# setup -------------------------------------------------------------------


library(targets)
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)

tar_load(regions)
tar_load(s)
tar_load(stn)
tar_load(h)
tar_load(h_0)
tar_load(cones)
tar_load(watersheds)
tar_load(nests)
tar_load(nest_likelihood)
tar_load(DEM)
tar_load(sea)
tar_load(cost_catchments)
tar_load(mean_max_mamu_cc)

dem <- merge(DEM, sea)
headings <- h_0

labels <- read.csv("GIS/region_labels.csv")
labels <- labels |> st_as_sf(coords = c("X", "Y"), crs = 3005)
labels <- labels[labels$label %in% regions$region, ]

# Pick a site for site-specific figures
site <- "Aaltanhash"

# Filter all inputs down to the appropriate site
w <- watersheds[watersheds$site == site, ]
cone <- cones[cones$site == site, ]
s_stn <- stn[stn$site == site, ]
rad_stn <- s[s$site == site, "geometry"] |> unique()
headings <- headings[headings$site == site, ]
h <- h[h$site == site, ]
s_pt <- s[s$site == site, "geometry"] # pull IRL station coordinates

# Grab land polygons
ca <- rnaturalearth::ne_countries(scale = 10) |> # Load countries polygons
  dplyr::filter(name %in% c("Canada", "United States of America"))

ca <- st_transform(ca, 3005)

# Grab background imagery
bg_tiles <- maptiles::get_tiles(w, provider = "Esri.WorldTopoMap", crop = TRUE)




# bins map ----------------------------------------------------------------

DEM <- terra::mask(DEM, regions)

# Bounding box, add buffer
# (Same bounding box as hex bins map)
bbox <- st_bbox(cost_catchments)
bbox[1] <- bbox[1] - 70000
bbox[2] <- bbox[2] - 50000
bbox[3] <- bbox[3] + 110000
bbox[4] <- bbox[4] + 50000

## Create hillshade effect
slope <- terrain(DEM, "slope", unit = "radians")
aspect <- terrain(DEM, "aspect", unit = "radians")
hill <- shade(slope, aspect, 30, 270)
names(hill) <- "shades"

# Hillshading, but we need a palette
pal_greys <- hcl.colors(1000, "Grays")

# Create hexagons layer
hex <- merge(mean_max_mamu_cc, stn)
hex <- hex |>
  st_drop_geometry() |>
  dplyr::select(lat, lon, bootmean) |>
  sf::st_as_sf(coords = c("lon", "lat"), 
               crs = 4326) |>
  sf::st_transform(3005)
hex <- cbind(hex, sf::st_coordinates(hex))


map <- 
  ggplot() +
  geom_sf(data = ca,
          color = "#BCBCBC",
          size = 0.1,
          fill = "#d0d0d0") +
  geom_spatraster(data = DEM,
                  alpha = 0.5) +
  scale_fill_hypso_tint_c("moon_hypso",
                          na.value = NA,
                          guide = "none") +
  # scale_fill_distiller(na.value = NA,
  #                      guide = "none",
  #                      palette = "Greys") +
  ggnewscale::new_scale_fill() + 
  geom_spatraster(data = hill,
                  alpha = 0.3,
                  show.legend = FALSE) +
  scale_fill_gradientn(colors = pal_greys, 
                       na.value = NA) +
  ggnewscale::new_scale_fill() + 
  stat_summary_hex(data = hex,
                   aes(x = X,
                       y = Y,
                       z = bootmean),
                   alpha = 0.8) +
  rcartocolor::scale_fill_carto_c(palette = "BurgYl",
                                  name = "Raw MAMU count") +
  geom_sf(data = regions,
          color = "white",
          linewidth = 0.9,
          fill = NA) +
  geom_sf(data = regions,
          color = "#202020",
          fill = NA) +
  geom_sf_text (data = labels,
                aes(label = label,
                    fontface = "bold")) +
  coord_sf(xlim = bbox[c(1,3)],
           ylim = bbox[c(2,4)],
           expand = F) +
  theme_minimal() +
  guides(fill = guide_colourbar(barwidth = 15,
                                barheight = 0.5,
                                ticks = FALSE)) +
  theme(axis.title = element_blank(),
        legend.position = "bottom")


inset <- 
  ggplot() +
  geom_sf(data = ca,
          color = NA,
          fill = "grey65") +
  geom_sf(data = regions,
          color = "white",
          linewidth = 1.1,
          fill = NA) +
  geom_sf(data = regions,
          color = "#202020",
          fill = NA,
          linewidth = 0.2) +
  geom_sf(data = st_as_sfc(bbox),
          fill = NA, 
          colour = "black",
          linewidth = 0.6) +
  coord_sf(xlim = c(237827, 1407983),
           ylim = c(325414, 1770177),
           expand = FALSE) +
  theme_void() #+
#theme(panel.border = element_rect(fill = NA, color = "grey15"))

# It's gonna look all crooked in the plot preview
# What matters is how it looks after ggsave
cowplot::ggdraw() +
  cowplot::draw_plot(map) +
  cowplot::draw_plot(inset,
                     x = 0.7,
                     y = 0.27,
                     width = 0.25)

# Save map
ggsave("figs/fig-1.png",
       device = "png",
       dpi = 300,
       width = 6.5,
       height = 7,
       units = "in")



# nest density curve ------------------------------------------------------

ggplot(nests, aes(x = nest_dist_km)) +
  geom_density() + 
  labs(x = "Distance from the coast (km)",
       y = "Nest density") +
  theme_minimal()

ggsave("figs/fig-2.png",
       device = "png",
       dpi = 300,
       width = 6.5,
       height = 4,
       units = "in")

# radar histogram ---------------------------------------------------------

rh <- ggplot(data = headings) + 
  geom_histogram(aes(x = heading), breaks = seq(0, 360, by = 10)) + 
  geom_vline(xintercept = h[["mean"]],
                      color = "red") +
  geom_vline(xintercept = h[["lower"]],
                      color = "grey",
                      linetype = "dashed") +
  geom_vline(xintercept = h[["upper"]],
                      color = "grey",
                      linetype = "dashed") +
  coord_polar() +
  scale_x_continuous(
    limits = c(0, 360), 
    breaks = seq(0, 330, by = 30),
    expand = expansion(mult = c(0, 0)) # Keeps the circular ends properly connected
  ) +
  theme_minimal() +
  theme(axis.title = element_blank(),
                 axis.text.y = element_blank())





# least cost paths --------------------------------------------------------

# Dissolve watersheds together by site
w <- w |>
  dplyr::select(site) |>
  dplyr::group_by(site) |>
  dplyr::summarize()

nl <- crop(nest_likelihood, w, mask = TRUE) # crop nest likelihood raster to our watershed

# Prepare watershed DEM
tmp <- crop(dem, w, mask = TRUE)
tmp2 <- crop(dem, cone, mask = TRUE) # also extract anything in the cone path - e.g. water - so birds can cross over water areas in the cone's path
tmp <- merge(tmp, tmp2)

# Calculate resistance surface for the selected watersheds
surf <- movecost::mc_surface(dtm = tmp, funct = "e")

# Calculate the origin pixel of the surface raster
# (In many cases, the station coordinate does not actually
# exactly overlap the raster)
p <- as.points(tmp, values = FALSE) # extract raster centroids
p <- st_as_sf(p)
origin <- p[st_nearest_feature(s_stn, p), ] # pull closest raster centroid to our station

# Calculate accumulated cost surface from origin
acc <- movecost::mc_accum(surf, origin = origin)

# Now let's make 100 'pseudo nests' within the watershed, all
# falling within the most likely habitat areas, and calculate
# path lengths of traveling from the origin to each pseudo-nest
# Sample 10 'nests' from nlf
pseudo_nests <- spatSample(nl,
                                  size = 100, 
                                  method = "weights", # incorporate the raster probability into the sampling
                                  na.rm = TRUE,
                                  as.points = TRUE) |>
  st_as_sf()

# Calculate least-cost paths between origin and the pseudonests
paths <- movecost::mc_paths(surf, 
                            origin = origin,
                            destin = pseudo_nests)

# Extract out maximum cost of longest flight path up to 15/30 km
max_path <- ifelse(s_stn$region == "HG", 15000, 30000) # Set max flight path in km - 30 km for mainland/Vancouver Island, 15 km for Haida Gwaii
max_cost <- paths$paths[paths$paths$length <= units::as_units(max_path, "m"), ] |> 
  dplyr::pull(cost) |> 
  max()

paths_geom <- paths$paths
paths_geom$pass_fail <- paths_geom$cost < max_cost

# Draw cost catchment boundary
cat <- movecost::mc_boundary(surface = surf, 
                             origin = origin,
                             limit = max_cost)

out <- cat$boundaries

# main panel plot ---------------------------------------------------------

w_bbox <- st_bbox(w) |>
  st_as_sfc()

# Overwrites previous inset
inset <- ggplot() +
  geom_sf(data = ca#,
          # color = NA,
          # fill = "grey65"
  ) +
  geom_sf(data = regions,
          fill = NA,
          color = "grey90") +
  geom_sf(data = w_bbox,
          color = "red",
          fill = NA,
          size = 1) +
  coord_sf(xlim = c(504003, 1300323),
           ylim = c(365000, 1054113)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

# theme_void() +
# theme(panel.background = element_rect(fill = "grey98",
#                                                         color = NA),
#                panel.border = element_rect(fill = NA,
#                                                     color = "grey15"))

# overwrites previous `p`
p <- 
  ggplot() +
  # Background
  geom_spatraster_rgb(data = bg_tiles) +
  
  # Radar station (IRL point)
  geom_sf(data = s_pt,
                   pch = 10, 
                   aes(color = "Radar station")) +
  scale_color_manual(values = "black",
                              guide = guide_legend(order = 1)) +
  ggnewscale::new_scale_color() +
  
  # Nest likelihood raster
  geom_spatraster(data = nl,
                           aes(fill = layer),
                           show.legend = FALSE) +
  scale_fill_princess_c("aura", 
                        name = "Nest likelihood", 
                        limits = c(0,1)) + 
  ggnewscale::new_scale_fill() +
  
  # Watershed boundaries
  geom_sf(data = watersheds[watersheds$site == site,],
                   aes(color = "Selected watersheds"),
          fill = NA) +
  scale_color_manual(values = "red",
                              guide = guide_legend(order = 4)) +
  ggnewscale::new_scale_color() +
  
  # Flight paths
  #ggnewscale::new_scale_color() +
  geom_sf(data = paths_geom[paths_geom$pass_fail == FALSE, ],
                   #aes(alpha = pass_fail),
                   #aes(color = pass_fail),
                   size = 0.4,
                   alpha = 0.4,
                   linetype = "11",
                   color = "darkcyan") +
  #scale_alpha_discrete(c(0.5, 0.2)) +
  #scale_color_manual(values = c("grey85", "grey10")) +
  geom_sf(data = paths_geom[paths_geom$pass_fail == TRUE, ],
                   aes(color = "Flight paths"),
                   size = 0.4) +
  scale_color_manual(values = "cyan") +
  # Pseudo nests
  geom_sf(data = pseudo_nests, 
                   aes(fill = "Pseudo nest"),
                   pch = 21,  
                   color = "black", 
                   size = 2) +
  scale_fill_manual(values = "cyan", name = "Pseudo nests") +
  # Cost catchment boundaries
  ggnewscale::new_scale_color() +
  geom_sf(data = out,
                   aes(color = "Cost catchment boundary"),
          fill = NA, 
          size = 1) +
  scale_color_manual(values = "black") +
  
  # Radar station (entry point)
  geom_sf(data = s_stn,
                   aes(shape = "Radar entry point"),
                   color = "orange",
                   size = 3) +
  scale_shape_manual(values = 16,
                              guide = guide_legend(order = 2)) +
  
  # Radar cone
  ggnewscale::new_scale_fill() +
  geom_sf(data = cone,
                   aes(fill = "Radar entry cone"),
          color = "orange",
          alpha = 0.3) +
  scale_fill_manual(values = "orange",
                             guide = guide_legend(order = 3)) +
  
  # Plot label
  # annotate("text", 
  #                   label = "A",
  #                   x = 828000,
  #                   y = 929000,
  #                   size = 8) +
  
  # Misc aesthetics
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  #ggtitle(site) +
  theme_minimal() +
  theme(legend.title = element_blank(),
                 #legend.position = "bottom",
        legend.spacing.y = unit(-10, "pt"),
        legend.box.margin = margin(t = 160)
        ) +
  guides(colour = guide_legend(ncol = 2,
                                        nrow = 2,
                                        byrow = TRUE))


# combine -----------------------------------------------------------------

final <- 
  p + patchwork::inset_element(
  rh, # created in 'watershed selection map' code chunk
  top = 1,
  right = 1.5,
  bottom = 0.5,
  left = 0.99
) + patchwork::inset_element(
  inset,
  top = 1,
  right = 0.45,
  bottom = 0.75,
  left = 0.25
)
  
# x <- ggpubr::ggarrange(inset, rh, ncol = 1, labels = c("B", "C"))
# 
# final <- ggpubr::ggarrange(p, x, nrow = 1, 
#                   widths = c(2, 1),
#                   labels = c("A", NA), 
#                   #common.legend = TRUE, 
#                   #legend = "bottom"
#                   )

final

ggsave("figs/fig-3.png",
                final,
                bg = "white",
                #width = 10,
                #height = 6,
                dpi = 300)




# full catchments map -----------------------------------------------------


map <- ggplot() +
  geom_sf(data = ca,
          color = "#BCBCBC",
          size = 0.1,
          fill = "#d0d0d0") +
  geom_spatraster(data = DEM,
                  alpha = 0.5) +
  scale_fill_hypso_tint_c("moon_hypso",
                          na.value = NA,
                          guide = "none") +
  # scale_fill_distiller(na.value = NA,
  #                      guide = "none",
  #                      palette = "Greys") +
  ggnewscale::new_scale_fill() + 
  geom_spatraster(data = hill,
                  alpha = 0.3) +
  scale_fill_gradientn(colors = pal_greys, na.value = NA) +
  ggnewscale::new_scale_fill() + 
  geom_sf(data = cost_catchments, 
          aes(fill = site),
          color = "grey95",
          linewidth = 0.05,
          alpha = 0.9) +
  scale_fill_manual(values = rep(rcartocolor::carto_pal(n = 12, "Vivid")[1:11], 6)) +
  geom_sf(data = regions,
          color = "white",
          linewidth = 0.9,
          fill = NA) +
  geom_sf(data = regions,
          color = "#202020",
          fill = NA) +
  geom_sf_text (data = labels,
                aes(label = label,
                    fontface = "bold")) +
  coord_sf(xlim = bbox[c(1,3)],
           ylim = bbox[c(2,4)],
           expand = F) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        legend.position = "none")



inset <- 
  ggplot() +
  geom_sf(data = ca,
          color = NA,
          fill = "grey65") +
  geom_sf(data = regions,
          color = "white",
          linewidth = 1.1,
          fill = NA) +
  geom_sf(data = regions,
          color = "#202020",
          fill = NA,
          linewidth = 0.2) +
  geom_sf(data = st_as_sfc(bbox),
          fill = NA, 
          colour = "black",
          linewidth = 0.6) +
  coord_sf(xlim = c(237827, 1407983),
           ylim = c(325414, 1770177),
           expand = FALSE) +
  theme_void() #+
#theme(panel.border = element_rect(fill = NA, color = "grey15"))

# It's gonna look all crooked in the plot preview
# What matters is how it looks after ggsave
cowplot::ggdraw() +
  cowplot::draw_plot(map) +
  cowplot::draw_plot(inset,
            x = 0.7,
            y = 0.27,
            width = 0.25)

# Save map
ggsave("figs/fig-4.png",
       device = "png",
       dpi = 300,
       width = 6.5,
       height = 7,
       units = "in")

  