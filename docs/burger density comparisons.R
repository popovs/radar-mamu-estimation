#' Comparing our derived density estimates for each catchment 
#' vs Burger (2002) "Conservation assessment of marbled murrelets 
#' in British Columbia: a review of the biology, populations, 
#' habitat associations, and conservation" estimates; Table A3-1

# Load up our densities
library(targets)
tar_load(cc_density)

# Calculate the suitable habitat within our catchments, exluding
# any overlapping areas (Burger catchments do NOT overlap, so this
# is already accounted for).
library(sf)
tar_load(cost_catchments)
tar_load(suitable_habitat)

cc_density <- merge(cost_catchments, cc_density)

# Drop pre-existing regional habitat area - we will recalculate it
cc_density <- dplyr::select(cc_density, -reg_area_ha)

# Load up Burger 2002 densities
b2002_density <- read.csv("data/Burger - 2002 - Table A3-1.csv", skip = 1)

# Merge the two together
d <- merge(cc_density, b2002_density, 
           by.x = "site", by.y = "station_popov")

# Label the points Burger considered to be outliers
d$outlier_yn <- grepl("\\*\\*", d$station_burger)

# Extract suitable habitat per Burger region groupings
burger_regions <- d |>
  dplyr::select(study_area) |>
  dplyr::group_by(study_area) |>
  dplyr::summarise()

plot(burger_regions, border = NA) # Popov study catchments, merged & grouped by Burger (2002) regions

# Intersect Burger regions w suitable habitat, then
# recalculate the total habitat area within each B. region.
ixn <- st_intersection(burger_regions, suitable_habitat)
ixn <- ixn |> 
  dplyr::group_by(study_area) |>
  dplyr::summarise()

plot(ixn, border = NA) # Habitat area within the merged/grouped catchments

ixn$reg_area_ha <- st_area(ixn) |> 
  units::set_units("ha")
ixn <- st_drop_geometry(ixn)

# Alright, now merge the Burger regional area into our
# dataset
d <- merge(d, ixn, by = "study_area")
d <- units::drop_units(d)
d <- st_drop_geometry(d)

# Out of curiosity, compare the two
d |>
  dplyr::group_by(study_area) |>
  dplyr::summarise(sum_ha_w_overlap = sum(sh_area_ha), # if we simply summed up area within catchments together
                   actual_ha = mean(reg_area_ha), # vs merging in GIS to remove overlaps
                   diff = sum_ha_w_overlap - actual_ha)

## 01 PLOTS ----

library(ggplot2)

## N MAMU PLOT ----
ggplot(d, aes(x = bootmean, y = mamu_count, color = outlier_yn)) +
  geom_point() +
  geom_linerange(aes(xmin = boot_min, xmax = boot_max), alpha = 0.4) +
  scale_color_manual(values = c("black", "red")) +
  geom_abline(slope = 1, intercept = 0) +
  ggrepel::geom_text_repel(data = d[d$boot_max > 1500, ],
                           aes(label = site),
                           nudge_x = 500,
                           nudge_y = 100) +
  coord_equal() +
  labs(x = "Popov MAMU count",
       y = "Burger MAMU count",
       caption = "In RED are points Burger (2002) considered outliers
       Error bars show 95% bootstrapped CI around mean station count") +
  theme(legend.position = "none")

## DENSITY PLOTS ----
ggplot(d, aes(x = density, y = density_a, color = outlier_yn)) +
  geom_point() +
  geom_linerange(aes(xmin = density_lwr, xmax = density_upr), alpha = 0.4) +
  scale_color_manual(values = c("black", "red")) +
  # ggrepel::geom_text_repel(aes(label = site),
  #                          min.segment.length = 0,
  #                          max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal() +
  labs(x = "Popov density",
       y = "Burger density (All mature + old growth)",
       caption = "In RED are points Burger (2002) considered outliers
       Error bars show 95% bootstrapped CI density bounds") +
  theme(legend.position = "none")

ggplot(d, aes(x = density, y = density_b, color = outlier_yn)) +
  geom_point() +
  geom_linerange(aes(xmin = density_lwr, xmax = density_upr), alpha = 0.4) +
  scale_color_manual(values = c("black", "red")) +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal() +
  labs(x = "Popov density",
       y = "Burger density (Most likely habitat)",
       caption = "In RED are points Burger (2002) considered outliers
       Error bars show 95% bootstrapped CI density bounds") +
  theme(legend.position = "none")

ggplot(d, aes(x = density, y = density_c, color = outlier_yn)) +
  geom_point() +
  geom_linerange(aes(xmin = density_lwr, xmax = density_upr), alpha = 0.4) +
  scale_color_manual(values = c("black", "red")) +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal() +
  labs(x = "Popov density",
       y = "Burger density (Highest realistic)",
       caption = "In RED are points Burger (2002) considered outliers
       Error bars show 95% bootstrapped CI density bounds") +
  theme(legend.position = "none")

## CATCHMENT AREA PLOT ----
ggplot(d, aes(x = area_ha, y = watershed_area_ha, color = outlier_yn)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal() +
  labs(x = "Popov catchment (ha)",
       y = "Burger catchment (ha)",
       caption = "In RED are points Burger (2002) considered outliers") +
  theme(legend.position = "none")

## HABITAT AREA PLOTS ----
ggplot(d, aes(x = sh_area_ha, y = habitat_ha_a, color = outlier_yn)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal() +
  labs(x = "Popov habitat (ha)",
       y = "Burger habitat (ha) (All mature + old growth)",
       caption = "In RED are points Burger (2002) considered outliers") +
  theme(legend.position = "none")

ggplot(d, aes(x = sh_area_ha, y = habitat_ha_b, color = outlier_yn)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal() +
  labs(x = "Popov habitat (ha)",
       y = "Burger habitat (ha) (Most likely habitat)",
       caption = "In RED are points Burger (2002) considered outliers") +
  theme(legend.position = "none")

ggplot(d, aes(x = sh_area_ha, y = habitat_ha_c, color = outlier_yn)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal() +
  labs(x = "Popov habitat (ha)",
       y = "Burger habitat (ha) (Highest realistic)",
       caption = "In RED are points Burger (2002) considered outliers") +
  theme(legend.position = "none")

## N SURVEYS PLOT ----
ggplot(d, aes(x = N, y = n_years_surveyed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal() +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  scale_y_continuous(breaks = seq(0, 10, 1)) +
  labs(x = "Popov N years per catchment",
       y = "Burger N years per catchment") +
  theme(legend.position = "none")


## 02 SUMMARY STATS ----

##### SUM METHOD 1 ----

# OUTLIERS INCLUDED

d |>
  dplyr::group_by(study_area) |>
  dplyr::summarise(popov_d = round(sum(bootmean) / mean(reg_area_ha), 3), # note using sh_area_ha here, to not overlap our catchment areas
                   popov_d_lwr = round(sum(boot_min, na.rm = TRUE) / mean(reg_area_ha), 3),
                   popov_d_upr = round(sum(boot_max, na.rm = TRUE) / mean(reg_area_ha), 3),
                   burger_a_d = round(sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_a, na.rm = TRUE), 3),
                   burger_b_d = round(sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_b, na.rm = TRUE), 3),
                   burger_c_d = round(sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_c, na.rm = TRUE), 3),
                   popov_n = sum(N),
                   burger_n = sum(n_years_surveyed)) |>
  dplyr::mutate(popov_d = paste0(popov_d, " [", popov_d_lwr, "-", popov_d_upr, "]")) |>
  dplyr::mutate(burger_d = paste0(burger_b_d, " [", burger_a_d, "-", burger_c_d, "]")) |>
  dplyr::select(study_area, popov_d, burger_d, popov_n, burger_n) |>
  knitr::kable()

d |>
  dplyr::group_by(study_area) |>
  dplyr::summarise(popov_d = sum(bootmean) / mean(reg_area_ha),
                   popov_d_lwr = sum(boot_min, na.rm = TRUE) / mean(reg_area_ha),
                   popov_d_upr = sum(boot_max, na.rm = TRUE) / mean(reg_area_ha),
                   burger_a_d = sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_a, na.rm = TRUE),
                   burger_b_d = sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_b, na.rm = TRUE),
                   burger_c_d = sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_c, na.rm = TRUE),) |>
  ggplot(aes(x = popov_d, y = burger_b_d)) +
  geom_point() +
  geom_linerange(aes(xmin = popov_d_lwr, xmax = popov_d_upr),
                 alpha = 0.3) +
  geom_linerange(aes(ymin = burger_a_d, ymax = burger_c_d),
                 alpha = 0.3) +
  geom_abline(slope = 1) +
  ggrepel::geom_text_repel(aes(label = study_area),
                           min.segment.length = 0) +
  scale_x_continuous(breaks = seq(0, 0.5, 0.01)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.01)) +
  coord_fixed() +
  labs(x = "Popov density estimates",
       y = "Burger density estimates",
       caption = "Error bars in X dimension show 95% CI Popov density est.
       Error bars in Y dimension show Burger min/max density est based on definition of habitat used.
       Outliers included.")
  
# OUTLIERS EXCLUDED

d |>
  dplyr::filter(outlier_yn == FALSE) |>
  dplyr::group_by(study_area) |>
  dplyr::summarise(popov_d = round(sum(bootmean) / mean(reg_area_ha), 3), # note using sh_area_ha here, to not overlap our catchment areas
                   popov_d_lwr = round(sum(boot_min, na.rm = TRUE) / mean(reg_area_ha), 3),
                   popov_d_upr = round(sum(boot_max, na.rm = TRUE) / mean(reg_area_ha), 3),
                   burger_a_d = round(sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_a, na.rm = TRUE), 3),
                   burger_b_d = round(sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_b, na.rm = TRUE), 3),
                   burger_c_d = round(sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_c, na.rm = TRUE), 3),
                   popov_n = sum(N),
                   burger_n = sum(n_years_surveyed)) |>
  dplyr::mutate(popov_d = paste0(popov_d, " [", popov_d_lwr, "-", popov_d_upr, "]")) |>
  dplyr::mutate(burger_d = paste0(burger_b_d, " [", burger_a_d, "-", burger_c_d, "]")) |>
  dplyr::select(study_area, popov_d, burger_d, popov_n, burger_n) |>
  knitr::kable()

d |>
  dplyr::filter(outlier_yn == FALSE) |>
  dplyr::group_by(study_area) |>
  dplyr::summarise(popov_d = sum(bootmean) / mean(reg_area_ha),
                   popov_d_lwr = sum(boot_min, na.rm = TRUE) / mean(reg_area_ha),
                   popov_d_upr = sum(boot_max, na.rm = TRUE) / mean(reg_area_ha),
                   burger_a_d = sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_a, na.rm = TRUE),
                   burger_b_d = sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_b, na.rm = TRUE),
                   burger_c_d = sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_c, na.rm = TRUE),) |>
  ggplot(aes(x = popov_d, y = burger_b_d)) +
  geom_point() +
  geom_linerange(aes(xmin = popov_d_lwr, xmax = popov_d_upr),
                 alpha = 0.3) +
  geom_linerange(aes(ymin = burger_a_d, ymax = burger_c_d),
                 alpha = 0.3) +
  geom_abline(slope = 1) +
  ggrepel::geom_text_repel(aes(label = study_area),
                           min.segment.length = 0) +
  scale_x_continuous(breaks = seq(0, 0.5, 0.01)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.01)) +
  coord_fixed() +
  labs(x = "Popov density estimates",
       y = "Burger density estimates",
       caption = "Error bars in X dimension show 95% CI Popov density est.
       Error bars in Y dimension show Burger min/max density est based on definition of habitat used.
       Outliers excluded.")



#### MEAN METHOD 2 ----

# This is incorrect. It essentially gives equal weight to all catchment
# densities, regardless of how large or small it is.

# OUTLIERS INCLUDED
d |>
  na.omit() |>
  dplyr::group_by(study_area) |>
  dplyr::summarise(popov_d = round(mean(density, na.rm = T), 3),
                   popov_d_lwr = round(mean(density_lwr, na.rm = T), 3),
                   popov_d_upr = round(mean(density_upr, na.rm = T), 3),
                   burger_a_d = round(mean(density_a, na.rm = T), 3),
                   burger_b_d = round(mean(density_b, na.rm = T), 3),
                   burger_c_d = round(mean(density_c, na.rm = T), 3),
                   popov_n = sum(N),
                   burger_n = sum(n_years_surveyed)) |>
  dplyr::mutate(popov_d = paste0(popov_d, " [", popov_d_lwr, "-", popov_d_upr, "]")) |>
  dplyr::mutate(burger_d = paste0(burger_b_d, " [", burger_a_d, "-", burger_c_d, "]")) |>
  dplyr::select(study_area, popov_d, burger_d, popov_n, burger_n) |>
  knitr::kable()


d |>
  na.omit() |>
  dplyr::group_by(study_area) |>
  dplyr::summarise(popov_d = mean(density, na.rm = T),
                   popov_d_lwr = mean(density_lwr, na.rm = T),
                   popov_d_upr = mean(density_upr, na.rm = T),
                   burger_a_d = mean(density_a, na.rm = T),
                   burger_b_d = mean(density_b, na.rm = T),
                   burger_c_d = mean(density_c, na.rm = T)) |>
  ggplot(aes(x = popov_d, y = burger_b_d)) +
  geom_point() +
  geom_linerange(aes(xmin = popov_d_lwr, xmax = popov_d_upr),
                 alpha = 0.3) +
  geom_linerange(aes(ymin = burger_a_d, ymax = burger_c_d),
                 alpha = 0.3) +
  geom_abline(slope = 1) +
  ggrepel::geom_text_repel(aes(label = study_area),
                           min.segment.length = 0) +
  scale_x_continuous(limits = c(0, 0.23), 
                     breaks = seq(0, 0.5, 0.01)) +
  scale_y_continuous(limits = c(0, 0.1),
                     breaks = seq(0, 0.5, 0.01)) +
  coord_fixed() +
  labs(x = "Popov density estimates",
       y = "Burger density estimates",
       caption = "Error bars in X dimension show 95% CI Popov density est.
       Error bars in Y dimension show Burger min/max density est based on definition of habitat used.
       Outliers included.")

# OUTLIERS EXCLUDED
d |>
  na.omit() |>
  dplyr::filter(outlier_yn == FALSE) |>
  dplyr::group_by(study_area) |>
  dplyr::summarise(popov_d = mean(density, na.rm = T),
                   popov_d_lwr = mean(density_lwr, na.rm = T),
                   popov_d_upr = mean(density_upr, na.rm = T),
                   burger_a_d = mean(density_a, na.rm = T),
                   burger_b_d = mean(density_b, na.rm = T),
                   burger_c_d = mean(density_c, na.rm = T)) |>
  ggplot(aes(x = popov_d, y = burger_b_d)) +
  geom_point() +
  geom_linerange(aes(xmin = popov_d_lwr, xmax = popov_d_upr),
                 alpha = 0.3) +
  geom_linerange(aes(ymin = burger_a_d, ymax = burger_c_d),
                 alpha = 0.3) +
  geom_abline(slope = 1) +
  ggrepel::geom_text_repel(aes(label = study_area),
                           min.segment.length = 0) +
  scale_x_continuous(limits = c(0, 0.23), 
                     breaks = seq(0, 0.5, 0.01)) +
  scale_y_continuous(limits = c(0, 0.1),
                     breaks = seq(0, 0.5, 0.01)) +
  coord_fixed() +
  labs(x = "Popov density estimates",
       y = "Burger density estimates",
       caption = "Error bars in X dimension show 95% CI Popov density est.
       Error bars in Y dimension show Burger min/max density est based on definition of habitat used.
       Outliers excluded.")


# Regional, all Burger data minus outliers. Does it match Table A3-2?

# SUM method - the better method to use.
b2002_density |>
  dplyr::mutate(outlier_yn = grepl("\\*\\*", station_burger)) |>
  dplyr::filter(outlier_yn == FALSE) |>
  dplyr::group_by(study_area) |>
  dplyr::summarise(burger_a_d = sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_a, na.rm = TRUE),
                   burger_b_d = sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_b, na.rm = TRUE),
                   burger_c_d = sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_c, na.rm = TRUE),)

# MEAN method - this is the one that lines up w what's reported in original pub.
b2002_density |>
  dplyr::mutate(outlier_yn = grepl("\\*\\*", station_burger)) |>
  dplyr::filter(outlier_yn == FALSE) |>
  dplyr::group_by(study_area) |>
  dplyr::summarise(N = dplyr::n(),
                   burger_a_d = mean(density_a, na.rm = TRUE),
                   burger_b_d = mean(density_b, na.rm = TRUE),
                   burger_c_d = mean(density_c, na.rm = TRUE))


# POOLED ----

# Across the whole study area, not the regions
# Includes only sites that were in both studies

total_area <- sum(ixn$reg_area_ha)

# SUM method
d |>
  dplyr::summarise(popov_d = round(sum(bootmean) / total_area, 3),
                   popov_d_lwr = round(sum(boot_min, na.rm = TRUE) / total_area, 3),
                   popov_d_upr = round(sum(boot_max, na.rm = TRUE) / total_area, 3),
                   burger_a_d = round(sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_a, na.rm = TRUE), 3),
                   burger_b_d = round(sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_b, na.rm = TRUE), 3),
                   burger_c_d = round(sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_c, na.rm = TRUE), 3),
                   popov_n = sum(N),
                   burger_n = sum(n_years_surveyed)) |>
  dplyr::mutate(popov_d = paste0(popov_d, " [", popov_d_lwr, "-", popov_d_upr, "]")) |>
  dplyr::mutate(burger_d = paste0(burger_b_d, " [", burger_a_d, "-", burger_c_d, "]")) |>
  dplyr::select(popov_d, burger_d, popov_n, burger_n) |>
  knitr::kable()


# MEAN method
d |>
  dplyr::summarise(popov_d = round(mean(density, na.rm = T), 3),
                   popov_d_lwr = round(mean(density_lwr, na.rm = T), 3),
                   popov_d_upr = round(mean(density_upr, na.rm = T), 3),
                   burger_a_d = round(mean(density_a, na.rm = T), 3),
                   burger_b_d = round(mean(density_b, na.rm = T), 3),
                   burger_c_d = round(mean(density_c, na.rm = T), 3),
                   popov_n = sum(N),
                   burger_n = sum(n_years_surveyed)) |>
  dplyr::mutate(popov_d = paste0(popov_d, " [", popov_d_lwr, "-", popov_d_upr, "]")) |>
  dplyr::mutate(burger_d = paste0(burger_b_d, " [", burger_a_d, "-", burger_c_d, "]")) |>
  dplyr::select(popov_d, burger_d, popov_n, burger_n) |>
  knitr::kable()



# All Burger (2002) - including those not in our study
# SUM method
b2002_density |>
  dplyr::summarise(burger_a_d = round(sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_a, na.rm = TRUE), 3),
                   burger_b_d = round(sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_b, na.rm = TRUE), 3),
                   burger_c_d = round(sum(mamu_count, na.rm = TRUE) / sum(habitat_ha_c, na.rm = TRUE), 3),
                   burger_n = sum(n_years_surveyed)) |>
  dplyr::mutate(burger_d = paste0(burger_b_d, " [", burger_a_d, "-", burger_c_d, "]")) |>
  dplyr::select(burger_d, burger_n) |>
  knitr::kable()

# MEAN method
b2002_density |>
  dplyr::mutate(outlier_yn = grepl("\\*\\*", station_burger)) |>
  #dplyr::filter(outlier_yn == FALSE) |>
  dplyr::summarise(burger_a_d = round(mean(density_a, na.rm = TRUE), 3),
                   burger_b_d = round(mean(density_b, na.rm = TRUE), 3),
                   burger_c_d = round(mean(density_c, na.rm = TRUE), 3),
                   burger_n = sum(n_years_surveyed)) |>
  dplyr::mutate(burger_d = paste0(burger_b_d, " [", burger_a_d, "-", burger_c_d, "]")) |>
  dplyr::select(burger_d, burger_n) |>
  knitr::kable()

# All Popov
tar_read(total_density)


## 03 CATCHMENT VS HABITAT AREA ----

plot(b2002_density$watershed_area_ha, b2002_density$habitat_ha_a)

b2002_density |>
  dplyr::select(study_area, watershed_area_ha, habitat_ha_a, habitat_ha_b, habitat_ha_c) |>
  tidyr::pivot_longer(cols = c(habitat_ha_a:habitat_ha_c),
                      names_to = "habitat_method",
                      values_to = "habitat_ha") |>
  ggplot(aes(x = watershed_area_ha,
               y = habitat_ha,
               color = habitat_method)) +
  geom_point() +
  geom_smooth(method = "lm") +
  khroma::scale_color_okabeito()


plot(cc_density$area_ha, cc_density$sh_area_ha)

## 04 CATCHMENT VS MAMU ----

plot(b2002_density$watershed_area_ha, b2002_density$mamu_count)

b2002_density |>
  dplyr::select(study_area, watershed_area_ha, mamu_count) |>
  ggplot(aes(x = watershed_area_ha,
             y = mamu_count,
             color = study_area)) +
  geom_point() +
  geom_smooth(method = "lm") +
  khroma::scale_color_okabeito()

d |>
  dplyr::select(study_area, area_ha, bootmean) |>
  ggplot(aes(x = area_ha,
             y = bootmean,
             color = study_area)) +
  geom_point() +
  geom_smooth(method = "lm") +
  khroma::scale_color_okabeito()


cc_density |>
  units::drop_units() |>
  dplyr::select(region, area_ha, bootmean) |>
  ggplot(aes(x = area_ha,
             y = bootmean,
             color = region)) +
  geom_point() +
  geom_smooth(method = "lm") +
  khroma::scale_color_okabeito()
