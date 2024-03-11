# 3. CALCULATE MAMU DENSITY AND POPULATION

# Following the creation of the catchments and standardization
# of the radar survey data, we can now calculate the density
# of birds within each catchment. From there, we will group 
# each catchment by conservation region (based on the 
# conservation region that the survey station was based in)
# and calculate the mean density of birds per conservation
# region. We will then multiply this mean density of birds
# per region times the amount of available habitat area to 
# come up with a population estimate across all of BC.

# 01 SETUP ----------------------------------------------------------------

# This script assumes you have run scripts 1 and 2, and that
# the predicted population + GIS files are either loaded into
# the present R environment OR can be loaded in from the 'data'
# directory.

# TODO: change this to load w final file name/update above 
# description as needed once scripts are finalized.
# Note: you need to have 'cons_reg' loaded up from 01-prepare-survey-data
# TODO: rename `c30_2` var/file as needed 
catchments <- sf::st_read('data/radar_derived_catchments_v2.gpkg')
c30_2 <- terra::rast("../Habitat/30km_from_shore_OVER_mountains_below_elev_cutoffs.tiff")

# TODO: update this description w file name
# Chop up c30_2 (aka 'mamu containment zone' aka 'mcz') by region
# and count up the number of raster cells in c30_2 using exact_extract()
cons_reg$mcz_cell_count <- exactextractr::exact_extract(c30_2, cons_reg, 'count')

# Calculate number of hectares of habitat are within each conservation region
# Cell resolution is 500m x 500m, aka 250,000 m2, aka 25 ha per cell
cons_reg$mcz_ha <- cons_reg$mcz_cell_count * 25

# 01-1 Merge prediction datasets with catchments ----

# Merge in `p` objects with catchments
p <- merge(p, catchments, by = "site", all.x = TRUE)
p_allyears <- merge(p_allyears, catchments, by = "site", all.x = TRUE)
p2022 <- merge(p2022, catchments, by = "site", all.x = TRUE)

# Check out the habitat area-MAMU numbers relationship
# Notably, increase in habitat area correlates to more
# birds for all regions except Vancouver Island.
library(units)
ggplot(p2022, aes(x = area_ha, y = pred_fit, color = region)) +
  geom_point() +
  geom_smooth(method = "lm") +
  #facet_wrap(~region) +
  labs(x = "Area (ha)",
       y = "2022 MAMU population estimate")


# 02 CALCULATE DENSITY ----------------------------------------------------

# Calculate the average density of birds per catchment by region
# I.e., mean number birds and mean area of catchments (in hectares)
# for each region

# First, calculate the density of birds in each catchment.
p$density <- p$pred_fit / p$area_ha
p$density_lwr <- p$ci_lwr / p$area_ha
p$density_upr <- p$ci_upr / p$area_ha

p_allyears$density <- p_allyears$pred_fit / p_allyears$area_ha
p_allyears$density_lwr <- p_allyears$ci_lwr / p_allyears$area_ha
p_allyears$density_upr <- p_allyears$ci_upr / p_allyears$area_ha

p2022$density <- p2022$pred_fit / p2022$area_ha
p2022$density_lwr <- p2022$ci_lwr / p2022$area_ha
p2022$density_upr <- p2022$ci_upr / p2022$area_ha

# Next, calculate the mean density per region (for datasets with 
# multiple years, calculate mean density per region by year).
library(dplyr)

density <- p %>% 
  filter(!is.na(area_ha)) %>%
  group_by(year, region) %>%
  summarize(density = mean(density),
            density_lwr = mean(density_lwr),
            density_upr = mean(density_upr))

density_allyears <- p_allyears %>% 
  filter(!is.na(area_ha)) %>%
  group_by(year, region) %>%
  summarize(density = mean(density),
            density_lwr = mean(density_lwr),
            density_upr = mean(density_upr))

density2022 <- p2022 %>% 
  filter(!is.na(area_ha)) %>%
  group_by(region) %>% 
  summarize(density = mean(density),
            density_lwr = mean(density_lwr),
            density_upr = mean(density_upr))

# Merge conservation region with density
density <- merge(density, cons_reg, by = "region", all.x = TRUE)
density_allyears <- merge(density_allyears, cons_reg, by = "region", all.x = TRUE)
density2022 <- merge(density2022, cons_reg, by = "region", all.x = TRUE)

# Finally, extrapolate the density by region to the total area
density$mamu <- density$density * density$mcz_ha
density$mamu_lwr <- density$density_lwr * density$mcz_ha
density$mamu_upr <- density$density_upr * density$mcz_ha

density_allyears$mamu <- density_allyears$density * density_allyears$mcz_ha
density_allyears$mamu_lwr <- density_allyears$density_lwr * density_allyears$mcz_ha
density_allyears$mamu_upr <- density_allyears$density_upr * density_allyears$mcz_ha

density2022$mamu <- density2022$density * density2022$mcz_ha
density2022$mamu_lwr <- density2022$density_lwr * density2022$mcz_ha
density2022$mamu_upr <- density2022$density_upr * density2022$mcz_ha

# Examine the data
options(scipen = 999)
head(density)
density2022

colSums(density[density$year == 1996, c("mamu", "mamu_lwr", "mamu_upr")])
colSums(density[density$year == 2005, c("mamu", "mamu_lwr", "mamu_upr")])
colSums(density[density$year == 2015, c("mamu", "mamu_lwr", "mamu_upr")])
colSums(density[density$year == 2022, c("mamu", "mamu_lwr", "mamu_upr")])
colSums(density2022[,c("mamu", "mamu_lwr", "mamu_upr")]) # TODO: why the heck are these different

density %>%
  filter(year == 2022) %>%
  select(region, mamu, mamu_lwr, mamu_upr, mcz_ha, density) %>%
  group_by(region) %>%
  summarize(mamu = round(sum(mamu)),
            mamu_lwr = round(sum(mamu_lwr)),
            mamu_upr = round(sum(mamu_upr)),
            area_ha = mean(mcz_ha) / 1000000,
            density = mean(density) * 1000) %>%
  arrange(desc(mamu))  #%>%
  #select(-region) %>%
  #colSums()

density2022 %>%
  select(region, mamu, mamu_lwr, mamu_upr, mcz_ha, density) %>%
  mutate(mamu = round(mamu),
         mamu_lwr = round(mamu_lwr),
         mamu_upr = round(mamu_upr),
         mcz_ha = mcz_ha / 1000000,
         density = density * 1000) %>%
  arrange(desc(mamu)) #%>%
  #select(-region) %>%
  #colSums() # note ignore the density column here, you need to calc that by hand
  
# Plot it
density %>%
  select(year, mamu, mamu_lwr, mamu_upr) %>%
  mutate(year_n = as.numeric(year) + 1995) %>%
  dplyr::group_by(year_n) %>%
  summarise(mamu = sum(mamu),
            mamu_lwr = sum(mamu_lwr),
            mamu_upr = sum(mamu_upr)) %>%
  ggplot(aes(x = year_n, y = mamu)) +
  geom_ribbon(aes(ymin = mamu_lwr,
                  ymax = mamu_upr),
              fill = "lightgrey",
              alpha = 0.3) +
  geom_point() +
  geom_line() +
  xlab("Year") +
  ylab("BC murrelet population") +
  theme_minimal()
