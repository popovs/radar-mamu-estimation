# 2. STANDARDIZE MAMU RADAR COUNTS

# In this script, MAMU counts from radar stations will be
# standardized following a generalized linear mixed model
# approach. The goal is to standardize radar counts to 
# eliminate any variation in counts that may be caused by
# differing seasons/day of year and/or suboptimal radar
# setups (smaller radar radius or lower radar tilt).


# 01 SETUP ----------------------------------------------------------------

# 01-1 Read in 'surveys' ----

# The first step of loading and cleaning the data is
# handled by the MAMU::process_radar_data() function. 

#library(devtools)
#devtools::install_github("popovs/MAMU")
library(MAMU)

# Load the survey data
# See `?MAMU::process_radar_data` for a description of how data is processed
# See `View(MAMU::process_radar_data)` for the data cleaning code "under the hood"
surveys <- MAMU::process_radar_data("data/ECCC_FLNR_MAMU-RadarData-20240307.xlsx")

# 01-2 Read in 'conservation regions' ----

# Read in the conservation regions data
library(sf)
cons_reg <- st_read("GIS/cons_reg.shp")
cons_reg <- st_transform(cons_reg, 3005) # Set BC Albers projection
cons_reg <- cons_reg[,1] # Drop everything but first column
names(cons_reg) <- c("region", "geometry") # Rename cols


# 01-3 Prepare `s` dataframe ----

# Model dataframe will be s
# Response var is 'mamuinpd' (i.e., count of incoming MAMU pre-dawn)
s <- data.frame(y = surveys$mamuinpd)

# Clean up any vars as needed
s$observer <- as.factor(surveys$observer)
s[["observer"]][is.na(s$observer)] <- "Unknown" # 15 records have NA observer. Give them value of "Unknown" to match other unknown observer records.

# Create survey spatial object
s$lat <- surveys$lat
s$lon <- surveys$lon
s <- s[!is.na(s$lon), ]
s <- st_as_sf(s, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Of potential interest - label Metro Vancouver points
# (including points in the interior where a bird would have
# to fly through Metro Vancouver)
s$yvr <- ifelse(s$lat < 49.7282 & s$lon > -123.2937,
                "Metro Vancouver",
                "Other")

# Now transform CRS
s <- st_transform(s, 3005) # Set BC Albers projection

# Add cleaned site name
s$site <- as.factor(surveys$new_name)

# Make other factors of interest
s$loc <- as.factor(surveys$loc)
s$year <- as.factor(surveys$year)

# Re-scale numerics
s$s_year <- scale(as.numeric(s$year))[,1]

s$tilt <- as.numeric(surveys$tilt)
s$s_tilt <- scale(s$tilt)[,1]

s$radius <- as.numeric(surveys$radius)
s$s_radius <- scale(s$radius)[,1]

s$doy <- as.numeric(surveys$doy)
s$s_doy <- scale(s$doy)[,1]

# Chuck in rescaled lat/lon for good measure
s$s_lat <- scale(s$lat)[,1]
s$s_lon <- scale(s$lon)[,1]

# Keep only complete cases
s <- s[complete.cases(st_drop_geometry(s)),]

# Add conservation region to each record
# Note that it's important to do this at the end, as the intersection
# will change order of the records around a little bit!
s <- st_intersection(s, cons_reg)

# Check work
plot(s[, "region"])


# Add various metrics of survey effort per catchment
library(dplyr)
library(ggplot2)

# how many times was a site sampled across ALL years?
s <- s %>% group_by(site) %>% mutate(total_effort = n())
# how many times was a site sampled within a year?
s <- s %>%
  group_by(site, year) %>%
  mutate(within_year_effort = n())
# how many times was a site sampled up until that year?
s <- s %>%
  arrange(site, year) %>%
  group_by(site) %>%
  mutate(count = 1,
         rolling_effort = cumsum(count)) %>%
  select(-count)

ggplot(s, aes(x = total_effort, y = y, color = region)) +
  geom_point(show.legend = FALSE) +
  geom_smooth(method = "lm") +
  ggtitle("MAMU vs Total count surveys across all years") +
  xlab("N surveys (total_effort)") +
  ylab("Count of MAMU") +
  scale_y_continuous(trans = "log10") +
  annotation_logticks(sides = "l") +
  facet_wrap(~region) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(s, aes(x = within_year_effort, y = y, color = region)) +
  geom_point(show.legend = FALSE) +
  geom_smooth(method = "lm") +
  ggtitle("MAMU vs Count of surveys within a year") +
  xlab("N surveys (within_year_effort)") +
  ylab("Count of MAMU") +
  scale_y_continuous(trans = "log10") +
  annotation_logticks(sides = "l") +
  facet_wrap(~region) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(s, aes(x = rolling_effort, y = y, color = region)) +
  geom_point(show.legend = FALSE) +
  geom_smooth(method = "lm") +
  ggtitle("MAMU vs Rolling count of surveys") +
  xlab("N surveys (rolling_effort)") +
  ylab("Count of MAMU") +
  scale_y_continuous(trans = "log10") +
  annotation_logticks(sides = "l") +
  facet_wrap(~region) +
  theme_minimal() +
  theme(legend.position = "none")

# Check to see if some years have unusually high or low sampling effort
# For example, we know that ~2012-2015 there weren't as many surveys done,
# and this may impact our results down the line
ggplot(s, aes(x = as.numeric(as.character(year)), y = within_year_effort, color = region)) +
  geom_point(show.legend = FALSE) +
  #geom_smooth(method = "lm") +
  ggtitle("Survey effort vs Year") +
  xlab("Year") +
  ylab("N surveys (within_year_effort)") +
  facet_wrap(~region) +
  theme_minimal() +
  theme(legend.position = "none")

# Cut out sites/surveys with only one survey over the entire study period
s <- s[s$total_effort > 1,]

# Cut out sites/surveys that were only ever surveyed in one year
# of the entire study period (can't get yearly trends from these)
efforty <- aggregate(year ~ site, s, function(x) length(unique(x)))
efforty <- efforty[efforty$year == 1, ] # this is now a list of sites that were only ever sampled within a single year
#View(s[s$site %in% efforty$site, ]) # double check that the surveys were only within 1 year per site
s <- s[!(s$site %in% efforty$site), ] # chuck em
rm(efforty)


# 01-4 Create `stn` object ----

# Create object of mean coordinate per survey station
# Also include count column to get n surveys at each
# unique site
stn <- s %>%
  st_transform(crs = 4326) %>%
  select(site, region, loc, yvr) %>%
  group_by(site) %>%
  aggregate(.,
            by = list(.$site),
            function(x) x = x[1]) %>%
  st_centroid() %>%
  select(-Group.1) %>%
  mutate("region" = region,
         "lat" = st_coordinates(.)[,2],
         "lon" = st_coordinates(.)[,1],
         "s_lat" = scale(lat)[,1],
         "s_lon" = scale(lon)[,1]) %>%
  select(site, region, loc, yvr, lat, lon, s_lat, s_lon)


# 01-5 Summary statistics ----

summary(s$total_effort)
summary(s$y)
psych::describeBy(s$y, s$region, mat = TRUE, quant = c(0.25, 0.75))
sum(s$y == 0)
sum(s$y == 0) / length(s$y) * 100



# 02 BUILD MODEL ----------------------------------------------------------

# Here, we will actually build our negative binomial model using the glmmTMB
# library. Several other (failed) candidate models are included below for
# reference, though they are commented out.

library(glmmTMB)

# WITH controlling for total sampling effort
# We're assuming a log relationship with sampling effort. At some point, 
# no matter how much you sample, you are not going to find more birds.
# Note: adding year:region interaction term makes the residuals a bit funky
m <- glmmTMB(y ~ s_doy + I(s_doy^2) + s_year + region + lat + lon
             + s_tilt + s_radius + (1|observer) + (s_year|site) + log(total_effort)
             + s_doy:region  + I(s_doy^2):region + (1|year),
             data = s,
             ziformula = ~1, # YES zero-inflation
             family = nbinom1)

# WITHOUT controlling for total sampling effort
m2 <- glmmTMB(y ~ s_doy + I(s_doy^2) + s_year + region + lat + lon
              + s_tilt + s_radius + (1|observer) + (s_year|site)
              + s_doy:region  + I(s_doy^2):region + (1|year),
              data = s,
              ziformula = ~1, # YES zero-inflation
              family = nbinom1)

# Sampling effort as an offset
# m3 <- glmmTMB(y ~ s_doy + I(s_doy^2) + s_year * region + lat + lon
#               + s_tilt + s_radius + (1|observer) + (s_year|site)
#               + s_doy:region  + I(s_doy^2):region,
#               data = s,
#               ziformula = ~1, # YES zero-inflation
#               offset = log(total_effort),
#               family = nbinom1)

# Sampling effort as RE
# This is not allowed technically, continuous var RE's are a no-no
# m3 <- glmmTMB(y ~ s_doy + I(s_doy^2) + s_year * region + lat + lon
#               + s_tilt + s_radius + (1|observer) + (s_year|site) + (1|total_effort)
#               + s_doy:region  + I(s_doy^2):region,
#               data = s,
#               ziformula = ~1, # YES zero-inflation
#               family = nbinom1)

# Sampling effort varies by region
# m4 <- glmmTMB(y ~ s_doy + I(s_doy^2) + s_year * region + lat + lon
#               + s_tilt + s_radius + (1|observer) + (s_year|site) + (region|total_effort)
#               + s_doy:region  + I(s_doy^2):region,
#               data = s,
#               ziformula = ~1, # YES zero-inflation
#               family = nbinom1)

# Sampling effort varies by site
# Not enough site-effort combos for this to work
# m5 <- glmmTMB(y ~ s_doy + I(s_doy^2) + s_year * region + lat + lon
#              + s_tilt + s_radius + (1|observer) + (s_year|site) + (site|total_effort)
#              + s_doy:region  + I(s_doy^2):region,
#              data = s,
#              ziformula = ~0, # NO zero-inflation
#              family = nbinom1)

# Sampling effort is controlled for by year
# m6 <- glmmTMB(y ~ s_doy + I(s_doy^2) + s_year * region + lat + lon
#               + s_tilt + s_radius + (1|observer) + (s_year|site) + log(within_year_effort)
#               + s_doy:region  + I(s_doy^2):region,
#               data = s,
#               ziformula = ~1, # YES zero-inflation
#               family = nbinom1)

bbmle::AICtab(m, m2)
DHARMa::testResiduals(m)
DHARMa::testResiduals(m2)

# CHANGE HERE TO PREDICT VALUES FROM DIFFERENT MODEL
top_model <- m2

# Clean up `m*` objects
rm(list = ls()[grep("^m", ls())])

# 03 BUILD NEWDATA --------------------------------------------------------------

# Build prediction newdata

# In this section, we will build our dataset to feed into our candidate
# model for prediction purposes. Note we are not making ecological inferences
# from this model, but rather simply trying to create a predicted resposne
# dataset that controls for variation in sampling effort and radar setups
# across the years.

# 03-1 ALL YEARS, SINGLE MEANS ----
# Mean covariate values across ALL years assigned to the prediction dataset
# This is used to get decadal trends in bird counts, smoothing over yearly variation

# Note for DOY: mean DOY for ALL years by region.
doy_reg <- aggregate(doy ~ region, s, FUN = function(x) round(mean(x)))
doy_reg <- merge(doy_reg, unique(st_drop_geometry(s[s$doy %in% doy_reg$doy, c("doy", "s_doy")])), by = "doy")

# We want to keep the same scaled years btwn original dataset in newdata one
s_years <- aggregate(s_year ~ year, s, mean) 

p_allyears <- expand.grid(site = unique(s$site),
                          year = unique(s$year))
names(p_allyears) <- c("site", "year")
p_allyears <- plyr::join(p_allyears, unique(st_drop_geometry(s[,c("site", "region")])), by = "site", match = "first") # note it's crucial to keep the same factor levels as in original dataset
p_allyears$observer <- s$observer[1] # Standardize our observer to Bernard
p_allyears <- merge(p_allyears, doy_reg, by = "region")
p_allyears <- merge(p_allyears, s_years, by = "year")
p_allyears$tilt <- max(s$tilt) # max tilt value so we don't 'undercount' birds
p_allyears$s_tilt <- max(s$s_tilt)
p_allyears$radius <- max(s$radius) # max radius value so we don't 'undercount' birds
p_allyears$s_radius <- max(s$s_radius)
p_allyears$total_effort <- median(s$total_effort) # median count (sample size) of # of surveys per site to standardize surveying effort

# Grab mean lat/long for each survey from the stn df
p_allyears <- merge(p_allyears, st_drop_geometry(stn[,names(stn) != "region"]), by = "site")


# 03-2 ALL YEARS, MEAN VARIABLES PER YEAR ----
# Mean covariate values PER YEAR assigned to each year of the prediction dataset
# This is used to get an accurate estimate for yearly bird population counts, including
# yearly variation in counts

# Mean variables by year and region
year_reg_vars <- s %>%
  sf::st_drop_geometry() %>%
  group_by(region, year) %>%
  summarise(doy = round(mean(doy)),
            total_effort = median(s$total_effort)) # we want to standardize effort across all records - we don't want to see any variation caused by survey effort across years
# Merge in the original s_year and s_doy values
year_reg_vars <- merge(year_reg_vars, unique(st_drop_geometry(s[s$doy %in% year_reg_vars$doy, c("doy", "s_doy")])), by = "doy")
year_reg_vars <- merge(year_reg_vars, unique(st_drop_geometry(s[s$year %in% year_reg_vars$year, c("year", "s_year")])), by = "year")

# Mean variables by year across all regions 
# Some regions aren't present in certain years, so we'll use the overall 
# yearly means to fill in those gaps).
yearly_vars <- s %>%
  sf::st_drop_geometry() %>%
  group_by(year) %>%
  summarise(doy = round(mean(doy)),
            total_effort = median(s$total_effort)) # we want to standardize effort across all records - we don't want to see any variation caused by survey effort across years
# Merge in the original s_year and s_doy values
yearly_vars <- merge(yearly_vars, unique(st_drop_geometry(s[s$doy %in% yearly_vars$doy, c("doy", "s_doy")])), by = "doy")
yearly_vars <- merge(yearly_vars, unique(st_drop_geometry(s[s$year %in% yearly_vars$year, c("year", "s_year")])), by = "year")


# Unfortunately, not all regions were present in each year.
# So for regions that were not present in a given year, we'll
# just fill in the yearly means in there.
yr_region <- expand.grid(region = unique(s$region),
                         year = unique(s$year))
p_x <- merge(year_reg_vars, yr_region, by = c("year", "region"), all.y = TRUE)

p_x1 <- p_x[which(!is.na(p_x$doy)),]
p_x2 <- p_x[which(is.na(p_x$doy)),c("year", "region")] # pull out all NA value year-region combos. For these we'll use the yearly mean overall
p_x2 <- merge(p_x2, yearly_vars, by = "year")

p_x <- rbind(p_x1, p_x2) # 216 rows
rm(p_x1, p_x2)

# We'll now merge p_x (which contains our mean variables either by region/year or just
# year) with p, our matrix of all possible site-year combos.
p <- expand.grid(site = unique(s$site),
                 year = unique(s$year))
names(p) <- c("site", "year")
p <- merge(p, unique(sf::st_drop_geometry(s[s,c("region", "site")])), by = "site")
p <- merge(p, p_x, by = c("region", "year"))
rm(p_x)
p$tilt <- max(s$tilt) # max tilt value so we don't 'undercount' birds
p$s_tilt <- max(s$s_tilt)
p$radius <- max(s$radius) # max radius value so we don't 'undercount' birds
p$s_radius <- max(s$s_radius)
p$observer <- s$observer[1] # Standardize our observer to Bernard
# Merge in station info
p <- merge(p, st_drop_geometry(stn[, names(stn) != "region"]), by = "site")

rm(year_reg_vars, yearly_vars, yr_region)


# 03-3 2022 ONLY ----

# Using doy_reg above 
# Not all regions are present in 2022 data

p2022 <- unique(st_drop_geometry(s[,c("region", "site")]))
p2022$observer <- s$observer[1] # Standardize our observer to Bernard
p2022$year <- unique(s[["year"]][s$year == 2022]) # do this to keep the same factor level
p2022$s_year <- unique(s[["s_year"]][s$year == 2022])
p2022 <- merge(p2022, doy_reg, by = "region")
p2022$tilt <- max(s$tilt) # max tilt value so we don't 'undercount' birds
p2022$s_tilt <- max(s$s_tilt)
p2022$radius <- max(s$radius) # max radius value so we don't 'undercount' birds
p2022$s_radius <- max(s$s_radius)
p2022$total_effort <- median(s$total_effort)
p2022 <- merge(p2022, st_drop_geometry(stn[, names(stn) != "region"]), by = "site")



# 04 PREDICT --------------------------------------------------------------

# 04-1 Define the prediction function ----

# Return predictions of model m, including confidence intervals.
# See Ben Bolker's glmmTMB FAQ for code methodology:
# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#glmmtmb

# if newdata is not supplied, it will just use original
# model `m` dataset.
glmmTMB_preds <- function(m, newdata) {
  if (missing(newdata)) {
    # Design matrix of fixed effects
    mm <- model.matrix(delete.response(terms(m)), m$frame)
  } else {
    # Design matrix of fixed effects
    mm <- model.matrix(delete.response(terms(m)), newdata)
  }
  
  # Linear predictor; back transformed w inverse link
  # function (in this case, `exp()` for negative binomial model)
  # `%*%` runs matrix multiplication
  preds <- as.data.frame(exp(drop(mm %*% fixef(m)[["cond"]]))) # should be the exact same as `predict(type = "response")` method
  preds$var <- diag(mm %*% vcov(m)[["cond"]] %*% t(mm))
  preds$SE <- sqrt(preds$var)
  preds$SE2 <- sqrt(preds$var + sigma(m)^2)
  
  preds$pred_fit <- round(unname(preds[[1]]))
  preds$ci_lwr <- round(exp(log(preds$pred_fit) - (1.959964 * preds$SE)))
  preds$ci_upr <- round(exp(log(preds$pred_fit) + (1.959964 * preds$SE)))
  
  preds <- as.data.frame(preds)
  preds <- preds[,c("pred_fit", "ci_lwr", "ci_upr")]
  return(preds)
}



# 04-2 Predicted vs observed ----

# Compare original dataset observed values to predicted
# Note that we don't expect it to match - we are smoothing out
# fluctuations caused by unoptimal radar setups etc., after all,
# so we expect some of the predicted to be higher than observed.
s <- cbind(s, glmmTMB_preds(m = top_model))
ggplot(s, aes(x = pred_fit, y = y)) + geom_point() + geom_abline() + coord_fixed()


# 04-3 Prediction datasets ----

# Keep in mind these are predictions for our sampled
# catchments only - not extrapolated to the whole province

# All years 
# Mean covariate values across ALL years applied to entire prediction dataset
p_allyears <- cbind(p_allyears, glmmTMB_preds(m = top_model, newdata = p_allyears))

# Yearly variation
# Mean covariate values PER YEAR applied to each year in the prediction dataset
p <- cbind(p, glmmTMB_preds(m = top_model, newdata = p))

p %>%
  group_by(year) %>%
  summarise(count = sum(pred_fit)) %>%
  print(n = 27)

# 2022 only
p2022 <- cbind(p2022, glmmTMB_preds(m = top_model, newdata = p2022))

# Note predictions for p[p$year == 2022,] will NOT equal p2022.
# The mean DOY and year is different there. So, it'll
# be close, but not an exact match.
sum(p2022$pred_fit) 
colSums(p2022[,c("pred_fit", "ci_lwr", "ci_upr")])

p2022 %>%
  dplyr::group_by(region) %>%
  summarise(count = sum(pred_fit),
            lwr = sum(ci_lwr),
            upr = sum(ci_upr))


# 05 CLEAN UP -------------------------------------------------------------

rm(doy_reg, s_years)
