##### ----------------------------------------------------------------------------------------------------------
##### TITLE:        DATA ANALYSIS AND STATISTICAL CODE FOR: CONFORMITY OF SPECIES-AREA RELATIONSHIPS IN ATOLLS
##### AUTHORS:      SEBASTIAN STEIBL, SIMON STEIGER, LUIS VALENTE, JAMES C. RUSSELL
##### LAST EDITED:  26 Mar 2024
##### ----------------------------------------------------------------------------------------------------------

## Description:
#   This script performs the statistical analysis for Steibl et al. (2025) Rainfall increases conformity and strength of
#   species-area relationships in atolls. It loads and pre-processes the data, executes the analysis, and generates tables
#   and figures for the manuscript.

## Reproducibility:
#   This project uses the 'renv' package to manage package dependencies and ensure a reproducible computational
#   environment. To set up the environment on your machine, please run: renv::restore()
#   This command reads the 'renv.lock' file and installs the exact versions of packages used in this analysis.

# Data:
#   The data used for this analysis is stored in the 'data/' directory.
#   We use the 'here' package to construct file paths relative to the project root,
#   ensuring that data files are correctly located regardless of the working directory.


####################################
### 1. ----- Load packages ----- ###
####################################

library("tidyverse")
library("ggh4x")
library("here")
library("patchwork")
library("brms")
library("coda")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("raster")

################################
### ----- 2. Load data ----- ###
################################

dat <- read.csv("data/SAR_species_matrix.csv")
isl <- read.csv("data/SAR_env-data_islets.csv")
atoll <- read.csv("data/SAR_env-data_atolls.csv")


###########################################
### ----- 3. Format, prepare data ----- ###
###########################################

# convert longitudes for 'sf' format
atoll$long <- ifelse(atoll$long < 0, atoll$long+360, atoll$long)

# Turn atoll and islet identifiers into factors
dat <- dat %>% mutate(across(1:2, as.factor))
isl <- isl %>% mutate(across(1:2, as.factor))
atoll <- atoll %>% mutate(across(1, as.factor))

# Turn pres-abs matrix into numerical values
dat <- dat %>% mutate(across(!c(1:2), as.numeric))

# Replace NA (:= encode absence) with zeros
dat <- dat %>% mutate(across(where(is.numeric), ~ replace_na(., 0)))

# Remove introduced species (:= encoded as '2' in original data frame)
dat <- dat %>% mutate(across(where(is.numeric), ~ replace(., . == 2, 0)))

# Convert area units from km² to ha
isl$area_sqkm <- isl$area_sqkm*100
atoll$total_atoll_area_sqkm <- atoll$total_atoll_area_sqkm*100

# Calculate total land area per atoll and merge to atoll-level env data frame
isl.area <- isl %>% group_by(atoll) %>% summarise(land_area = sum(area_sqkm)) %>% as.data.frame()
atoll <- atoll %>% right_join(., isl.area, by = "atoll")
rm(isl.area) # no longer needed, keep tidy

# Calculate number of islets per atoll and merge to atoll-level env data frame
n_isl <- isl %>% group_by(atoll) %>% count(atoll) %>% as.data.frame() %>% rename(n_islands = n)
atoll <- atoll %>% right_join(., n_isl, by = "atoll")
rm(n_isl) # no longer needed, keep tidy

# Calculate isolation of atolls based on minimum distance of either to nearest high island or nearest continental landmass
atoll <- atoll %>% mutate(isolation = pmin(distance_high_island_km, distance_continent_km))

########################################################
### ----- 4. Create map plot of studied atolls ----- ###
########################################################

# Load rainfall raster file
rainf_df <- read.csv("data/rainfall_raster.csv")

# Load world map data
mapWorld <- ne_countries(scale = "medium", returnclass = "sf")
mapCoast <- ne_coastline(scale = "medium", returnclass = "sf")

# Define a custom CRS for Pacific-centered map (rotated by 180°)
crs_pacific <- "+proj=longlat +lon_wrap=180 +datum=WGS84"

pacify <- function(map, lon_0=-180) {
  map %>% 
    st_break_antimeridian(lon_0 = lon_0) %>%
    st_transform(crs = crs_pacific)
}

ne_fac <- function(f) {
  f(scale = "medium", returnclass = "sf") %>%
    st_set_crs(4326) %>%
    pacify()
}

# Re-project the world map to the Pacific-centered CRS
pacific_world <- ne_fac(ne_countries)
pacific_coast <- ne_fac(ne_coastline)

# create map layer
map <- ggplot() + 
  geom_tile(data = rainf_df, aes(x = lon, y = -lat, fill = log(rainfall * 365), colour = log(rainfall * 365))) +
  scale_fill_gradientn(colours = c("#fbe725", "#3dbc74", "#306b84","black"),
                         values = c(1.0, 0.8, 0.6, 0),
                         name = "log(rainfall)") +
  scale_colour_gradientn(colours = c("#fbe725", "#3dbc74", "#306b84","black"),
                       values = c(1.0, 0.8, 0.6, 0),
                       name = "log(rainfall)") +
  geom_sf(data = pacific_world, fill = "#EBE3D5", color = NA) +
  geom_sf(data = pacific_coast, linewidth = 0.15) +
  coord_sf(xlim = c(46.5, 305),
           ylim = c(-30, 30), expand = TRUE) +
  theme_classic() +
  geom_segment(aes(x = 34.8, xend = 34.8, y = -30.1, yend = 30.1), linewidth = 0.8) +
  geom_segment(aes(x = 59.9, xend = 300.1, y = -33, yend = -33), linewidth = 0.8) +
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(size = 7),
    axis.ticks = element_line(linewidth = 0.8),
    title = element_text(size = 8),
    legend.position = c(0.935, 0.5),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, hjust = 0.5),
    legend.key.height = unit(12, "pt"),
    legend.background = element_rect(linewidth = 0.5, linetype = "solid", colour = "black"))

# Extract coordinates for 'sf' format
spat.atoll <- atoll %>% st_as_sf(., coords = c("long", "lat"), crs = 4326)
coords.atoll <- spat.atoll %>% st_coordinates(extract())

# Create plot
p1 <- map + 
  geom_point(data = spat.atoll,
             aes(x = coords.atoll[,1], y = coords.atoll[,2]), 
             size = 1.1, shape = 21, alpha = 0.9, stroke = 1.9, colour = "#C96868")
p1

# Export as SVG vector file
#ggsave(p1, filename = "fig01_atoll_map.svg", dpi = 300, width = 158.5, height = 70, units = "mm")

################################################
### ----- 5. Calculate alpha diversity ----- ###
################################################

# ... at the islet level
isl.ric <- dat %>% 
              group_by(islet, atoll) %>%
              summarise(ric = rowSums(across(where(is.numeric)))) %>%
              as.data.frame()
isl <- isl %>% right_join(., isl.ric, by = c("islet", "atoll"))
rm(isl.ric) # no longer needed, keep tidy

# ... at the atoll level
atoll.ric <- dat %>% 
              group_by(atoll) %>% 
              summarise(ric = sum(map_lgl(across(where(is.numeric)), ~ any(. > 0)))) %>% 
              as.data.frame()
atoll <- atoll %>% right_join(., atoll.ric, by = "atoll")
rm(atoll.ric) # no longer needed, keep tidy

# Because no species-based information on floral community of Ontong Java Atoll was reported, the above formula for calculating
# atoll-level species does not apply for this atoll. The reported total number of indigenous plants in the atoll (58 species) must
# therefore be entered automatically here for this atoll
# see: Bayliss-Smith, T. (1973) Section II: Environment, chapter 2.7,  p. 47. University of Cambridge
atoll[which(atoll$atoll == "Ontong_Java"), "ric"] <- 58

# Drop n=1 island where zero native plants are occurring
isl <- isl %>% filter(ric > 0)

# prepare exportable raw data frame for Bayesian analysis in Stan and Turing.jl languages
SAR_analysis <- isl %>% right_join(., atoll[,c("atoll",
                                               "annual_precipitation_mm",
                                               "hurricanes_50km",
                                               "isolation",
                                               "ric")], by = "atoll")
# export
write.csv(SAR_analysis, "SAR_analysis.csv", row.names = FALSE)
rm(SAR_analysis) # no longer needed keep tidy

# Add atoll-level environmental information to isl data frame
isl <- isl %>% left_join(., atoll[,c("atoll",
                                     "annual_precipitation_mm",
                                     "isolation",
                                     "hurricanes_50km")], by = "atoll")

# z-standardize and log-transform predictors for downstream modelling
isl$z.area <- (log(isl$area_sqkm) - mean(log(isl$area_sqkm))) / sd(log(isl$area_sqkm))
isl$z.rainfall <- (log(isl$annual_precipitation_mm) - mean(log(isl$annual_precipitation_mm))) / sd(log(isl$annual_precipitation_mm))
isl$z.cyclones <- (log(isl$hurricanes_50km+1) - mean(log(isl$hurricanes_50km+1))) / sd(log(isl$hurricanes_50km+1))



p <- ggplot(dat = isl, aes(x = area_sqkm, y = ric)) +
  geom_point(shape = 21, size = 1.5, stroke = 1.2, alpha = 0.8) +
  geom_smooth(method = "lm", fill = "#E3E3E3", color = "#C96868", linewidth = 1) +
  scale_x_continuous(trans = "log10",
                     name = "Island area [ha]") +
  scale_y_continuous(trans = "log10",
                     name = "Species richness") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  theme_classic() +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 8)) +
  facet_wrap(.~atoll)

#ggsave(p, filename = "atoll-wise_SAR.svg", dpi = 300, width = 159, height = 90, unit = "mm")

##########################################################################################################
### ----- 6. Assess effect of environmental drivers on slope of Species-Area Relationships (SAR) ----- ###
##########################################################################################################

# Define model
mod_SAR.slope <- bf(log(ric) ~ 1 + z.area * (z.rainfall + z.cyclones) + (1 | atoll))

# Run model
SAR.slope.model <- brm(formula = mod_SAR.slope,
                       data = isl,
                       family = gaussian(),
                       chains = 4,
                       warmup = 2500,
                       iter = 5000,
                       save_pars = save_pars(all = TRUE))
# Evaluate model
plot(SAR.slope.model)
pp_check(SAR.slope.model, type = "dens_overlay") + theme_classic()
gelman.diag(as.mcmc(SAR.slope.model)[,c(1:6)], multivariate = FALSE)
geweke.diag(as.mcmc(SAR.slope.model)[,c(1:6)])

# Leave-One-Out (LOO) cross-validation 
loo(SAR.slope.model)

# Explore model output
summary(SAR.slope.model)
bayes_R2(SAR.slope.model)

# Test hypotheses that rainfall and/or cyclones influence the slope of SAR within atolls
hypothesis(SAR.slope.model, "z.area:z.rainfall > 0")
hypothesis(SAR.slope.model, "z.area:z.cyclones < 0")

conditional_effects(SAR.slope.model)

#################################################################################################################
### ----- 7. Assess effect of environmental drivers on variability in  Species-Area Relationships (SAR) ----- ###
#################################################################################################################

### Bayesian Hierarchical Model

# Define model
mod_SAR.sigma <- bf(log(ric) ~ 1 + z.area + (1 | atoll),
                    sigma ~ 1 + z.area + z.rainfall + z.cyclones)

# Run model
SAR.sigma.model <- brm(formula = mod_SAR.sigma,
                     data = isl,
                     family = gaussian(),
                     chains = 4,
                     warmup = 2500,
                     iter = 5000,
                     save_pars = save_pars(all = TRUE))

# Evaluate model
plot(SAR.sigma.model)
pp_check(SAR.sigma.model, type = "dens_overlay") + theme_classic()
gelman.diag(as.mcmc(SAR.sigma.model)[,c(1:9)], multivariate = FALSE)
geweke.diag(as.mcmc(SAR.sigma.model)[,c(1:9)])

# Leave-one-out cross-validation
loo(SAR.sigma.model)

# Explore model output
summary(SAR.sigma.model)
bayes_R2(SAR.sigma.model)

# Test hypotheses on effect of biogeographic and environmental drivers on SAR predictability
hypothesis(SAR.sigma.model, "sigma_z.area < 0")
hypothesis(SAR.sigma.model, "sigma_z.rainfall < 0")
hypothesis(SAR.sigma.model, "sigma_z.cyclones > 0")

# Extract slope estimates
mu.slope.interv.area <- fixef(SAR.sigma.model)["z.area",c("Estimate", "Q2.5", "Q97.5")]
sigma.slope.interv.area <- fixef(SAR.sigma.model)["sigma_z.area",c("Estimate", "Q2.5", "Q97.5")]
sigma.slope.interv.rainfall <- fixef(SAR.sigma.model)["sigma_z.rainfall",c("Estimate", "Q2.5", "Q97.5")]
sigma.slope.interv.cyclones <- fixef(SAR.sigma.model)["sigma_z.cyclones",c("Estimate", "Q2.5", "Q97.5")]

# Obtain conditional effects
c_eff_mu <- conditional_effects(SAR.sigma.model, dpar = "mu")
c_eff_sigma <- conditional_effects(SAR.sigma.model, dpar = "sigma")

# Back-transform predictors (un-z-standardise)
c_mu_area <- c_eff_mu[["z.area"]]
c_mu_area$log_area <- (c_mu_area$z.area * sd(log(isl$area_sqkm)) + mean(log(isl$area_sqkm)))

c_sigma_area <- c_eff_sigma[["z.area"]]
c_sigma_area$log_area <- (c_sigma_area$z.area * sd(log(isl$area_sqkm))) + mean(log(isl$area_sqkm))

c_sigma_rainfall <- c_eff_sigma[["z.rainfall"]]
c_sigma_rainfall$log_rainfall <- (c_sigma_rainfall$z.rainfall * sd(log(isl$annual_precipitation_mm))) + mean(log(isl$annual_precipitation_mm))

c_sigma_cyclones <- c_eff_sigma[["z.cyclones"]]
c_sigma_cyclones$log_cyclones <- (c_sigma_cyclones$z.cyclones * sd(log(isl$hurricanes_50km+1))) + mean(log(isl$hurricanes_50km+1))

# Plot conditional effects
p2.1 <- ggplot(c_mu_area, aes(x = log_area, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#E3E3E3") +
  geom_point(inherit.aes = FALSE,
             data = isl, aes(x = log(area_sqkm), y = log(ric)),
             shape = 21, size = 1.5, stroke = 1.2, alpha = 0.6) +
  geom_line(color = "#C96868", linewidth = 1, alpha = 0.9) +
  xlab("log(island area [ha])") +
  ylab("log(species richness)") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  scale_x_continuous(limits = c(-4.7, 6.6),
                     breaks = seq(-4.5, 6.5, 2.2)) +
  scale_y_continuous(limits = c(0, 4.6),
                     breaks = seq(0, 4.5, 0.5)) +
  ggtitle("A. Island-level SAR of atolls") +
  annotate(geom = "text",
           x = -4.5, y = 4.5,
           label = paste0("β = ", round(mu.slope.interv.area[1], 3), " [", round(mu.slope.interv.area[2],3), "; ", round(mu.slope.interv.area[3], 3), "]",
                          "\nPP(β > 0) = ", hypothesis(SAR.sigma.model, "z.area > 0")$hypothesis$Post.Prob*100, "%"),
           hjust = 0,
           vjust = 1,
           size = 2,
           family = "serif") +
  theme_classic() +
  theme(axis.text = element_text(size = 7, family = "serif"),
        axis.title = element_text(size = 8, family = "serif"),
        plot.title = element_text(size = 9, family = "serif"))  

p2.2 <- ggplot(c_sigma_area, aes(x = log_area, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#E3E3E3") +
  geom_line(color = "#C96868", linewidth = 1) +
  xlab("log(island area [ha])") +
  ylab("σ(SAR)") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  scale_x_continuous(limits = c(-4.7, 6.5),
                     breaks = seq(-4.5, 6.5, 2.2)) +
  scale_y_continuous(limits = c(0.18, 1),
                     breaks = seq(0.25, 1, 0.25)) +
  ggtitle("B. Area effect") +
  annotate(geom = "text",
           x = 6.5, y = 1,
           label = paste0("β = ", round(sigma.slope.interv.area[1], 3), " [", round(sigma.slope.interv.area[2],3), "; ", round(sigma.slope.interv.area[3], 3), "]",
                          "\nPP(β < 0) = ", hypothesis(SAR.sigma.model, "sigma_z.area < 0")$hypothesis$Post.Prob*100, "%"),
           hjust = 1,
           vjust = 1,
           size = 2,
           family = "serif") +
  theme_classic() +
  theme(axis.text = element_text(size = 7, family = "serif"),
        axis.title = element_text(size = 8, family = "serif"),
        plot.title = element_text(size = 9, family = "serif"))

p2.3 <- ggplot(c_sigma_rainfall, aes(x = log_rainfall, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#E3E3E3") +
  geom_line(color = "#C96868", linewidth = 1) +
  xlab("log(rainfall [mm])") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  scale_x_continuous(limits = c(6, 8.2),
                     breaks = seq(6, 8, 0.5)) +
  scale_y_continuous(limits = c(0.18, 1),
                     breaks = seq(0.25, 1, 0.25)) +
  ggtitle("C. Rainfall effect") +
  annotate(geom = "text",
           x = 8, y = 1,
           label = paste0("β = ", round(sigma.slope.interv.rainfall[1], 3), " [", round(sigma.slope.interv.rainfall[2],3), "; ", round(sigma.slope.interv.rainfall[3], 3), "]",
                          "\nPP(β < 0) = ", hypothesis(SAR.sigma.model, "sigma_z.rainfall < 0")$hypothesis$Post.Prob*100, "%"),
           hjust = 1,
           vjust = 1,
           size = 2,
           family = "serif") +
  theme_classic() +
  theme(axis.text = element_text(size = 7, family = "serif"),
        axis.title.x = element_text(size = 8, family = "serif"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 9, family = "serif"))

p2.4 <- ggplot(c_sigma_cyclones, aes(x = log_cyclones, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#E3E3E3") +
  geom_line(color = "#C96868", linewidth = 1) +
  xlab("log(cyclone frequency)") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  scale_x_continuous(limits = c(0, 2.64),
                     breaks = seq(0, 2.5, 0.5)) +
  scale_y_continuous(limits = c(0.18, 1),
                     breaks = seq(0.25, 1, 0.25)) +
  ggtitle("D. Cyclone effect") +
  annotate(geom = "text",
           x = 2.5, y = 1,
           label = paste0("β = ", round(sigma.slope.interv.cyclones[1], 3), " [", round(sigma.slope.interv.cyclones[2],3), "; ", round(sigma.slope.interv.cyclones[3], 3), "]",
                          "\nPP(β > 0) = ", hypothesis(SAR.sigma.model, "sigma_z.cyclones > 0")$hypothesis$Post.Prob*100, "%"),
           hjust = 1,
           vjust = 1,
           size = 2,
           family = "serif") +
  theme_classic() +
  theme(axis.text = element_text(size = 7, family = "serif"),
        axis.title.x = element_text(size = 8, family = "serif"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 9, family = "serif"))

# Stitch plot panels together
p2 <- p2.1 / (p2.2 + p2.3 + p2.4)
p2

# Save and export as SVG vector file
#ggsave(p2, filename = "fig02_island_SAR.svg", dpi = 300, width = 158.5, height = 140, units = "mm")

#########################################################################
### ----- 8. Evaluate atoll-level SAR and environmental drivers ----- ###
#########################################################################

# z-standardize and log-transform predictors
atoll$z.area <- (log(atoll$total_atoll_area_sqkm) - mean(log(atoll$total_atoll_area_sqkm))) / sd(log(atoll$total_atoll_area_sqkm))
atoll$z.isolation <- (log(atoll$isolation) - mean(log(atoll$isolation))) / sd(log(atoll$isolation))
atoll$z.rainfall <- (log(atoll$annual_precipitation_mm) - mean(log(atoll$annual_precipitation_mm))) / sd(log(atoll$annual_precipitation_mm))
atoll$z.cyclones <- (log(atoll$hurricanes_50km+1) - mean(log(atoll$hurricanes_50km+1))) / sd(log(atoll$hurricanes_50km+1))

# define model
mod_atolldiv <- bf(log(ric) ~ 1 + z.area + z.isolation + z.rainfall + z.cyclones)

# Run model
SAR.atoll.model <- brm(formula = mod_atolldiv,
                        data = atoll,
                        family = gaussian(),
                        chains = 4,
                        warmup = 2500,
                        iter = 5000,
                        save_pars = save_pars(all = TRUE))

# Check model
plot(SAR.atoll.model)
pp_check(SAR.atoll.model, type = "dens_overlay") + theme_classic()
gelman.diag(as.mcmc(SAR.atoll.model)[,c(1:7)], multivariate = FALSE)
geweke.diag(as.mcmc(SAR.atoll.model)[,c(1:7)])

# Leave-one-out (LOO) cross-validation
loo(SAR.atoll.model, moment_match = TRUE)

# Explore 
summary(SAR.atoll.model)
bayes_R2(SAR.atoll.model)

# Obtain posterior probabilities for effects of biogeographic and environmental drivers of atoll-level species richness
hypothesis(SAR.atoll.model, "z.area > 0")
hypothesis(SAR.atoll.model, "z.isolation < 0")
hypothesis(SAR.atoll.model, "z.rainfall > 0")
hypothesis(SAR.atoll.model, "z.cyclones > 0")

# Extract conditional effects for plotting
c_eff_atoll <- conditional_effects(SAR.atoll.model)

# Back-transform predictors (un-z-standardise)
c_atoll_area <- c_eff_atoll[["z.area"]]
c_atoll_area$log_area <- (c_atoll_area$z.area * sd(log(atoll$total_atoll_area_sqkm)) + mean(log(atoll$total_atoll_area_sqkm)))

c_atoll_isolation <- c_eff_atoll[["z.isolation"]]
c_atoll_isolation$log_isolation <- (c_atoll_isolation$z.isolation + sd(log(atoll$isolation)) + mean(log(atoll$isolation)))

c_atoll_rainfall <- c_eff_atoll[["z.rainfall"]]
c_atoll_rainfall$log_rainfall <- (c_atoll_rainfall$z.rainfall * sd(log(atoll$annual_precipitation_mm))) + mean(log(atoll$annual_precipitation_mm))

c_atoll_cyclones <- c_eff_atoll[["z.cyclones"]]
c_atoll_cyclones$log_cyclones <- (c_atoll_cyclones$z.cyclones * sd(log(atoll$hurricanes_50km+1))) + mean(log(atoll$hurricanes_50km+1))

# Extract slope estimates
slope.interv.area <- fixef(SAR.atoll.model)["z.area",c("Estimate", "Q2.5", "Q97.5")]
slope.interv.isolation <- fixef(SAR.atoll.model)["z.isolation", c("Estimate", "Q2.5", "Q97.5")]
slope.interv.rainfall <- fixef(SAR.atoll.model)["z.rainfall",c("Estimate", "Q2.5", "Q97.5")]
slope.interv.cyclones <- fixef(SAR.atoll.model)["z.cyclones",c("Estimate", "Q2.5", "Q97.5")]

# Plot correlations
p3.1 <- ggplot(c_atoll_area, aes(x = log_area, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#E3E3E3") +
  geom_line(color = "#C96868", linewidth = 1) +
  geom_point(inherit.aes = FALSE,
             data = atoll, aes(x = log(total_atoll_area_sqkm), y = log(ric)),
             shape = 21, size = 1.5, stroke = 1.2, alpha = 0.6) +
  xlab("log(total atoll land area [ha])") +
  ylab("log(species richness)") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  scale_x_continuous(limits = c(6, 12),
                     breaks = seq(6, 12, 2)) +
  scale_y_continuous(limits = c(1.3, 4.6),
                     breaks = seq(1.5, 4.5, 0.5)) +
  ggtitle("A. Atoll-level SAR") +
  annotate(geom = "text",
           x = 12, y = 1.5,
           label = paste0("β = ", round(slope.interv.area[1], 3), " [", round(slope.interv.area[2],3), "; ", round(slope.interv.area[3], 3), "]",
                          "\nPP(β > 0) = ", hypothesis(SAR.atoll.model, "z.area > 0")$hypothesis$Post.Prob*100, "%"),
           hjust = 1,
           vjust = 0,
           size = 2,
           family = "serif") +
  theme_classic() +
  theme(axis.text = element_text(size = 7, family = "serif"),
        axis.title = element_text(size = 8, family = "serif"),
        plot.title = element_text(size = 9, family = "serif"))

p3.2 <- ggplot(c_atoll_isolation, aes(x = log_isolation, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#E3E3E3") +
  geom_line(color = "#C96868", linewidth = 1) +
  geom_point(inherit.aes = FALSE,
             data = atoll, aes(x = log(isolation), y = log(ric)),
             shape = 21, size = 1.5, stroke = 1.2, alpha = 0.6) +
  xlab("log(distance nearest landmass [km])") +
  ylab("log(species richness)") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  scale_x_continuous(limits = c(4, 8.5),
                     breaks = seq(4, 8.5, 0.5)) +
  scale_y_continuous(limits = c(1.3, 4.6),
                     breaks = seq(1.5, 4.5, 0.5)) +
  ggtitle("B. Isolation effect") +
  annotate(geom = "text",
           x = 8, y = 1.5,
           label = paste0("β = ", round(slope.interv.isolation[1], 3), " [", round(slope.interv.isolation[2],3), "; ", round(slope.interv.isolation[3], 3), "]",
                          "\nPP(β < 0) = ", hypothesis(SAR.atoll.model, "z.isolation < 0")$hypothesis$Post.Prob*100, "%"),
           hjust = 1,
           vjust = 0,
           size = 2,
           family = "serif") +
  theme_classic() +
  theme(axis.text = element_text(size = 7, family = "serif"),
        axis.title = element_text(size = 8, family = "serif"),
        plot.title = element_text(size = 9, family = "serif"))

p3.3 <- ggplot(c_atoll_rainfall, aes(x = log_rainfall, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#E3E3E3") +
  geom_line(color = "#C96868", linewidth = 1) +
  geom_point(inherit.aes = FALSE,
             data = atoll, aes(x = log(annual_precipitation_mm), y = log(ric)),
             shape = 21, size = 1.5, stroke = 1.2, alpha = 0.6) +
  xlab("log(annual rainfall [mm])") +
  ylab("log(species richness)") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  scale_x_continuous(limits = c(6, 8.2),
                     breaks = seq(6, 8, 0.5)) +
  scale_y_continuous(limits = c(1.3, 4.6),
                     breaks = seq(1.5, 4.5, 0.5)) +
  ggtitle("B. Rainfall effect") +
  annotate(geom = "text",
           x = 8, y = 1.5,
           label = paste0("β = ", round(slope.interv.rainfall[1], 3), " [", round(slope.interv.rainfall[2],3), "; ", round(slope.interv.rainfall[3], 3), "]",
                          "\nPP(β > 0) = ", hypothesis(SAR.atoll.model, "z.rainfall > 0")$hypothesis$Post.Prob*100, "%"),
           hjust = 1,
           vjust = 0,
           size = 2,
           family = "serif") +
  theme_classic() +
  theme(axis.text = element_text(size = 7, family = "serif"),
        axis.title.x = element_text(size = 8, family = "serif"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 9, family = "serif"))

p3.4 <- ggplot(c_atoll_cyclones, aes(x = log_cyclones, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#E3E3E3") +
  geom_line(color = "#C96868", linewidth = 1) +
  geom_point(inherit.aes = FALSE,
             data = atoll, aes(x = log(hurricanes_50km+1), y = log(ric)),
             shape = 21, size = 1.5, stroke = 1.2, alpha = 0.6, position = position_jitter(width = 0.01)) +
  xlab("log(cyclone frequency)") +
  ylab("log(species richness)") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  scale_x_continuous(limits = c(-0.1, 2.7),
                     breaks = seq(0, 2.5, 0.5)) +
  scale_y_continuous(limits = c(1.3, 4.6),
                     breaks = seq(1.5, 4.5, 0.5)) +
  ggtitle("C. Cyclone effect") +
  annotate(geom = "text",
           x = 2.5, y = 1.5,
           label = paste0("β = ", round(slope.interv.cyclones[1], 3), " [", round(slope.interv.cyclones[2],3), "; ", round(slope.interv.cyclones[3], 3), "]",
                          "\nPP(β > 0) = ", hypothesis(SAR.atoll.model, "z.cyclones > 0")$hypothesis$Post.Prob*100, "%"),
           hjust = 1,
           vjust = 0,
           size = 2,
           family = "serif") +
  theme_classic() +
  theme(axis.text = element_text(size = 7, family = "serif"),
        axis.title.x = element_text(size = 8, family = "serif"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 9, family = "serif"))

# Stitch plot panels together
p3 <- p3.1 + p3.3 + p3.4
p3

# Export as SVG vector file
#ggsave(p3, filename = "fig03_atoll_SAR.svg", dpi = 300, width = 158.5, height = 70, units = "mm")
