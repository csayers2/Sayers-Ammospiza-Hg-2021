
# This script is very sensitive to the amount of memory space, so let's start fresh here
rm(list = ls())

# To test for correlations between tidal marsh sparrow THg concentrations and landscape characteristics surrounding
# sampling points, we attempted to define a statistically relevant spatial extent, or “scale of effect,” for the
# independent variables. To accomplish this, we compared linear mixed models of identical structure, but with
# spatial attributes extracted within different buffer radii (Jackson and Fahrig 2012).

#------------------------------------------ LOADING THE DATA --------------------------------------------------------
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(raster)
library(sp)
library(AICcmodavg)
library(lmerTest)
library(MuMIn)
library(lme4)
library(usdm)
library(ggResidpanel)

HgSamples <- read_excel("Sayers_SHARP_Hg_Samples.xlsx")

# Creating new column with natural-log transformed THg concentrations + date
HgSamples <- HgSamples %>%
  mutate(lHgConcentration = log(HgConcentration)) %>%
  mutate(Date = make_date(Year, Month, Day),
         JulianDate = yday(Date))

#------------------------------------------- MANIPULATING RASTER DATA ------------------------------------------

# Making a spatial data frame from coordinates
coords <- cbind(HgSamples$Longitude, HgSamples$Latitude)
points_spdf <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
plot(points_spdf)

# Importing contiguous USA NLCD raster for 2016 (this is not available yet using FedData::get_nlcd)
# NLCD 2016 Land Cover (CONUS) - https://www.mrlc.gov/data
# It's very imoportant that the .img stays within the original NLCD folder - I think it communicates with metadata, etc.
# We cannot put the raster within the GitHub repository because the files are too large
# You can download the raster separately, save it to your desktop, and call it independently
# file.choose() # this function helps you find the exact file pathway
nlcd_raster <- raster("/Users/chrissayers/Documents/Senior Thesis Stuff/NLCD_2016_Land_Cover_L48_20190424/NLCD_2016_Land_Cover_L48_20190424.img")
crs(nlcd_raster) # Albers equal area projection, WGS84 datum, we will need to reproject spdf to match this
plot(nlcd_raster)

# Reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# Now there is agreement between raster crs and coordinates crs

# Creating a spatial extent to crop the NLCD raster and to make sure all points can have an X km buffer
crop_extent <- raster(ncol = 3, nrow = 3)
extent(crop_extent) <- extent(points_spdf) * 1.5
crs(crop_extent) <- crs(points_spdf)

# Cropping raster so that it is easier to extract from -- takes a bit of time
nlcd_raster <- crop(nlcd_raster, crop_extent)

# Plotting to make sure it is all there
plot(nlcd_raster)
points(points_spdf)

#--------- Mercury Deposition Network Wet Deposition --------------------------
# Read in the MDN wet deposition raster info (units in µg/m2)
library(raster)
library(sp)

wetdep <- raster("Hg_dep_2017/Hg_dep_2017.tif")
res(wetdep) # raster has a resolution of 2338x2338

crs(wetdep) # Albers equal area projection, NAD83 datum, we will need to reproject spdf to match this
plot(wetdep)

# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

#------------------------------------------- 5KM BUFFER MODEL ------------------------------------------
# reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# now there is agreement between raster crs and coordinates crs

# Extracting NLCD values
lc <- raster::extract(x = nlcd_raster, y = points_spdf, buffer = 5000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc <- as.data.frame(lc) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc)

# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.lc) <- c("ID", "Open Water", "Developed, Open Space",
                       "Developed, Low Intensity", "Developed, Medium Intensity",
                       "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                       "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                       "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                       "Emergent Herbaceous Wetlands")
head(prop.lc)

# lumping landcover types
Hglandcover <- data.frame(cbind(Water = rowSums(prop.lc[,c("Open Water")]),
                                Developed = rowSums(prop.lc[,c("Developed, Open Space", 
                                                               "Developed, Low Intensity",
                                                               "Developed, Medium Intensity", 
                                                               "Developed, High Intensity")]),
                                Barren = prop.lc$"Barren Land",
                                Forest = rowSums(prop.lc[,c("Deciduous Forest",
                                                            "Evergreen Forest",
                                                            "Mixed Forest", "Shrub/Scrub")]),
                                Cultivated = rowSums(prop.lc[,c("Herbaceous", "Pasture/Hay",
                                                                "Cultivated Crops")]), 
                                Wetlands = rowSums(prop.lc[,c("Woody Wetlands",
                                                              "Emergent Herbaceous Wetlands")])))

# just checking to see if all landcover categories sum to 1 -- THEY DO!
prop.lc <- prop.lc %>%
  mutate(Sum = rowSums(prop.lc[,c("Open Water", "Developed, Open Space", "Developed, Low Intensity",
                                  "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
                                  "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                                  "Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                  "Emergent Herbaceous Wetlands")]))
head(prop.lc$Sum)

#--------- Mercury Deposition Network Wet Deposition --------------------------
# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 5000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within Xkm of each marsh centroid
avg.wetdep <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T))

# Making a single data frame for the model
HgData5km <- data.frame(cbind(HgSamples, Hglandcover, avg.wetdep))

# This is going to help extraction speeds increase
rm(lc, prop.lc, Hglandcover, wetdep.ex, avg.wetdep)

# Making the models ------------------------------
landmod5 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
                 data = HgData5km, REML = F)
depmod5 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID),
                data = HgData5km, REML = F)
bothmod5 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest)
                 + scale(AvgWetDep) + (1 | SiteID), data = HgData5km, REML = F)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

# Checking model assumptions -------------------------------------
# Checking for multicollinearity
A = data.frame(HgData5km$Developed, HgData5km$Forest, HgData5km$AvgWetDep)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3) # none have collinearity problems, VIF < 3

# Checking for homogeneity of variance & normality of residuals
resid_panel(landmod5, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(depmod5, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(bothmod5, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(landmod5)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(landmod5)$SiteID[,1]) # we are normal

ggqqplot(ranef(depmod5)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(depmod5)$SiteID[,1]) # we are normal

ggqqplot(ranef(bothmod5)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(bothmod5)$SiteID[,1]) # we are normal


#------------------------------------------- 6KM BUFFER MODEL ------------------------------------------
# reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# now there is agreement between raster crs and coordinates crs

# Extracting NLCD values
lc <- raster::extract(x = nlcd_raster, y = points_spdf, buffer = 6000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc <- as.data.frame(lc) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc)

# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.lc) <- c("ID", "Unclassified", "Open Water", "Developed, Open Space",
                       "Developed, Low Intensity", "Developed, Medium Intensity",
                       "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                       "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                       "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                       "Emergent Herbaceous Wetlands")
head(prop.lc)

# lumping landcover types
Hglandcover <- data.frame(cbind(Water = rowSums(prop.lc[,c("Unclassified", "Open Water")]),
                                Developed = rowSums(prop.lc[,c("Developed, Open Space", 
                                                               "Developed, Low Intensity",
                                                               "Developed, Medium Intensity", 
                                                               "Developed, High Intensity")]),
                                Barren = prop.lc$"Barren Land",
                                Forest = rowSums(prop.lc[,c("Deciduous Forest",
                                                            "Evergreen Forest",
                                                            "Mixed Forest", "Shrub/Scrub")]),
                                Cultivated = rowSums(prop.lc[,c("Herbaceous", "Pasture/Hay",
                                                                "Cultivated Crops")]), 
                                Wetlands = rowSums(prop.lc[,c("Woody Wetlands",
                                                              "Emergent Herbaceous Wetlands")])))

# just checking to see if all landcover categories sum to 1 -- THEY DO!
prop.lc <- prop.lc %>%
  mutate(Sum = rowSums(prop.lc[,c("Unclassified","Open Water", "Developed, Open Space", "Developed, Low Intensity",
                                  "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
                                  "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                                  "Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                  "Emergent Herbaceous Wetlands")]))
head(prop.lc$Sum)

#--------- Mercury Deposition Network Wet Deposition --------------------------
# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 6000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within Xkm of each marsh centroid
avg.wetdep <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T))

# Making a single data frame for the model
HgData6km <- data.frame(cbind(HgSamples, Hglandcover, avg.wetdep))

# This is going to help extraction speeds increase
rm(lc, prop.lc, Hglandcover, wetdep.ex, avg.wetdep)

# Making the models ------------------------------
landmod6 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
                 data = HgData6km, REML = F)
depmod6 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID),
                data = HgData6km, REML = F)
bothmod6 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest)
                 + scale(AvgWetDep) + (1 | SiteID), data = HgData6km, REML = F)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

# Checking model assumptions -------------------------------------
# Checking for multicollinearity
A = data.frame(HgData6km$Developed, HgData6km$Forest, HgData6km$AvgWetDep)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3) # none have collinearity problems, VIF < 3

# Checking for homogeneity of variance & normality of residuals
resid_panel(landmod6, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(depmod6, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(bothmod6, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(landmod6)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(landmod6)$SiteID[,1]) # we are normal

ggqqplot(ranef(depmod6)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(depmod6)$SiteID[,1]) # we are normal

ggqqplot(ranef(bothmod6)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(bothmod6)$SiteID[,1]) # we are normal


#------------------------------------------- 7KM BUFFER MODEL ------------------------------------------
# reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# now there is agreement between raster crs and coordinates crs

# Extracting NLCD values
lc <- raster::extract(x = nlcd_raster, y = points_spdf, buffer = 7000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc <- as.data.frame(lc) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc)

# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.lc) <- c("ID", "Unclassified", "Open Water", "Developed, Open Space",
                       "Developed, Low Intensity", "Developed, Medium Intensity",
                       "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                       "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                       "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                       "Emergent Herbaceous Wetlands")
head(prop.lc)

# lumping landcover types
Hglandcover <- data.frame(cbind(Water = rowSums(prop.lc[,c("Unclassified", "Open Water")]),
                                Developed = rowSums(prop.lc[,c("Developed, Open Space", 
                                                               "Developed, Low Intensity",
                                                               "Developed, Medium Intensity", 
                                                               "Developed, High Intensity")]),
                                Barren = prop.lc$"Barren Land",
                                Forest = rowSums(prop.lc[,c("Deciduous Forest",
                                                            "Evergreen Forest",
                                                            "Mixed Forest", "Shrub/Scrub")]),
                                Cultivated = rowSums(prop.lc[,c("Herbaceous", "Pasture/Hay",
                                                                "Cultivated Crops")]), 
                                Wetlands = rowSums(prop.lc[,c("Woody Wetlands",
                                                              "Emergent Herbaceous Wetlands")])))

# just checking to see if all landcover categories sum to 1 -- THEY DO!
prop.lc <- prop.lc %>%
  mutate(Sum = rowSums(prop.lc[,c("Unclassified","Open Water", "Developed, Open Space", "Developed, Low Intensity",
                                  "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
                                  "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                                  "Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                  "Emergent Herbaceous Wetlands")]))
head(prop.lc$Sum)

#--------- Mercury Deposition Network Wet Deposition --------------------------
# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 7000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within Xkm of each marsh centroid
avg.wetdep <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T))

# Making a single data frame for the model
HgData7km <- data.frame(cbind(HgSamples, Hglandcover, avg.wetdep))

# This is going to help extraction speeds increase
rm(lc, prop.lc, Hglandcover, wetdep.ex, avg.wetdep)

# Making the models ------------------------------
landmod7 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
                 data = HgData7km, REML = F)
depmod7 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID),
                data = HgData7km, REML = F)
bothmod7 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest)
                 + scale(AvgWetDep) + (1 | SiteID), data = HgData7km, REML = F)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

# Checking model assumptions -------------------------------------
# Checking for multicollinearity
A = data.frame(HgData7km$Developed, HgData7km$Forest, HgData7km$AvgWetDep)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3) # none have collinearity problems, VIF < 3

# Checking for homogeneity of variance & normality of residuals
resid_panel(landmod7, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(depmod7, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(bothmod7, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(landmod7)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(landmod7)$SiteID[,1]) # we are normal

ggqqplot(ranef(depmod7)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(depmod7)$SiteID[,1]) # we are normal

ggqqplot(ranef(bothmod7)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(bothmod7)$SiteID[,1]) # we are normal


#------------------------------------------- 8KM BUFFER MODEL ------------------------------------------
# reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# now there is agreement between raster crs and coordinates crs

# Extracting NLCD values
lc <- raster::extract(x = nlcd_raster, y = points_spdf, buffer = 8000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc <- as.data.frame(lc) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc)

# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.lc) <- c("ID", "Unclassified", "Open Water", "Developed, Open Space",
                       "Developed, Low Intensity", "Developed, Medium Intensity",
                       "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                       "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                       "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                       "Emergent Herbaceous Wetlands")
head(prop.lc)

# lumping landcover types
Hglandcover <- data.frame(cbind(Water = rowSums(prop.lc[,c("Unclassified", "Open Water")]),
                                Developed = rowSums(prop.lc[,c("Developed, Open Space", 
                                                               "Developed, Low Intensity",
                                                               "Developed, Medium Intensity", 
                                                               "Developed, High Intensity")]),
                                Barren = prop.lc$"Barren Land",
                                Forest = rowSums(prop.lc[,c("Deciduous Forest",
                                                            "Evergreen Forest",
                                                            "Mixed Forest", "Shrub/Scrub")]),
                                Cultivated = rowSums(prop.lc[,c("Herbaceous", "Pasture/Hay",
                                                                "Cultivated Crops")]), 
                                Wetlands = rowSums(prop.lc[,c("Woody Wetlands",
                                                              "Emergent Herbaceous Wetlands")])))

# just checking to see if all landcover categories sum to 1 -- THEY DO!
prop.lc <- prop.lc %>%
  mutate(Sum = rowSums(prop.lc[,c("Unclassified","Open Water", "Developed, Open Space", "Developed, Low Intensity",
                                  "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
                                  "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                                  "Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                  "Emergent Herbaceous Wetlands")]))
head(prop.lc$Sum)

#--------- Mercury Deposition Network Wet Deposition --------------------------
# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 8000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within Xkm of each marsh centroid
avg.wetdep <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T))

# Making a single data frame for the model
HgData8km <- data.frame(cbind(HgSamples, Hglandcover, avg.wetdep))

# This is going to help extraction speeds increase
rm(lc, prop.lc, Hglandcover, wetdep.ex, avg.wetdep)

# Making the models ------------------------------
landmod8 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
                 data = HgData8km, REML = F)
depmod8 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID),
                data = HgData8km, REML = F)
bothmod8 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest)
                 + scale(AvgWetDep) + (1 | SiteID), data = HgData8km, REML = F)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

# Checking model assumptions -------------------------------------
# Checking for multicollinearity
A = data.frame(HgData8km$Developed, HgData8km$Forest, HgData8km$AvgWetDep)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3) # none have collinearity problems, VIF < 3

# Checking for homogeneity of variance & normality of residuals
resid_panel(landmod8, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(depmod8, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(bothmod8, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(landmod8)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(landmod8)$SiteID[,1]) # we are normal

ggqqplot(ranef(depmod8)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(depmod8)$SiteID[,1]) # we are normal

ggqqplot(ranef(bothmod8)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(bothmod8)$SiteID[,1]) # we are normal


#------------------------------------------- 9KM BUFFER MODEL ------------------------------------------
# reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# now there is agreement between raster crs and coordinates crs

# Extracting NLCD values
lc <- raster::extract(x = nlcd_raster, y = points_spdf, buffer = 9000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc <- as.data.frame(lc) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc)

# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.lc) <- c("ID", "Unclassified", "Open Water", "Developed, Open Space",
                       "Developed, Low Intensity", "Developed, Medium Intensity",
                       "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                       "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                       "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                       "Emergent Herbaceous Wetlands")
head(prop.lc)

# lumping landcover types
Hglandcover <- data.frame(cbind(Water = rowSums(prop.lc[,c("Unclassified", "Open Water")]),
                                Developed = rowSums(prop.lc[,c("Developed, Open Space", 
                                                               "Developed, Low Intensity",
                                                               "Developed, Medium Intensity", 
                                                               "Developed, High Intensity")]),
                                Barren = prop.lc$"Barren Land",
                                Forest = rowSums(prop.lc[,c("Deciduous Forest",
                                                            "Evergreen Forest",
                                                            "Mixed Forest", "Shrub/Scrub")]),
                                Cultivated = rowSums(prop.lc[,c("Herbaceous", "Pasture/Hay",
                                                                "Cultivated Crops")]), 
                                Wetlands = rowSums(prop.lc[,c("Woody Wetlands",
                                                              "Emergent Herbaceous Wetlands")])))

# just checking to see if all landcover categories sum to 1 -- THEY DO!
prop.lc <- prop.lc %>%
  mutate(Sum = rowSums(prop.lc[,c("Unclassified","Open Water", "Developed, Open Space", "Developed, Low Intensity",
                                  "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
                                  "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                                  "Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                  "Emergent Herbaceous Wetlands")]))
head(prop.lc$Sum)

#--------- Mercury Deposition Network Wet Deposition --------------------------
# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 9000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within Xkm of each marsh centroid
avg.wetdep <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T))

# Making a single data frame for the model
HgData9km <- data.frame(cbind(HgSamples, Hglandcover, avg.wetdep))

# This is going to help extraction speeds increase
rm(lc, prop.lc, Hglandcover, wetdep.ex, avg.wetdep)

# Making the models ------------------------------
landmod9 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
                 data = HgData9km, REML = F)
depmod9 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID),
                data = HgData9km, REML = F)
bothmod9 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest)
                 + scale(AvgWetDep) + (1 | SiteID), data = HgData9km, REML = F)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

# Checking model assumptions -------------------------------------
# Checking for multicollinearity
A = data.frame(HgData9km$Developed, HgData9km$Forest, HgData9km$AvgWetDep)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3) # none have collinearity problems, VIF < 3

# Checking for homogeneity of variance & normality of residuals
resid_panel(landmod9, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(depmod9, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(bothmod9, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(landmod9)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(landmod9)$SiteID[,1]) # we are normal

ggqqplot(ranef(depmod9)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(depmod9)$SiteID[,1]) # we are normal

ggqqplot(ranef(bothmod9)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(bothmod9)$SiteID[,1]) # we are normal


#------------------------------------------- 10KM BUFFER MODEL ------------------------------------------
# reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# now there is agreement between raster crs and coordinates crs

# Extracting NLCD values
lc <- raster::extract(x = nlcd_raster, y = points_spdf, buffer = 10000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc <- as.data.frame(lc) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc)

# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.lc) <- c("ID", "Unclassified", "Open Water", "Developed, Open Space",
                       "Developed, Low Intensity", "Developed, Medium Intensity",
                       "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                       "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                       "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                       "Emergent Herbaceous Wetlands")
head(prop.lc)

# lumping landcover types
Hglandcover <- data.frame(cbind(Water = rowSums(prop.lc[,c("Unclassified", "Open Water")]),
                                Developed = rowSums(prop.lc[,c("Developed, Open Space", 
                                                               "Developed, Low Intensity",
                                                               "Developed, Medium Intensity", 
                                                               "Developed, High Intensity")]),
                                Barren = prop.lc$"Barren Land",
                                Forest = rowSums(prop.lc[,c("Deciduous Forest",
                                                            "Evergreen Forest",
                                                            "Mixed Forest", "Shrub/Scrub")]),
                                Cultivated = rowSums(prop.lc[,c("Herbaceous", "Pasture/Hay",
                                                                "Cultivated Crops")]), 
                                Wetlands = rowSums(prop.lc[,c("Woody Wetlands",
                                                              "Emergent Herbaceous Wetlands")])))

# just checking to see if all landcover categories sum to 1 -- THEY DO!
prop.lc <- prop.lc %>%
  mutate(Sum = rowSums(prop.lc[,c("Unclassified","Open Water", "Developed, Open Space", "Developed, Low Intensity",
                                  "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
                                  "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                                  "Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                  "Emergent Herbaceous Wetlands")]))
head(prop.lc$Sum)

#--------- Mercury Deposition Network Wet Deposition --------------------------
# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 10000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within Xkm of each marsh centroid
avg.wetdep <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T))

# Making a single data frame for the model
HgData10km <- data.frame(cbind(HgSamples, Hglandcover, avg.wetdep))

# This is going to help extraction speeds increase
rm(lc, prop.lc, Hglandcover, wetdep.ex, avg.wetdep)

# Making the models ------------------------------
landmod10 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
                 data = HgData10km, REML = F)
depmod10 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID),
                data = HgData10km, REML = F)
bothmod10 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest)
                 + scale(AvgWetDep) + (1 | SiteID), data = HgData10km, REML = F)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

# Checking model assumptions -------------------------------------
# Checking for multicollinearity
A = data.frame(HgData10km$Developed, HgData10km$Forest, HgData10km$AvgWetDep)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3) # none have collinearity problems, VIF < 3

# Checking for homogeneity of variance & normality of residuals
resid_panel(landmod10, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(depmod10, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(bothmod10, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(landmod10)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(landmod10)$SiteID[,1]) # we are normal

ggqqplot(ranef(depmod10)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(depmod10)$SiteID[,1]) # we are normal

ggqqplot(ranef(bothmod10)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(bothmod10)$SiteID[,1]) # we are normal


#------------------------------------------- 15KM BUFFER MODEL ------------------------------------------
# reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# now there is agreement between raster crs and coordinates crs

# Extracting NLCD values
lc <- raster::extract(x = nlcd_raster, y = points_spdf, buffer = 15000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc <- as.data.frame(lc) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc)

# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.lc) <- c("ID", "Unclassified", "Open Water", "Developed, Open Space",
                       "Developed, Low Intensity", "Developed, Medium Intensity",
                       "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                       "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                       "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                       "Emergent Herbaceous Wetlands")
head(prop.lc)

# lumping landcover types
Hglandcover <- data.frame(cbind(Water = rowSums(prop.lc[,c("Unclassified", "Open Water")]),
                                Developed = rowSums(prop.lc[,c("Developed, Open Space", 
                                                               "Developed, Low Intensity",
                                                               "Developed, Medium Intensity", 
                                                               "Developed, High Intensity")]),
                                Barren = prop.lc$"Barren Land",
                                Forest = rowSums(prop.lc[,c("Deciduous Forest",
                                                            "Evergreen Forest",
                                                            "Mixed Forest", "Shrub/Scrub")]),
                                Cultivated = rowSums(prop.lc[,c("Herbaceous", "Pasture/Hay",
                                                                "Cultivated Crops")]), 
                                Wetlands = rowSums(prop.lc[,c("Woody Wetlands",
                                                              "Emergent Herbaceous Wetlands")])))

# just checking to see if all landcover categories sum to 1 -- THEY DO!
prop.lc <- prop.lc %>%
  mutate(Sum = rowSums(prop.lc[,c("Unclassified","Open Water", "Developed, Open Space", "Developed, Low Intensity",
                                  "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
                                  "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                                  "Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                  "Emergent Herbaceous Wetlands")]))
head(prop.lc$Sum)

#--------- Mercury Deposition Network Wet Deposition --------------------------
# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 15000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within Xkm of each marsh centroid
avg.wetdep <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T))

# Making a single data frame for the model
HgData15km <- data.frame(cbind(HgSamples, Hglandcover, avg.wetdep))

# This is going to help extraction speeds increase
rm(lc, prop.lc, Hglandcover, wetdep.ex, avg.wetdep)

# Making the models ------------------------------
landmod15 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
                 data = HgData15km, REML = F)
depmod15 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID),
                data = HgData15km, REML = F)
bothmod15 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest)
                 + scale(AvgWetDep) + (1 | SiteID), data = HgData15km, REML = F)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

# Checking model assumptions -------------------------------------
# Checking for multicollinearity
A = data.frame(HgData15km$Developed, HgData15km$Forest, HgData15km$AvgWetDep)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3) # none have collinearity problems, VIF < 3

# Checking for homogeneity of variance & normality of residuals
resid_panel(landmod15, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(depmod15, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(bothmod15, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(landmod15)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(landmod15)$SiteID[,1]) # we are normal

ggqqplot(ranef(depmod15)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(depmod15)$SiteID[,1]) # we are normal

ggqqplot(ranef(bothmod15)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(bothmod15)$SiteID[,1]) # we are normal


#------------------------------------------- 20KM BUFFER MODEL ------------------------------------------
# reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# now there is agreement between raster crs and coordinates crs

# Extracting NLCD values
lc <- raster::extract(x = nlcd_raster, y = points_spdf, buffer = 20000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc <- as.data.frame(lc) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc)

# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.lc) <- c("ID", "Unclassified", "Open Water", "Developed, Open Space",
                       "Developed, Low Intensity", "Developed, Medium Intensity",
                       "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                       "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                       "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                       "Emergent Herbaceous Wetlands")
head(prop.lc)

# lumping landcover types
Hglandcover <- data.frame(cbind(Water = rowSums(prop.lc[,c("Unclassified", "Open Water")]),
                                Developed = rowSums(prop.lc[,c("Developed, Open Space", 
                                                               "Developed, Low Intensity",
                                                               "Developed, Medium Intensity", 
                                                               "Developed, High Intensity")]),
                                Barren = prop.lc$"Barren Land",
                                Forest = rowSums(prop.lc[,c("Deciduous Forest",
                                                            "Evergreen Forest",
                                                            "Mixed Forest", "Shrub/Scrub")]),
                                Cultivated = rowSums(prop.lc[,c("Herbaceous", "Pasture/Hay",
                                                                "Cultivated Crops")]), 
                                Wetlands = rowSums(prop.lc[,c("Woody Wetlands",
                                                              "Emergent Herbaceous Wetlands")])))

# just checking to see if all landcover categories sum to 1 -- THEY DO!
prop.lc <- prop.lc %>%
  mutate(Sum = rowSums(prop.lc[,c("Unclassified","Open Water", "Developed, Open Space", "Developed, Low Intensity",
                                  "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
                                  "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                                  "Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                  "Emergent Herbaceous Wetlands")]))
head(prop.lc$Sum)

#--------- Mercury Deposition Network Wet Deposition --------------------------
# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 20000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within Xkm of each marsh centroid
avg.wetdep <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T))

# Making a single data frame for the model
HgData20km <- data.frame(cbind(HgSamples, Hglandcover, avg.wetdep))

# This is going to help extraction speeds increase
rm(lc, prop.lc, Hglandcover, wetdep.ex, avg.wetdep)

# Making the models ------------------------------
landmod20 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
                 data = HgData20km, REML = F)
depmod20 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID),
                data = HgData20km, REML = F)
bothmod20 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest)
                 + scale(AvgWetDep) + (1 | SiteID), data = HgData20km, REML = F)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

# Checking model assumptions -------------------------------------
# Checking for multicollinearity
A = data.frame(HgData20km$Developed, HgData20km$Forest, HgData20km$AvgWetDep)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3) # none have collinearity problems, VIF < 3

# Checking for homogeneity of variance & normality of residuals
resid_panel(landmod20, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(depmod20, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(bothmod20, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(landmod20)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(landmod20)$SiteID[,1]) # we are normal

ggqqplot(ranef(depmod20)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(depmod20)$SiteID[,1]) # we are normal

ggqqplot(ranef(bothmod20)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(bothmod20)$SiteID[,1]) # we are normal


#------------------------------------------- 25KM BUFFER MODEL ------------------------------------------
# reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# now there is agreement between raster crs and coordinates crs

# Extracting NLCD values
lc <- raster::extract(x = nlcd_raster, y = points_spdf, buffer = 25000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc <- as.data.frame(lc) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc)

# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.lc) <- c("ID", "Unclassified", "Open Water", "Developed, Open Space",
                       "Developed, Low Intensity", "Developed, Medium Intensity",
                       "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                       "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                       "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                       "Emergent Herbaceous Wetlands")
head(prop.lc)

# lumping landcover types
Hglandcover <- data.frame(cbind(Water = rowSums(prop.lc[,c("Unclassified", "Open Water")]),
                                Developed = rowSums(prop.lc[,c("Developed, Open Space", 
                                                               "Developed, Low Intensity",
                                                               "Developed, Medium Intensity", 
                                                               "Developed, High Intensity")]),
                                Barren = prop.lc$"Barren Land",
                                Forest = rowSums(prop.lc[,c("Deciduous Forest",
                                                            "Evergreen Forest",
                                                            "Mixed Forest", "Shrub/Scrub")]),
                                Cultivated = rowSums(prop.lc[,c("Herbaceous", "Pasture/Hay",
                                                                "Cultivated Crops")]), 
                                Wetlands = rowSums(prop.lc[,c("Woody Wetlands",
                                                              "Emergent Herbaceous Wetlands")])))

# just checking to see if all landcover categories sum to 1 -- THEY DO!
prop.lc <- prop.lc %>%
  mutate(Sum = rowSums(prop.lc[,c("Unclassified","Open Water", "Developed, Open Space", "Developed, Low Intensity",
                                  "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
                                  "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                                  "Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                  "Emergent Herbaceous Wetlands")]))
head(prop.lc$Sum)

#--------- Mercury Deposition Network Wet Deposition --------------------------
# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 25000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within Xkm of each marsh centroid
avg.wetdep <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T))

# Making a single data frame for the model
HgData25km <- data.frame(cbind(HgSamples, Hglandcover, avg.wetdep))

# This is going to help extraction speeds increase
rm(lc, prop.lc, Hglandcover, wetdep.ex, avg.wetdep)

# Making the models ------------------------------
landmod25 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
                 data = HgData25km, REML = F)
depmod25 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID),
                data = HgData25km, REML = F)
bothmod25 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest)
                 + scale(AvgWetDep) + (1 | SiteID), data = HgData25km, REML = F)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

# Checking model assumptions -------------------------------------
# Checking for multicollinearity
A = data.frame(HgData25km$Developed, HgData25km$Forest, HgData25km$AvgWetDep)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3) # none have collinearity problems, VIF < 3

# Checking for homogeneity of variance & normality of residuals
resid_panel(landmod25, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(depmod25, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(bothmod25, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(landmod25)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(landmod25)$SiteID[,1]) # we are normal

ggqqplot(ranef(depmod25)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(depmod25)$SiteID[,1]) # we are normal

ggqqplot(ranef(bothmod25)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(bothmod25)$SiteID[,1]) # we are normal


#------------------------------------------- 30KM BUFFER MODEL ------------------------------------------
library(sp)
library(raster)
# reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# now there is agreement between raster crs and coordinates crs

# splitting data set so extraction can run
points_spdf1 <- points_spdf[1:20,]
points_spdf2 <- points_spdf[21:40,]
points_spdf3 <- points_spdf[41:60,]
points_spdf4 <- points_spdf[61:80,]
points_spdf5 <- points_spdf[81:100,]
points_spdf6 <- points_spdf[101:120,]

# Extracting NLCD values and turning into dataframes
lc1 <- raster::extract(x = nlcd_raster, y = points_spdf1, buffer = 30000, df = T)
lc2 <- raster::extract(x = nlcd_raster, y = points_spdf2, buffer = 30000, df = T)
lc3 <- raster::extract(x = nlcd_raster, y = points_spdf3, buffer = 30000, df = T)
lc4 <- raster::extract(x = nlcd_raster, y = points_spdf4, buffer = 30000, df = T)
lc5 <- raster::extract(x = nlcd_raster, y = points_spdf5, buffer = 30000, df = T)
lc6 <- raster::extract(x = nlcd_raster, y = points_spdf6, buffer = 30000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc1 <- as.data.frame(lc1) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc1)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc2 <- as.data.frame(lc2) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc2)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc3 <- as.data.frame(lc3) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc3)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc4 <- as.data.frame(lc4) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc4)


# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc5 <- as.data.frame(lc5) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc5)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc6 <- as.data.frame(lc6) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc6)

# joining data frames together
prop.lc  <- prop.lc1 %>% 
  full_join(prop.lc2,  by.y = "ID") %>% 
  full_join(prop.lc3,  by.y = "ID") %>%
  full_join(prop.lc4,  by.y = "ID") %>% 
  full_join(prop.lc5,  by.y = "ID") %>% 
  full_join(prop.lc6,  by.y = "ID")
# renaming ID's
prop.lc[1:120, "ID"] <- 1:120

# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.lc) <- c("ID", "Unclassified", "Open Water", "Developed, Open Space",
                       "Developed, Low Intensity", "Developed, Medium Intensity",
                       "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                       "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                       "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                       "Emergent Herbaceous Wetlands")
head(prop.lc)

# lumping landcover types
Hglandcover <- data.frame(cbind(Water = rowSums(prop.lc[,c("Unclassified", "Open Water")]),
                                Developed = rowSums(prop.lc[,c("Developed, Open Space", 
                                                               "Developed, Low Intensity",
                                                               "Developed, Medium Intensity", 
                                                               "Developed, High Intensity")]),
                                Barren = prop.lc$"Barren Land",
                                Forest = rowSums(prop.lc[,c("Deciduous Forest",
                                                            "Evergreen Forest",
                                                            "Mixed Forest", "Shrub/Scrub")]),
                                Cultivated = rowSums(prop.lc[,c("Herbaceous", "Pasture/Hay",
                                                                "Cultivated Crops")]), 
                                Wetlands = rowSums(prop.lc[,c("Woody Wetlands",
                                                              "Emergent Herbaceous Wetlands")])))

# just checking to see if all landcover categories sum to 1 -- THEY DO!
prop.lc <- prop.lc %>%
  mutate(Sum = rowSums(prop.lc[,c("Unclassified","Open Water", "Developed, Open Space", "Developed, Low Intensity",
                                  "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
                                  "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                                  "Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                  "Emergent Herbaceous Wetlands")]))
head(prop.lc$Sum)

#--------- Mercury Deposition Network Wet Deposition --------------------------
# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 30000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within Xkm of each marsh centroid
avg.wetdep <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T))

# Making a single data frame for the model
HgData30km <- data.frame(cbind(HgSamples, Hglandcover, avg.wetdep))

# This is going to help extraction speeds increase
rm(prop.lc, prop.lc1, prop.lc2, prop.lc3, prop.lc4, prop.lc5,
   prop.lc6, Hglandcover, wetdep.ex, avg.wetdep)

# Making the models ------------------------------
landmod30 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
                  data = HgData30km, REML = F)
depmod30 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID),
                 data = HgData30km, REML = F)
bothmod30 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest)
                  + scale(AvgWetDep) + (1 | SiteID), data = HgData30km, REML = F)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

# Checking model assumptions -------------------------------------
# Checking for multicollinearity
A = data.frame(HgData30km$Developed, HgData30km$Forest, HgData30km$AvgWetDep)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3) # none have collinearity problems, VIF < 3

# Checking for homogeneity of variance & normality of residuals
resid_panel(landmod30, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(depmod30, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(bothmod30, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(landmod30)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(landmod30)$SiteID[,1]) # we are normal

ggqqplot(ranef(depmod30)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(depmod30)$SiteID[,1]) # we are normal

ggqqplot(ranef(bothmod30)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(bothmod30)$SiteID[,1]) # we are normal

#------------------------------------------- 35KM BUFFER MODEL ------------------------------------------
# reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# now there is agreement between raster crs and coordinates crs

# splitting data set so extraction can run
points_spdf1 <- points_spdf[1:20,]
points_spdf2 <- points_spdf[21:40,]
points_spdf3 <- points_spdf[41:60,]
points_spdf4 <- points_spdf[61:80,]
points_spdf5 <- points_spdf[81:100,]
points_spdf6 <- points_spdf[101:120,]

# Extracting NLCD values and turning into dataframes
lc1 <- raster::extract(x = nlcd_raster, y = points_spdf1, buffer = 35000, df = T)
lc2 <- raster::extract(x = nlcd_raster, y = points_spdf2, buffer = 35000, df = T)
lc3 <- raster::extract(x = nlcd_raster, y = points_spdf3, buffer = 35000, df = T)
lc4 <- raster::extract(x = nlcd_raster, y = points_spdf4, buffer = 35000, df = T)
lc5 <- raster::extract(x = nlcd_raster, y = points_spdf5, buffer = 35000, df = T)
lc6 <- raster::extract(x = nlcd_raster, y = points_spdf6, buffer = 35000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc1 <- as.data.frame(lc1) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc1)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc2 <- as.data.frame(lc2) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc2)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc3 <- as.data.frame(lc3) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc3)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc4 <- as.data.frame(lc4) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc4)


# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc5 <- as.data.frame(lc5) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc5)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc6 <- as.data.frame(lc6) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc6)

# joining data frames together
prop.lc  <- prop.lc1 %>% 
  full_join(prop.lc2,  by.y = "ID") %>% 
  full_join(prop.lc3,  by.y = "ID") %>%
  full_join(prop.lc4,  by.y = "ID") %>% 
  full_join(prop.lc5,  by.y = "ID") %>% 
  full_join(prop.lc6,  by.y = "ID")
# renaming ID's
prop.lc[1:120, "ID"] <- 1:120

# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.lc) <- c("ID", "Unclassified", "Open Water", "Developed, Open Space",
                       "Developed, Low Intensity", "Developed, Medium Intensity",
                       "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                       "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                       "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                       "Emergent Herbaceous Wetlands")
head(prop.lc)

# lumping landcover types
Hglandcover <- data.frame(cbind(Water = rowSums(prop.lc[,c("Unclassified", "Open Water")]),
                                Developed = rowSums(prop.lc[,c("Developed, Open Space", 
                                                               "Developed, Low Intensity",
                                                               "Developed, Medium Intensity", 
                                                               "Developed, High Intensity")]),
                                Barren = prop.lc$"Barren Land",
                                Forest = rowSums(prop.lc[,c("Deciduous Forest",
                                                            "Evergreen Forest",
                                                            "Mixed Forest", "Shrub/Scrub")]),
                                Cultivated = rowSums(prop.lc[,c("Herbaceous", "Pasture/Hay",
                                                                "Cultivated Crops")]), 
                                Wetlands = rowSums(prop.lc[,c("Woody Wetlands",
                                                              "Emergent Herbaceous Wetlands")])))

# just checking to see if all landcover categories sum to 1 -- THEY DO!
prop.lc <- prop.lc %>%
  mutate(Sum = rowSums(prop.lc[,c("Unclassified","Open Water", "Developed, Open Space", "Developed, Low Intensity",
                                  "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
                                  "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                                  "Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                  "Emergent Herbaceous Wetlands")]))
head(prop.lc$Sum)

#--------- Mercury Deposition Network Wet Deposition --------------------------
# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 35000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within Xkm of each marsh centroid
avg.wetdep <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T))

# Making a single data frame for the model
HgData35km <- data.frame(cbind(HgSamples, Hglandcover, avg.wetdep))

# This is going to help extraction speeds increase
rm(lc, lc1, lc2, lc3, lc4, lc5, lc6, prop.lc, prop.lc1, prop.lc2, prop.lc3, prop.lc4, prop.lc5,
   prop.lc6, Hglandcover, wetdep.ex, avg.wetdep)

# Making the models ------------------------------
landmod35 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
                  data = HgData35km, REML = F)
depmod35 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID),
                 data = HgData35km, REML = F)
bothmod35 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest)
                  + scale(AvgWetDep) + (1 | SiteID), data = HgData35km, REML = F)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

# Checking model assumptions -------------------------------------
# Checking for multicollinearity
A = data.frame(HgData35km$Developed, HgData35km$Forest, HgData35km$AvgWetDep)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3) # none have collinearity problems, VIF < 3

# Checking for homogeneity of variance & normality of residuals
resid_panel(landmod35, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(depmod35, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(bothmod35, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(landmod35)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(landmod35)$SiteID[,1]) # we are normal

ggqqplot(ranef(depmod35)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(depmod35)$SiteID[,1]) # we are normal

ggqqplot(ranef(bothmod35)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(bothmod35)$SiteID[,1]) # we are normal


#------------------------------------------- 40KM BUFFER MODEL ------------------------------------------
# reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# now there is agreement between raster crs and coordinates crs

# splitting data set so extraction can run
points_spdf1 <- points_spdf[1:20,]
points_spdf2 <- points_spdf[21:40,]
points_spdf3 <- points_spdf[41:60,]
points_spdf4 <- points_spdf[61:80,]
points_spdf5 <- points_spdf[81:100,]
points_spdf6 <- points_spdf[101:120,]

# Extracting NLCD values and turning into dataframes
lc1 <- raster::extract(x = nlcd_raster, y = points_spdf1, buffer = 40000, df = T)
lc2 <- raster::extract(x = nlcd_raster, y = points_spdf2, buffer = 40000, df = T)
lc3 <- raster::extract(x = nlcd_raster, y = points_spdf3, buffer = 40000, df = T)
lc4 <- raster::extract(x = nlcd_raster, y = points_spdf4, buffer = 40000, df = T)
lc5 <- raster::extract(x = nlcd_raster, y = points_spdf5, buffer = 40000, df = T)
lc6 <- raster::extract(x = nlcd_raster, y = points_spdf6, buffer = 40000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc1 <- as.data.frame(lc1) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc1)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc2 <- as.data.frame(lc2) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc2)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc3 <- as.data.frame(lc3) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc3)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc4 <- as.data.frame(lc4) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc4)


# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc5 <- as.data.frame(lc5) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc5)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc6 <- as.data.frame(lc6) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc6)

# joining data frames together
prop.lc  <- prop.lc1 %>% 
  full_join(prop.lc2,  by.y = "ID") %>% 
  full_join(prop.lc3,  by.y = "ID") %>%
  full_join(prop.lc4,  by.y = "ID") %>% 
  full_join(prop.lc5,  by.y = "ID") %>% 
  full_join(prop.lc6,  by.y = "ID")
# renaming ID's
prop.lc[1:120, "ID"] <- 1:120


# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.lc) <- c("ID", "Unclassified", "Open Water", "Developed, Open Space",
                       "Developed, Low Intensity", "Developed, Medium Intensity",
                       "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                       "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                       "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                       "Emergent Herbaceous Wetlands")
head(prop.lc)

# lumping landcover types
Hglandcover <- data.frame(cbind(Water = rowSums(prop.lc[,c("Unclassified", "Open Water")]),
                                Developed = rowSums(prop.lc[,c("Developed, Open Space", 
                                                               "Developed, Low Intensity",
                                                               "Developed, Medium Intensity", 
                                                               "Developed, High Intensity")]),
                                Barren = prop.lc$"Barren Land",
                                Forest = rowSums(prop.lc[,c("Deciduous Forest",
                                                            "Evergreen Forest",
                                                            "Mixed Forest", "Shrub/Scrub")]),
                                Cultivated = rowSums(prop.lc[,c("Herbaceous", "Pasture/Hay",
                                                                "Cultivated Crops")]), 
                                Wetlands = rowSums(prop.lc[,c("Woody Wetlands",
                                                              "Emergent Herbaceous Wetlands")])))

# just checking to see if all landcover categories sum to 1 -- THEY DO!
prop.lc <- prop.lc %>%
  mutate(Sum = rowSums(prop.lc[,c("Unclassified","Open Water", "Developed, Open Space", "Developed, Low Intensity",
                                  "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
                                  "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                                  "Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                  "Emergent Herbaceous Wetlands")]))
head(prop.lc$Sum)

##--------- Mercury Deposition Network Wet Deposition --------------------------
# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 40000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within Xkm of each marsh centroid
avg.wetdep <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T))

# Making a single data frame for the model
HgData40km <- data.frame(cbind(HgSamples, Hglandcover, avg.wetdep))

# This is going to help extraction speeds increase
rm(lc, lc1, lc2, lc3, lc4, lc5, lc6, prop.lc, prop.lc1, prop.lc2, prop.lc3, prop.lc4, prop.lc5,
   prop.lc6, Hglandcover, wetdep.ex, avg.wetdep)

# Making the models ------------------------------
landmod40 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
                  data = HgData40km, REML = F)
depmod40 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID),
                 data = HgData40km, REML = F)
bothmod40 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest)
                  + scale(AvgWetDep) + (1 | SiteID), data = HgData40km, REML = F)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

# Checking model assumptions -------------------------------------
# Checking for multicollinearity
A = data.frame(HgData40km$Developed, HgData40km$Forest, HgData40km$AvgWetDep)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3) # none have collinearity problems, VIF < 3

# Checking for homogeneity of variance & normality of residuals
resid_panel(landmod40, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(depmod40, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(bothmod40, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(landmod40)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(landmod40)$SiteID[,1]) # we are normal

ggqqplot(ranef(depmod40)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(depmod40)$SiteID[,1]) # we are normal

ggqqplot(ranef(bothmod40)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(bothmod40)$SiteID[,1]) # we are normal


#------------------------------------------- 45KM BUFFER MODEL ------------------------------------------
# reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# now there is agreement between raster crs and coordinates crs

# splitting data set so extraction can run
points_spdf1 <- points_spdf[1:20,]
points_spdf2 <- points_spdf[21:40,]
points_spdf3 <- points_spdf[41:60,]
points_spdf4 <- points_spdf[61:80,]
points_spdf5 <- points_spdf[81:100,]
points_spdf6 <- points_spdf[101:120,]

# Extracting NLCD values and turning into dataframes
lc1 <- raster::extract(x = nlcd_raster, y = points_spdf1, buffer = 45000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc1 <- as.data.frame(lc1) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc1)

rm(lc1)

lc2 <- raster::extract(x = nlcd_raster, y = points_spdf2, buffer = 45000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc2 <- as.data.frame(lc2) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc2)

rm(lc2)

lc3 <- raster::extract(x = nlcd_raster, y = points_spdf3, buffer = 45000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc3 <- as.data.frame(lc3) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc3)

rm(lc3)

lc4 <- raster::extract(x = nlcd_raster, y = points_spdf4, buffer = 45000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc4 <- as.data.frame(lc4) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc4)

rm(lc4)

lc5 <- raster::extract(x = nlcd_raster, y = points_spdf5, buffer = 45000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc5 <- as.data.frame(lc5) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc5)

rm(lc5)

lc6 <- raster::extract(x = nlcd_raster, y = points_spdf6, buffer = 45000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc6 <- as.data.frame(lc6) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc6)

rm(lc6)

# joining data frames together
prop.lc  <- prop.lc1 %>% 
  full_join(prop.lc2,  by.y = "ID") %>% 
  full_join(prop.lc3,  by.y = "ID") %>%
  full_join(prop.lc4,  by.y = "ID") %>% 
  full_join(prop.lc5,  by.y = "ID") %>% 
  full_join(prop.lc6,  by.y = "ID")
# renaming ID's
prop.lc[1:120, "ID"] <- 1:120

# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.lc) <- c("ID", "Unclassified", "Open Water", "Developed, Open Space",
                       "Developed, Low Intensity", "Developed, Medium Intensity",
                       "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                       "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                       "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                       "Emergent Herbaceous Wetlands")
head(prop.lc)

# lumping landcover types
Hglandcover <- data.frame(cbind(Water = rowSums(prop.lc[,c("Unclassified", "Open Water")]),
                                Developed = rowSums(prop.lc[,c("Developed, Open Space", 
                                                               "Developed, Low Intensity",
                                                               "Developed, Medium Intensity", 
                                                               "Developed, High Intensity")]),
                                Barren = prop.lc$"Barren Land",
                                Forest = rowSums(prop.lc[,c("Deciduous Forest",
                                                            "Evergreen Forest",
                                                            "Mixed Forest", "Shrub/Scrub")]),
                                Cultivated = rowSums(prop.lc[,c("Herbaceous", "Pasture/Hay",
                                                                "Cultivated Crops")]), 
                                Wetlands = rowSums(prop.lc[,c("Woody Wetlands",
                                                              "Emergent Herbaceous Wetlands")])))

# just checking to see if all landcover categories sum to 1 -- THEY DO!
prop.lc <- prop.lc %>%
  mutate(Sum = rowSums(prop.lc[,c("Unclassified","Open Water", "Developed, Open Space", "Developed, Low Intensity",
                                  "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
                                  "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                                  "Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                  "Emergent Herbaceous Wetlands")]))
head(prop.lc$Sum)

#--------- Mercury Deposition Network Wet Deposition --------------------------
# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 45000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within Xkm of each marsh centroid
avg.wetdep <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T))

# Making a single data frame for the model
HgData45km <- data.frame(cbind(HgSamples, Hglandcover, avg.wetdep))

# This is going to help extraction speeds increase
rm(prop.lc, prop.lc1, prop.lc2, prop.lc3, prop.lc4, prop.lc5,
   prop.lc6, Hglandcover, wetdep.ex, avg.wetdep)

# Making the models ------------------------------
landmod45 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
                  data = HgData45km, REML = F)
depmod45 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID),
                 data = HgData45km, REML = F)
bothmod45 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest)
                  + scale(AvgWetDep) + (1 | SiteID), data = HgData45km, REML = F)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

# Checking model assumptions -------------------------------------
# Checking for multicollinearity
A = data.frame(HgData45km$Developed, HgData45km$Forest, HgData45km$AvgWetDep)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3) # none have collinearity problems, VIF < 3

# Checking for homogeneity of variance & normality of residuals
resid_panel(landmod45, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(depmod45, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(bothmod45, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(landmod45)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(landmod45)$SiteID[,1]) # we are normal

ggqqplot(ranef(depmod45)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(depmod45)$SiteID[,1]) # we are normal

ggqqplot(ranef(bothmod45)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(bothmod45)$SiteID[,1]) # we are normal

#------------------------------------------- 50KM BUFFER MODEL ------------------------------------------
# reproject spdf to NLCD crs so we can extract the data
points_spdf <- spTransform(points_spdf, crs(nlcd_raster))
# now there is agreement between raster crs and coordinates crs

# splitting data set so extraction can run
points_spdf1 <- points_spdf[1:20,]
points_spdf2 <- points_spdf[21:40,]
points_spdf3 <- points_spdf[41:60,]
points_spdf4 <- points_spdf[61:80,]
points_spdf5 <- points_spdf[81:100,]
points_spdf6 <- points_spdf[101:120,]

# Extracting NLCD values and turning into dataframes
lc1 <- raster::extract(x = nlcd_raster, y = points_spdf1, buffer = 50000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc1 <- as.data.frame(lc1) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc1)

rm(lc1)

lc2 <- raster::extract(x = nlcd_raster, y = points_spdf2, buffer = 50000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc2 <- as.data.frame(lc2) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc2)

rm(lc2)

lc3 <- raster::extract(x = nlcd_raster, y = points_spdf3, buffer = 50000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc3 <- as.data.frame(lc3) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc3)

rm(lc3)

lc4 <- raster::extract(x = nlcd_raster, y = points_spdf4, buffer = 50000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc4 <- as.data.frame(lc4) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc4)

rm(lc4)

lc5 <- raster::extract(x = nlcd_raster, y = points_spdf5, buffer = 50000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc5 <- as.data.frame(lc5) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc5)

rm(lc5)

lc6 <- raster::extract(x = nlcd_raster, y = points_spdf6, buffer = 50000, df = T)

# Calculating percent landcover within Xkm radius of each marsh centroid and creating a data frame
prop.lc6 <- as.data.frame(lc6) %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%        # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_type), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_type, pland)                 # convert to long format
head(prop.lc6)

rm(lc6)

# joining data frames together
prop.lc  <- prop.lc1 %>% 
  full_join(prop.lc2,  by.y = "ID") %>% 
  full_join(prop.lc3,  by.y = "ID") %>%
  full_join(prop.lc4,  by.y = "ID") %>% 
  full_join(prop.lc5,  by.y = "ID") %>% 
  full_join(prop.lc6,  by.y = "ID")
# renaming ID's
prop.lc[1:120, "ID"] <- 1:120

# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.lc) <- c("ID", "Unclassified", "Open Water", "Developed, Open Space",
                       "Developed, Low Intensity", "Developed, Medium Intensity",
                       "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                       "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                       "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                       "Emergent Herbaceous Wetlands")
head(prop.lc)

# lumping landcover types
Hglandcover <- data.frame(cbind(Water = rowSums(prop.lc[,c("Unclassified", "Open Water")]),
                                Developed = rowSums(prop.lc[,c("Developed, Open Space", 
                                                               "Developed, Low Intensity",
                                                               "Developed, Medium Intensity", 
                                                               "Developed, High Intensity")]),
                                Barren = prop.lc$"Barren Land",
                                Forest = rowSums(prop.lc[,c("Deciduous Forest",
                                                            "Evergreen Forest",
                                                            "Mixed Forest", "Shrub/Scrub")]),
                                Cultivated = rowSums(prop.lc[,c("Herbaceous", "Pasture/Hay",
                                                                "Cultivated Crops")]), 
                                Wetlands = rowSums(prop.lc[,c("Woody Wetlands",
                                                              "Emergent Herbaceous Wetlands")])))

# just checking to see if all landcover categories sum to 1 -- THEY DO!
prop.lc <- prop.lc %>%
  mutate(Sum = rowSums(prop.lc[,c("Unclassified","Open Water", "Developed, Open Space", "Developed, Low Intensity",
                                  "Developed, Medium Intensity", "Developed, High Intensity", "Barren Land",
                                  "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                                  "Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                  "Emergent Herbaceous Wetlands")]))
head(prop.lc$Sum)

#--------- Mercury Deposition Network Wet Deposition --------------------------
# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 50000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within Xkm of each marsh centroid
avg.wetdep <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T))

# Making a single data frame for the model
HgData50km <- data.frame(cbind(HgSamples, Hglandcover, avg.wetdep))

# This is going to help extraction speeds increase
rm(prop.lc, prop.lc1, prop.lc2, prop.lc3, prop.lc4, prop.lc5,
   prop.lc6, Hglandcover, wetdep.ex, avg.wetdep)

# Making the models ------------------------------
landmod50 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
                  data = HgData50km, REML = F)
depmod50 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID),
                 data = HgData50km, REML = F)
bothmod50 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest)
                  + scale(AvgWetDep) + (1 | SiteID), data = HgData50km, REML = F)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

# Checking model assumptions -------------------------------------
# Checking for multicollinearity
A = data.frame(HgData50km$Developed, HgData50km$Forest, HgData50km$AvgWetDep)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3) # none have collinearity problems, VIF < 3

# Checking for homogeneity of variance & normality of residuals
resid_panel(landmod50, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(depmod50, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for homogeneity of variance & normality of residuals
resid_panel(bothmod50, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
# residual plots look fairly homogeneous
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(landmod50)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(landmod50)$SiteID[,1]) # we are normal

ggqqplot(ranef(depmod50)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(depmod50)$SiteID[,1]) # we are normal

ggqqplot(ranef(bothmod50)$SiteID[,1])  # looks fairly normal
shapiro.test(ranef(bothmod50)$SiteID[,1]) # we are normal


#------------------------------------------- MODEL SET SUMMARY-----------------------------------------------------
# creating candidate model set
wetdepcandset = list(depmod5, depmod6, depmod7, depmod8, bothmod9, depmod10,
                     depmod15, depmod20, depmod25, depmod30, depmod35, depmod40,
                     depmod45, depmod50)
depmodnames = c("depmod5", "depmod6", "depmod7", "depmod8", "depmod9", "depmod10",
                "depmod15", "depmod20", "depmod25", "depmod30", "depmod35",
                "depmod40", "depmod45", "depmod50")
depmodelsetsummary <- as.data.frame(aictab(cand.set = wetdepcandset, modnames = depmodnames))


landcandset = list(landmod5, landmod6, landmod7, landmod8, bothmod9, landmod10,
                   landmod15, landmod20, landmod25, landmod30, landmod35, landmod40,
                   landmod45, landmod50)
landmodnames = c("landmod5", "landmod6", "landmod7", "landmod8", "landmod9", "landmod10",
                 "landmod15", "landmod20", "landmod25", "landmod30", "landmod35",
                 "landmod40", "landmod45", "landmod50")
landmodelsetsummary <- as.data.frame(aictab(cand.set = landcandset, modnames = landmodnames))


# The model used in this comparison contained independent variables that would reasonably be influenced
# by this radius selection: developed and forested land cover, in addition to Hg wet deposition,
# within a 1–50 km radius around each sampling point. We defined the scale of effect to be
# the buffer radius that resulted in the lowest model AICc value.

# It's important to note that regardless of which full model was used (landmod or bothmod) 
# or whether "Species" was included as a fixed effect, the AIC hierarchy was always the same, and 
# 30 km was always the best scale of effect. We did not consider the model structure of "depmod" to
# make biological sense in this situation without the influence of land cover.

bothcandset = list(bothmod5, bothmod6, bothmod7, bothmod8, bothmod9, bothmod10,
                   bothmod15, bothmod20, bothmod25, bothmod30, bothmod35, bothmod40,
                   bothmod45, bothmod50)
bothmodnames = c("bothmod5", "bothmod6", "bothmod7", "bothmod8", "bothmod9", "bothmod10",
                 "bothmod15", "bothmod20", "bothmod25", "bothmod30", "bothmod35", 
                 "bothmod40", "bothmod45", "bothmod50")
bothmodelsetsummary <- as.data.frame(aictab(cand.set = bothcandset, modnames = bothmodnames))

bothmodelsetsummary[, "Radius"] <- c(30, 35, 25, 40, 20, 5, 45, 6, 50, 7, 8,
                                     9, 10, 15)

# 30km Radius is the "winner" but is statistically equivalent to 20-45km models since ∆AIC<2
bothmodelsetsummary %>% 
  ggplot() + # use ∆AIC to interpret easier
  #   geom_rect(data = bothmodelsetsummary, inherit.aes = FALSE,
  #             aes(xmin = 20, xmax = 45, ymin = -Inf, ymax = 2), 
  #             fill = 'gray', alpha = 0.2) +
  #  geom_smooth(data = bothmodelsetsummary, mapping = aes(x = Radius, y = Delta_AICc),
  #              color = "#548235", se = F, size = 3) +
  #geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2), fill = "gray80") +
  #geom_smooth(mapping = aes(x = Radius, y = Delta_AICc),
  #            color = "black", se = F) +
  geom_line(mapping = aes(x = Radius, y = Delta_AICc), color = "gray60", linetype = "dashed") +
  geom_point(mapping = aes(x = Radius, y = Delta_AICc)) +
  theme_classic() +
  labs(x = "Radius (km)", y = "∆AICc") +
  # ggtitle("ln(THg) ~ % Developed + % Forest + Annual Precipitation + (1 | Site)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"),
        axis.title.x = element_text(size = 14, face = "bold", color = "black"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"))

ggsave("Sayers_AICPlot.jpg", dpi = 750, width = 174, height = 116, units = "mm")
