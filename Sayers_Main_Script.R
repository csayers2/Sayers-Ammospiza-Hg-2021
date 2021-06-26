
#--------------------------------------- LOADING & PARTITIONING THE DATA ------------------------------------------
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(lubridate)

HgSamples <- read_excel("Sayers_SHARP_Hg_Samples.xlsx")

# Creating SALS only data set
SALSData <- HgSamples %>% 
  filter(Species == "SALS") %>%
  group_by(SiteID) %>%
  mutate(MeanSALS = mean(HgConcentration),
         SD.SALS = sd(HgConcentration))
# computing cumulative statistics for SALS
summary(SALSData$HgConcentration) 
sd(SALSData$HgConcentration)

# Creating SESP only data set
SESPData <- HgSamples %>% 
  filter(Species == "SESP") %>%
  group_by(SiteID) %>%
  mutate(MeanSESP = mean(HgConcentration),
         SD.SESP = sd(HgConcentration))
# computing cumulative statistics for SESP
summary(SESPData$HgConcentration)
sd(SESPData$HgConcentration)

# Calculating cumulative species means per site
SiteData <- HgSamples %>%
  group_by(SiteID) %>%
  mutate(MeanTMS = mean(HgConcentration),
         SD.TMS = sd(HgConcentration))

# Location summary statistics (just change the Site ID or State name to look up a region)
LocationSALSData <- filter(SALSData, State %in% c("NJ"))
sd(LocationSALSData$HgConcentration)
summary(LocationSALSData$HgConcentration)

LocationSALSData <- filter(SALSData, SiteID %in% c("31427_HHO153a"))
sd(LocationSALSData$HgConcentration)
summary(LocationSALSData$HgConcentration)

# How many birds are above the toxicity reference values published in Jackson et al. (2011)?
# Change filters to answer different questions
count(HgSamples, Species == "SALS" & HgConcentration >= 0.7)
count(HgSamples, Species == "SALS" & HgConcentration >= 0.7 & HgConcentration < 1.2)
count(HgSamples, Species == "SALS" & HgConcentration >= 1.2)

count(HgSamples, Species == "SESP" & HgConcentration >= 0.7)
count(HgSamples, Species == "SESP" & HgConcentration >= 0.7 & HgConcentration < 1.2)
count(HgSamples, Species == "SESP" & HgConcentration >= 1.2)

count(HgSamples, HgConcentration >= 0.7)

count(HgSamples, Species == "SALS" & State == "NJ")
count(HgSamples, Species == "SALS" & State == "NJ" & HgConcentration >= 0.7)

count(HgSamples, Species == "SESP" & State == "NJ")
count(HgSamples, Species == "SESP" & State == "NJ" & HgConcentration >= 0.7)


#------------------------------------------- DATA EXPLORATION --------------------------------------------------
# much of the strategy below is from Zurr et al. (2010)  https://doi.org/10.1111/j.2041-210X.2009.00001.x

# OUTLIERS & NORMALITY OF RESPONSE VARIABLE ----------

#checking cumulative (non-species specific) data (n = 120)
summary(HgSamples$HgConcentration) # mean = 0.4919
sd(HgSamples$HgConcentration) # sd = 0.3021768
ggdensity(HgSamples$HgConcentration, xlab = "THg Concentration (µg/g ww)") # some high outliers 
ggqqplot(HgSamples$HgConcentration, ylab = "THg Concentration (µg/g ww)") # tails stray from normal
shapiro.test(HgSamples$HgConcentration) # W = 0.88055, p-value = 2.242e-08, not normal
# log-transformation may be necessary for linear models

# Creating new column with natural-log transformed THg concentrations + date
HgSamples <- HgSamples %>%
  mutate(lHgConcentration = log(HgConcentration)) %>%
  mutate(Date = make_date(Year, Month, Day),
         JulianDate = yday(Date))
ggdensity(HgSamples$lHgConcentration, xlab = "ln(THg) Concentration (µg/g ww)")
ggqqplot(HgSamples$lHgConcentration) # much better than before
shapiro.test(HgSamples$lHgConcentration) # W = 0.9787, p-value = 0.05391, just barely normal!

# Boxplot & Cleavland dotplot of raw data by species
boxplot(HgSamples$HgConcentration ~ HgSamples$Species)
# we have 4 high outliers (2 SALS, 2 SESP) in the center of the range
ggplot(HgSamples, aes(x = HgConcentration, y = reorder(SiteID, Latitude), color = Species)) +
  geom_point() +
  labs(x = "THg Concentration (µg/g ww)", y = "Site ID from S to N")

# Cleavland dotplot of natural-log transformed data, no apparent outliers anymore
ggplot(HgSamples, aes(x = lHgConcentration, y = reorder(SiteID, Latitude), color = Species)) +
  geom_point() +
  labs(x = "ln(THg) Concentration (µg/g ww)", y = "Site ID from S to N")

# Violin plot of species outliers, SALS is slightly higher
ggplot(data = HgSamples,
       mapping = aes(x = Species, y = HgConcentration, fill = Species)) +
  geom_violin(size = 1) +
  geom_boxplot(width = 0.1, fill="white", size = 1) +
  stat_summary(fun = mean, geom = "point", fill = "black", shape = 8, size = 3) +
  scale_fill_manual(values=c("#808080","#C0C0C0")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Species", y = "THg Concentration (µg/g ww)") +
  scale_x_discrete(labels = c("Saltmarsh\nSparrow", "Seaside\nSparrow")) +
  geom_hline(yintercept = c(0.7), size = 1, linetype = "dashed", color = "black") + 
  # EC10 specified by Jackson et al. (2011)
  scale_y_continuous(limits = c(0.0, 1.6))

#------------------------------------------- MANIPULATING RASTER DATA ------------------------------------------
library(raster)
library(sp)

# Making a spatial data frame from Hg sample coordinates
coords <- cbind(HgSamples$Longitude, HgSamples$Latitude)
points_spdf <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
plot(points_spdf)

# Importing contiguous USA NLCD raster for 2016 (this is not available yet using FedData::get_nlcd)
# NLCD 2016 Land Cover (CONUS) - https://www.mrlc.gov/data
# It's very imoportant that the .img stays within the original NLCD folder - I think it communicates with metadata, etc.
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

#----------------------------------------- EXTRACTING RASTERS @ 30KM BUFFER ------------------------------------------
# splitting the data set so the extraction can run
points_spdf1 <- points_spdf[1:20,]
points_spdf2 <- points_spdf[21:40,]
points_spdf3 <- points_spdf[41:60,]
points_spdf4 <- points_spdf[61:80,]
points_spdf5 <- points_spdf[81:100,]
points_spdf6 <- points_spdf[101:120,]

# I chose a 30 km buffer radius based on the results of RadiusSelectionWetDep.R
# We deemed 30 km to be the most appropriate scale of effect, since the 30 km model produced the lowest AIC score

# Extracting NLCD values and turning into dataframes
lc1 <- raster::extract(x = nlcd_raster, y = points_spdf1, buffer = 30000, df = T)
# Calculating percent landcover within 30km radius of each marsh centroid and creating a data frame
prop.lc1 <- as.data.frame(lc1) %>%
  setNames(c("ID", "lc_class")) %>%        # rename for ease
  group_by(ID, lc_class) %>%               # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_class, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_class), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_class, pland)                 # convert to long format
head(prop.lc1)
rm(lc1)

lc2 <- raster::extract(x = nlcd_raster, y = points_spdf2, buffer = 30000, df = T)
# Calculating percent landcover within 30km radius of each marsh centroid and creating a data frame
prop.lc2 <- as.data.frame(lc2) %>%
  setNames(c("ID", "lc_class")) %>%        # rename for ease
  group_by(ID, lc_class) %>%               # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_class, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_class), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_class, pland)                 # convert to long format
head(prop.lc2)
rm(lc2)

lc3 <- raster::extract(x = nlcd_raster, y = points_spdf3, buffer = 30000, df = T)
# Calculating percent landcover within 30km radius of each marsh centroid and creating a data frame
prop.lc3 <- as.data.frame(lc3) %>%
  setNames(c("ID", "lc_class")) %>%        # rename for ease
  group_by(ID, lc_class) %>%               # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_class, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_class), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_class, pland)                 # convert to long format
head(prop.lc3)
rm(lc3)

lc4 <- raster::extract(x = nlcd_raster, y = points_spdf4, buffer = 30000, df = T)
# Calculating percent landcover within 30km radius of each marsh centroid and creating a data frame
prop.lc4 <- as.data.frame(lc4) %>%
  setNames(c("ID", "lc_class")) %>%        # rename for ease
  group_by(ID, lc_class) %>%               # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_class, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_class), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_class, pland)                 # convert to long format
head(prop.lc4)
rm(lc4)

lc5 <- raster::extract(x = nlcd_raster, y = points_spdf5, buffer = 30000, df = T)
# Calculating percent landcover within 30km radius of each marsh centroid and creating a data frame
prop.lc5 <- as.data.frame(lc5) %>%
  setNames(c("ID", "lc_class")) %>%        # rename for ease
  group_by(ID, lc_class) %>%               # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_class, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_class), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_class, pland)                 # convert to long format
head(prop.lc5)
rm(lc5)

lc6 <- raster::extract(x = nlcd_raster, y = points_spdf6, buffer = 30000, df = T)
# Calculating percent landcover within 30km radius of each marsh centroid and creating a data frame
prop.lc6 <- as.data.frame(lc6) %>%
  setNames(c("ID", "lc_class")) %>%        # rename for ease
  group_by(ID, lc_class) %>%               # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_class, pland) %>%   # keep only these vars
  tidyr::complete(ID, nesting(lc_class), 
                  fill = list(pland = 0)) %>%   # fill in implicit landcover 0s
  tidyr::spread(lc_class, pland)                 # convert to long format
head(prop.lc6)
rm(lc6)

# joining data frames together
prop.landcover30km <- prop.lc1 %>% 
  full_join(prop.lc2) %>% 
  full_join(prop.lc3) %>%
  full_join(prop.lc4) %>% 
  full_join(prop.lc5) %>% 
  full_join(prop.lc6)
# renaming ID's
prop.landcover30km[1:120, "ID"] <- 1:120

# Renaming NLCD data frame columns so they are easier to interpret
colnames(prop.landcover30km) <- c("ID", "Unclassified", "Open Water", "Developed, Open Space",
                                  "Developed, Low Intensity", "Developed, Medium Intensity",
                                  "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                                  "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                                  "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                                  "Emergent Herbaceous Wetlands")

# Checking to see if all landcover categories sum to 1 -- THEY DO!
prop.landcover30km <- prop.landcover30km %>%
  mutate(Sum = rowSums(prop.landcover30km[,c("Unclassified", "Open Water", "Developed, Open Space",
                                             "Developed, Low Intensity", "Developed, Medium Intensity",
                                             "Developed, High Intensity", "Barren Land", "Deciduous Forest",
                                             "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                                             "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                                             "Emergent Herbaceous Wetlands")]))
head(prop.landcover30km$Sum)

# Lumping Unclassified and Water cover together
prop.landcover30km <- prop.landcover30km %>%
  mutate(Water = rowSums(prop.landcover30km[,c("Unclassified", "Open Water")]))
head(prop.landcover30km)


#--------- Mercury Deposition Network Wet Deposition --------------------------
# Read in the MDN wet deposition raster info (units in µg/m2)
library(raster)
library(sp)

wetdep <- raster("Hg_dep_2017/Hg_dep_2017.tif")
crs(wetdep) # Albers equal area projection, NAD83 datum, we will need to reproject spdf to match this
plot(wetdep)

# Reproject spdf to MDN crs
points_spdf <- spTransform(points_spdf, crs(wetdep))
# Now there is agreement between raster crs and coordinates crs

# Extracting deposition values
wetdep.ex <- raster::extract(x = wetdep, y = points_spdf, buffer = 30000, df = T)
wetdep.ex <- as.data.frame(wetdep.ex)

# Takes average wet deposition in 2017 within 30km of each marsh centroid
avg.wetdep30km <- wetdep.ex %>%
  setNames(c("ID", "Hg_dep_2017")) %>%
  group_by(ID) %>%
  summarise(AvgWetDep = mean(Hg_dep_2017, na.rm = T)) 

#----------------------------------------- DATA EXPLORATION CONTINUED ------------------------------------------

# OUTLIERS OF PREDICTOR VARIABLES ----------

# Plotting distribution of raw landscape types over the sampling range
library(reshape2)

# unlumped, we have some outlier sites which could be problematic for the models
# most NLCD categories have a median coverage of <5% of the 30km buffer — we should lump relevant ones
lcdf <- melt(data = prop.landcover30km, id.vars = "ID",
             measure.vars = c("Water", "Developed, Open Space", 
                              "Developed, Low Intensity","Developed, Medium Intensity",
                              "Developed, High Intensity", "Barren Land","Deciduous Forest",
                              "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceous",
                              "Pasture/Hay", "Cultivated Crops", "Woody Wetlands",
                              "Emergent Herbaceous Wetlands"))

ggplot(data = lcdf, mapping = aes(x = value, y = variable, fill = variable)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", fill = "black", shape = 8, size = 3) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  labs(x = "Coverage within 30 km Radius", y = "NLCD Land Cover Class") +
  scale_x_continuous(labels = scales::percent) +
  geom_vline(xintercept = c(0.05), linetype = "dashed", color = "black")

# lumping relevant landcover types
Hglandcover30km <- prop.landcover30km %>%
  mutate(Water = rowSums(prop.landcover30km[,c("Unclassified", "Open Water")]),
         Developed = rowSums(prop.landcover30km[,c("Developed, Open Space", 
                                                   "Developed, Low Intensity",
                                                   "Developed, Medium Intensity", 
                                                   "Developed, High Intensity")]),
         Barren = prop.landcover30km$"Barren Land",
         Forest = rowSums(prop.landcover30km[,c("Deciduous Forest", "Evergreen Forest",
                                                "Mixed Forest", "Shrub/Scrub")]),
         Cultivated = rowSums(prop.landcover30km[,c("Herbaceous", "Pasture/Hay", "Cultivated Crops")]), 
         Wetlands = rowSums(prop.landcover30km[,c("Woody Wetlands","Emergent Herbaceous Wetlands")]))

# plotting distribution of lumped landscape types over the sampling range
# lumped, 6 "developed" outliers could be problematic in the models
lumpedlcdf <- melt(data = Hglandcover30km, id.vars = "ID", measure.vars = c("Water", "Developed", "Forest",
                                                                            "Cultivated", "Wetlands"))

ggplot(data = lumpedlcdf, mapping = aes(x = variable, y = value, fill = variable)) +
  geom_violin(size = 1) +
  geom_boxplot(width = 0.1, fill="white", size = 1) +
  stat_summary(fun = mean, geom = "point", fill = "black", shape = 8, size = 3) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold")) +
  labs(x = "Landcover Type", y = "% Coverage within 30km") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = c(0.05), linetype = "dashed", color = "black")

# Cleavland dotplot of final predictor variables
# 6 "developed" outliers could be problematic in the models
contx <- data.frame(cbind(HgSamples, Hglandcover30km, avg.wetdep30km)) %>%
  melt(id.vars = c("ID", "SiteID", "Latitude", "Species", "HgConcentration"),
       measure.vars = c("Water", "Developed", "Forest", "Cultivated", "Wetlands",
                        "AvgWetDep", "JulianDate"))
ggplot(contx, aes(x = value, y = reorder(SiteID, Latitude))) +
  geom_point() +
  facet_wrap(~ variable, scales = "free", nrow = 2) + 
  labs(y = "Site ID from S to N")

# Relationships between X & Y
ggplot(contx, aes(x = value, y = HgConcentration)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ variable, scales = "free", nrow = 2) + 
  labs(y = "THg Concentration (µg/g ww)")

#--------------------------------------- ANALYSIS + LINEAR MODELS ------------------------------------------
# one-way ANOVA to formally test species differences (F1,118 = 5.371, P = 0.0222)
# response variable must be log-transformed to satisfy ANOVA assumptions
# we should definitely include Species in our full model
anova <- aov(log(HgConcentration) ~ Species, data = HgSamples)
summary(anova)

#Making a single data frame for the model
HgData30km <- data.frame(cbind(HgSamples, Hglandcover30km, avg.wetdep30km))

#CHECKING MULTICOLLINEARITY OF PREDICTORS ---------------
library(ggpubr)
library(MuMIn)
library(lme4)
library(usdm)

# wittle the predictor variables down from this list, testing for multicollinearity among
# lumped and unlumped variables

#X = data.frame(HgData30km$Developed..Open.Space, HgData30km$Developed..Low.Intensity,
#               HgData30km$Developed..Medium.Intensity, HgData30km$Developed..High.Intensity,
#               HgData30km$Deciduous.Forest, HgData30km$Evergreen.Forest, HgData30km$Mixed.Forest,
#               HgData30km$Shrub.Scrub, HgData30km$Woody.Wetlands,
#               HgData30km$Emergent.Herbaceous.Wetlands, HgData30km$AvgWetDep,
#               HgData30km$Latitude, HgData30km$Longitude, HgData30km$Weight, HgData30km$JulianDate)

# We initially included the proportion of wetland land cover, latitude, and longitude as fixed effects
# However, independent variables were highly correlated at larger spatial scales, as indicated by variation
# inflation factors (VIF > 3) and correlation coefficients (R2 > 0.7). After removing these factors from
# the full model, there was no evidence of multicollinearity at any of the spatial scales we evaluated.

A = data.frame(HgData30km$Developed,
               HgData30km$Forest,
               HgData30km$AvgWetDep,
               HgData30km$JulianDate)
vif(A)
vifcor(A, th = 0.7)
vifstep(A, th = 3)

fullmodel <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + scale(AvgWetDep) 
                  + Species + scale(JulianDate) + (1 | SiteID), data = HgData30km, REML = F)

# CHECKING MODEL ASSUMPTIONS -------------------------------------
# Checking for homogeneity of variance & normality of residuals
mean(residuals(fullmodel)) # very very close to 0

library(ggResidpanel)
resid_panel(fullmodel, plots = "all", type = NA, bins = 30,
            smoother = T, qqline = T, qqbands = T, scale = 1,
            theme = "bw", axis.text.size = 10, title.text.size = 12,
            title.opt = TRUE, nrow = NULL)
shapiro.test(residuals(fullmodel)) # not normal
# residual plots looks great
# normality plot has a few tail stragglers, but the rest looks good

# Checking for normality of random effects
ggqqplot(ranef(fullmodel)$SiteID[,1]) # few tail stragglers, but the rest looks good
shapiro.test(ranef(fullmodel)$SiteID[,1]) # we are normal

# Checking for autocorrelation/independence
library(lawstat)
acf(HgData30km$HgConcentration) # raw data is autocorrelated
acf(residuals(fullmodel)) # random effects variable corrects for this
runs.test(residuals(fullmodel)) # we do not have autocorrelated data

#CANDIDATE MODEL SET (32 MODELS) -----------------------------------------------
library(AICcmodavg)
library(lmerTest)
fullmodel <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + scale(AvgWetDep) 
                  + Species + scale(JulianDate) + (1 | SiteID), data = HgData30km, REML = F)

mod1 <- lmer(log(HgConcentration) ~ scale(Developed) + (1 | SiteID), data = HgData30km, REML = F)
mod2 <- lmer(log(HgConcentration) ~ scale(Forest) + (1 | SiteID), data = HgData30km, REML = F)
mod3 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + (1 | SiteID), data = HgData30km, REML = F)
mod4 <- lmer(log(HgConcentration) ~ Species + (1 | SiteID), data = HgData30km, REML = F)
mod5 <- lmer(log(HgConcentration) ~ scale(JulianDate) + (1 | SiteID), data = HgData30km, REML = F)
mod6 <- lmer(log(HgConcentration) ~ (1 | SiteID), data = HgData30km, REML = F)

mod7 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + (1 | SiteID),
             data = HgData30km, REML = F)
mod8 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(AvgWetDep) + (1 | SiteID),
             data = HgData30km, REML = F)
mod9 <- lmer(log(HgConcentration) ~ scale(Developed) + Species + (1 | SiteID),
             data = HgData30km, REML = F)
mod10 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(JulianDate) + (1 | SiteID),
              data = HgData30km, REML = F)
mod11 <- lmer(log(HgConcentration) ~ scale(Forest) + scale(AvgWetDep) + (1 | SiteID),
              data = HgData30km, REML = F)
mod12 <- lmer(log(HgConcentration) ~ scale(Forest) + Species + (1 | SiteID),
              data = HgData30km, REML = F)
mod13 <- lmer(log(HgConcentration) ~ scale(Forest) + scale(JulianDate) + (1 | SiteID),
              data = HgData30km, REML = F)
mod14 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + Species + (1 | SiteID),
              data = HgData30km, REML = F)
mod15 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + scale(JulianDate) + (1 | SiteID),
              data = HgData30km, REML = F)
mod16 <- lmer(log(HgConcentration) ~ Species + scale(JulianDate) + (1 | SiteID),
              data = HgData30km, REML = F)

mod17 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + scale(AvgWetDep) 
              + (1 | SiteID), data = HgData30km, REML = F)
mod18 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + Species
              + (1 | SiteID), data = HgData30km, REML = F)
mod19 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + scale(JulianDate)
              + (1 | SiteID), data = HgData30km, REML = F)
mod20 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(AvgWetDep) + Species
              + (1 | SiteID), data = HgData30km, REML = F)
mod21 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(AvgWetDep) + scale(JulianDate)
              + (1 | SiteID), data = HgData30km, REML = F)
mod22 <- lmer(log(HgConcentration) ~ scale(Developed) + Species + scale(JulianDate)
              + (1 | SiteID), data = HgData30km, REML = F)
mod23 <- lmer(log(HgConcentration) ~ scale(Forest) + scale(AvgWetDep) + Species
              + (1 | SiteID), data = HgData30km, REML = F)
mod24 <- lmer(log(HgConcentration) ~ scale(Forest) + scale(AvgWetDep) + scale(JulianDate)
              + (1 | SiteID), data = HgData30km, REML = F)
mod25 <- lmer(log(HgConcentration) ~ scale(Forest) + Species + scale(JulianDate)
              + (1 | SiteID), data = HgData30km, REML = F)
mod26 <- lmer(log(HgConcentration) ~ scale(AvgWetDep) + Species + scale(JulianDate)
              + (1 | SiteID), data = HgData30km, REML = F)

mod27 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + scale(AvgWetDep) + Species
              + (1 | SiteID), data = HgData30km, REML = F)
mod28 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + scale(AvgWetDep) + scale(JulianDate)
              + (1 | SiteID), data = HgData30km, REML = F)
mod29 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(Forest) + Species + scale(JulianDate)
              + (1 | SiteID), data = HgData30km, REML = F)
mod30 <- lmer(log(HgConcentration) ~ scale(Developed) + scale(AvgWetDep) + Species + scale(JulianDate)
              + (1 | SiteID), data = HgData30km, REML = F)
mod31 <- lmer(log(HgConcentration) ~ scale(Forest) + scale(AvgWetDep) + Species + scale(JulianDate)
              + (1 | SiteID), data = HgData30km, REML = F)

# Creating candidate model set
candset = list(fullmodel, mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13,
               mod14, mod15, mod16, mod17, mod18, mod19, mod20, mod21, mod22, mod23,
               mod24, mod25, mod26, mod27, mod28, mod29, mod30, mod31)
modnames = c("fullmod", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9", "mod10",
             "mod11", "mod12", "mod13", "mod14", "mod15", "mod16", "mod17", "mod18", "mod19", "mod20", 
             "mod21", "mod22", "mod23", "mod24", "mod25", "mod26", "mod27", "mod28", "mod29", "mod30",
             "mod31")
modelsetsummary <- as.data.frame(aictab(cand.set = candset, modnames = modnames))
View(modelsetsummary)

# Model averaging for entire model set
modelavg <- model.avg(candset)
summary(modelavg)

confint(modelavg) # unconditional 95% CI
MuMIn::importance(modelavg)
