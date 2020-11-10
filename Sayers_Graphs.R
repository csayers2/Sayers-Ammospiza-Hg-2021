
# This script requires objects from Sayers_SeniorThesis.R

#------------------------------------------ GRAPHS ------------------------------------------

# Ecotoxicology allows free color art, so let's make all the figures more engaging

# Google satellite imagery as a background
library(ggplot2)
library(ggmap)
library(sp)
library(raster)
# register_google(key = "Enter your key here", write = TRUE)

# For google map, you have to give the center of the window you desire
# Possibility for the map type argument: terrain / satellite / roadmap / hybrid
# get the map info
satellitemap <- get_map(location = c(lon = -72.5, lat = 41), maptype = "hybrid", source = "google", zoom = 6)

SALSggmap <- ggmap(satellitemap) +
  stat_summary_2d(data = SALSData, mapping = aes(x = Longitude, y = Latitude, z = MeanSALS), fun = mean,
                  binwidth = c(.25,.25)) + 
  coord_quickmap(xlim = c(-77.5, -67.5), ylim = c(37, 45)) +
  scale_fill_gradient(name = "Mean Blood THg\n(µg/g ww)", high = "red", low = "yellow") +
  labs(x = "Longitude", y = "Latitude") + 
  ggtitle("Saltmarsh Sparrow") +
  theme(plot.title = element_text(size = 32, face = "bold", color = "black", hjust = 0.5),
        axis.title.x = element_text(size = 32, face = "bold"),
        axis.title.y = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 24)) +
  scale_x_continuous(limits = c(-77.5, -67.5)) +
  scale_y_continuous(limits = c(37, 45))

SESPggmap <- ggmap(satellitemap) +
  stat_summary_2d(data = SESPData, mapping = aes(x = Longitude, y = Latitude, z = MeanSESP), fun = mean,
                  binwidth = c(.25,.25)) + 
  coord_quickmap(xlim = c(-77.5, -67.5), ylim = c(37, 45)) +
  scale_fill_gradient(name = "Mean Blood THg\n(µg/g ww)", high = "red", low = "yellow") +
  labs(x = "Longitude", y = "Latitude") + 
  ggtitle("Seaside Sparrow") +
  theme(plot.title = element_text(size = 32, face = "bold", color = "black", hjust = 0.5),
        axis.title.x = element_text(size = 32, face = "bold"),
        axis.title.y = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 24)) +
  scale_x_continuous(limits = c(-77.5, -67.5)) +
  scale_y_continuous(limits = c(37, 45))

ggarrange(SALSggmap, SESPggmap, heights = c(1, 1), nrow = 1, ncol = 2,
          common.legend = TRUE, legend = "right")

#violin plot
ggplot(data = HgSamples, mapping = aes(x = Species, y = HgConcentration, fill = Species)) +
  geom_violin(size = 1) +
  geom_boxplot(width = 0.1, fill = "white", size = 1) +
  stat_summary(fun = mean, geom = "point", fill = "black", shape = 8, size = 5) +
  scale_fill_manual(values = c("#E69F00", "#548235")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 32, face = "bold", color = "black"),
        axis.text.y = element_text(size = 32, color = "black")) +
  labs(x = "Species", y = "THg Concentration (µg/g ww)") +
  scale_x_discrete(labels = c("Saltmarsh Sparrow", "Seaside Sparrow")) +
  geom_hline(yintercept = c(0.7, 1.2), size = 1, linetype = "dashed", color = "red") + 
  scale_y_continuous(limits = c(0.0, 1.6))






























