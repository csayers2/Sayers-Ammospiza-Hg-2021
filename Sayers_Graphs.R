# This script requires objects from Sayers_SeniorThesis.R

#------------------------------------------ GRAPHS ------------------------------------------
library(dplyr)
library(ggpubr)
library(ggplot2)
library(mapproj)

# Ecotoxicology allows free color art, so let's make all the figures more engaging

#figure 1 - map
states <- map_data("state")

# Dots to show individual samples

#Saltmarsh
SALSData <- SALSData %>%
  arrange(HgConcentration) # this is so that the points will be plotted big -> small

SALSmap.man <- ggplot() + 
  geom_polygon(data = states, aes(x = long, y = lat, group = group), color = "black", size = 0.2, fill = "gray") +
  geom_point(data = SALSData, mapping = aes(x = Longitude, y = Latitude, color = HgConcentration, size = HgConcentration)) +
  scale_color_gradient(high = "red", low = "yellow", breaks = c(0.4, 0.8, 1.2)) +
  scale_size(range = c(5,1)) + # setting the point size
  labs(x = "Longitude", y = "Latitude", color = "Blood THg\n(µg/g ww)") + 
  # ggtitle("Saltmarsh Sparrow") +
  coord_quickmap(xlim = c(-110, -60), ylim = c(0, 30)) +
  theme_classic() +
  coord_quickmap(xlim = c(-77, -68.5), ylim = c(37.5, 44.5)) +
  theme(plot.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12)) +
  guides(size = F) + # removing size legend
  annotate("text", label = "a)", x = -76.75, y = 44.5, size = 5, color = "black", face = "bold")

#Seaside
SESPData <- SESPData %>%
  arrange(HgConcentration) # this is so that the points will be plotted big -> small

SESPmap.man <- ggplot() + 
  geom_polygon(data = states, aes(x = long, y = lat, group = group), color = "black", size = 0.2, fill = "gray") +
  geom_point(data = SESPData, mapping = aes(x = Longitude, y = Latitude, color = HgConcentration, size = HgConcentration)) +
  scale_color_gradient(high = "red", low = "yellow", breaks = c(0.4, 0.8, 1.2)) +
  scale_size(range = c(5,1)) + # setting the point size
  labs(x = "Longitude", y = "Latitude", color = "Blood THg\n(µg/g ww)") + 
  # ggtitle("Seaside Sparrow") +
  coord_quickmap(xlim = c(-110, -60), ylim = c(0, 30)) +
  theme_classic() +
  coord_quickmap(xlim = c(-77, -68.5), ylim = c(37.5, 44.5)) +
  theme(plot.title = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12)) +
  guides(size = F) + # removing size legend
  annotate("text", label = "b)", x = -76.75, y = 44.5, size = 5, color = "black")

ggarrange(SALSmap.man, SESPmap.man, heights = c(1, 1), nrow = 1, ncol = 2,
          common.legend = TRUE, legend = "right")

ggsave("Sayers_DotMap.jpg", dpi = 750, width = 174, height = 116, units = "mm")



#violin plot
ggplot(data = HgSamples, mapping = aes(x = Species, y = HgConcentration, fill = Species)) +
  geom_violin(size = 0.75) +
  geom_boxplot(width = 0.1, fill = "white", size = .75) +
  stat_summary(fun = mean, geom = "point", fill = "black", shape = 8, size = 3) +
  scale_fill_manual(values = c("#E69F00", "#548235")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        axis.text.y = element_text(size = 12, color = "black")) +
  labs(x = "Species", y = "THg Concentration (µg/g ww)") +
  scale_x_discrete(labels = c("Saltmarsh Sparrow", "Seaside Sparrow")) +
  geom_hline(yintercept = c(0.7, 1.2), size = .75, linetype = "dashed", color = c("darkorange", "red3")) + 
  scale_y_continuous(limits = c(0.0, 1.6))

ggsave("Sayers_ViolinPlot.jpg", dpi = 750, width = 174, height = 116, units = "mm")
