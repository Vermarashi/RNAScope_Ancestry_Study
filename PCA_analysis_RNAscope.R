#################################################################################
#################################################################################
# Title: RNAScope-Ancestry PCA Analysis						
# Author: Rashi Verma								
# Date: 2025-10-02								
# Description:									
#   - Perform principal component analysis using shared high-quality SNPs	
#   - Evaluate genetic ancestry clustering of MECA participants			
#   - Compare MECA samples with 1000 Genomes reference populations		
#										
# Input:									
#   - Quality-controlled genotype data						
#										
# Output:									
#   - PCA coordinates								
#   - PCA plots									
#   - Ancestry clustering results						
#################################################################################
#################################################################################

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read in the eigenvec file
data <- read.csv("merged_pca.eigenvec", sep = " ", header = FALSE)
colnames(data) <- c("FID", "IID", paste0("PC", 1:(ncol(data) - 2)))

# Define population categories (for ancestry)
data$Population <- recode(data$FID, 
                          "YRI" = "African", "GWD" = "African", "ESN" = "African", "LWK" = "African", "MSL" = "African",
                          "GBR" = "European", "IBS" = "European", "CEU" = "European", "FIN" = "European", "TSI" = "European", 
                          "ACB" = "Admixed", "ASW" = "Admixed", 
                          "PEL" = "American", "AFR" = "Patients")

# Define colors for populations
colors = c("European" = "palegreen", "African" = "#283593", "Admixed" = "lightskyblue", 
           "American" = "yellow", "Patients" = "hotpink")

# Define population categories (for sub-ancestry)
data$Population <- recode(data$FID, 
                          "YRI" = "Yoruba_West", "GWD" = "Senegambian_West", "ESN" = "Esan_West", "LWK" = "Luhya_East", "MSL" = "Mende_West",
                          "ACB" = "Admixed", "ASW" = "Admixed", 
                          "AFR" = "Patients")

colors = c("Yoruba_West" = "palegreen", "Senegambian_West" = "#283593", "Esan_West" = "#00F5FF", "Mende_West" = "lightskyblue",
           "Luhya_East" = "yellow", "Patients" = "hotpink", "Admixed" = "orange")

tiff(file="PCA_1_3.tiff", unit= "in", res = 600, width = 10, height = 6)
ggplot(data, aes(x = PC1, y = PC2, color = Population, fill = Population)) +
  geom_point(data = data[data$Population != "Patients", ], 
             size = 1, stroke = 0.5) +
  geom_point(data = data[data$Population == "American", ], 
             size = 1, shape = 21, color = "darkgray", fill = "yellow", show.legend = FALSE) +
  geom_point(data = data[data$Population == "Patients", ], 
             size = 1, stroke = 0.5, shape = 17) +
  # Customizing plot
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors, guide = "none") +  # Disable duplicate legend for fill
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PCA Plot of PC1 vs PC2") +
  theme_bw() +  # Apply white background theme
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_line(color = "darkgrey"),  # Minor grid lines in black
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 12),  # Increase the size of the axis numbers
    legend.title = element_text(size = 6), #24
    legend.text = element_text(size = 6), #24
    legend.position = "bottom") + guides(color = guide_legend(override.aes = list(size = 5)))

dev.off()

