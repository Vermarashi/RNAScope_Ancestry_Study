#################################################################################
#################################################################################
# Title: ADMIXTURE-Based Ancestry Benchmarking			
# Author: Rashi Verma						
# Date: 2026-08-02						
# Description:							
#   - Estimate ancestry proportions using ADMIXTURE		
#   - Benchmark ADMIXTURE estimates against RNAScope-Ancestry	
#   - Generate ancestry proportion estimates for comparison	
#								
# Input:							
#   - Quality-controlled genotype data				
#								
# Output:							
#   - ADMIXTURE ancestry proportions				
#   - Benchmarking results					
#################################################################################
#################################################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pophelper")
library(pophelper)
library(ggplot2)
library(ggthemes) 


# Load the K=3 Q-file
## Ensure your working directory contains the file
q_data <- read.table("merged_relaxed.3.Q", header=FALSE)
colnames(q_data) <- c("Q1", "Q2", "Q3")

# Add your labels from your labels object
labels <- read.table("update_ids_m.txt", header=FALSE)
q_data$Population <- labels$V4

# Run PCA
## We perform PCA on the ancestry proportions
pca_res <- prcomp(q_data[, 1:3], scale. = TRUE)

## Prepare data for plotting
# Create a new dataframe with the PC1 and PC2 coordinates
pca_df <- data.frame(PC1 = pca_res$x[,1], 
                     PC2 = pca_res$x[,2], 
                     Population = q_data$Population)


pub_theme <- theme_classic() +
  theme(
    text = element_text(size = 14, family = "sans"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    # Removing both major and minor grid lines
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
  )

my_colors <- c(
  "AFRICAN" = "darkblue",   
  "EUROPEAN" = "lightgreen",  
  "ADMIXED" = "#377EB8",   
  "AMERICAN" = "yellow",
  "Query1" = "hotpink",    
  "Query2" = "#FF7F00"     
)

outlines <- c(
  "AFRICAN"  = "NA", 
  "EUROPEAN" = "NA", 
  "ADMIXED"  = "NA", 
  "AMERICAN" = "black",        
  "Query1"   = "NA", 
  "Query2"   = "NA"
)

## Generating the plot with uniform round points (shape 21)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Population, fill = Population)) +
  # shape = 21 allows for both color (border) and fill (inside)
  geom_point(size = 1.0, alpha = 0.7, stroke = 0.2, shape = 21) + 
  scale_color_manual(values = outlines) + 
  scale_fill_manual(values = my_colors) + # Ensuring fill matches the border color
  pub_theme +
  labs(
    title = "PCA of Ancestral Components (K=4)",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 1), "%)")
  )