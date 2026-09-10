#################################################################################
#################################################################################
# Title: Ancestry-Associated Gene Expression Analysis
# Author: Rashi Verma
# Date: 2025-10-08
# Description:
#   - Correlate ancestry estimates with gene expression
#   - Identify genes associated with African ancestry components
#   - Correct for multiple testing
#
# Input:
#   - Gene expression matrix
#   - Ancestry estimates
#
# Output:
#   - Ancestry-associated genes
#   - Correlation statistics
#   - Multiple-testing adjusted results
#################################################################################
#################################################################################

# Load required libraries
library(dplyr)
library(ggplot2)
library(patchwork)
library(edgeR)
library(limma)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Step 1: Load the ancestry proportion and expression data
ancestry <- read.csv("visit1_subancestry.csv", header = TRUE)
cpm <- read.csv("logcpm.csv")

## Check and remove duplicates in first column
duplicates <- duplicated(cpm[, 1])
cpm <- cpm[!duplicated(cpm[, 1]), ]

## Remove "X" prefix from column names (if present)
colnames(cpm) <- sub("^X", "", colnames(cpm))

## Match ancestry samples with CPM samples
matched_samples <- intersect(colnames(cpm)[-1], ancestry$Individual)

## Subset CPM and ancestry data
colnames(cpm)[1] <- "gene_id"
cpm <- cpm[, c("gene_id", matched_samples)]
ancestry <- ancestry %>% filter(Individual %in% matched_samples)

## Transpose CPM for correlation analysis
rownames(cpm) <- cpm$gene_id
cpm <- cpm[ , -1]
cpm <- t(cpm)

## Align ancestry with CPM
rownames(ancestry) <- ancestry$Individual
ancestry <- ancestry[rownames(cpm), ]
stopifnot(all(rownames(expr)==rownames(ancestry)))

# Step 2: correlation for sub-ancestry
## Initialize results dataframe
correlation_results <- data.frame(
  gene_id = colnames(cpm), 
  WestAfrican_cor = NA, WestAfrican_p = NA, 
  EastAfrican_cor = NA, EastAfrican_p = NA
)

## Calculate Spearman correlation and p-values
for (gene in colnames(cpm)) {
  test_west <- cor.test(cpm[, gene], ancestry$WestAfrican, method = "spearman")
  test_east <- cor.test(cpm[, gene], ancestry$EastAfrican, method = "spearman")
  
  correlation_results[correlation_results$gene_id == gene, "WestAfrican_cor"] <- test_west$estimate
  correlation_results[correlation_results$gene_id == gene, "WestAfrican_p"] <- test_west$p.value
  correlation_results[correlation_results$gene_id == gene, "EastAfrican_cor"] <- test_east$estimate
  correlation_results[correlation_results$gene_id == gene, "EastAfrican_p"] <- test_east$p.value
}

write.csv(correlation_results, "correlation_results.csv", row.names = FALSE)

## Adjust p-values using BH method
correlation_results <- correlation_results %>%
  mutate(
    WestAfrican_adj_p = p.adjust(WestAfrican_p, method = "BH"),
    EastAfrican_adj_p = p.adjust(EastAfrican_p, method = "BH")
  )

## POSITIVELY correlated genes
significant_west_pos <- correlation_results %>%
  filter(WestAfrican_adj_p < 0.05) %>%
  arrange(desc(WestAfrican_cor)) %>%
  head(50)

significant_east_pos <- correlation_results %>%
  filter(EastAfrican_adj_p < 0.05) %>%
  arrange(desc(EastAfrican_cor)) %>%
  head(50)

## Plot positive West African
plot_west_pos_dot <- ggplot(significant_west_pos, aes(x = reorder(gene_id, WestAfrican_cor), y = WestAfrican_cor, color = WestAfrican_cor,  size = -log10(WestAfrican_adj_p))) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "#FFE4B5", high = "#FF8C00") +
  coord_flip() +
  labs(
    title = "Positively correlated (West African)",
    x = "Genes", 
    y = "Correlation",
    color = "Correlation",
    size = "-log10(FDR)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.2, face = "bold"),
    axis.text.y = element_text(size = 4, face = "bold"),
    axis.text.x = element_text(size = 4, face = "bold"),
    axis.title.y = element_text(size = 7),
    axis.title.x = element_text(size = 7)
  )

## Plot positive East African
plot_east_pos_dot <- ggplot(significant_east_pos, aes(x = reorder(gene_id, EastAfrican_cor), y = EastAfrican_cor, color = EastAfrican_cor,  size = -log10(EastAfrican_adj_p))) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "#ADD8E6", high = "#00008B") +
  coord_flip() +
  labs(
    title = "Positively correlated (East African)",
    x = "Genes", 
    y = "Correlation",
    color = "Correlation",
    size = "-log10(FDR)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.2, face = "bold"),
    axis.text.y = element_text(size = 4, face = "bold"),
    axis.text.x = element_text(size = 4, face = "bold"),
    axis.title.y = element_text(size = 7),
    axis.title.x = element_text(size = 7)
  )

## NEGATIVELY correlated genes
significant_west_neg <- correlation_results %>%
  filter(WestAfrican_adj_p < 0.05) %>%
  arrange(WestAfrican_cor) %>%
  head(50)

significant_east_neg <- correlation_results %>%
  filter(EastAfrican_adj_p < 0.05) %>%
  arrange(EastAfrican_cor) %>%
  head(50)

## Plot negative West African
plot_west_neg_dot <- ggplot(significant_west_neg, aes(x = reorder(gene_id, WestAfrican_cor), y = WestAfrican_cor, color = WestAfrican_cor,  size = -log10(WestAfrican_adj_p))) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "#FFE4B5", high = "#FF8C00") +
  coord_flip() +
  labs(
    title = "Negetively correlated (West African)",
    x = "Genes", 
    y = "Correlation",
    color = "Correlation",
    size = "-log10(Adj.p.val)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.2, face = "bold"),
    axis.text.y = element_text(size = 4, face = "bold"),
    axis.text.x = element_text(size = 4, face = "bold"),
    axis.title.y = element_text(size = 7),
    axis.title.x = element_text(size = 7)
  )

## Plot negative East African
plot_east_neg_dot <- ggplot(significant_east_neg, aes(x = reorder(gene_id, EastAfrican_cor), y = EastAfrican_cor, color = EastAfrican_cor,  size = -log10(EastAfrican_adj_p))) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "#ADD8E6", high = "#00008B") +
  coord_flip() +
  labs(
    title = "Negetively correlated (East African)",
    x = "Genes", 
    y = "Correlation",
    color = "Correlation",
    size = "-log10(FDR)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.2, face = "bold"),
    axis.text.y = element_text(size = 4, face = "bold"),
    axis.text.x = element_text(size = 4, face = "bold"),
    axis.title.y = element_text(size = 7),
    axis.title.x = element_text(size = 7)
  )

## Combine all plots
combined_positive_plot <- plot_west_pos_dot / plot_east_pos_dot + 
  plot_annotation(title = "Top 50 Positively Correlated Genes with African Ancestry Fractions",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold")))

combined_negative_plot <- plot_west_neg_dot / plot_east_neg_dot + 
  plot_annotation(title = "Top 50 Negatively Correlated Genes with African Ancestry Fractions",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold")))

## Print plots
print(combined_positive_plot)
print(combined_negative_plot)

## Combine all four plots in a 2x2 grid
all_plots <- (plot_west_pos_dot | plot_east_pos_dot) / 
  (plot_west_neg_dot | plot_east_neg_dot) 

## Print combined plot
print(all_plots)

tiff(file="correlation.tiff", unit= "in", res = 600, width = 8, height = 8)
pdf(file = "correlation.pdf", width = 8, height = 8)
(plot_west_pos_dot | plot_east_pos_dot) / 
  (plot_west_neg_dot | plot_east_neg_dot)
dev.off()


total_genes <- nrow(correlation_results)

# Number of positively correlated genes (significant and correlation > 0)
positive_west_count <- correlation_results %>%
  filter(WestAfrican_adj_p < 0.05 & WestAfrican_cor > 0) %>%
  nrow()

positive_east_count <- correlation_results %>%
  filter(EastAfrican_adj_p < 0.05 & EastAfrican_cor > 0) %>%
  nrow()

# Print results
cat("Total number of genes tested:", total_genes, "\n")
cat("Number of positively correlated genes with West African ancestry:", positive_west_count, "\n")
cat("Number of positively correlated genes with East African ancestry:", positive_east_count, "\n")