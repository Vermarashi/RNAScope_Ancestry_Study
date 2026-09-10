###########################################################################################
###########################################################################################
# Title: Multivariable West African Ancestry–Expression Analysis
# Author: Rashi Verma
# Date: 2026-07-25
# Description:
#   - Evaluate associations between West African ancestry and gene expression
#   - Adjust for demographic and technical covariates
#   - Identify independently associated genes
#
# Input:
#   - Gene expression matrix
#   - West African ancestry estimates
#   - Age, sex, and batch information
#
# Output:
#   - Multivariable association results
#   - Significant ancestry-associated genes
###########################################################################################
###########################################################################################

# Load Libraries
library(dplyr)
library(limma)
library(ggplot2)

# Load data
cpm <- read.csv("logcpm.csv", check.names = FALSE)                                   
ancestry <- read.csv("visit2_subancestry.csv", header = TRUE)  

# Match samples
cpm <- cpm[!duplicated(cpm[,1]),]
colnames(cpm) <- sub("^X", "", colnames(cpm))
ancestry$Individual <-as.character(ancestry$Individual)
matched_samples <- intersect(colnames(cpm)[-1], ancestry$Individual)
cat("Matched samples:", length(matched_samples),"\n")
colnames(cpm)[1] <- "gene_id"
cpm <- cpm[,c("gene_id", matched_samples)]
ancestry <- ancestry %>% filter(Individual %in% matched_samples)

# Expression matrix
rownames(cpm) <- cpm$gene_id
cpm <- cpm[, -1]
expr <- t(cpm)
expr <- as.matrix(expr)
rownames(ancestry) <-ancestry$Individual
ancestry <-ancestry[rownames(expr),]
stopifnot(all(rownames(expr) == rownames(ancestry)))

# Covariates
## Age
ancestry$Age <-as.numeric(ancestry$Age)
ancestry$Age[is.na(ancestry$Age)] <-median(ancestry$Age,na.rm = TRUE)

## Sex
ancestry$Sex <- factor(ancestry$male, levels = c(0,1), labels = c("Female","Male"))

## Batch
ancestry$Batch <-factor(ancestry$Batch)

# Remove missing values
keep <- complete.cases(ancestry[,c("WestAfrican","Age","Sex","Batch")])
ancestry <-ancestry[keep,]
expr <-expr[rownames(ancestry),]
stopifnot(all(rownames(expr)==rownames(ancestry)))

# Check West/East ancestry relationship
cat("West-East correlation:", cor(ancestry$WestAfrican, ancestry$EastAfrican, method="spearman"),"\n")
summary(ancestry$WestAfrican)

# Linear regression model
design <- model.matrix(~ WestAfrican + Age + Sex + Batch, data = ancestry)
colnames(design)

# limma model
fit <- lmFit(t(expr), design)
fit <- eBayes(fit)

# West African ancestry association
west_results <- topTable(fit, coef = "WestAfrican", number = Inf, adjust.method = "BH")
write.csv(west_results,"WestAfrican_linear_regression_results.csv", row.names = TRUE)

# Significant genes
west_sig <- west_results %>% filter(adj.P.Val < 0.05)
west_sig

# Direction of association
## Positive West association is equivalent to negative East association
## Negative West association is Equivalent to positive East association
west_positive <- west_sig %>% filter(logFC > 0)
west_negative <- west_sig %>% filter(logFC < 0)
dim(west_positive)
dim(west_negative)

# Top 50 genes plot
top50 <- west_results %>% filter(adj.P.Val < 0.05) %>% arrange(adj.P.Val) %>% head(50)
top50$Gene <-rownames(top50)

ggplot(top50,aes(x = logFC, y = reorder(Gene, logFC))) +
  geom_point(aes(size = -log10(adj.P.Val))) +
  theme_bw() +
  labs(
    x = "West African ancestry association (β)",
    y = "Gene",
    size = "-log10(FDR)",
    title = "Top genes associated with West African ancestry"
  )

# Save positive/negative gene lists
write.csv(west_positive, "WestAfrican_positive_genes.csv")
write.csv(west_negative, "WestAfrican_negative_genes.csv")