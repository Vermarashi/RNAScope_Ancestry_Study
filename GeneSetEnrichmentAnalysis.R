###########################################################################################
###########################################################################################
# Title: Gene Set Enrichment Analysis of Ancestry-Associated Genes
# Author: Rashi Verma
# Date: 2026-07-25
# Description:
#   - Perform gene-set enrichment analysis using ancestry-associated genes
#   - Identify enriched biological processes and molecular functions
#   - Generate enrichment plots and summary results
#
# Input:
#   - Ancestry-associated gene statistics
#
# Output:
#   - Enriched gene sets/pathways
#   - Enrichment statistics
#   - Visualization plots
###########################################################################################
###########################################################################################

# Run GSEA
## Load Libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

west_results <- read.csv("west_results", check.names = FALSE) 
west_results$Gene <- rownames(west_results)

gene_map <- bitr(west_results$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb =org.Hs.eg.db)
west_ranked <- merge(west_results, gene_map, by.x = "Gene", by.y = "SYMBOL")

## Remove duplicated Entrez IDs
west_ranked <- west_ranked[!duplicated(west_ranked$ENTREZID), ]
geneList <- west_ranked$logFC
names(geneList) <- west_ranked$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

# BP
gsea_west <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

#CC
gsea_west_CC <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

#MF
gsea_west_MF <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

p_BP <- dotplot(
  gsea_west_BP,
  showCategory = 15,
  title = "Biological Process"
) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    plot.title   = element_text(size = 11, face = "bold", color = "black")
  )

p_CC <- dotplot(
  gsea_west_CC,
  showCategory = 15,
  title = "Cellular Component"
) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    plot.title   = element_text(size = 11, face = "bold", color = "black")
  )

p_MF <- dotplot(
  gsea_west_MF,
  showCategory = 15,
  title = "Molecular Function"
) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 12, color = "black", ),
    axis.text.y  = element_text(size = 12, color = "black"),
    plot.title   = element_text(size = 11, face = "bold", color = "black")
  )

combined_GSEA <- p_BP | p_CC | p_MF
print(combined_GSEA)
