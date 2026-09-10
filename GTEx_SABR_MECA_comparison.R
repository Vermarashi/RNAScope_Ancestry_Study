###############################################################################
###############################################################################
# Title: Comparison of MECA Ancestry-Associated Genes with SABR and GTEx eQTLs
# Author: Rashi Verma
# Date: 2026-09-08
# Description:
#   - Load gene-level results from SABR eQTL, GTEx eQTL, and MECA analyses
#   - Identify significant genes from each dataset
#   - Standardize gene identifiers and map Ensembl IDs to gene symbols
#   - Calculate pairwise and three-way overlaps among SABR, GTEx, and MECA gene sets
#   - Generate a Venn diagram to visualize gene-set overlaps
#   - Export overlapping gene lists and an overlap summary table
#
# Input:
#   - SABR_GTEX_MECA.csv
#   - SABR eQTL results
#   - GTEx eQTL results
#   - MECA ancestry-associated gene results
#
# Significance thresholds:
#   - SABR eQTL: q-value < 0.05
#   - GTEx eQTL: q-value < 0.05
#   - MECA: adjusted P-value < 0.05
#
# Output:
#   - MECA_SABR_GTEx_Venn.tiff
#   - SABR_GTEx_overlap_genes.csv
#   - SABR_MECA_overlap_genes.csv
#   - GTEx_MECA_overlap_genes.csv
#   - SABR_GTEx_MECA_three_way_overlap_genes.csv
#   - SABR_GTEx_MECA_overlap_summary.csv
#
# Software:
#   - R
#   - readr
#   - dplyr
#   - ggplot2
#   - ggVennDiagram
#   - AnnotationDbi
#   - org.Hs.eg.db
#   - VennDiagram
###############################################################################
###############################################################################

# Load packages
library(readr)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(VennDiagram)
library(grid)

# Read CSV file
dat <- read_csv(
  "SABR_GTEX_MECA.csv",
  na = c("", "NA")
)


# Rename columns based on the CSV structure
# SABR = columns 1-3
# column 4 = blank
# GTEx = columns 5-7
# columns 8-9 = blank
# MECA = columns 10-11
names(dat)[1]  <- "SABR_gene_id"
names(dat)[2]  <- "SABR_gene_symbol"
names(dat)[3]  <- "SABR_qval"

names(dat)[5]  <- "GTEx_gene_id"
names(dat)[6]  <- "GTEx_gene_name"
names(dat)[7]  <- "GTEx_qval"

names(dat)[10] <- "MECA_gene_symbol"
names(dat)[11] <- "MECA_adjP"

# Keep only required columns
dat <- dat %>%
  dplyr::select(
    SABR_gene_id,
    SABR_gene_symbol,
    SABR_qval,
    GTEx_gene_id,
    GTEx_gene_name,
    GTEx_qval,
    MECA_gene_symbol,
    MECA_adjP
  )


# Remove the first row
# The first row contains:
# gene_id / gene_symbol / qval / gene_id / gene_name /
# qval / gene_symbol / adj.P.Val
dat <- dat %>%
  dplyr::filter(SABR_gene_id != "gene_id")

# Make sure p-values are numeric
dat <- dat %>%
  dplyr::mutate(
    SABR_qval = as.numeric(SABR_qval),
    GTEx_qval = as.numeric(GTEx_qval),
    MECA_adjP = as.numeric(MECA_adjP)
  )

# Remove Ensembl version numbers
# Example:
# ENSG00000187608.10
# becomes
# ENSG00000187608
dat <- dat %>%
  dplyr::mutate(
    SABR_gene_id_clean =
      sub("\\..*$", "", SABR_gene_id),
    
    GTEx_gene_id_clean =
      sub("\\..*$", "", GTEx_gene_id)
  )

#  Select significant SABR eQTL genes
# q-value < 0.05
sabr_sig <- dat %>%
  dplyr::filter(
    !is.na(SABR_qval),
    SABR_qval < 0.05,
    !is.na(SABR_gene_id_clean),
    SABR_gene_id_clean != ""
  ) %>%
  dplyr::select(
    SABR_gene_id_clean,
    SABR_gene_symbol,
    SABR_qval
  ) %>%
  dplyr::distinct(
    SABR_gene_id_clean,
    .keep_all = TRUE
  )

# Select significant GTEx eGenes
# q-value < 0.05
gtex_sig <- dat %>%
  dplyr::filter(
    !is.na(GTEx_qval),
    GTEx_qval < 0.05,
    !is.na(GTEx_gene_id_clean),
    GTEx_gene_id_clean != ""
  ) %>%
  dplyr::select(
    GTEx_gene_id_clean,
    GTEx_gene_name,
    GTEx_qval
  ) %>%
  dplyr::distinct(
    GTEx_gene_id_clean,
    .keep_all = TRUE
  )

# Select significant MECA genes
# adjusted P-value < 0.05
meca_sig <- dat %>%
  dplyr::filter(
    !is.na(MECA_adjP),
    MECA_adjP < 0.05,
    !is.na(MECA_gene_symbol),
    MECA_gene_symbol != ""
  ) %>%
  dplyr::select(
    MECA_gene_symbol,
    MECA_adjP
  ) %>%
  dplyr::distinct(
    MECA_gene_symbol,
    .keep_all = TRUE
  )

# Map SABR Ensembl IDs to gene symbols
sabr_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(sabr_sig$SABR_gene_id_clean),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
)

sabr_map <- sabr_map %>%
  dplyr::filter(
    !is.na(ENSEMBL),
    !is.na(SYMBOL)
  ) %>%
  dplyr::distinct(
    ENSEMBL,
    .keep_all = TRUE
  )


# Add mapped SABR symbols
sabr_sig <- sabr_sig %>%
  dplyr::left_join(
    sabr_map,
    by = c(
      "SABR_gene_id_clean" = "ENSEMBL"
    )
  ) %>%
  dplyr::mutate(
    SABR_symbol_final = dplyr::coalesce(
      SYMBOL,
      SABR_gene_symbol
    )
  )

# Map GTEx Ensembl IDs to gene symbols
gtex_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(gtex_sig$GTEx_gene_id_clean),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
)

gtex_map <- gtex_map %>%
  dplyr::filter(
    !is.na(ENSEMBL),
    !is.na(SYMBOL)
  ) %>%
  dplyr::distinct(
    ENSEMBL,
    .keep_all = TRUE
  )

# Add mapped GTEx symbols
gtex_sig <- gtex_sig %>%
  dplyr::left_join(
    gtex_map,
    by = c(
      "GTEx_gene_id_clean" = "ENSEMBL"
    )
  )

# Create final gene sets
# Use gene SYMBOLS for the Venn diagram.
# This allows comparison with MECA, which is in symbols.
sabr_genes <- sabr_sig %>%
  dplyr::pull(SABR_symbol_final) %>%
  na.omit() %>%
  toupper() %>%
  trimws() %>%
  unique()

gtex_genes <- gtex_sig %>%
  dplyr::pull(SYMBOL) %>%
  na.omit() %>%
  toupper() %>%
  trimws() %>%
  unique()

meca_genes <- meca_sig %>%
  dplyr::pull(MECA_gene_symbol) %>%
  na.omit() %>%
  toupper() %>%
  trimws() %>%
  unique()

# Remove empty gene names
sabr_genes <- sabr_genes[sabr_genes != ""]
gtex_genes <- gtex_genes[gtex_genes != ""]
meca_genes <- meca_genes[meca_genes != ""]

# Calculate pairwise overlaps
sabr_gtex <- intersect(
  sabr_genes,
  gtex_genes
)

sabr_meca <- intersect(
  sabr_genes,
  meca_genes
)

gtex_meca <- intersect(
  gtex_genes,
  meca_genes
)

# Calculate three-way overlap
three_way <- Reduce(
  intersect,
  list(
    sabr_genes,
    gtex_genes,
    meca_genes
  )
)

print(three_way)

# Create Venn diagram
gene_sets <- list(
  SABR_eQTL = sabr_genes,
  GTEx_eQTL = gtex_genes,
  MECA      = meca_genes
)

venn_plot <- venn.diagram(
  x = gene_sets,
  
  category.names = c(
    "SABR eQTL",
    "GTEx eQTL",
    "MECA"
  ),
  
  filename = NULL,
  
  fill = c(
    "#C6DBEF",   # light blue
    "#CDECCF",   # light green
    "#F4C6D7"    # light pink
  ),
  
  alpha = 0.55,
  
  col = c(
    "#6B9EC4",
    "#75B57D",
    "#C17C99"
  ),
  
  lwd = 1.5,
  
  # Numbers inside Venn
  cex = 12 / 12,
  
  # Category labels
  cat.cex = 12 / 12,
  
  fontfamily = "Arial",
  cat.fontfamily = "Arial",
  
  cat.col = "black",
  
  cat.pos = c(-20, 20, 180),
  
  cat.dist = c(0.05, 0.05, 0.05),
  
  margin = 0.1
)

grid.newpage()
grid.draw(venn_plot)

# Save Venn diagram as high-resolution TIFF
tiff(filename = "MECA_SABR_GTEx_Venn.tiff", width = 8, height = 8, units = "in", res = 600, compression = "lzw")
print(venn_plot)
dev.off()

# Save pairwise overlap genes
write.csv(data.frame(gene_symbol = sabr_gtex), "SABR_GTEx_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(gene_symbol = sabr_meca), "SABR_MECA_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(gene_symbol = gtex_meca), "GTEx_MECA_overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(gene_symbol = three_way), "SABR_GTEx_MECA_three_way_overlap_genes.csv", row.names = FALSE)

# Create overlap summary table
overlap_summary <- data.frame(
  Comparison = c(
    "SABR eQTL ∩ GTEx eQTL",
    "SABR eQTL ∩ MECA",
    "GTEx eQTL ∩ MECA",
    "SABR eQTL ∩ GTEx eQTL ∩ MECA"
  ),
  
  Number_of_genes = c(
    length(sabr_gtex),
    length(sabr_meca),
    length(gtex_meca),
    length(three_way)
  )
)

print(overlap_summary)

write.csv(overlap_summary,"SABR_GTEx_MECA_overlap_summary.csv", row.names = FALSE)
