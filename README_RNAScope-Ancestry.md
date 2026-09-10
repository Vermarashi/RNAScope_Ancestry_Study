###############################################################################
###############################################################################
# Title: RNAScope-Ancestry: MECA RNA-seq Ancestry Inference Workflow
# Author: Rashi Verma
# Date: 2026-09-08
# Description:
#   - Infer genetic ancestry from short-read RNA-seq data
#   - Process RNA-seq reads and identify variants from expressed regions
#   - Perform SNP quality control, filtering, and linkage disequilibrium pruning
#   - Infer ancestry using PCA and Rye
#   - Benchmark ancestry estimates using ADMIXTURE
#   - Evaluate concordance between ancestry inference approaches
#   - Quantify RNA-seq gene expression
#   - Identify ancestry-associated gene expression patterns
#   - Perform gene-set enrichment analysis
#
# Study:
#   - Morehouse-Emory Cardiovascular Center for Health Equity (MECA)
#
# Input:
#   - Short-read RNA-seq data
#   - GRCh38 reference genome
#   - 1000 Genomes Project reference populations
#
# Output:
#   - Quality-controlled SNP datasets
#   - PCA-based ancestry estimates
#   - Rye ancestry estimates
#   - ADMIXTURE ancestry estimates
#   - Ancestry-associated gene-expression results
#   - Gene-set enrichment results
###############################################################################
###############################################################################

## Overview

This repository contains the computational workflow used to infer genetic ancestry from short-read RNA-seq data from the Morehouse-Emory Cardiovascular Center for Health Equity (MECA) study and to evaluate ancestry-associated gene expression.

The workflow includes:

1. Extraction of unmapped RNA-seq reads
2. FASTQ generation and adapter/quality trimming
3. STAR/Bowtie2 alignment
4. BAM processing and quality control
5. GATK RNA-seq variant calling
6. Variant processing and ancestry-specific SNP selection
7. PCA-based ancestry inference
8. ADMIXTURE-based benchmarking
9. Comparison of ancestry estimates
10. RNA-seq expression quantification
11. Ancestry–gene-expression correlation analysis
12. Gene-set enrichment analysis

---

## Repository structure

```text
RNAScope-Ancestry/
│
├── README.md
│
├── scripts/
│   │
│   ├── 01_preprocessing/
│   │   ├── 00_EXTRACT_UNMAPPED.sh
│   │   ├── 01_Bam2Fastq.sh
│   │   └── 02_TrimGalore.sh
│   │
│   ├── 02_alignment/
│   │   └── 03A_STAR-Bowtie2_GRCH.sh
│   │
│   ├── 03_bam_processing/
│   │   ├── 04_SAMTOOLS.sh
│   │   └── 05_Samstats.sh
│   │
│   ├── 04_variant_calling/
│   │   ├── 05A_GATK-1.sh
│   │   ├── 05B_GATK-2.sh
│   │   ├── 05C_GATK-3.sh
│   │   ├── 05D_GATK-4.sh
│   │   └── 05E_GATK-5.sh
│   │
│   ├── 05_ancestry/
│   │   ├── 06A_meca_ancestry_workflow.txt
│   │   ├── 06B_meca_sub_ancestry_workflow.txt
│   │   ├── 06C_common_variants.py
│   │   ├── 06D_update_famID.py
│   │   └── 06E_update_ids.txt
│   │
│   ├── 06_pca_admixture/
│   │   ├── 07A_PCA_analysis_RNAScope.R
│   │   ├── 07B_PCA_ADMIXTURE.R
│   │   └── 08_comparison.R
│   │
│   ├── 07_expression/
│   │   ├── 09A_Stringtie2countMatrix_1.sh
│   │   └── 09B_prepDE.py
│   │
│   └── 08_expression_ancestry/
│       ├── 10A_correlation.R
│       ├── 10B_west_covar_corr.R
│       └── 10C_GeneSetEnrichmentAnalysis.R
│
└── results/
    ├── ancestry/
    ├── pca/
    ├── admixture/
    └── enrichment/
```

> The exact folder names can be changed, but the order above reflects the analytical logic of the workflow rather than simply the original script numbering.

---

# Workflow

## Step 1. RNA-seq preprocessing

### `00_EXTRACT_UNMAPPED.sh`

Extracts unmapped reads from the initial alignment/BAM files for downstream processing.

### `01_Bam2Fastq.sh`

Converts BAM files to FASTQ files for subsequent read processing.

### `02_TrimGalore.sh`

Performs adapter and quality trimming and generates trimmed FASTQ files. Trim Galore is also used to assess read quality through FastQC.

---

## Step 2. RNA-seq alignment

### `03A_STAR-Bowtie2_GRCH.sh`

Aligns trimmed RNA-seq reads to the GRCh38 reference genome using STAR/Bowtie2.

The output is a coordinate-sorted BAM file for each sample.

Conceptually:

```text
R1 + R2 FASTQ
      ↓
STAR/Bowtie2
      ↓
sample BAM
```

---

## Step 3. BAM processing and QC

### `04_SAMTOOLS.sh`

Performs BAM-level processing such as merging, sorting, and indexing as required by the alignment workflow.

### `05_Samstats.sh`

Generates alignment/chromosome-level statistics using `samtools idxstats`.

This script is primarily a QC/statistics step and is not a core ancestry-inference step.

---

# Step 4. RNA-seq variant calling

The GATK scripts follow the RNA-seq variant-calling workflow.

### `05A_GATK-1.sh`

Marks PCR/optical duplicates using GATK MarkDuplicates.

```text
Aligned BAM
    ↓
MarkDuplicates
    ↓
marked-duplicates BAM
```

### `05B_GATK-2.sh`

Processes the RNA-seq BAM using:

- SplitNCigarReads
- BaseRecalibrator
- ApplyBQSR

The final output is a BQSR-applied BAM suitable for variant calling.

```text
marked BAM
    ↓
SplitNCigarReads
    ↓
BaseRecalibrator
    ↓
ApplyBQSR
    ↓
BQSR-applied BAM
```

### `05CC_GATK-33.sh`

Runs GATK HaplotypeCaller in GVCF mode to generate variant information for individual samples.

```text
one sample BAM
      ↓
HaplotypeCaller -ERC GVCF
      ↓
one sample GVCF
```

### `05D_GATK-4.sh`

Performs cohort-level genotyping using:

- GenomicsDBImport
- GenotypeGVCFs

Conceptually:

```text
Sample1.g.vcf
Sample2.g.vcf
Sample3.g.vcf
       ...
SampleN.g.vcf
       ↓
GenomicsDBImport
       ↓
GenotypeGVCFs
       ↓
one cohort-level VCF
```

The cohort VCF contains genotype information for all included MECA samples.

### `05E_GATK-5.sh`

Applies GATK variant-level filtration to the cohort VCF.

The resulting filtered VCF is used for downstream SNP quality control and ancestry analysis.

---

# Step 5. Ancestry SNP processing

### `06A_meca_ancestry_workflow.txt`

Main MECA ancestry-processing workflow.

This stage takes the filtered variant data and prepares the SNP dataset for ancestry inference.

### `06B_meca_sub_ancestry_workflow.txt`

Performs the sub-ancestry analysis.

### `06C_common_variants.py`

Identifies common/shared variants required for comparison between MECA RNA-seq variants and reference populations. This file needed while running 06A_meca_ancestry_workflow.txt and 06B_meca_sub_ancestry_workflow.txt for identifying common variants.

### `06D_update_famID.py`

Updates PLINK family/sample identifiers as required for downstream analysis. This file needed while running 06A_meca_ancestry_workflow.txt and 06B_meca_sub_ancestry_workflow.txt for updating fam ID. This step is optional. 

### `06E_update_ids.txt`

Contains commands for updating sample identifiers. This file s needed for updating the famIDs.

---

# Step 6. PCA and ancestry inference

### `07A_PCA_analysis_RNAScope.R`

Performs PCA-based ancestry analysis using the RNAScope-Ancestry workflow.

The PCA results are used to evaluate clustering of MECA samples relative to reference populations.

### `07B_PCA_ADMIXTURE.R`

Performs ADMIXTURE-based ancestry estimation for benchmarking against RNAScope-Ancestry.

### `08_comparison.R`

Compares ancestry estimates obtained using the different approaches.

This analysis is used to evaluate concordance between RNAScope-Ancestry and ADMIXTURE-based estimates.

---

# Step 7. RNA-seq expression quantification

### `09A_Stringtie2countMatrix_1.sh`

Processes StringTie transcript-level results to generate expression/count information.

### `09B_prepDE.py`

Prepares the StringTie results for downstream differential-expression/expression analyses.

---

# Step 8. Ancestry-associated gene expression

### `10A_correlation.R`

Tests associations between ancestry estimates and gene expression.

### `10B_west_covar_corr.R`

Performs the multivariable analysis of West African ancestry and gene expression while accounting for the specified demographic/technical covariates.

### `10C_GeneSetEnrichmentAnalysis.R`

Performs gene-set/pathway enrichment analysis of ancestry-associated genes.

---

# Recommended manuscript-level workflow

For the manuscript, the computational workflow can be summarized as:

```text
RNA-seq reads
     ↓
Unmapped read extraction
     ↓
BAM → FASTQ
     ↓
Trim Galore
     ↓
STAR/Bowtie2 alignment
     ↓
Coordinate-sorted BAM
     ↓
MarkDuplicates
     ↓
SplitNCigarReads
     ↓
Base Quality Score Recalibration
     ↓
HaplotypeCaller
     ↓
Per-sample GVCFs
     ↓
Joint genotyping
     ↓
Cohort VCF
     ↓
Variant filtration
     ↓
Common/shared SNP selection
     ↓
MAF / missingness / HWE / LD filtering
     ↓
PCA + Rye ancestry inference
     ↓
ADMIXTURE benchmarking
     ↓
Ancestry estimates
     ↓
Validation using MAGE samples (DNA samples present in 1KGP
     ↓
Gene-expression association
     ↓
Gene-set enrichment
```

---

# Reproducibility

Before publication, add a `software_versions.txt` file containing the versions of all major software used, including where applicable:

- Trim Galore v 0.6.10
- FastQC v v0.12.1
- STAR v 2.7.11b
- Bowtie2 v 2.5.2
- SAMtools v 1.22.1
- GATK v 4.6.2.0
- PLINK v 90b7.2
- Rye v 0.1
- ADMIXTURE v 1.3.0
- StringTie v 2.2.1
- Python v 3.12.3
- R v 4.4.2

---

# Citation

If you use this workflow or RNAScope-Ancestry, please cite the associated manuscript:

> Verma R, et al. RNAScope-Ancestry: ancestry inference from RNA-seq data. [Manuscript citation to be added after publication.]

