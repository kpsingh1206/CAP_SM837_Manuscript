Single-cell RNA-seq analysis of CAP and SM837 in A375 melanoma cells
This repository contains the R analysis script accompanying a study of direct and indirect cold atmospheric plasma (CAP) treatment, the chromone SM837, and their combinations in A375 melanoma cells. The manuscript describes metabolic-activity experiments and single-cell transcriptomic profiling 48 hours after treatment.
Study overview
The study compares untreated control cells, direct CAP, indirect CAP using plasma-treated medium (PTM), SM837 alone, and SM837 combined with direct or indirect CAP. The manuscript reports CAP IC50 exposure times of 6 seconds for direct CAP and 13 seconds for PTM, and an SM837 IC50 of 83 µM. The scRNA-seq experiment used the IC50 CAP exposures, with or without 20 µM SM837. Combination-treatment synergy was assessed using the Bliss independence model.
Repository contents
- [Code_28_08_2026.R](Code_28_08_2026.R) — end-to-end R workflow for quality control, metadata harmonization, Seurat integration and clustering, cell-cycle analysis, cluster composition, marker analysis, pseudobulk differential expression, pathway enrichment, and correlation analyses.
- README.md — project overview and instructions.
The analysis script creates these output folders under the project output directory: 01_QC, 02_Objects, 03_Integration, 04_Figures, 05_Composition, 06_CellCycle, 07_ClusterMarkers, 08_DEG_Pseudobulk, 09_DEG_Pseudobulk_ByCluster, 10_Pathway, 11_Correlation, and 12_Logs.
Conditions included in the analysis
The script analyzes:
- DMEM control
- Direct CAP
- Indirect CAP (PTM)
- SM837 alone
- SM837 plus direct CAP
- SM837 plus indirect CAP
The two “CAP from below” conditions (CAP_von_unten and SM837+CAP_von_unten) are excluded from all manuscript analyses by the current script.
Input data
The script starts from three raw Seurat objects and does not use previously cleaned or integrated objects:
- _1_Lea1_Seurat.rds
- Lea2_Seurat.rds
- Lea3_Seurat.rds
Each object must contain an RNA assay and Sample_Name metadata. These data files are not included in this public repository. Update the SOURCE_DIR setting near the beginning of the script to point to the directory containing the three input files. The script creates the output directory manuscript_27.08.2026 within that location.
The manuscript draft states that the datasets are available from the corresponding author on reasonable request; it does not provide a public data accession.
Analysis and statistical notes
The script uses 2,000 variable features, 30 integration dimensions, 30 clustering dimensions, a clustering resolution of 0.5, and random seed 123. The integrated assay is used for dimensional reduction, clustering, and UMAP. RNA counts are used for expression analyses.
The primary treatment differential-expression analysis is replicate-aware pseudobulk DESeq2 with the model:
~ replicate_id + condition_id
The planned contrasts compare each treatment with the DMEM control and compare the two combination treatments. Differential-expression results are filtered at FDR < 0.05 and absolute log2 fold change > 0.25. The script treats Lea1, Lea2, and Lea3 as experimental runs. It notes that SM837+CAP_indirekt is absent from Lea2, so comparisons involving this condition use two runs and are exploratory.
Running the analysis
Use an R installation with the required packages available. The script checks for packages and stops if any are missing; it does not install them automatically. Install CRAN packages with install.packages() and Bioconductor packages with BiocManager::install() as appropriate.
Required packages:
Seurat, SeuratObject, dplyr, ggplot2, pheatmap,
DESeq2, clusterProfiler, ReactomePA, org.Hs.eg.db, enrichplot
After updating SOURCE_DIR and confirming that the three input files are present, run the script from R or RStudio:
source("Code_28_08_2026.R")
The script writes package-version information and analysis outputs to the project output directory. Review the generated logs and QC summaries before interpreting downstream results.
Citation and license
The supplied manuscript is a draft and does not contain a finalized title or publication citation. Add the final citation here when available. No reuse license is currently specified for this repository.
