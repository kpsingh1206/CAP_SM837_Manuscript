# CAP and SM837 single-cell RNA-seq analysis in A375 melanoma cells

This repository contains the R workflow and selected analysis outputs for a study of cold atmospheric plasma (CAP), the chromone compound SM837, and their combination in A375 melanoma cells. The workflow analyzes single-cell RNA-seq data collected 48 hours after treatment. The associated study also evaluates metabolic activity and combination effects; synergy is assessed with the Bliss independence model.

## Contents

- [Study design](#study-design)
- [Repository structure](#repository-structure)
- [Data requirements](#data-requirements)
- [Software and installation](#software-and-installation)
- [Run the analysis](#run-the-analysis)
- [Analysis settings and interpretation](#analysis-settings-and-interpretation)
- [Outputs](#outputs)
- [Data availability, citation, and license](#data-availability-citation-and-license)

## Study design

The current analysis includes these treatment groups:

| Group | Description |
|---|---|
| `Kontrolle_DMEM` | DMEM control |
| `CAP_direkt` | Direct CAP treatment |
| `CAP_indirekt` | Indirect CAP treatment using plasma-treated medium (PTM) |
| `SM837_alleine` | SM837 alone |
| `SM837+CAP_direkt` | SM837 combined with direct CAP |
| `SM837+CAP_indirekt` | SM837 combined with indirect CAP |

The two CAP-from-below groups, `CAP_von_unten` and `SM837+CAP_von_unten`, are intentionally excluded from the manuscript analysis by the current script.

## Repository structure

| Path | Contents |
|---|---|
| [`scripts/Code_28_08_2026.R`](scripts/Code_28_08_2026.R) | Main R analysis pipeline |
| [`data/QC/`](data/QC/) | Selected pre-QC violin plots and count-versus-feature plots for Lea1–Lea3 |
| [`data/Integration/`](data/Integration/) | Selected integration summaries, PCA elbow plot, cluster counts, and feature information |
| `README.md` | Project overview and instructions |

The repository contains selected QC and integration artifacts. The full analysis outputs are written to the output directory configured by the script.

## Data requirements

The pipeline starts from three raw Seurat objects; it does not use previously cleaned or integrated objects as input:

- `_1_Lea1_Seurat.rds`
- `Lea2_Seurat.rds`
- `Lea3_Seurat.rds`

Each object must contain an RNA assay and a `Sample_Name` metadata field. These RDS files are not included in the repository.

At the top of the script, update `SOURCE_DIR` to the local directory containing all three files. The script defines `PROJECT_DIR` as a `manuscript_27.08.2026` subdirectory inside `SOURCE_DIR`, then creates the analysis output folders there. The path currently in the script is machine-specific.

## Software and installation

Use an R installation compatible with Seurat 5 and the required Bioconductor packages. The script checks for the following packages and stops if any are missing; it does not install packages automatically:

- CRAN: `Seurat`, `SeuratObject`, `dplyr`, `ggplot2`, `pheatmap`
- Bioconductor: `DESeq2`, `clusterProfiler`, `ReactomePA`, `org.Hs.eg.db`, `enrichplot`
- `BiocManager` is required to check the Bioconductor environment.

Install the packages in R:

```r
install.packages(c(
  "Seurat", "SeuratObject", "dplyr", "ggplot2", "pheatmap", "BiocManager"
))

BiocManager::install(c(
  "DESeq2", "clusterProfiler", "ReactomePA", "org.Hs.eg.db", "enrichplot"
))
```

The script checks that Seurat is version 5 or newer and records the R, Bioconductor, and required package versions in `package_versions_27.08.2026.csv` in the output directory.

## Run the analysis

1. Place the three RDS input files in one directory.
2. Edit `SOURCE_DIR` near the beginning of `scripts/Code_28_08_2026.R` so it points to that directory.
3. Install the required packages.
4. Open the repository as the R working directory and run:

```r
source("scripts/Code_28_08_2026.R")
```

Review the log files and QC summaries in the generated output before interpreting downstream results. The script creates its output under `SOURCE_DIR/manuscript_27.08.2026`, rather than writing the complete run into the repository.

## Analysis settings and interpretation

- The fixed random seed is `123`.
- The workflow uses 2,000 variable features, 30 dimensions for integration and clustering, and clustering resolution `0.5`.
- The integrated assay is used for dimensional reduction, clustering, and UMAP. RNA counts are used for expression analyses.
- The primary treatment differential-expression analysis is replicate-aware pseudobulk DESeq2 with the design `~ replicate_id + condition_id`.
- The script treats Lea1, Lea2, and Lea3 as experimental runs. It defines treatment-versus-DMEM and combination contrasts; differential-expression results use FDR < 0.05 and absolute log2 fold change > 0.25.
- `SM837+CAP_indirekt` is absent from Lea2. Comparisons involving this group therefore use two runs and should be treated as exploratory.
- The pipeline also performs QC, metadata harmonization, Seurat integration and clustering, cell-cycle analysis, cluster-composition analysis, cluster-marker analysis, pathway enrichment, and correlation analyses.

Cells from the same run are not independent biological replicates. Interpret cell-level plots as descriptive and use the replicate-aware pseudobulk results for treatment-level inference.

## Outputs

The script creates these folders under `PROJECT_DIR`:

| Folder | Main contents |
|---|---|
| `01_QC` | Quality-control summaries and plots |
| `02_Objects` | Processed Seurat objects |
| `03_Integration` | Integration and dimensional-reduction results |
| `04_Figures` | UMAP and other figures |
| `05_Composition` | Cluster-composition summaries |
| `06_CellCycle` | Cell-cycle analysis |
| `07_ClusterMarkers` | Cluster-marker results |
| `08_DEG_Pseudobulk` | Primary pseudobulk differential-expression results |
| `09_DEG_Pseudobulk_ByCluster` | Cluster-stratified pseudobulk results |
| `10_Pathway` | Pathway-enrichment results |
| `11_Correlation` | Correlation analyses |
| `12_Logs` | Run logs and diagnostics |

## Data availability, citation, and license

The three raw Seurat objects are not distributed in this repository. The project description states that the data are available from the corresponding author on reasonable request; update this statement if a public accession becomes available.

A final publication citation is not yet specified here. Add the finalized title, authors, journal, year, and DOI when available.

No `LICENSE` file is currently present. Until a license is added, this repository does not grant permission to reuse or redistribute its code or data.

