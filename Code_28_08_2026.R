# ==============================================================================
# LEA scRNA-seq MANUSCRIPT MASTER PIPELINE
# Date: 27.08.2026
# Output folder:
# F:/UMR_data_scRNA/seurat_all/manuscript_27.08.2026
#
# RAW INPUT FILES
#   _1_Lea1_Seurat.rds
#   Lea2_Seurat.rds
#   Lea3_Seurat.rds
#
# This script starts only from the three raw Seurat objects.
# Existing cleaned/integrated RDS files are NOT used as inputs.
#
# MAIN PRINCIPLE
# - Integrated assay: PCA, clustering, UMAP.
# - RNA counts: marker discovery and expression analyses.
# - Primary treatment DEG: replicate-aware pseudobulk DESeq2.
#
# STATISTICAL ASSUMPTION
# Lea1, Lea2 and Lea3 are treated as replicate experimental runs.
# Primary DEG model: ~ replicate_id + condition_id
# ==============================================================================


# ==============================================================================
# 0. PATHS
# ==============================================================================

SOURCE_DIR  <- "F:/UMR_data_scRNA/seurat_all"
PROJECT_DIR <- file.path(SOURCE_DIR, "manuscript_27.08.2026")

RAW_LEA1 <- file.path(SOURCE_DIR, "_1_Lea1_Seurat.rds")
RAW_LEA2 <- file.path(SOURCE_DIR, "Lea2_Seurat.rds")
RAW_LEA3 <- file.path(SOURCE_DIR, "Lea3_Seurat.rds")

stopifnot(dir.exists(SOURCE_DIR))
stopifnot(file.exists(RAW_LEA1))
stopifnot(file.exists(RAW_LEA2))
stopifnot(file.exists(RAW_LEA3))

dir.create(PROJECT_DIR, recursive = TRUE, showWarnings = FALSE)

DIR_QC          <- file.path(PROJECT_DIR, "01_QC")
DIR_OBJECTS     <- file.path(PROJECT_DIR, "02_Objects")
DIR_INTEGRATION <- file.path(PROJECT_DIR, "03_Integration")
DIR_FIGURES     <- file.path(PROJECT_DIR, "04_Figures")
DIR_COMPOSITION <- file.path(PROJECT_DIR, "05_Composition")
DIR_CELL_CYCLE  <- file.path(PROJECT_DIR, "06_CellCycle")
DIR_MARKERS     <- file.path(PROJECT_DIR, "07_ClusterMarkers")
DIR_DEG_PB      <- file.path(PROJECT_DIR, "08_DEG_Pseudobulk")
DIR_DEG_CLUSTER <- file.path(PROJECT_DIR, "09_DEG_Pseudobulk_ByCluster")
DIR_PATHWAY     <- file.path(PROJECT_DIR, "10_Pathway")
DIR_CORRELATION <- file.path(PROJECT_DIR, "11_Correlation")
DIR_LOGS        <- file.path(PROJECT_DIR, "12_Logs")

for (d in c(
  DIR_QC, DIR_OBJECTS, DIR_INTEGRATION, DIR_FIGURES,
  DIR_COMPOSITION, DIR_CELL_CYCLE, DIR_MARKERS,
  DIR_DEG_PB, DIR_DEG_CLUSTER, DIR_PATHWAY,
  DIR_CORRELATION, DIR_LOGS
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("SOURCE_DIR :", SOURCE_DIR, "\n")
cat("PROJECT_DIR:", PROJECT_DIR, "\n")


# ==============================================================================
# 1. PACKAGES AND REPRODUCIBILITY CHECK
# ==============================================================================

# ------------------------------------------------------------------------------
# 1.1 Check BiocManager
# ------------------------------------------------------------------------------
# BiocManager is required for checking the Bioconductor environment.
# Package installation should NOT be performed automatically inside the
# manuscript analysis script.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  stop(
    "BiocManager is not installed. Install it once using:\n",
    'install.packages("BiocManager")'
  )
}


# ------------------------------------------------------------------------------
# 1.2 Packages required for the complete manuscript analysis
# ------------------------------------------------------------------------------

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "ggplot2",
  "pheatmap",
  "DESeq2",
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db",
  "enrichplot"
)


# ------------------------------------------------------------------------------
# 1.3 Check whether all required packages are installed
# ------------------------------------------------------------------------------

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  
  stop(
    paste0(
      "\nThe following required package(s) are not installed:\n",
      paste(missing_packages, collapse = ", "),
      "\n\nInstall the missing packages before running the manuscript analysis."
    )
  )
}


# ------------------------------------------------------------------------------
# 1.4 Load packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  
  library(Seurat)
  library(SeuratObject)
  
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  
  library(DESeq2)
  
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(enrichplot)
  
})

cat(
  "\n====================================================\n",
  "All manuscript-analysis packages loaded successfully.\n",
  "====================================================\n"
)


# ------------------------------------------------------------------------------
# 1.5 Set random seed
# ------------------------------------------------------------------------------

# A fixed random seed improves reproducibility for stochastic analyses such as
# PCA initialization, clustering, UMAP, and related procedures.

RANDOM_SEED <- 123

set.seed(RANDOM_SEED)

cat(
  "Random seed:",
  RANDOM_SEED,
  "\n"
)


# ------------------------------------------------------------------------------
# 1.6 Record software versions
# ------------------------------------------------------------------------------

cat("\nR version:\n")
print(R.version.string)

cat("\nBioconductor version:\n")
print(BiocManager::version())

cat("\nImportant package versions:\n")

package_versions <- sapply(
  required_packages,
  function(pkg) {
    as.character(packageVersion(pkg))
  }
)

print(package_versions)


# ------------------------------------------------------------------------------
# 1.7 Basic version checks
# ------------------------------------------------------------------------------

# This pipeline was developed for Seurat v5.
if (packageVersion("Seurat") < "5.0.0") {
  stop(
    "Seurat version 5.0.0 or newer is required for this manuscript pipeline."
  )
}

# DESeq2 must be available for the primary pseudobulk DEG analysis.
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop(
    "DESeq2 is required for pseudobulk differential-expression analysis."
  )
}


# ------------------------------------------------------------------------------
# 1.8 Display library locations
# ------------------------------------------------------------------------------

cat("\nR library paths:\n")
print(.libPaths())


# ------------------------------------------------------------------------------
# END PACKAGE SETUP
# ------------------------------------------------------------------------------

cat("\nPackage setup completed successfully.\n\n")


# ==============================================================================
# SAVE SOFTWARE VERSION INFORMATION
# ==============================================================================

package_version_table <- data.frame(
  Package = names(package_versions),
  Version = unname(package_versions),
  stringsAsFactors = FALSE
)

write.csv(
  package_version_table,
  file.path(
    PROJECT_DIR,
    "package_versions_27.08.2026.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 2. BIOLOGICAL LABELS AND ANALYSIS PARAMETERS
# ==============================================================================


# ------------------------------------------------------------------------------
# 2.1 Biological treatment order
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# 2.1 Treatments excluded from ALL manuscript analyses
# ------------------------------------------------------------------------------

EXCLUDED_TREATMENTS <- c(
  "CAP_von_unten",
  "SM837+CAP_von_unten"
)


# Final biological treatment order used in all downstream analyses.
# The two "von_unten" conditions are intentionally excluded.

desired_order <- c(
  "Kontrolle_DMEM",
  "CAP_direkt",
  "CAP_indirekt",
  "SM837_alleine",
  "SM837+CAP_direkt",
  "SM837+CAP_indirekt"
)


# Confirm excluded treatments are not accidentally present in desired_order.

stopifnot(
  length(
    intersect(
      EXCLUDED_TREATMENTS,
      desired_order
    )
  ) == 0
)

# ------------------------------------------------------------------------------
# 2.2 Model-safe treatment IDs
# ------------------------------------------------------------------------------

condition_map <- c(
  "Kontrolle_DMEM"     = "Kontrolle_DMEM",
  "CAP_direkt"         = "CAP_direkt",
  "CAP_indirekt"       = "CAP_indirekt",
  "SM837_alleine"      = "SM837_alleine",
  "SM837+CAP_direkt"   = "SM837_CAP_direkt",
  "SM837+CAP_indirekt" = "SM837_CAP_indirekt"
)


safe_condition_order <- unname(
  condition_map[
    desired_order
  ]
)


stopifnot(
  length(safe_condition_order) ==
    length(desired_order)
)

stopifnot(
  !anyNA(safe_condition_order)
)

stopifnot(
  !anyDuplicated(safe_condition_order)
)

# ------------------------------------------------------------------------------
# 2.3 Validate treatment mapping
# ------------------------------------------------------------------------------

# Every biological treatment must have a corresponding model-safe ID.

stopifnot(
  length(safe_condition_order) ==
    length(desired_order)
)

stopifnot(
  !any(is.na(safe_condition_order))
)

# Model-safe IDs must also be unique.

stopifnot(
  !anyDuplicated(safe_condition_order)
)


# ------------------------------------------------------------------------------
# 2.4 Quality-control thresholds
# ------------------------------------------------------------------------------

# These thresholds are retained from the previous analysis.
#
# Their final suitability must be confirmed using the QC distributions
# generated separately for Lea1, Lea2 and Lea3.

QC_MIN_FEATURES <- 1000
QC_MAX_FEATURES <- 7000
QC_MAX_MT       <- 25


# Validate QC parameters.

stopifnot(
  QC_MIN_FEATURES > 0
)

stopifnot(
  QC_MAX_FEATURES >
    QC_MIN_FEATURES
)

stopifnot(
  QC_MAX_MT > 0,
  QC_MAX_MT <= 100
)


# ------------------------------------------------------------------------------
# 2.5 Seurat integration and clustering parameters
# ------------------------------------------------------------------------------

# Number of highly variable genes selected per dataset.

N_VARIABLE_FEATURES <- 2000


# Number of dimensions used for anchor integration.

N_INTEGRATION_DIMS <- 30


# Number of principal components calculated.
#
# Fifty PCs are calculated so that dimensionality can be examined using
# an elbow plot before interpreting downstream results.

N_PCA <- 50


# Number of PCs used for neighbor graph construction, clustering and UMAP.
#
# Initially retained as 30 for consistency with the previous analysis.
# This value will be checked against the PCA elbow plot.

N_CLUSTER_DIMS <- 30


# Seurat graph-based clustering resolution.

CLUSTER_RESOLUTION <- 0.5


# Create dimension vectors once so that the same values are used consistently.

integration_dims <- 1:N_INTEGRATION_DIMS
cluster_dims     <- 1:N_CLUSTER_DIMS


# Validate dimensionality parameters.

stopifnot(
  N_INTEGRATION_DIMS <= N_PCA
)

stopifnot(
  N_CLUSTER_DIMS <= N_PCA
)

stopifnot(
  CLUSTER_RESOLUTION > 0
)


# ------------------------------------------------------------------------------
# 2.6 Differential-expression thresholds
# ------------------------------------------------------------------------------

# Benjamini-Hochberg adjusted p-value / FDR threshold.

FDR_CUTOFF <- 0.05


# Minimum absolute log2 fold-change required for reporting a DEG.

LOG2FC_CUTOFF <- 0.25


stopifnot(
  FDR_CUTOFF > 0,
  FDR_CUTOFF < 1
)

stopifnot(
  LOG2FC_CUTOFF >= 0
)


# ------------------------------------------------------------------------------
# 2.7 Pseudobulk safeguards
# ------------------------------------------------------------------------------

# Minimum number of cells contributing to one
# replicate × treatment pseudobulk sample.

MIN_CELLS_PER_PB <- 20


# Lea1, Lea2 and Lea3 are expected to provide the three experimental
# replicates. For primary manuscript inference, all three should be
# represented in each treatment arm.

MIN_REPLICATES_PER_ARM <- 3


# Gene prefilter:
#
# A gene must have at least MIN_GENE_COUNT counts in at least
# MIN_SAMPLES_WITH_COUNT pseudobulk samples before DESeq2 fitting.

MIN_GENE_COUNT <- 10

MIN_SAMPLES_WITH_COUNT <- 3


stopifnot(
  MIN_CELLS_PER_PB > 0
)

stopifnot(
  MIN_REPLICATES_PER_ARM >= 2
)

stopifnot(
  MIN_GENE_COUNT >= 0
)

stopifnot(
  MIN_SAMPLES_WITH_COUNT >= 1
)


# ------------------------------------------------------------------------------
# 2.8 Analysis switches
# ------------------------------------------------------------------------------

# Primary replicate-aware pseudobulk DEG within individual Seurat clusters.

RUN_CLUSTER_PSEUDOBULK <- TRUE


# Exploratory cell-level Seurat FindMarkers analysis.
#
# FALSE by default because cell-level tests are not used as the primary
# replicate-level inferential analysis.

RUN_EXPLORATORY_CELL_LEVEL_DEG <- FALSE


# GO / KEGG / Reactome analysis for whole-population pseudobulk DEG.

RUN_PATHWAY_ALL_CELLS <- TRUE


# Cluster-specific pathway analysis can generate many output files.
# Keep disabled initially and enable selectively after reviewing
# cluster-specific DEG results.

RUN_PATHWAY_BY_CLUSTER <- FALSE


# ------------------------------------------------------------------------------
# 2.9 Display parameter summary
# ------------------------------------------------------------------------------

cat("\n====================================================\n")
cat("BIOLOGICAL AND ANALYSIS PARAMETERS\n")
cat("====================================================\n")

cat("\nTreatment order:\n")
print(desired_order)

cat("\nModel-safe condition IDs:\n")
print(safe_condition_order)

cat("\nQC thresholds:\n")
cat(
  "nFeature_RNA >", QC_MIN_FEATURES, "\n",
  "nFeature_RNA <", QC_MAX_FEATURES, "\n",
  "percent.mt     <", QC_MAX_MT, "%\n"
)

cat("\nSeurat parameters:\n")
cat(
  "Variable features :", N_VARIABLE_FEATURES, "\n",
  "Integration dims  :", N_INTEGRATION_DIMS, "\n",
  "PCA calculated    :", N_PCA, "\n",
  "Clustering dims   :", N_CLUSTER_DIMS, "\n",
  "Resolution        :", CLUSTER_RESOLUTION, "\n"
)

cat("\nDEG parameters:\n")
cat(
  "FDR cutoff          :", FDR_CUTOFF, "\n",
  "|log2FC| cutoff     :", LOG2FC_CUTOFF, "\n"
)

cat("\nPseudobulk parameters:\n")
cat(
  "Minimum cells/PB       :", MIN_CELLS_PER_PB, "\n",
  "Minimum replicates/arm :", MIN_REPLICATES_PER_ARM, "\n",
  "Minimum gene count     :", MIN_GENE_COUNT, "\n",
  "Minimum samples        :", MIN_SAMPLES_WITH_COUNT, "\n"
)

cat("\nSection 2 completed successfully.\n\n")


# ==============================================================================
# 3. PLANNED DIFFERENTIAL-EXPRESSION COMPARISONS
# ==============================================================================


# ------------------------------------------------------------------------------
# 3.1 Define biological comparisons
# ------------------------------------------------------------------------------

comparisons <- data.frame(
  
  test_label = c(
    "CAP_direkt",
    "CAP_indirekt",
    "SM837_alleine",
    "SM837+CAP_direkt",
    "SM837+CAP_indirekt",
    "SM837+CAP_direkt"
  ),
  
  reference_label = c(
    rep(
      "Kontrolle_DMEM",
      5
    ),
    "SM837+CAP_indirekt"
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# 3.2 Validate biological labels
# ------------------------------------------------------------------------------

invalid_test_labels <- setdiff(
  unique(comparisons$test_label),
  desired_order
)

if (length(invalid_test_labels) > 0) {
  
  stop(
    paste0(
      "Invalid test treatment label(s): ",
      paste(
        invalid_test_labels,
        collapse = ", "
      )
    )
  )
}


invalid_reference_labels <- setdiff(
  unique(comparisons$reference_label),
  desired_order
)

if (length(invalid_reference_labels) > 0) {
  
  stop(
    paste0(
      "Invalid reference treatment label(s): ",
      paste(
        invalid_reference_labels,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------------------------
# 3.3 Create model-safe treatment IDs
# ------------------------------------------------------------------------------

comparisons$test_id <- unname(
  condition_map[
    comparisons$test_label
  ]
)

comparisons$reference_id <- unname(
  condition_map[
    comparisons$reference_label
  ]
)


if (
  anyNA(comparisons$test_id) ||
  anyNA(comparisons$reference_id)
) {
  
  stop(
    "At least one planned comparison could not be mapped ",
    "to a model-safe condition ID."
  )
}


# ------------------------------------------------------------------------------
# 3.4 Create comparison names
# ------------------------------------------------------------------------------

comparisons$comparison_name <- paste0(
  comparisons$test_label,
  "_vs_",
  comparisons$reference_label
)


comparisons$comparison_id <- paste0(
  comparisons$test_id,
  "_vs_",
  comparisons$reference_id
)


# ------------------------------------------------------------------------------
# 3.5 Validate comparison structure
# ------------------------------------------------------------------------------

same_condition <- (
  comparisons$test_id ==
    comparisons$reference_id
)

if (any(same_condition)) {
  
  stop(
    "Invalid comparison detected: ",
    "test and reference condition are identical."
  )
}


if (anyDuplicated(comparisons$comparison_id)) {
  
  duplicated_comparisons <- comparisons$comparison_id[
    duplicated(comparisons$comparison_id)
  ]
  
  stop(
    paste0(
      "Duplicated comparison(s) detected: ",
      paste(
        duplicated_comparisons,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------------------------
# 3.6 Assign comparison type
# ------------------------------------------------------------------------------

comparisons$comparison_type <- ifelse(
  comparisons$reference_label ==
    "Kontrolle_DMEM",
  "Treatment_vs_Control",
  "Treatment_vs_Treatment"
)


# ------------------------------------------------------------------------------
# 3.7 Classify comparisons according to replication
# ------------------------------------------------------------------------------

# SM837+CAP_indirekt is absent from Lea2.
#
# Therefore comparisons involving this condition contain only
# Lea1 and Lea3 and are classified as exploratory n = 2 analyses.

comparisons$analysis_class <- ifelse(
  
  comparisons$test_label ==
    "SM837+CAP_indirekt" |
    comparisons$reference_label ==
    "SM837+CAP_indirekt",
  
  "Exploratory_n2",
  
  "Primary_n3"
)


# ------------------------------------------------------------------------------
# 3.8 Assign comparison number
# ------------------------------------------------------------------------------

comparisons$comparison_number <- seq_len(
  nrow(comparisons)
)


# ------------------------------------------------------------------------------
# 3.9 Reorder columns
# ------------------------------------------------------------------------------

comparisons <- comparisons %>%
  dplyr::select(
    comparison_number,
    comparison_name,
    comparison_id,
    comparison_type,
    analysis_class,
    test_label,
    reference_label,
    test_id,
    reference_id
  )


# ------------------------------------------------------------------------------
# 3.10 Display planned comparisons
# ------------------------------------------------------------------------------

cat(
  "\n====================================================\n"
)

cat(
  "PLANNED DIFFERENTIAL-EXPRESSION COMPARISONS\n"
)

cat(
  "====================================================\n\n"
)

print(
  comparisons
)

cat(
  "\nTotal planned comparisons:",
  nrow(comparisons),
  "\n"
)


# ------------------------------------------------------------------------------
# 3.11 Save planned comparisons
# ------------------------------------------------------------------------------

if (!dir.exists(DIR_LOGS)) {
  
  dir.create(
    DIR_LOGS,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


write.csv(
  comparisons,
  file.path(
    DIR_LOGS,
    "planned_comparisons.csv"
  ),
  row.names = FALSE
)


cat(
  "\nPlanned comparisons saved to:\n",
  file.path(
    DIR_LOGS,
    "planned_comparisons.csv"
  ),
  "\n"
)


cat(
  "\nSection 3 completed successfully.\n\n"
)


# ==============================================================================
# 4. HELPER FUNCTIONS — FINAL UPDATED VERSION
# ==============================================================================


# ==============================================================================
# 4.1 SAFE FILE / DIRECTORY NAME
# ==============================================================================

safe_name <- function(x) {
  
  x <- as.character(x)
  
  gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    x
  )
}


# ==============================================================================
# 4.2 SAFE PDF SAVING
# ==============================================================================

save_pdf <- function(
    plot_object,
    file,
    width = 10,
    height = 7
) {
  
  dir.create(
    dirname(file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  grDevices::pdf(
    file = file,
    width = width,
    height = height,
    onefile = TRUE
  )
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  print(plot_object)
  
  invisible(file)
}


# ==============================================================================
# 4.3 CLEAN AND VALIDATE SAMPLE METADATA
# ==============================================================================

clean_metadata <- function(
    obj,
    dataset_id
) {
  
  cat(
    "\n====================================================\n",
    "CLEANING METADATA: ", dataset_id, "\n",
    "====================================================\n",
    sep = ""
  )
  
  
  # ---------------------------------------------------------------------------
  # 1. Confirm Seurat object
  # ---------------------------------------------------------------------------
  
  if (!inherits(obj, "Seurat")) {
    
    stop(
      dataset_id,
      " is not a Seurat object."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 2. Check RNA assay safely
  # ---------------------------------------------------------------------------
  
  assay_names <- names(
    obj@assays
  )
  
  
  if (
    is.null(assay_names) ||
    length(assay_names) == 0
  ) {
    
    stop(
      dataset_id,
      ": no assays were found."
    )
  }
  
  
  if (!any(assay_names == "RNA")) {
    
    stop(
      dataset_id,
      ": RNA assay was not found. Available assays: ",
      paste(
        assay_names,
        collapse = ", "
      )
    )
  }
  
  
  SeuratObject::DefaultAssay(obj) <- "RNA"
  
  
  cat(
    dataset_id,
    " assay(s): ",
    paste(
      assay_names,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
  
  
  # ---------------------------------------------------------------------------
  # 3. Check Sample_Name metadata
  # ---------------------------------------------------------------------------
  
  metadata_names <- colnames(
    obj@meta.data
  )
  
  
  if (!any(metadata_names == "Sample_Name")) {
    
    stop(
      dataset_id,
      ": Sample_Name metadata column was not found."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 4. Confirm cells are present
  # ---------------------------------------------------------------------------
  
  n_before <- ncol(obj)
  
  
  if (n_before == 0) {
    
    stop(
      dataset_id,
      ": object contains zero cells."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 5. Add explicit dataset / replicate metadata
  # ---------------------------------------------------------------------------
  
  obj$dataset_id <- rep(
    dataset_id,
    ncol(obj)
  )
  
  
  obj$replicate_id <- rep(
    dataset_id,
    ncol(obj)
  )
  
  
  # ---------------------------------------------------------------------------
  # 6. Convert original Sample_Name to character
  # ---------------------------------------------------------------------------
  
  sample_name_raw <- as.character(
    obj@meta.data$Sample_Name
  )
  
  
  if (
    length(sample_name_raw) !=
    ncol(obj)
  ) {
    
    stop(
      dataset_id,
      ": Sample_Name length does not equal number of cells."
    )
  }
  
  
  if (anyNA(sample_name_raw)) {
    
    stop(
      dataset_id,
      ": NA detected in original Sample_Name."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 7. Preserve ORIGINAL labels
  # ---------------------------------------------------------------------------
  
  obj$Sample_Name_original <-
    sample_name_raw
  
  
  # ---------------------------------------------------------------------------
  # 8. Save RAW Sample_Name distribution
  # ---------------------------------------------------------------------------
  
  raw_sample_table <- as.data.frame(
    table(
      Sample_Name = sample_name_raw,
      useNA = "ifany"
    )
  )
  
  
  write.csv(
    raw_sample_table,
    file.path(
      DIR_QC,
      paste0(
        dataset_id,
        "_SampleName_RAW.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ---------------------------------------------------------------------------
  # 9. Count technical exclusions
  # ---------------------------------------------------------------------------
  
  n_multiplet <- sum(
    sample_name_raw ==
      "Multiplet"
  )
  
  
  n_undetermined <- sum(
    sample_name_raw ==
      "Undetermined"
  )
  
  
  # ---------------------------------------------------------------------------
  # 10. Count study-design exclusions
  # ---------------------------------------------------------------------------
  
  # >>> UPDATED <<<
  # These two treatment groups are excluded from ALL downstream analyses.
  
  n_cap_von_unten <- sum(
    sample_name_raw ==
      "CAP_von_unten"
  )
  
  
  n_sm837_cap_von_unten <- sum(
    sample_name_raw ==
      "SM837+CAP_von_unten"
  )
  
  
  cat(
    dataset_id,
    ": Multiplet = ",
    n_multiplet,
    "; Undetermined = ",
    n_undetermined,
    "; CAP_von_unten = ",
    n_cap_von_unten,
    "; SM837+CAP_von_unten = ",
    n_sm837_cap_von_unten,
    "\n",
    sep = ""
  )
  
  
  # ---------------------------------------------------------------------------
  # 11. Define excluded cells
  # ---------------------------------------------------------------------------
  
  technical_exclusion <- (
    sample_name_raw ==
      "Multiplet" |
      sample_name_raw ==
      "Undetermined"
  )
  
  
  # >>> UPDATED <<<
  
  design_exclusion <- (
    sample_name_raw %in%
      EXCLUDED_TREATMENTS
  )
  
  
  keep_cells <- !(
    technical_exclusion |
      design_exclusion
  )
  
  
  # ---------------------------------------------------------------------------
  # 12. Identify cells to retain
  # ---------------------------------------------------------------------------
  
  cells_to_keep <- colnames(obj)[
    keep_cells
  ]
  
  
  if (length(cells_to_keep) == 0) {
    
    stop(
      dataset_id,
      ": sample cleaning would remove every cell."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 13. Remove excluded cells
  # ---------------------------------------------------------------------------
  
  # >>> UPDATED / IMPORTANT <<<
  
  obj <- subset(
    obj,
    cells = cells_to_keep
  )
  
  
  # ---------------------------------------------------------------------------
  # 14. Recover labels after subsetting
  # ---------------------------------------------------------------------------
  
  sample_name_clean <- as.character(
    obj@meta.data$Sample_Name_original
  )
  
  
  # ---------------------------------------------------------------------------
  # 15. Correct known Lea3 typo
  # ---------------------------------------------------------------------------
  
  sample_name_clean[
    sample_name_clean ==
      "SM1837+CAP_indirekt"
  ] <- "SM837+CAP_indirekt"
  
  
  # ---------------------------------------------------------------------------
  # 16. Standardize SM837-alone label
  # ---------------------------------------------------------------------------
  
  sample_name_clean[
    sample_name_clean ==
      "SM837"
  ] <- "SM837_alleine"
  
  
  # ---------------------------------------------------------------------------
  # 17. Check NA values
  # ---------------------------------------------------------------------------
  
  if (anyNA(sample_name_clean)) {
    
    stop(
      dataset_id,
      ": NA detected after Sample_Name cleaning."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 18. Confirm technical exclusions are absent
  # ---------------------------------------------------------------------------
  
  if (
    any(
      sample_name_clean ==
      "Multiplet"
    ) ||
    any(
      sample_name_clean ==
      "Undetermined"
    )
  ) {
    
    stop(
      dataset_id,
      ": Multiplet or Undetermined remains after cleaning."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 19. Confirm excluded treatments are absent
  # ---------------------------------------------------------------------------
  
  excluded_remaining <- intersect(
    unique(
      sample_name_clean
    ),
    EXCLUDED_TREATMENTS
  )
  
  
  if (
    length(
      excluded_remaining
    ) > 0
  ) {
    
    stop(
      dataset_id,
      ": excluded treatment group(s) remain after cleaning: ",
      paste(
        excluded_remaining,
        collapse = ", "
      )
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 20. Confirm known typo is gone
  # ---------------------------------------------------------------------------
  
  if (
    any(
      grepl(
        "SM1837",
        sample_name_clean,
        fixed = TRUE
      )
    )
  ) {
    
    stop(
      dataset_id,
      ": SM1837 typo remains after cleaning."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 21. Validate biological treatment labels
  # ---------------------------------------------------------------------------
  
  observed_labels <- sort(
    unique(
      sample_name_clean
    )
  )
  
  
  unexpected_labels <- setdiff(
    observed_labels,
    desired_order
  )
  
  
  if (
    length(
      unexpected_labels
    ) > 0
  ) {
    
    stop(
      dataset_id,
      ": unexpected Sample_Name value(s): ",
      paste(
        unexpected_labels,
        collapse = ", "
      )
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 22. Report missing treatments
  # ---------------------------------------------------------------------------
  
  missing_labels <- setdiff(
    desired_order,
    observed_labels
  )
  
  
  if (
    length(
      missing_labels
    ) > 0
  ) {
    
    warning(
      dataset_id,
      " is missing treatment(s): ",
      paste(
        missing_labels,
        collapse = ", "
      )
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 23. Create ordered biological factor
  # ---------------------------------------------------------------------------
  
  obj$Sample_Name <- factor(
    sample_name_clean,
    levels = desired_order
  )
  
  
  if (anyNA(obj$Sample_Name)) {
    
    stop(
      dataset_id,
      ": Sample_Name factor conversion generated NA values."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 24. Create model-safe condition ID
  # ---------------------------------------------------------------------------
  
  condition_id_character <- unname(
    condition_map[
      as.character(
        obj$Sample_Name
      )
    ]
  )
  
  
  if (
    anyNA(
      condition_id_character
    )
  ) {
    
    stop(
      dataset_id,
      ": condition mapping generated NA values."
    )
  }
  
  
  obj$condition_id <- factor(
    condition_id_character,
    levels = safe_condition_order
  )
  
  
  if (
    anyNA(
      obj$condition_id
    )
  ) {
    
    stop(
      dataset_id,
      ": condition_id factor conversion generated NA values."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 25. Final cell number
  # ---------------------------------------------------------------------------
  
  n_after <- ncol(obj)
  
  
  # ---------------------------------------------------------------------------
  # 26. Validate number of removed cells
  # ---------------------------------------------------------------------------
  
  expected_removed <-
    n_multiplet +
    n_undetermined +
    n_cap_von_unten +
    n_sm837_cap_von_unten
  
  
  observed_removed <-
    n_before -
    n_after
  
  
  if (
    observed_removed !=
    expected_removed
  ) {
    
    stop(
      dataset_id,
      ": cell-removal audit mismatch. ",
      "Expected removed = ",
      expected_removed,
      "; observed removed = ",
      observed_removed,
      "."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 27. Save cleaning audit
  # ---------------------------------------------------------------------------
  
  cleaning_audit <- data.frame(
    
    dataset_id =
      dataset_id,
    
    cells_raw =
      n_before,
    
    multiplet_removed =
      n_multiplet,
    
    undetermined_removed =
      n_undetermined,
    
    technical_cells_removed =
      n_multiplet +
      n_undetermined,
    
    CAP_von_unten_removed =
      n_cap_von_unten,
    
    SM837_CAP_von_unten_removed =
      n_sm837_cap_von_unten,
    
    study_design_cells_removed =
      n_cap_von_unten +
      n_sm837_cap_von_unten,
    
    cells_after_sample_cleaning =
      n_after,
    
    total_cells_removed =
      observed_removed,
    
    percent_removed = round(
      100 *
        observed_removed /
        n_before,
      3
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  write.csv(
    cleaning_audit,
    file.path(
      DIR_QC,
      paste0(
        dataset_id,
        "_sample_cleaning_audit.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ---------------------------------------------------------------------------
  # 28. Save cleaned treatment distribution
  # ---------------------------------------------------------------------------
  
  cleaned_sample_table <- as.data.frame(
    table(
      Sample_Name =
        obj$Sample_Name,
      useNA = "ifany"
    )
  )
  
  
  write.csv(
    cleaned_sample_table,
    file.path(
      DIR_QC,
      paste0(
        dataset_id,
        "_SampleName_CLEAN.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ---------------------------------------------------------------------------
  # 29. Sample_Tag × Sample_Name check
  # ---------------------------------------------------------------------------
  
  if (
    any(
      colnames(
        obj@meta.data
      ) ==
      "Sample_Tag"
    )
  ) {
    
    sample_tag_table <- as.data.frame(
      table(
        
        Sample_Tag =
          as.character(
            obj@meta.data$Sample_Tag
          ),
        
        Sample_Name =
          as.character(
            obj$Sample_Name
          ),
        
        useNA = "ifany"
      )
    )
    
    
    write.csv(
      sample_tag_table,
      file.path(
        DIR_QC,
        paste0(
          dataset_id,
          "_SampleTag_by_SampleName.csv"
        )
      ),
      row.names = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 30. Display final summary
  # ---------------------------------------------------------------------------
  
  cat(
    "\n",
    dataset_id,
    " metadata cleaning completed successfully.\n",
    "Cells before cleaning : ",
    n_before,
    "\n",
    "Technical removals    : ",
    n_multiplet +
      n_undetermined,
    "\n",
    "Design exclusions     : ",
    n_cap_von_unten +
      n_sm837_cap_von_unten,
    "\n",
    "Cells after cleaning  : ",
    n_after,
    "\n",
    "Total cells removed   : ",
    observed_removed,
    "\n",
    sep = ""
  )
  
  
  cat(
    "\nFinal treatment distribution:\n"
  )
  
  
  print(
    table(
      obj$Sample_Name,
      useNA = "ifany"
    )
  )
  
  
  return(obj)
}


# ==============================================================================
# 4.4 CALCULATE MITOCHONDRIAL RNA FRACTION
# ==============================================================================

# >>> UPDATED <<<
# This replaces the old version containing:
#
# if (!"RNA" %in% Assays(obj))
#
# That old assay check caused the:
# "'match' requires vector arguments"
# error in your R environment.

add_percent_mt <- function(
    obj,
    dataset_id
) {
  
  # ---------------------------------------------------------------------------
  # 1. Confirm Seurat object
  # ---------------------------------------------------------------------------
  
  if (!inherits(obj, "Seurat")) {
    
    stop(
      dataset_id,
      ": object is not a Seurat object."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 2. Confirm RNA assay safely
  # ---------------------------------------------------------------------------
  
  assay_names <- names(
    obj@assays
  )
  
  
  if (
    is.null(assay_names) ||
    length(assay_names) == 0
  ) {
    
    stop(
      dataset_id,
      ": no assays were found."
    )
  }
  
  
  if (
    !any(
      assay_names ==
      "RNA"
    )
  ) {
    
    stop(
      dataset_id,
      ": RNA assay missing."
    )
  }
  
  
  SeuratObject::DefaultAssay(obj) <- "RNA"
  
  
  # ---------------------------------------------------------------------------
  # 3. Get feature names from RNA assay
  # ---------------------------------------------------------------------------
  
  feature_names <- rownames(
    obj[["RNA"]]
  )
  
  
  if (
    is.null(feature_names) ||
    length(feature_names) == 0
  ) {
    
    stop(
      dataset_id,
      ": RNA assay contains no feature names."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 4. Identify mitochondrial genes
  # ---------------------------------------------------------------------------
  
  mt_genes <- grep(
    "^MT-",
    feature_names,
    value = TRUE
  )
  
  
  cat(
    dataset_id,
    " - mitochondrial genes found: ",
    length(mt_genes),
    "\n",
    sep = ""
  )
  
  
  if (
    length(mt_genes) == 0
  ) {
    
    stop(
      dataset_id,
      ": no genes match '^MT-'. ",
      "Check whether rownames are HGNC symbols or Ensembl IDs."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 5. Calculate mitochondrial percentage
  # ---------------------------------------------------------------------------
  
  obj[["percent.mt"]] <-
    Seurat::PercentageFeatureSet(
      object = obj,
      assay = "RNA",
      features = mt_genes
    )
  
  
  # ---------------------------------------------------------------------------
  # 6. Validate mitochondrial percentage
  # ---------------------------------------------------------------------------
  
  if (
    !any(
      colnames(
        obj@meta.data
      ) ==
      "percent.mt"
    )
  ) {
    
    stop(
      dataset_id,
      ": percent.mt was not created."
    )
  }
  
  
  if (
    anyNA(
      obj$percent.mt
    )
  ) {
    
    stop(
      dataset_id,
      ": NA values detected in percent.mt."
    )
  }
  
  
  if (
    any(
      obj$percent.mt < 0 |
      obj$percent.mt > 100
    )
  ) {
    
    stop(
      dataset_id,
      ": invalid percent.mt values detected."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 7. Display mitochondrial summary
  # ---------------------------------------------------------------------------
  
  cat(
    dataset_id,
    " - percent.mt range: ",
    round(
      min(
        obj$percent.mt
      ),
      3
    ),
    " to ",
    round(
      max(
        obj$percent.mt
      ),
      3
    ),
    "%\n",
    sep = ""
  )
  
  
  cat(
    dataset_id,
    " - median percent.mt: ",
    round(
      median(
        obj$percent.mt
      ),
      3
    ),
    "%\n",
    sep = ""
  )
  
  
  return(obj)
}


# ==============================================================================
# 4.5 QC SUMMARY
# ==============================================================================

qc_summary <- function(
    obj,
    dataset_id,
    stage
) {
  
  required_qc_columns <- c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  )
  
  
  missing_qc_columns <- setdiff(
    required_qc_columns,
    colnames(
      obj@meta.data
    )
  )
  
  
  if (
    length(
      missing_qc_columns
    ) > 0
  ) {
    
    stop(
      dataset_id,
      ": missing QC column(s): ",
      paste(
        missing_qc_columns,
        collapse = ", "
      )
    )
  }
  
  
  data.frame(
    
    dataset_id =
      dataset_id,
    
    stage =
      stage,
    
    cells =
      ncol(obj),
    
    median_nFeature_RNA =
      median(
        obj$nFeature_RNA,
        na.rm = TRUE
      ),
    
    mean_nFeature_RNA =
      mean(
        obj$nFeature_RNA,
        na.rm = TRUE
      ),
    
    q05_nFeature_RNA =
      unname(
        stats::quantile(
          obj$nFeature_RNA,
          probs = 0.05,
          na.rm = TRUE
        )
      ),
    
    q95_nFeature_RNA =
      unname(
        stats::quantile(
          obj$nFeature_RNA,
          probs = 0.95,
          na.rm = TRUE
        )
      ),
    
    median_nCount_RNA =
      median(
        obj$nCount_RNA,
        na.rm = TRUE
      ),
    
    mean_nCount_RNA =
      mean(
        obj$nCount_RNA,
        na.rm = TRUE
      ),
    
    q05_nCount_RNA =
      unname(
        stats::quantile(
          obj$nCount_RNA,
          probs = 0.05,
          na.rm = TRUE
        )
      ),
    
    q95_nCount_RNA =
      unname(
        stats::quantile(
          obj$nCount_RNA,
          probs = 0.95,
          na.rm = TRUE
        )
      ),
    
    median_percent_mt =
      median(
        obj$percent.mt,
        na.rm = TRUE
      ),
    
    mean_percent_mt =
      mean(
        obj$percent.mt,
        na.rm = TRUE
      ),
    
    q95_percent_mt =
      unname(
        stats::quantile(
          obj$percent.mt,
          probs = 0.95,
          na.rm = TRUE
        )
      ),
    
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 4.6 SAVE QC PLOTS
# ==============================================================================

save_qc_plots <- function(
    obj,
    dataset_id,
    stage
) {
  
  # ---------------------------------------------------------------------------
  # Violin plots
  # ---------------------------------------------------------------------------
  
  vln <- Seurat::VlnPlot(
    obj,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
    group.by = "Sample_Name",
    pt.size = 0,
    ncol = 3
  )
  
  
  if (
    requireNamespace(
      "patchwork",
      quietly = TRUE
    )
  ) {
    
    vln <- vln +
      patchwork::plot_annotation(
        title = paste(
          dataset_id,
          stage
        )
      )
  }
  
  
  save_pdf(
    vln,
    file.path(
      DIR_QC,
      paste0(
        dataset_id,
        "_",
        stage,
        "_QC_violin.pdf"
      )
    ),
    width = 14,
    height = 8
  )
  
  
  # ---------------------------------------------------------------------------
  # nCount vs nFeature scatter
  # ---------------------------------------------------------------------------
  
  scatter <- Seurat::FeatureScatter(
    obj,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    group.by = "Sample_Name"
  )
  
  
  save_pdf(
    scatter,
    file.path(
      DIR_QC,
      paste0(
        dataset_id,
        "_",
        stage,
        "_nCount_vs_nFeature.pdf"
      )
    ),
    width = 9,
    height = 7
  )
  
  
  invisible(
    list(
      violin = vln,
      scatter = scatter
    )
  )
}


# ==============================================================================
# 4.7 QUALITY-CONTROL FILTERING
# ==============================================================================

filter_qc <- function(
    obj,
    dataset_id
) {
  
  n_before <- ncol(obj)
  
  
  if (
    n_before == 0
  ) {
    
    stop(
      dataset_id,
      ": zero cells available before QC filtering."
    )
  }
  
  
  conditions_before <- unique(
    as.character(
      obj$Sample_Name
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Apply QC filters
  # ---------------------------------------------------------------------------
  
  filtered <- subset(
    obj,
    subset =
      nFeature_RNA >
      QC_MIN_FEATURES &
      nFeature_RNA <
      QC_MAX_FEATURES &
      percent.mt <
      QC_MAX_MT
  )
  
  
  n_after <- ncol(
    filtered
  )
  
  
  if (
    n_after == 0
  ) {
    
    stop(
      dataset_id,
      ": QC filtering removed every cell."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Detect complete treatment loss
  # ---------------------------------------------------------------------------
  
  conditions_after <- unique(
    as.character(
      filtered$Sample_Name
    )
  )
  
  
  lost_conditions <- setdiff(
    conditions_before,
    conditions_after
  )
  
  
  if (
    length(
      lost_conditions
    ) > 0
  ) {
    
    warning(
      dataset_id,
      ": QC filtering removed all cells from treatment(s): ",
      paste(
        lost_conditions,
        collapse = ", "
      )
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Restore factor levels
  # ---------------------------------------------------------------------------
  
  filtered$Sample_Name <- factor(
    as.character(
      filtered$Sample_Name
    ),
    levels = desired_order
  )
  
  
  filtered$condition_id <- factor(
    as.character(
      filtered$condition_id
    ),
    levels = safe_condition_order
  )
  
  
  # >>> UPDATED <<<
  # Confirm factor reconstruction did not create NAs.
  
  if (
    anyNA(
      filtered$Sample_Name
    ) ||
    anyNA(
      filtered$condition_id
    )
  ) {
    
    stop(
      dataset_id,
      ": factor reconstruction after QC generated NA values."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Save QC filtering audit
  # ---------------------------------------------------------------------------
  
  qc_audit <- data.frame(
    
    dataset_id =
      dataset_id,
    
    cells_before_QC =
      n_before,
    
    cells_after_QC =
      n_after,
    
    cells_removed_QC =
      n_before -
      n_after,
    
    percent_removed_QC =
      round(
        100 *
          (n_before -
             n_after) /
          n_before,
        3
      ),
    
    stringsAsFactors = FALSE
  )
  
  
  write.csv(
    qc_audit,
    file.path(
      DIR_QC,
      paste0(
        dataset_id,
        "_QC_filter_audit.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  condition_count_table <- as.data.frame(
    table(
      Sample_Name =
        filtered$Sample_Name
    )
  )
  
  
  write.csv(
    condition_count_table,
    file.path(
      DIR_QC,
      paste0(
        dataset_id,
        "_cells_by_condition_AFTER_QC.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  cat(
    dataset_id,
    " cells before QC: ",
    n_before,
    "; after QC: ",
    n_after,
    "\n",
    sep = ""
  )
  
  
  return(filtered)
}


# ==============================================================================
# 4.8 JOIN SEURAT v5 RNA LAYERS ONLY WHEN NECESSARY
# ==============================================================================

join_rna_layers_if_needed <- function(
    obj
) {
  
  # >>> UPDATED <<<
  # Avoid old Assays(obj) membership check.
  
  assay_names <- names(
    obj@assays
  )
  
  
  if (
    !any(
      assay_names ==
      "RNA"
    )
  ) {
    
    stop(
      "RNA assay is missing."
    )
  }
  
  
  SeuratObject::DefaultAssay(obj) <- "RNA"
  
  
  # Legacy Seurat assay does not need JoinLayers().
  
  if (
    !inherits(
      obj[["RNA"]],
      "Assay5"
    )
  ) {
    
    return(obj)
  }
  
  
  layer_names <- SeuratObject::Layers(
    obj[["RNA"]]
  )
  
  
  split_count_layers <- grep(
    "^counts\\.",
    layer_names,
    value = TRUE
  )
  
  
  split_data_layers <- grep(
    "^data\\.",
    layer_names,
    value = TRUE
  )
  
  
  needs_joining <- (
    length(
      split_count_layers
    ) > 0 ||
      length(
        split_data_layers
      ) > 0
  )
  
  
  if (
    needs_joining
  ) {
    
    cat(
      "\nJoining split RNA layers:\n",
      paste(
        layer_names,
        collapse = ", "
      ),
      "\n"
    )
    
    
    obj <- SeuratObject::JoinLayers(
      object = obj,
      assay = "RNA"
    )
  }
  
  
  return(obj)
}


# ==============================================================================
# 4.9 BUILD PSEUDOBULK COUNT MATRIX
# ==============================================================================

make_pseudobulk <- function(
    obj,
    group_columns,
    min_cells =
      MIN_CELLS_PER_PB
) {
  
  if (
    min_cells < 1
  ) {
    
    stop(
      "min_cells must be >= 1."
    )
  }
  
  
  obj <- join_rna_layers_if_needed(
    obj
  )
  
  
  SeuratObject::DefaultAssay(obj) <- "RNA"
  
  
  metadata <- obj@meta.data
  
  
  missing_columns <- setdiff(
    group_columns,
    colnames(
      metadata
    )
  )
  
  
  if (
    length(
      missing_columns
    ) > 0
  ) {
    
    stop(
      "Missing pseudobulk metadata column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  
  group_data <- metadata[
    ,
    group_columns,
    drop = FALSE
  ]
  
  
  if (
    anyNA(
      group_data
    )
  ) {
    
    stop(
      "NA detected in pseudobulk grouping metadata."
    )
  }
  
  
  group_data[] <- lapply(
    group_data,
    as.character
  )
  
  
  empty_group_values <- vapply(
    group_data,
    function(x) {
      any(
        !nzchar(
          trimws(x)
        )
      )
    },
    logical(1)
  )
  
  
  if (
    any(
      empty_group_values
    )
  ) {
    
    stop(
      "Empty values detected in pseudobulk grouping column(s): ",
      paste(
        names(
          empty_group_values
        )[
          empty_group_values
        ],
        collapse = ", "
      )
    )
  }
  
  
  pb_id <- do.call(
    paste,
    c(
      group_data,
      sep = "___"
    )
  )
  
  
  obj$pb_id <- pb_id
  
  
  cell_counts <- table(
    obj$pb_id
  )
  
  
  keep_ids <- names(
    cell_counts[
      cell_counts >=
        min_cells
    ]
  )
  
  
  if (
    length(
      keep_ids
    ) == 0
  ) {
    
    stop(
      "No pseudobulk group contains at least ",
      min_cells,
      " cells."
    )
  }
  
  
  pb_counts <- Seurat::AggregateExpression(
    object = obj,
    assays = "RNA",
    group.by = "pb_id",
    return.seurat = FALSE,
    verbose = FALSE
  )$RNA
  
  
  keep_ids <- intersect(
    colnames(
      pb_counts
    ),
    keep_ids
  )
  
  
  if (
    length(
      keep_ids
    ) == 0
  ) {
    
    stop(
      "No aggregated pseudobulk IDs matched minimum-cell groups."
    )
  }
  
  
  pb_counts <- pb_counts[
    ,
    keep_ids,
    drop = FALSE
  ]
  
  
  metadata_pb <- obj@meta.data
  
  
  pb_meta <- metadata_pb %>%
    dplyr::select(
      pb_id,
      dplyr::all_of(
        group_columns
      )
    ) %>%
    dplyr::distinct() %>%
    as.data.frame()
  
  
  if (
    anyDuplicated(
      pb_meta$pb_id
    )
  ) {
    
    stop(
      "One pseudobulk ID maps to more than one metadata combination."
    )
  }
  
  
  rownames(
    pb_meta
  ) <- pb_meta$pb_id
  
  
  if (
    !setequal(
      colnames(
        pb_counts
      ),
      rownames(
        pb_meta
      )
    )
  ) {
    
    stop(
      "Pseudobulk count columns and metadata IDs do not match."
    )
  }
  
  
  pb_meta <- pb_meta[
    colnames(
      pb_counts
    ),
    ,
    drop = FALSE
  ]
  
  
  stopifnot(
    identical(
      colnames(
        pb_counts
      ),
      rownames(
        pb_meta
      )
    )
  )
  
  
  pb_meta$n_cells <- as.integer(
    cell_counts[
      rownames(
        pb_meta
      )
    ]
  )
  
  
  if (
    anyNA(
      pb_meta$n_cells
    )
  ) {
    
    stop(
      "Unable to assign cell counts to pseudobulk metadata."
    )
  }
  
  
  pb_counts <- as.matrix(
    pb_counts
  )
  
  
  pb_counts <- round(
    pb_counts
  )
  
  
  if (
    anyNA(
      pb_counts
    )
  ) {
    
    stop(
      "NA detected in pseudobulk count matrix."
    )
  }
  
  
  if (
    any(
      pb_counts < 0
    )
  ) {
    
    stop(
      "Negative values detected in pseudobulk count matrix."
    )
  }
  
  
  storage.mode(
    pb_counts
  ) <- "integer"
  
  
  return(
    list(
      counts =
        pb_counts,
      metadata =
        pb_meta
    )
  )
}


# ==============================================================================
# 4.10 PREPARE AND FIT DESeq2 PSEUDOBULK MODEL
# ==============================================================================

prepare_deseq2 <- function(
    pb_counts,
    pb_meta
) {
  
  if (
    ncol(
      pb_counts
    ) < 2
  ) {
    
    stop(
      "Too few pseudobulk samples for DESeq2."
    )
  }
  
  
  if (
    !identical(
      colnames(
        pb_counts
      ),
      rownames(
        pb_meta
      )
    )
  ) {
    
    stop(
      "Pseudobulk count columns and metadata rows are not aligned."
    )
  }
  
  
  required_metadata <- c(
    "replicate_id",
    "condition_id"
  )
  
  
  missing_metadata <- setdiff(
    required_metadata,
    colnames(
      pb_meta
    )
  )
  
  
  if (
    length(
      missing_metadata
    ) > 0
  ) {
    
    stop(
      "Missing DESeq2 metadata column(s): ",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  }
  
  
  pb_meta$replicate_id <- droplevels(
    factor(
      as.character(
        pb_meta$replicate_id
      )
    )
  )
  
  
  if (
    length(
      levels(
        pb_meta$replicate_id
      )
    ) < 2
  ) {
    
    stop(
      "Fewer than two replicate levels are represented."
    )
  }
  
  
  pb_meta$condition_id <- factor(
    as.character(
      pb_meta$condition_id
    ),
    levels =
      safe_condition_order
  )
  
  
  if (
    anyNA(
      pb_meta$condition_id
    )
  ) {
    
    stop(
      "NA detected in DESeq2 condition_id."
    )
  }
  
  
  pb_meta$condition_id <- droplevels(
    pb_meta$condition_id
  )
  
  
  if (
    length(
      levels(
        pb_meta$condition_id
      )
    ) < 2
  ) {
    
    stop(
      "Fewer than two treatment conditions are represented."
    )
  }
  
  
  if (
    "Kontrolle_DMEM" %in%
    levels(
      pb_meta$condition_id
    )
  ) {
    
    pb_meta$condition_id <- stats::relevel(
      pb_meta$condition_id,
      ref = "Kontrolle_DMEM"
    )
  }
  
  
  keep_gene <- rowSums(
    pb_counts >=
      MIN_GENE_COUNT
  ) >=
    MIN_SAMPLES_WITH_COUNT
  
  
  pb_counts_filtered <- pb_counts[
    keep_gene,
    ,
    drop = FALSE
  ]
  
  
  if (
    nrow(
      pb_counts_filtered
    ) == 0
  ) {
    
    stop(
      "No genes passed the pseudobulk gene-count filter."
    )
  }
  
  
  cat(
    "\nPseudobulk genes before filtering: ",
    nrow(
      pb_counts
    ),
    "\n",
    "Pseudobulk genes after filtering : ",
    nrow(
      pb_counts_filtered
    ),
    "\n",
    sep = ""
  )
  
  
  model_matrix <- stats::model.matrix(
    ~ replicate_id +
      condition_id,
    data = pb_meta
  )
  
  
  if (
    qr(
      model_matrix
    )$rank <
    ncol(
      model_matrix
    )
  ) {
    
    stop(
      "DESeq2 design matrix is not full rank. ",
      "Inspect replicate x treatment coverage."
    )
  }
  
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData =
      pb_counts_filtered,
    colData =
      pb_meta,
    design =
      ~ replicate_id +
      condition_id
  )
  
  
  dds <- DESeq2::DESeq(
    dds,
    quiet = TRUE
  )
  
  
  return(dds)
}


# ==============================================================================
# 4.11 CHECK REPLICATION FOR A PLANNED CONTRAST
# ==============================================================================

has_replication <- function(
    pb_meta,
    test_id,
    reference_id,
    verbose = TRUE
) {
  
  condition_character <- as.character(
    pb_meta$condition_id
  )
  
  
  replicate_character <- as.character(
    pb_meta$replicate_id
  )
  
  
  test_replicates <- sort(
    unique(
      replicate_character[
        condition_character ==
          test_id
      ]
    )
  )
  
  
  reference_replicates <- sort(
    unique(
      replicate_character[
        condition_character ==
          reference_id
      ]
    )
  )
  
  
  test_replicates <- test_replicates[
    !is.na(
      test_replicates
    )
  ]
  
  
  reference_replicates <- reference_replicates[
    !is.na(
      reference_replicates
    )
  ]
  
  
  shared_replicates <- intersect(
    test_replicates,
    reference_replicates
  )
  
  
  if (
    verbose
  ) {
    
    cat(
      "\nReplication check\n",
      "-----------------\n",
      "Test condition      : ",
      test_id,
      "\n",
      "Test replicates     : ",
      ifelse(
        length(
          test_replicates
        ) > 0,
        paste(
          test_replicates,
          collapse = ", "
        ),
        "NONE"
      ),
      "\n",
      "Reference condition : ",
      reference_id,
      "\n",
      "Reference replicates: ",
      ifelse(
        length(
          reference_replicates
        ) > 0,
        paste(
          reference_replicates,
          collapse = ", "
        ),
        "NONE"
      ),
      "\n",
      "Shared replicates   : ",
      ifelse(
        length(
          shared_replicates
        ) > 0,
        paste(
          shared_replicates,
          collapse = ", "
        ),
        "NONE"
      ),
      "\n",
      sep = ""
    )
  }
  
  
  return(
    length(
      test_replicates
    ) >=
      MIN_REPLICATES_PER_ARM &&
      length(
        reference_replicates
      ) >=
      MIN_REPLICATES_PER_ARM &&
      length(
        shared_replicates
      ) >=
      MIN_REPLICATES_PER_ARM
  )
}


# ==============================================================================
# 4.12 EXTRACT ONE DESeq2 CONTRAST
# ==============================================================================

extract_contrast <- function(
    dds,
    test_id,
    reference_id,
    test_label,
    reference_label
) {
  
  if (
    identical(
      test_id,
      reference_id
    )
  ) {
    
    stop(
      "Test and reference conditions are identical."
    )
  }
  
  
  available_conditions <- levels(
    dds$condition_id
  )
  
  
  if (
    !test_id %in%
    available_conditions
  ) {
    
    stop(
      "Test condition is absent from DESeq2 model: ",
      test_id
    )
  }
  
  
  if (
    !reference_id %in%
    available_conditions
  ) {
    
    stop(
      "Reference condition is absent from DESeq2 model: ",
      reference_id
    )
  }
  
  
  result_object <- DESeq2::results(
    dds,
    contrast = c(
      "condition_id",
      test_id,
      reference_id
    ),
    alpha =
      FDR_CUTOFF
  )
  
  
  result_df <- as.data.frame(
    result_object
  )
  
  
  result_df$gene <- rownames(
    result_df
  )
  
  
  comparison_name <- paste0(
    test_label,
    "_vs_",
    reference_label
  )
  
  
  result_df$comparison <-
    comparison_name
  
  result_df$test_condition <-
    test_label
  
  result_df$reference_condition <-
    reference_label
  
  
  result_df <- result_df %>%
    dplyr::select(
      gene,
      comparison,
      test_condition,
      reference_condition,
      dplyr::everything()
    ) %>%
    dplyr::arrange(
      padj
    )
  
  
  significant <- result_df %>%
    dplyr::filter(
      !is.na(
        padj
      ),
      padj <
        FDR_CUTOFF,
      abs(
        log2FoldChange
      ) >
        LOG2FC_CUTOFF
    )
  
  
  up <- significant %>%
    dplyr::filter(
      log2FoldChange >
        0
    )
  
  
  down <- significant %>%
    dplyr::filter(
      log2FoldChange <
        0
    )
  
  
  return(
    list(
      all =
        result_df,
      significant =
        significant,
      up =
        up,
      down =
        down,
      comparison_name =
        comparison_name,
      test_label =
        test_label,
      reference_label =
        reference_label
    )
  )
}


# ==============================================================================
# 4.13 SYMBOL TO ENTREZ MAPPING
# ==============================================================================

map_symbols_to_entrez <- function(
    symbols
) {
  
  symbols <- unique(
    as.character(
      symbols
    )
  )
  
  
  symbols <- symbols[
    !is.na(
      symbols
    ) &
      nzchar(
        trimws(
          symbols
        )
      )
  ]
  
  
  if (
    length(
      symbols
    ) == 0
  ) {
    
    return(
      data.frame(
        SYMBOL =
          character(0),
        ENTREZID =
          character(0),
        stringsAsFactors =
          FALSE
      )
    )
  }
  
  
  mapped <- tryCatch(
    
    suppressMessages(
      clusterProfiler::bitr(
        symbols,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb =
          org.Hs.eg.db
      )
    ),
    
    error = function(e) {
      
      warning(
        "Gene-symbol mapping failed: ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  
  if (
    is.null(
      mapped
    ) ||
    nrow(
      mapped
    ) == 0
  ) {
    
    return(
      data.frame(
        SYMBOL =
          character(0),
        ENTREZID =
          character(0),
        stringsAsFactors =
          FALSE
      )
    )
  }
  
  
  mapped <- mapped[
    ,
    c(
      "SYMBOL",
      "ENTREZID"
    ),
    drop = FALSE
  ]
  
  
  return(
    unique(
      mapped
    )
  )
}


# ==============================================================================
# 4.14 EXPERIMENT-SPECIFIC PATHWAY ENRICHMENT
# ==============================================================================

run_enrichment <- function(
    all_de,
    selected_de,
    direction,
    comparison_name,
    parent_dir
) {
  
  direction <- toupper(
    as.character(
      direction
    )
  )
  
  
  if (
    !direction %in%
    c(
      "UP",
      "DOWN"
    )
  ) {
    
    stop(
      "direction must be either 'UP' or 'DOWN'."
    )
  }
  
  
  if (
    !"gene" %in%
    colnames(
      all_de
    )
  ) {
    
    stop(
      "all_de has no 'gene' column."
    )
  }
  
  
  if (
    !"gene" %in%
    colnames(
      selected_de
    )
  ) {
    
    stop(
      "selected_de has no 'gene' column."
    )
  }
  
  
  selected_symbols <- unique(
    as.character(
      selected_de$gene
    )
  )
  
  
  selected_symbols <- selected_symbols[
    !is.na(
      selected_symbols
    ) &
      nzchar(
        trimws(
          selected_symbols
        )
      )
  ]
  
  
  universe_symbols <- unique(
    as.character(
      all_de$gene
    )
  )
  
  
  universe_symbols <- universe_symbols[
    !is.na(
      universe_symbols
    ) &
      nzchar(
        trimws(
          universe_symbols
        )
      )
  ]
  
  
  if (
    length(
      selected_symbols
    ) < 10
  ) {
    
    message(
      comparison_name,
      " / ",
      direction,
      ": fewer than 10 selected genes; enrichment skipped."
    )
    
    return(
      invisible(
        NULL
      )
    )
  }
  
  
  if (
    length(
      universe_symbols
    ) < 10
  ) {
    
    message(
      comparison_name,
      " / ",
      direction,
      ": tested-gene universe is too small; enrichment skipped."
    )
    
    return(
      invisible(
        NULL
      )
    )
  }
  
  
  selected_map <- map_symbols_to_entrez(
    selected_symbols
  )
  
  
  universe_map <- map_symbols_to_entrez(
    universe_symbols
  )
  
  
  selected_entrez <- unique(
    selected_map$ENTREZID
  )
  
  
  universe_entrez <- unique(
    universe_map$ENTREZID
  )
  
  
  selected_entrez <- intersect(
    selected_entrez,
    universe_entrez
  )
  
  
  if (
    length(
      selected_entrez
    ) < 10
  ) {
    
    message(
      comparison_name,
      " / ",
      direction,
      ": fewer than 10 genes mapped to the tested Entrez universe; enrichment skipped."
    )
    
    return(
      invisible(
        NULL
      )
    )
  }
  
  
  if (
    length(
      universe_entrez
    ) < 10
  ) {
    
    message(
      comparison_name,
      " / ",
      direction,
      ": mapped enrichment universe is too small; enrichment skipped."
    )
    
    return(
      invisible(
        NULL
      )
    )
  }
  
  
  output_directory <- file.path(
    parent_dir,
    safe_name(
      comparison_name
    ),
    direction
  )
  
  
  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  mapped_selected_symbols <- unique(
    selected_map$SYMBOL
  )
  
  
  mapped_universe_symbols <- unique(
    universe_map$SYMBOL
  )
  
  
  mapping_summary <- data.frame(
    
    comparison =
      comparison_name,
    
    direction =
      direction,
    
    selected_symbols =
      length(
        selected_symbols
      ),
    
    mapped_selected_symbols =
      length(
        mapped_selected_symbols
      ),
    
    selected_entrez =
      length(
        selected_entrez
      ),
    
    selected_mapping_percent =
      round(
        100 *
          length(
            mapped_selected_symbols
          ) /
          length(
            selected_symbols
          ),
        2
      ),
    
    universe_symbols =
      length(
        universe_symbols
      ),
    
    mapped_universe_symbols =
      length(
        mapped_universe_symbols
      ),
    
    universe_entrez =
      length(
        universe_entrez
      ),
    
    universe_mapping_percent =
      round(
        100 *
          length(
            mapped_universe_symbols
          ) /
          length(
            universe_symbols
          ),
        2
      ),
    
    stringsAsFactors =
      FALSE
  )
  
  
  write.csv(
    mapping_summary,
    file.path(
      output_directory,
      "mapping_summary.csv"
    ),
    row.names = FALSE
  )
  
  
  unmapped_selected <- setdiff(
    selected_symbols,
    mapped_selected_symbols
  )
  
  
  write.csv(
    data.frame(
      gene =
        unmapped_selected
    ),
    file.path(
      output_directory,
      "unmapped_selected_genes.csv"
    ),
    row.names = FALSE
  )
  
  
  # ---------------------------------------------------------------------------
  # GO Biological Process
  # ---------------------------------------------------------------------------
  
  go_result <- tryCatch(
    
    clusterProfiler::enrichGO(
      gene =
        selected_entrez,
      universe =
        universe_entrez,
      OrgDb =
        org.Hs.eg.db,
      keyType =
        "ENTREZID",
      ont =
        "BP",
      pAdjustMethod =
        "BH",
      pvalueCutoff =
        FDR_CUTOFF,
      qvalueCutoff =
        FDR_CUTOFF,
      readable =
        TRUE
    ),
    
    error = function(e) {
      
      warning(
        comparison_name,
        " / ",
        direction,
        " GO enrichment failed: ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  
  if (
    !is.null(
      go_result
    )
  ) {
    
    go_df <- as.data.frame(
      go_result
    )
    
    
    write.csv(
      go_df,
      file.path(
        output_directory,
        "GO_BP_results.csv"
      ),
      row.names = FALSE
    )
    
    
    if (
      nrow(
        go_df
      ) > 0
    ) {
      
      go_plot <- enrichplot::dotplot(
        go_result,
        showCategory =
          min(
            20,
            nrow(
              go_df
            )
          )
      ) +
        ggplot2::ggtitle(
          paste(
            comparison_name,
            direction,
            "GO Biological Process"
          )
        )
      
      
      save_pdf(
        go_plot,
        file.path(
          output_directory,
          "GO_BP_dotplot.pdf"
        ),
        width = 10,
        height = 8
      )
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # KEGG
  # ---------------------------------------------------------------------------
  
  kegg_result <- tryCatch(
    
    clusterProfiler::enrichKEGG(
      gene =
        selected_entrez,
      universe =
        universe_entrez,
      organism =
        "hsa",
      pvalueCutoff =
        FDR_CUTOFF,
      pAdjustMethod =
        "BH",
      qvalueCutoff =
        FDR_CUTOFF
    ),
    
    error = function(e) {
      
      warning(
        comparison_name,
        " / ",
        direction,
        " KEGG enrichment failed: ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  
  if (
    !is.null(
      kegg_result
    )
  ) {
    
    kegg_df <- as.data.frame(
      kegg_result
    )
    
    
    write.csv(
      kegg_df,
      file.path(
        output_directory,
        "KEGG_results.csv"
      ),
      row.names = FALSE
    )
    
    
    if (
      nrow(
        kegg_df
      ) > 0
    ) {
      
      kegg_plot <- enrichplot::dotplot(
        kegg_result,
        showCategory =
          min(
            20,
            nrow(
              kegg_df
            )
          )
      ) +
        ggplot2::ggtitle(
          paste(
            comparison_name,
            direction,
            "KEGG"
          )
        )
      
      
      save_pdf(
        kegg_plot,
        file.path(
          output_directory,
          "KEGG_dotplot.pdf"
        ),
        width = 10,
        height = 8
      )
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # Reactome
  # ---------------------------------------------------------------------------
  
  reactome_result <- tryCatch(
    
    ReactomePA::enrichPathway(
      gene =
        selected_entrez,
      universe =
        universe_entrez,
      organism =
        "human",
      pvalueCutoff =
        FDR_CUTOFF,
      pAdjustMethod =
        "BH",
      qvalueCutoff =
        FDR_CUTOFF,
      readable =
        TRUE
    ),
    
    error = function(e) {
      
      warning(
        comparison_name,
        " / ",
        direction,
        " Reactome enrichment failed: ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  
  if (
    !is.null(
      reactome_result
    )
  ) {
    
    reactome_df <- as.data.frame(
      reactome_result
    )
    
    
    write.csv(
      reactome_df,
      file.path(
        output_directory,
        "Reactome_results.csv"
      ),
      row.names = FALSE
    )
    
    
    if (
      nrow(
        reactome_df
      ) > 0
    ) {
      
      reactome_plot <- enrichplot::dotplot(
        reactome_result,
        showCategory =
          min(
            20,
            nrow(
              reactome_df
            )
          )
      ) +
        ggplot2::ggtitle(
          paste(
            comparison_name,
            direction,
            "Reactome"
          )
        )
      
      
      save_pdf(
        reactome_plot,
        file.path(
          output_directory,
          "Reactome_dotplot.pdf"
        ),
        width = 10,
        height = 8
      )
    }
  }
  
  
  return(
    invisible(
      list(
        GO_BP =
          go_result,
        KEGG =
          kegg_result,
        Reactome =
          reactome_result,
        mapping_summary =
          mapping_summary
      )
    )
  )
}


# ==============================================================================
# END SECTION 4
# ==============================================================================

cat(
  "\n====================================================\n",
  "Section 4: All helper functions defined successfully.\n",
  "====================================================\n"
)


# ==============================================================================
# 5. LOAD THE THREE RAW OBJECTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 5.1 Confirm that all three raw input files exist
# ------------------------------------------------------------------------------

raw_files <- c(
  Lea1 = RAW_LEA1,
  Lea2 = RAW_LEA2,
  Lea3 = RAW_LEA3
)


missing_raw_files <- raw_files[
  !file.exists(raw_files)
]


if (length(missing_raw_files) > 0) {
  
  stop(
    paste0(
      "The following raw Seurat file(s) were not found:\n",
      paste(
        names(missing_raw_files),
        "=",
        missing_raw_files,
        collapse = "\n"
      )
    )
  )
}


cat(
  "\n====================================================\n",
  "RAW INPUT FILES FOUND SUCCESSFULLY\n",
  "====================================================\n"
)

print(raw_files)


# ------------------------------------------------------------------------------
# 5.2 Load the three raw Seurat objects
# ------------------------------------------------------------------------------

cat("\nLoading Lea1 raw object...\n")

lea1_raw <- readRDS(
  RAW_LEA1
)


cat("Loading Lea2 raw object...\n")

lea2_raw <- readRDS(
  RAW_LEA2
)


cat("Loading Lea3 raw object...\n")

lea3_raw <- readRDS(
  RAW_LEA3
)


# ------------------------------------------------------------------------------
# 5.3 Confirm that all loaded objects are Seurat objects
# ------------------------------------------------------------------------------

if (!inherits(lea1_raw, "Seurat")) {
  
  stop(
    "Lea1 raw file was loaded, but it is not a Seurat object."
  )
}


if (!inherits(lea2_raw, "Seurat")) {
  
  stop(
    "Lea2 raw file was loaded, but it is not a Seurat object."
  )
}


if (!inherits(lea3_raw, "Seurat")) {
  
  stop(
    "Lea3 raw file was loaded, but it is not a Seurat object."
  )
}


cat(
  "\nAll three files were successfully loaded as Seurat objects.\n"
)


# ------------------------------------------------------------------------------
# 5.4 Check basic object dimensions
# ------------------------------------------------------------------------------

basic_object_summary <- data.frame(
  
  dataset = c(
    "Lea1",
    "Lea2",
    "Lea3"
  ),
  
  cells = c(
    ncol(lea1_raw),
    ncol(lea2_raw),
    ncol(lea3_raw)
  ),
  
  features = c(
    nrow(lea1_raw),
    nrow(lea2_raw),
    nrow(lea3_raw)
  ),
  
  stringsAsFactors = FALSE
)


cat(
  "\n====================================================\n",
  "RAW OBJECT DIMENSIONS\n",
  "====================================================\n"
)

print(
  basic_object_summary
)


# ------------------------------------------------------------------------------
# 5.5 Confirm required RNA assay exists
# ------------------------------------------------------------------------------

rna_assay_check <- c(
  
  Lea1 = any(
    names(lea1_raw@assays) == "RNA"
  ),
  
  Lea2 = any(
    names(lea2_raw@assays) == "RNA"
  ),
  
  Lea3 = any(
    names(lea3_raw@assays) == "RNA"
  )
)


if (!all(rna_assay_check)) {
  
  stop(
    paste0(
      "RNA assay is missing from: ",
      paste(
        names(
          rna_assay_check[
            !rna_assay_check
          ]
        ),
        collapse = ", "
      )
    )
  )
}


cat(
  "\nRNA assay found in all three objects.\n"
)

# ------------------------------------------------------------------------------
# 5.6 Confirm required Sample_Name metadata exists
# ------------------------------------------------------------------------------

sample_name_check <- c(
  
  Lea1 = "Sample_Name" %in%
    colnames(lea1_raw[[]]),
  
  Lea2 = "Sample_Name" %in%
    colnames(lea2_raw[[]]),
  
  Lea3 = "Sample_Name" %in%
    colnames(lea3_raw[[]])
  
)


if (!all(sample_name_check)) {
  
  stop(
    paste0(
      "Sample_Name metadata is missing from: ",
      paste(
        names(
          sample_name_check[
            !sample_name_check
          ]
        ),
        collapse = ", "
      )
    )
  )
}


cat(
  "Sample_Name metadata found in all three objects.\n"
)


# ------------------------------------------------------------------------------
# 5.7 Create input-file manifest
# ------------------------------------------------------------------------------

file_information <- file.info(
  raw_files
)


input_manifest <- data.frame(
  
  dataset = names(raw_files),
  
  file = unname(raw_files),
  
  size_MB = round(
    file_information$size /
      (1024^2),
    2
  ),
  
  modified = file_information$mtime,
  
  cells = c(
    ncol(lea1_raw),
    ncol(lea2_raw),
    ncol(lea3_raw)
  ),
  
  features = c(
    nrow(lea1_raw),
    nrow(lea2_raw),
    nrow(lea3_raw)
  ),
  
  RNA_assay = unname(
    rna_assay_check
  ),
  
  Sample_Name_metadata = unname(
    sample_name_check
  ),
  
  stringsAsFactors = FALSE
)


# Ensure log directory exists.

dir.create(
  DIR_LOGS,
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  input_manifest,
  file.path(
    DIR_LOGS,
    "input_manifest.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 5.8 Display raw cell counts
# ------------------------------------------------------------------------------

cat(
  "\n====================================================\n",
  "RAW CELL COUNTS\n",
  "====================================================\n"
)

cat(
  "Lea1:",
  format(
    ncol(lea1_raw),
    big.mark = ","
  ),
  "cells\n"
)

cat(
  "Lea2:",
  format(
    ncol(lea2_raw),
    big.mark = ","
  ),
  "cells\n"
)

cat(
  "Lea3:",
  format(
    ncol(lea3_raw),
    big.mark = ","
  ),
  "cells\n"
)


cat(
  "\nTotal raw cells:",
  format(
    ncol(lea1_raw) +
      ncol(lea2_raw) +
      ncol(lea3_raw),
    big.mark = ","
  ),
  "\n"
)


# ------------------------------------------------------------------------------
# 5.9 Display raw Sample_Name distributions
# ------------------------------------------------------------------------------

cat(
  "\n====================================================\n",
  "RAW Sample_Name DISTRIBUTIONS\n",
  "====================================================\n"
)


cat("\nLea1:\n")

print(
  table(
    lea1_raw$Sample_Name,
    useNA = "ifany"
  )
)


cat("\nLea2:\n")

print(
  table(
    lea2_raw$Sample_Name,
    useNA = "ifany"
  )
)


cat("\nLea3:\n")

print(
  table(
    lea3_raw$Sample_Name,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------------------------
# END SECTION 5
# ------------------------------------------------------------------------------

cat(
  "\n====================================================\n",
  "Section 5: Raw objects loaded and validated successfully.\n",
  "====================================================\n"
)


# ==============================================================================
# 6. SAMPLE METADATA CLEANING
# ==============================================================================

lea1 <- clean_metadata(
  lea1_raw,
  "Lea1"
)

lea2 <- clean_metadata(
  lea2_raw,
  "Lea2"
)

lea3 <- clean_metadata(
  lea3_raw,
  "Lea3"
)

rm(
  lea1_raw,
  lea2_raw,
  lea3_raw
)

invisible(
  gc()
)


# ==============================================================================
# 7. MITOCHONDRIAL PERCENTAGE
# ==============================================================================


# ------------------------------------------------------------------------------
# 7.1 Calculate mitochondrial RNA percentage
# ------------------------------------------------------------------------------

lea1 <- add_percent_mt(
  lea1,
  "Lea1"
)

lea2 <- add_percent_mt(
  lea2,
  "Lea2"
)

lea3 <- add_percent_mt(
  lea3,
  "Lea3"
)


# ------------------------------------------------------------------------------
# 7.2 Confirm percent.mt was created
# ------------------------------------------------------------------------------

stopifnot(
  "percent.mt" %in% colnames(lea1@meta.data),
  "percent.mt" %in% colnames(lea2@meta.data),
  "percent.mt" %in% colnames(lea3@meta.data)
)


# ------------------------------------------------------------------------------
# 7.3 Confirm no missing or invalid mitochondrial percentages
# ------------------------------------------------------------------------------

stopifnot(
  !anyNA(lea1$percent.mt),
  !anyNA(lea2$percent.mt),
  !anyNA(lea3$percent.mt)
)

stopifnot(
  all(lea1$percent.mt >= 0 & lea1$percent.mt <= 100),
  all(lea2$percent.mt >= 0 & lea2$percent.mt <= 100),
  all(lea3$percent.mt >= 0 & lea3$percent.mt <= 100)
)


# ------------------------------------------------------------------------------
# 7.4 Display mitochondrial-percentage summaries
# ------------------------------------------------------------------------------

cat(
  "\n====================================================\n",
  "MITOCHONDRIAL RNA PERCENTAGE SUMMARY\n",
  "====================================================\n"
)

cat("\nLea1:\n")
print(
  summary(
    lea1$percent.mt
  )
)

cat("\nLea2:\n")
print(
  summary(
    lea2$percent.mt
  )
)

cat("\nLea3:\n")
print(
  summary(
    lea3$percent.mt
  )
)


cat(
  "\nSection 7 completed successfully.\n\n"
)

# ==============================================================================
# 8. QC BEFORE FILTERING
# ==============================================================================


# ------------------------------------------------------------------------------
# 8.1 Generate QC summaries before filtering
# ------------------------------------------------------------------------------

qc_before <- dplyr::bind_rows(
  qc_summary(
    lea1,
    "Lea1",
    "Before_QC"
  ),
  qc_summary(
    lea2,
    "Lea2",
    "Before_QC"
  ),
  qc_summary(
    lea3,
    "Lea3",
    "Before_QC"
  )
)


# ------------------------------------------------------------------------------
# 8.2 Validate QC summary
# ------------------------------------------------------------------------------

stopifnot(
  nrow(qc_before) == 3
)

stopifnot(
  all(
    qc_before$dataset_id %in%
      c(
        "Lea1",
        "Lea2",
        "Lea3"
      )
  )
)

stopifnot(
  all(
    qc_before$cells > 0
  )
)


# ------------------------------------------------------------------------------
# 8.3 Display QC summary
# ------------------------------------------------------------------------------

cat(
  "\n====================================================\n",
  "QC SUMMARY BEFORE FILTERING\n",
  "====================================================\n"
)

print(
  qc_before
)


# ------------------------------------------------------------------------------
# 8.4 Save QC summary
# ------------------------------------------------------------------------------

write.csv(
  qc_before,
  file.path(
    DIR_QC,
    "QC_summary_BEFORE.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 8.5 Generate QC plots before filtering
# ------------------------------------------------------------------------------

save_qc_plots(
  lea1,
  "Lea1",
  "BEFORE"
)

save_qc_plots(
  lea2,
  "Lea2",
  "BEFORE"
)

save_qc_plots(
  lea3,
  "Lea3",
  "BEFORE"
)


cat(
  "\nSection 8 completed successfully.\n\n"
)


# ==============================================================================
# 8.6 QC THRESHOLD IMPACT BEFORE FILTERING
# ==============================================================================

qc_threshold_impact <- function(
    obj,
    dataset_id
) {
  
  low_features <- (
    obj$nFeature_RNA <=
      QC_MIN_FEATURES
  )
  
  high_features <- (
    obj$nFeature_RNA >=
      QC_MAX_FEATURES
  )
  
  high_mt <- (
    obj$percent.mt >=
      QC_MAX_MT
  )
  
  fail_any <- (
    low_features |
      high_features |
      high_mt
  )
  
  
  data.frame(
    
    dataset_id =
      dataset_id,
    
    total_cells =
      ncol(obj),
    
    low_features =
      sum(low_features),
    
    high_features =
      sum(high_features),
    
    high_mitochondrial =
      sum(high_mt),
    
    fail_any_QC =
      sum(fail_any),
    
    retain_after_QC =
      sum(!fail_any),
    
    percent_removed =
      round(
        100 *
          sum(fail_any) /
          ncol(obj),
        2
      ),
    
    percent_retained =
      round(
        100 *
          sum(!fail_any) /
          ncol(obj),
        2
      ),
    
    stringsAsFactors = FALSE
  )
}


qc_threshold_impact_all <- dplyr::bind_rows(
  
  qc_threshold_impact(
    lea1,
    "Lea1"
  ),
  
  qc_threshold_impact(
    lea2,
    "Lea2"
  ),
  
  qc_threshold_impact(
    lea3,
    "Lea3"
  )
)


cat(
  "\n====================================================\n",
  "QC THRESHOLD IMPACT BEFORE FILTERING\n",
  "====================================================\n"
)

print(
  qc_threshold_impact_all
)


write.csv(
  qc_threshold_impact_all,
  file.path(
    DIR_QC,
    "QC_threshold_impact_BEFORE_filtering.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 8.7 QC THRESHOLD IMPACT BY TREATMENT
# ==============================================================================

qc_impact_by_treatment <- function(
    obj,
    dataset_id
) {
  
  qc_df <- data.frame(
    
    dataset_id = dataset_id,
    
    Sample_Name = as.character(
      obj$Sample_Name
    ),
    
    low_features = (
      obj$nFeature_RNA <=
        QC_MIN_FEATURES
    ),
    
    high_features = (
      obj$nFeature_RNA >=
        QC_MAX_FEATURES
    ),
    
    high_mt = (
      obj$percent.mt >=
        QC_MAX_MT
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  qc_df$fail_any <- (
    qc_df$low_features |
      qc_df$high_features |
      qc_df$high_mt
  )
  
  
  qc_summary_by_treatment <- qc_df %>%
    dplyr::group_by(
      dataset_id,
      Sample_Name
    ) %>%
    dplyr::summarise(
      
      total_cells =
        dplyr::n(),
      
      low_features =
        sum(low_features),
      
      high_features =
        sum(high_features),
      
      high_mt =
        sum(high_mt),
      
      fail_any_QC =
        sum(fail_any),
      
      retain_after_QC =
        sum(!fail_any),
      
      percent_removed = round(
        100 *
          sum(fail_any) /
          dplyr::n(),
        2
      ),
      
      percent_retained = round(
        100 *
          sum(!fail_any) /
          dplyr::n(),
        2
      ),
      
      .groups = "drop"
    )
  
  
  return(
    qc_summary_by_treatment
  )
}


qc_impact_by_treatment_all <- dplyr::bind_rows(
  
  qc_impact_by_treatment(
    lea1,
    "Lea1"
  ),
  
  qc_impact_by_treatment(
    lea2,
    "Lea2"
  ),
  
  qc_impact_by_treatment(
    lea3,
    "Lea3"
  )
)


cat(
  "\n====================================================\n",
  "QC IMPACT BY TREATMENT\n",
  "====================================================\n"
)

print(
  qc_impact_by_treatment_all,
  n = Inf
)


write.csv(
  qc_impact_by_treatment_all,
  file.path(
    DIR_QC,
    "QC_threshold_impact_BY_TREATMENT.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. QC FILTERING
# ==============================================================================


# ------------------------------------------------------------------------------
# 9.1 Apply QC filters
# ------------------------------------------------------------------------------

lea1 <- filter_qc(
  lea1,
  "Lea1"
)

lea2 <- filter_qc(
  lea2,
  "Lea2"
)

lea3 <- filter_qc(
  lea3,
  "Lea3"
)


# ------------------------------------------------------------------------------
# 9.2 Verify excluded treatments remain absent
# ------------------------------------------------------------------------------

stopifnot(
  !any(
    as.character(lea1$Sample_Name) %in%
      EXCLUDED_TREATMENTS
  ),
  !any(
    as.character(lea2$Sample_Name) %in%
      EXCLUDED_TREATMENTS
  ),
  !any(
    as.character(lea3$Sample_Name) %in%
      EXCLUDED_TREATMENTS
  )
)


# ------------------------------------------------------------------------------
# 9.3 Display cells retained after QC
# ------------------------------------------------------------------------------

cat(
  "\n====================================================\n",
  "CELL COUNTS AFTER QC FILTERING\n",
  "====================================================\n"
)


cat(
  "Lea1:",
  format(
    ncol(lea1),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Lea2:",
  format(
    ncol(lea2),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Lea3:",
  format(
    ncol(lea3),
    big.mark = ","
  ),
  "\n"
)


cat(
  "\nTotal cells after QC:",
  format(
    ncol(lea1) +
      ncol(lea2) +
      ncol(lea3),
    big.mark = ","
  ),
  "\n"
)


cat(
  "\nSection 9 completed successfully.\n\n"
)