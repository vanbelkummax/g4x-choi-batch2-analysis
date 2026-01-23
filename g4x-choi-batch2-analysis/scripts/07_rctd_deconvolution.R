#!/usr/bin/env Rscript
# RCTD Deconvolution for G4X Gastric Cancer Progression
# Using GSE134520 reference (44K cells, 18 cell types)
#
# Critical fixes applied:
# 1. Pin Python environment with reticulate
# 2. sparseMatrix with 1-based indexing (indices + 1L)
# 3. nUMI from FULL library BEFORE gene intersection
# 4. Weights reindexed to results_df with NA fill
# 5. Use layers['counts'] for G4X (X is log-transformed)

suppressPackageStartupMessages({
    library(spacexr)
    library(Matrix)
    library(reticulate)
    library(dplyr)
    library(jsonlite)
})

# CRITICAL: Pin Python environment with full path
conda_path <- file.path(Sys.getenv("HOME"), "miniforge3/envs/enact")
Sys.setenv(RETICULATE_PYTHON = file.path(conda_path, "bin/python"))
use_condaenv(conda_path, required = TRUE)

# ============================================================================
# Configuration
# ============================================================================
config <- list(
    g4x_path = "results/pilot/merged_pilot.h5ad",
    ref_path = "output/GSE134520_reference.h5ad",
    output_dir = "output",

    # RCTD parameters (relaxed for 337-gene panel)
    gene_cutoff = 0.000125,
    fc_cutoff = 0.3,
    gene_cutoff_reg = 2e-04,
    fc_cutoff_reg = 0.5,
    UMI_min = 1,
    CELL_MIN_INSTANCE = 3,

    # Processing
    max_cores = 12,
    doublet_mode = "doublet"  # or "full" or "multi"
)

cat("=== RCTD Deconvolution Pipeline ===\n")
cat(sprintf("G4X data: %s\n", config$g4x_path))
cat(sprintf("Reference: %s\n", config$ref_path))

# ============================================================================
# Helper Functions
# ============================================================================

#' Read h5ad file using anndata via reticulate
#' @param path Path to h5ad file
#' @param use_layer Layer name for counts (NULL = use X)
#' @return List with counts matrix, cell metadata, gene names
read_h5ad_sparse <- function(path, use_layer = NULL) {
    cat(sprintf("Reading: %s\n", path))

    ad <- import("anndata")
    np <- import("numpy")

    adata_py <- ad$read_h5ad(path)

    # Get counts matrix
    # Check if layer exists using Python's __contains__
    builtins <- import_builtins()
    has_layer <- FALSE
    if (!is.null(use_layer)) {
        layer_keys <- builtins$list(adata_py$layers$keys())
        has_layer <- use_layer %in% layer_keys
    }

    if (has_layer) {
        cat(sprintf("  Using layer: %s\n", use_layer))
        counts_py <- adata_py$layers[[use_layer]]
    } else {
        cat("  Using X matrix\n")
        counts_py <- adata_py$X
    }

    # Check if sparse
    scipy_sparse <- import("scipy.sparse")
    is_sparse <- scipy_sparse$issparse(counts_py)

    if (is_sparse) {
        # Convert to CSC for efficient column operations
        counts_csc <- scipy_sparse$csc_matrix(counts_py)

        # Extract sparse matrix components
        data <- as.numeric(counts_csc$data)
        indices <- as.integer(counts_csc$indices)
        indptr <- as.integer(counts_csc$indptr)
        shape <- as.integer(counts_csc$shape)

        # FIXED: sparseMatrix with 1-based indexing
        counts <- sparseMatrix(
            i = indices + 1L,  # Python 0-based → R 1-based
            p = indptr,
            x = data,
            dims = shape,
            repr = "C"
        )
    } else {
        # Dense matrix
        counts <- as.matrix(np$array(counts_py))
        counts <- Matrix(counts, sparse = TRUE)
    }

    # Get gene and cell names
    genes <- as.character(np$array(adata_py$var_names))
    cells <- as.character(np$array(adata_py$obs_names))

    rownames(counts) <- cells
    colnames(counts) <- genes

    # Transpose: cells x genes → genes x cells (RCTD expects genes x cells)
    counts <- t(counts)

    # Get cell metadata
    obs_df <- as.data.frame(adata_py$obs)
    rownames(obs_df) <- cells

    cat(sprintf("  Shape: %d genes x %d cells\n", nrow(counts), ncol(counts)))

    list(
        counts = counts,
        obs = obs_df,
        genes = genes,
        cells = cells
    )
}

# ============================================================================
# Load Data
# ============================================================================

# Load G4X data (use 'counts' layer - X is log-transformed!)
cat("\n--- Loading G4X Data ---\n")
g4x <- read_h5ad_sparse(config$g4x_path, use_layer = "counts")

# Load reference
cat("\n--- Loading Reference Data ---\n")
ref <- read_h5ad_sparse(config$ref_path, use_layer = NULL)

# ============================================================================
# Prepare Reference
# ============================================================================

cat("\n--- Preparing Reference ---\n")

# Get cell types
if ("cell_type" %in% colnames(ref$obs)) {
    cell_types <- factor(ref$obs$cell_type)
} else {
    stop("Reference must have 'cell_type' column in obs")
}
names(cell_types) <- ref$cells

cat(sprintf("Cell types in reference:\n"))
print(table(cell_types))

# CRITICAL: nUMI from FULL counts BEFORE gene intersection
nUMI_ref_full <- colSums(ref$counts)
cat(sprintf("Reference nUMI range: %.0f - %.0f\n", min(nUMI_ref_full), max(nUMI_ref_full)))

# ============================================================================
# Gene Intersection
# ============================================================================

cat("\n--- Computing Gene Overlap ---\n")

g4x_genes <- rownames(g4x$counts)
ref_genes <- rownames(ref$counts)

common_genes <- intersect(g4x_genes, ref_genes)
cat(sprintf("G4X genes: %d\n", length(g4x_genes)))
cat(sprintf("Reference genes: %d\n", length(ref_genes)))
cat(sprintf("Common genes: %d\n", length(common_genes)))

if (length(common_genes) < 100) {
    stop("Too few common genes for reliable deconvolution!")
}

# Save gene overlap stats
gene_stats <- list(
    n_g4x = length(g4x_genes),
    n_ref = length(ref_genes),
    n_common = length(common_genes),
    common_genes = common_genes,
    g4x_only = setdiff(g4x_genes, ref_genes),
    ref_only = setdiff(ref_genes, g4x_genes)[1:100]  # First 100
)
write_json(gene_stats, file.path(config$output_dir, "gene_overlap_stats.json"),
           auto_unbox = TRUE, pretty = TRUE)

# Subset to common genes
ref_counts_common <- ref$counts[common_genes, ]
g4x_counts_common <- g4x$counts[common_genes, ]

# ============================================================================
# Create RCTD Reference Object
# ============================================================================

cat("\n--- Creating RCTD Reference ---\n")

# CRITICAL: Use nUMI from FULL counts (before subsetting)
reference <- Reference(
    counts = ref_counts_common,
    cell_types = cell_types,
    nUMI = nUMI_ref_full,
    min_UMI = config$UMI_min,
    require_int = FALSE
)

cat(sprintf("Reference created: %d cell types, %d cells\n",
            length(unique(cell_types)), ncol(ref_counts_common)))

# Save reference checkpoint
saveRDS(reference, file.path(config$output_dir, "rctd_reference.rds"))
cat("Reference saved to: output/rctd_reference.rds\n")

# ============================================================================
# Create RCTD Query Object (SpatialRNA)
# ============================================================================

cat("\n--- Creating SpatialRNA Query ---\n")

# Get spatial coordinates
if ("spatial" %in% names(g4x$obs)) {
    # This won't work directly - need to handle differently
    coords_x <- g4x$obs$cell_x
    coords_y <- g4x$obs$cell_y
} else if ("cell_x" %in% colnames(g4x$obs)) {
    coords_x <- g4x$obs$cell_x
    coords_y <- g4x$obs$cell_y
} else {
    stop("No spatial coordinates found in G4X data!")
}

coords <- data.frame(
    x = as.numeric(coords_x),
    y = as.numeric(coords_y),
    row.names = g4x$cells
)

# nUMI for query
nUMI_g4x <- colSums(g4x$counts)

# Create SpatialRNA object
query <- SpatialRNA(
    coords = coords,
    counts = g4x_counts_common,
    nUMI = nUMI_g4x,
    require_int = FALSE
)

cat(sprintf("Query created: %d cells, %d genes\n",
            ncol(g4x_counts_common), nrow(g4x_counts_common)))

# ============================================================================
# Run RCTD
# ============================================================================

cat("\n--- Running RCTD Deconvolution ---\n")
cat(sprintf("Parameters:\n"))
cat(sprintf("  gene_cutoff: %.6f\n", config$gene_cutoff))
cat(sprintf("  fc_cutoff: %.2f\n", config$fc_cutoff))
cat(sprintf("  CELL_MIN_INSTANCE: %d\n", config$CELL_MIN_INSTANCE))
cat(sprintf("  doublet_mode: %s\n", config$doublet_mode))
cat(sprintf("  max_cores: %d\n", config$max_cores))

# Create RCTD object
rctd <- create.RCTD(
    spatialRNA = query,
    reference = reference,
    max_cores = config$max_cores,
    gene_cutoff = config$gene_cutoff,
    fc_cutoff = config$fc_cutoff,
    gene_cutoff_reg = config$gene_cutoff_reg,
    fc_cutoff_reg = config$fc_cutoff_reg,
    UMI_min = config$UMI_min,
    CELL_MIN_INSTANCE = config$CELL_MIN_INSTANCE
)

# Run RCTD
start_time <- Sys.time()
rctd <- run.RCTD(rctd, doublet_mode = config$doublet_mode)
end_time <- Sys.time()

cat(sprintf("RCTD completed in: %.1f minutes\n",
            as.numeric(difftime(end_time, start_time, units = "mins"))))

# Save RCTD object
saveRDS(rctd, file.path(config$output_dir, "rctd_results.rds"))
cat("RCTD results saved to: output/rctd_results.rds\n")

# ============================================================================
# Extract Results
# ============================================================================

cat("\n--- Extracting Results ---\n")

# Get results
results <- rctd@results

# Get all cell IDs from query
all_cells <- g4x$cells

# Initialize results data frame
results_df <- data.frame(
    cell_id = all_cells,
    stringsAsFactors = FALSE
)

# Get cell type assignments
if (config$doublet_mode == "doublet") {
    # Doublet mode: first_type, second_type, spot_class
    results_df$cell_type_rctd <- sapply(all_cells, function(cell) {
        if (cell %in% names(results$results_df$first_type)) {
            as.character(results$results_df[cell, "first_type"])
        } else {
            NA_character_
        }
    })

    results_df$second_type <- sapply(all_cells, function(cell) {
        if (cell %in% names(results$results_df$second_type)) {
            as.character(results$results_df[cell, "second_type"])
        } else {
            NA_character_
        }
    })

    results_df$spot_class <- sapply(all_cells, function(cell) {
        if (cell %in% rownames(results$results_df)) {
            as.character(results$results_df[cell, "spot_class"])
        } else {
            NA_character_
        }
    })
} else {
    # Full mode: just get first type
    results_df$cell_type_rctd <- sapply(all_cells, function(cell) {
        if (cell %in% names(results$results_df$first_type)) {
            as.character(results$results_df[cell, "first_type"])
        } else {
            NA_character_
        }
    })
}

# Get weights matrix
if (!is.null(results$weights)) {
    weights <- as.matrix(results$weights)

    # FIXED: Reindex weights to match results_df$cell_id
    weights_aligned <- matrix(
        NA_real_,
        nrow = length(all_cells),
        ncol = ncol(weights),
        dimnames = list(all_cells, colnames(weights))
    )

    # Fill in available weights
    common_cells <- intersect(all_cells, rownames(weights))
    weights_aligned[common_cells, ] <- weights[common_cells, ]

    # Replace NA with 0 for cells without weights
    weights_aligned[is.na(weights_aligned)] <- 0

    weights <- weights_aligned
} else {
    weights <- NULL
}

# ============================================================================
# Save Results
# ============================================================================

cat("\n--- Saving Results ---\n")

# Save annotations
write.csv(results_df, file.path(config$output_dir, "g4x_rctd_annotations.csv"),
          row.names = FALSE)
cat("Annotations saved to: output/g4x_rctd_annotations.csv\n")

# Save weights
if (!is.null(weights)) {
    weights_df <- as.data.frame(weights)
    weights_df$cell_id <- rownames(weights)
    write.csv(weights_df, file.path(config$output_dir, "g4x_rctd_weights.csv"),
              row.names = FALSE)
    cat("Weights saved to: output/g4x_rctd_weights.csv\n")
}

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n=== RCTD Summary ===\n")

# Cell type distribution
cat("\nCell type distribution:\n")
print(table(results_df$cell_type_rctd, useNA = "ifany"))

# Spot class distribution (doublet mode)
if ("spot_class" %in% colnames(results_df)) {
    cat("\nSpot class distribution:\n")
    print(table(results_df$spot_class, useNA = "ifany"))
}

# Success rate
n_total <- nrow(results_df)
n_assigned <- sum(!is.na(results_df$cell_type_rctd))
cat(sprintf("\nAssignment rate: %d / %d (%.1f%%)\n",
            n_assigned, n_total, 100 * n_assigned / n_total))

# Rejection rate
if ("spot_class" %in% colnames(results_df)) {
    n_reject <- sum(results_df$spot_class == "reject", na.rm = TRUE)
    cat(sprintf("Rejection rate: %d / %d (%.1f%%)\n",
                n_reject, n_total, 100 * n_reject / n_total))
}

# Per-sample summary
if ("sample_id" %in% colnames(g4x$obs)) {
    cat("\nPer-sample cell type distribution:\n")
    results_df$sample_id <- g4x$obs[results_df$cell_id, "sample_id"]
    print(table(results_df$sample_id, results_df$cell_type_rctd))
}

cat("\n=== RCTD Pipeline Complete ===\n")
