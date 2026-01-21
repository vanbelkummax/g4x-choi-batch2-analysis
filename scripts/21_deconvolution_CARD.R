#!/usr/bin/env Rscript
#' =============================================================================
#' Script 21: CARD Deconvolution (Spatial-Aware)
#' =============================================================================
#' CARD uses Conditional Autoregressive (CAR) spatial prior
#' - Neighboring spots influence each other's composition estimates
#' - More accurate than reference-free methods for spatial data
#' - Handles PDAC's complex TME with dense stroma

suppressPackageStartupMessages({
    library(CARD)
    library(Seurat)
    library(Matrix)
    library(dplyr)
    library(ggplot2)
})

#' -----------------------------------------------------------------------------
#' CONFIGURATION
#' -----------------------------------------------------------------------------
DATA_DIR <- "/home/user/spatial-hackathon-2026/outputs/adata/polymathic"
OUTPUT_DIR <- "/home/user/spatial-hackathon-2026/outputs/deconvolution"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Sample metadata
METADATA <- data.frame(
    sample_id = c("YP12A", "YP14A", "YP16A", "YP17A", "YP13A", "YP15A", "YP18A", "YP20A"),
    response = c("R", "R", "R", "R", "NR", "NR", "NR", "NR"),
    stringsAsFactors = FALSE
)

#' -----------------------------------------------------------------------------
#' HELPER FUNCTIONS
#' -----------------------------------------------------------------------------

load_spatial_data <- function(h5ad_path) {
    #' Load spatial data from h5ad and convert to CARD format
    #' Returns: list with counts matrix and coordinates

    # Use anndata2ri or manual conversion
    # For now, we'll read from Seurat objects if available

    # Try loading from Seurat object path
    seurat_path <- gsub("_polymathic.h5ad$", "_seurat.rds", h5ad_path)

    if (file.exists(seurat_path)) {
        message(sprintf("Loading Seurat object: %s", basename(seurat_path)))
        obj <- readRDS(seurat_path)

        counts <- GetAssayData(obj, slot = "counts")

        # Get spatial coordinates
        if ("spatial" %in% names(obj@images)) {
            coords <- GetTissueCoordinates(obj)
        } else {
            # Try to extract from metadata
            coords <- obj@meta.data[, c("array_row", "array_col")]
            colnames(coords) <- c("x", "y")
        }

        return(list(
            counts = counts,
            coords = as.data.frame(coords)
        ))
    }

    stop(sprintf("Seurat object not found: %s", seurat_path))
}

load_reference_data <- function(ref_path = NULL) {
    #' Load scRNA-seq reference for deconvolution
    #' Returns: list with counts matrix and cell type annotations

    # Use PDAC reference from Peng et al. or similar
    # For hackathon, we'll use a curated reference

    if (is.null(ref_path)) {
        # Default: use PDAC scRNA-seq reference from Polymath KB
        ref_path <- "/home/user/data/hackathon/references/pdac_scrna_reference.rds"
    }

    if (!file.exists(ref_path)) {
        message("Reference not found. Creating synthetic reference from markers...")
        return(create_synthetic_reference())
    }

    ref <- readRDS(ref_path)

    return(list(
        counts = GetAssayData(ref, slot = "counts"),
        meta = ref@meta.data,
        cell_type_col = "cell_type"
    ))
}

create_synthetic_reference <- function() {
    #' Create a synthetic reference based on marker genes
    #' Used when no external reference is available

    # PDAC cell type markers (from literature)
    markers <- list(
        "Ductal_cancer" = c("KRT19", "KRT7", "EPCAM", "MUC1", "CEACAM5"),
        "Acinar" = c("PRSS1", "CELA3A", "CPA1", "AMY2A", "PNLIP"),
        "Ductal_normal" = c("KRT19", "KRT7", "CFTR", "AQP1"),
        "CAF_myCAF" = c("ACTA2", "TAGLN", "MYH11", "ACTG2"),
        "CAF_iCAF" = c("IL6", "CXCL1", "CXCL12", "LIF", "PDGFRA"),
        "Endothelial" = c("PECAM1", "CDH5", "VWF", "CLDN5"),
        "Macrophage" = c("CD68", "CD163", "CSF1R", "MRC1"),
        "T_cell" = c("CD3D", "CD3E", "CD4", "CD8A"),
        "B_cell" = c("CD19", "MS4A1", "CD79A", "PAX5"),
        "Stellate" = c("RGS5", "PDGFRB", "MCAM", "CSPG4")
    )

    message("Note: Using marker-based pseudo-reference. Consider using real scRNA-seq reference.")

    return(list(
        markers = markers,
        synthetic = TRUE
    ))
}

run_card_deconvolution <- function(sp_counts, sp_coords, sc_counts, sc_meta,
                                    cell_type_col = "cell_type") {
    #' Run CARD deconvolution with spatial prior
    #'
    #' @param sp_counts Spatial counts matrix (genes x spots)
    #' @param sp_coords Spatial coordinates (x, y)
    #' @param sc_counts scRNA-seq counts matrix (genes x cells)
    #' @param sc_meta Cell metadata with cell type annotations
    #' @param cell_type_col Column name for cell types
    #' @return CARD object with deconvolution results

    message("Creating CARD object...")

    # Ensure coordinates are data frame with correct format
    sp_coords <- as.data.frame(sp_coords)
    colnames(sp_coords) <- c("x", "y")
    rownames(sp_coords) <- colnames(sp_counts)

    # Create CARD object
    CARD_obj <- createCARDObject(
        sc_count = sc_counts,
        sc_meta = sc_meta,
        spatial_count = sp_counts,
        spatial_location = sp_coords,
        ct.varname = cell_type_col,
        ct.select = NULL,  # Use all cell types
        sample.varname = "sample",
        minCountGene = 100,
        minCountSpot = 5
    )

    message("Running CARD deconvolution...")
    CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)

    return(CARD_obj)
}

visualize_card_results <- function(CARD_obj, sample_name, output_dir) {
    #' Visualize CARD deconvolution results

    message(sprintf("Generating visualizations for %s...", sample_name))

    # Get proportions
    props <- CARD_obj@Proportion_CARD

    # 1. Pie chart of overall composition
    avg_props <- colMeans(props)

    pdf(file.path(output_dir, sprintf("%s_card_composition.pdf", sample_name)),
        width = 10, height = 8)

    pie_df <- data.frame(
        cell_type = names(avg_props),
        proportion = avg_props
    ) %>%
        filter(proportion > 0.01) %>%
        arrange(desc(proportion))

    p1 <- ggplot(pie_df, aes(x = "", y = proportion, fill = cell_type)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y") +
        theme_minimal() +
        labs(title = sprintf("%s - Cell Type Composition (CARD)", sample_name),
             fill = "Cell Type") +
        theme(axis.text = element_blank())
    print(p1)

    # 2. Spatial visualization of top cell types
    coords <- CARD_obj@spatial_location

    top_types <- names(sort(avg_props, decreasing = TRUE)[1:4])

    for (ct in top_types) {
        spatial_df <- data.frame(
            x = coords$x,
            y = coords$y,
            proportion = props[, ct]
        )

        p <- ggplot(spatial_df, aes(x = x, y = y, color = proportion)) +
            geom_point(size = 1.5) +
            scale_color_viridis_c(option = "plasma") +
            theme_minimal() +
            coord_fixed() +
            labs(title = sprintf("%s - %s Proportion", sample_name, ct),
                 color = "Proportion") +
            theme(axis.title = element_blank())
        print(p)
    }

    dev.off()

    message(sprintf("Saved: %s_card_composition.pdf", sample_name))
}

#' -----------------------------------------------------------------------------
#' MAIN
#' -----------------------------------------------------------------------------

main <- function() {
    message("=" |> rep(60) |> paste(collapse = ""))
    message("CARD Deconvolution for PDAC Spatial Transcriptomics")
    message("=" |> rep(60) |> paste(collapse = ""))

    # Load reference data
    ref <- load_reference_data()

    if (isTRUE(ref$synthetic)) {
        message("\nWARNING: Using synthetic reference. Results will be approximate.")
        message("For best results, provide PDAC scRNA-seq reference.")
        message("Suggested: Peng et al. 2019 (GSE111672)")

        # Skip full deconvolution with synthetic reference
        # Instead, run marker-based analysis

        results <- list()

        for (i in 1:nrow(METADATA)) {
            sample_id <- METADATA$sample_id[i]
            response <- METADATA$response[i]

            h5ad_path <- file.path(DATA_DIR, sprintf("%s_polymathic.h5ad", sample_id))

            if (!file.exists(h5ad_path)) {
                message(sprintf("Skipping %s (file not found)", sample_id))
                next
            }

            message(sprintf("\nProcessing %s (%s)...", sample_id, response))

            # For now, just record that we need proper reference
            results[[sample_id]] <- list(
                sample_id = sample_id,
                response = response,
                status = "needs_reference"
            )
        }

        # Save status
        saveRDS(results, file.path(OUTPUT_DIR, "card_status_needs_reference.rds"))

        message("\n" |> paste(rep("=", 60), collapse = ""))
        message("ACTION REQUIRED:")
        message("1. Download PDAC scRNA-seq reference (Peng et al. 2019)")
        message("2. Save to: /home/user/data/hackathon/references/pdac_scrna_reference.rds")
        message("3. Re-run this script")
        message(rep("=", 60) |> paste(collapse = ""))

        return(invisible(results))
    }

    # Process each sample
    all_results <- list()

    for (i in 1:nrow(METADATA)) {
        sample_id <- METADATA$sample_id[i]
        response <- METADATA$response[i]

        h5ad_path <- file.path(DATA_DIR, sprintf("%s_polymathic.h5ad", sample_id))

        if (!file.exists(h5ad_path)) {
            message(sprintf("Skipping %s (file not found)", sample_id))
            next
        }

        message(sprintf("\n[%d/%d] Processing %s (%s)...",
                        i, nrow(METADATA), sample_id, response))

        tryCatch({
            # Load spatial data
            sp_data <- load_spatial_data(h5ad_path)

            # Run CARD
            CARD_obj <- run_card_deconvolution(
                sp_counts = sp_data$counts,
                sp_coords = sp_data$coords,
                sc_counts = ref$counts,
                sc_meta = ref$meta,
                cell_type_col = ref$cell_type_col
            )

            # Visualize
            visualize_card_results(CARD_obj, sample_id, OUTPUT_DIR)

            # Save results
            saveRDS(CARD_obj,
                    file.path(OUTPUT_DIR, sprintf("%s_card_results.rds", sample_id)))

            # Store proportions
            all_results[[sample_id]] <- list(
                proportions = CARD_obj@Proportion_CARD,
                response = response,
                status = "success"
            )

            message(sprintf("✓ %s completed", sample_id))

        }, error = function(e) {
            message(sprintf("✗ %s failed: %s", sample_id, e$message))
            all_results[[sample_id]] <- list(
                response = response,
                status = "failed",
                error = e$message
            )
        })
    }

    # Save combined results
    saveRDS(all_results, file.path(OUTPUT_DIR, "card_all_results.rds"))

    # Summary statistics
    message("\n" |> paste(rep("=", 60), collapse = ""))
    message("SUMMARY")
    message(rep("=", 60) |> paste(collapse = ""))

    success <- sum(sapply(all_results, function(x) x$status == "success"))
    message(sprintf("Processed: %d/%d samples", success, nrow(METADATA)))

    if (success > 0) {
        # Compare R vs NR
        R_samples <- names(all_results)[sapply(all_results, function(x)
            x$response == "R" && x$status == "success")]
        NR_samples <- names(all_results)[sapply(all_results, function(x)
            x$response == "NR" && x$status == "success")]

        if (length(R_samples) > 0 && length(NR_samples) > 0) {
            message("\nCell type proportions (mean):")

            R_props <- do.call(rbind, lapply(R_samples, function(s)
                colMeans(all_results[[s]]$proportions)))
            NR_props <- do.call(rbind, lapply(NR_samples, function(s)
                colMeans(all_results[[s]]$proportions)))

            comparison <- data.frame(
                cell_type = colnames(R_props),
                R_mean = colMeans(R_props),
                NR_mean = colMeans(NR_props)
            ) %>%
                mutate(
                    diff = R_mean - NR_mean,
                    enriched_in = ifelse(diff > 0, "R", "NR")
                ) %>%
                arrange(desc(abs(diff)))

            print(comparison)

            write.csv(comparison,
                      file.path(OUTPUT_DIR, "card_r_vs_nr_comparison.csv"),
                      row.names = FALSE)
        }
    }

    message("\nResults saved to: ", OUTPUT_DIR)
    return(invisible(all_results))
}

# Run if called directly
if (!interactive()) {
    main()
}
