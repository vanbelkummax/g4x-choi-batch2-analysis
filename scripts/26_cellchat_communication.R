#!/usr/bin/env Rscript
#' =============================================================================
#' Script 26: CellChat Cell-Cell Communication Analysis
#' =============================================================================
#' Analyzes ligand-receptor interactions between cell types in PDAC samples
#' Compares communication patterns between Responders and Non-Responders

suppressPackageStartupMessages({
    library(CellChat)
    library(Seurat)
    library(SeuratDisk)
    library(patchwork)
    library(ggplot2)
    library(dplyr)
})

# Configuration
DATA_DIR <- "/home/user/spatial-hackathon-2026/outputs/adata/annotated"
OUTPUT_DIR <- "/home/user/spatial-hackathon-2026/outputs/cellchat"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Sample metadata
METADATA <- data.frame(
    sample_id = c("YP12A", "YP15A", "YP12C", "YP15C", "YP03A", "YP03C", "YP04C"),
    response = c("R", "R", "R", "R", "NR", "NR", "NR"),
    timepoint = c("Pre", "Pre", "Post", "Post", "Pre", "Post", "Post"),
    stringsAsFactors = FALSE
)

#' Convert h5ad to Seurat object
h5ad_to_seurat <- function(h5ad_path) {
    message(sprintf("Converting %s to Seurat...", basename(h5ad_path)))

    # Convert via intermediate h5seurat
    h5seurat_path <- gsub("\\.h5ad$", ".h5seurat", h5ad_path)

    tryCatch({
        Convert(h5ad_path, dest = "h5seurat", overwrite = TRUE)
        obj <- LoadH5Seurat(h5seurat_path)
        return(obj)
    }, error = function(e) {
        message(sprintf("Conversion failed: %s", e$message))
        return(NULL)
    })
}

#' Run CellChat on a single sample
run_cellchat_single <- function(seurat_obj, sample_name, cell_type_col = "cell_type") {
    message(sprintf("\nRunning CellChat for %s...", sample_name))

    # Get expression data
    data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

    # Get cell type labels
    if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
        # Try alternatives
        alt_cols <- c("celltype", "CellType", "cluster", "leiden")
        for (col in alt_cols) {
            if (col %in% colnames(seurat_obj@meta.data)) {
                cell_type_col <- col
                break
            }
        }
    }

    labels <- seurat_obj@meta.data[[cell_type_col]]

    if (is.null(labels)) {
        message("  No cell type labels found, skipping...")
        return(NULL)
    }

    meta <- data.frame(labels = labels, row.names = colnames(data.input))

    # Create CellChat object
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

    # Set ligand-receptor database (human)
    CellChatDB <- CellChatDB.human
    CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
    cellchat@DB <- CellChatDB.use

    # Preprocessing
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    # Inference
    message("  Computing communication probability...")
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)

    return(cellchat)
}

#' Compare CellChat results between groups
compare_cellchat <- function(cellchat_list, group_labels) {
    message("\nComparing CellChat results...")

    # Merge for comparison
    cellchat_merged <- mergeCellChat(cellchat_list, add.names = group_labels)

    return(cellchat_merged)
}

#' Main analysis
main <- function() {
    message("=" |> rep(60) |> paste(collapse = ""))
    message("CellChat Cell-Cell Communication Analysis")
    message("=" |> rep(60) |> paste(collapse = ""))

    # Process samples
    cellchat_R <- list()
    cellchat_NR <- list()

    for (i in 1:nrow(METADATA)) {
        sample_id <- METADATA$sample_id[i]
        response <- METADATA$response[i]

        # Find h5ad file
        h5ad_path <- file.path(DATA_DIR, sprintf("%s_annotated.h5ad", sample_id))

        if (!file.exists(h5ad_path)) {
            message(sprintf("File not found: %s, skipping...", basename(h5ad_path)))
            next
        }

        # Convert to Seurat
        seurat_obj <- h5ad_to_seurat(h5ad_path)

        if (is.null(seurat_obj)) {
            next
        }

        # Run CellChat
        cellchat <- tryCatch({
            run_cellchat_single(seurat_obj, sample_id)
        }, error = function(e) {
            message(sprintf("  CellChat failed: %s", e$message))
            NULL
        })

        if (!is.null(cellchat)) {
            if (response == "R") {
                cellchat_R[[sample_id]] <- cellchat
            } else {
                cellchat_NR[[sample_id]] <- cellchat
            }
        }
    }

    # Save individual results
    message("\nSaving individual CellChat results...")

    all_cellchat <- c(cellchat_R, cellchat_NR)

    if (length(all_cellchat) > 0) {
        saveRDS(all_cellchat, file.path(OUTPUT_DIR, "cellchat_all_samples.rds"))

        # Extract summary statistics
        results_summary <- data.frame()

        for (name in names(all_cellchat)) {
            cc <- all_cellchat[[name]]

            # Number of interactions
            n_interactions <- nrow(cc@LR$LRsig)

            # Top pathways
            if (!is.null(cc@netP$pathways)) {
                top_pathways <- head(names(sort(cc@netP$pathways, decreasing = TRUE)), 5)
            } else {
                top_pathways <- NA
            }

            row <- data.frame(
                sample = name,
                response = ifelse(name %in% names(cellchat_R), "R", "NR"),
                n_interactions = n_interactions,
                top_pathway_1 = ifelse(length(top_pathways) > 0, top_pathways[1], NA),
                top_pathway_2 = ifelse(length(top_pathways) > 1, top_pathways[2], NA),
                top_pathway_3 = ifelse(length(top_pathways) > 2, top_pathways[3], NA)
            )
            results_summary <- rbind(results_summary, row)
        }

        write.csv(results_summary, file.path(OUTPUT_DIR, "cellchat_summary.csv"),
                  row.names = FALSE)

        message("\nResults saved to:")
        message(sprintf("  %s", file.path(OUTPUT_DIR, "cellchat_all_samples.rds")))
        message(sprintf("  %s", file.path(OUTPUT_DIR, "cellchat_summary.csv")))

        # Print summary
        message("\n" |> paste(rep("=", 60), collapse = ""))
        message("SUMMARY")
        message(rep("=", 60) |> paste(collapse = ""))
        print(results_summary)
    } else {
        message("\nNo CellChat results generated.")
    }

    return(all_cellchat)
}

# Run
if (!interactive()) {
    main()
}
