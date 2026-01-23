#!/usr/bin/env Rscript
# Quick extraction of RCTD results from saved RDS

library(Matrix)

cat("=== Extracting RCTD Results ===\n")

# Load RCTD results
rctd <- readRDS("output/rctd_results.rds")

# Get results_df directly
results_df <- rctd@results$results_df

cat(sprintf("Total cells in results: %d\n", nrow(results_df)))
cat(sprintf("Columns: %s\n", paste(names(results_df), collapse=", ")))

# Create output dataframe with cell_id
out_df <- data.frame(
    cell_id = rownames(results_df),
    cell_type_rctd = as.character(results_df$first_type),
    second_type = as.character(results_df$second_type),
    spot_class = as.character(results_df$spot_class),
    stringsAsFactors = FALSE
)

# Cell type distribution
cat("\n--- Cell Type Distribution ---\n")
print(table(out_df$cell_type_rctd, useNA = "ifany"))

# Spot class distribution
cat("\n--- Spot Class Distribution ---\n")
print(table(out_df$spot_class, useNA = "ifany"))

# Save annotations
write.csv(out_df, "output/g4x_rctd_annotations.csv", row.names = FALSE)
cat("\nAnnotations saved to: output/g4x_rctd_annotations.csv\n")

# Extract and save weights
if (!is.null(rctd@results$weights)) {
    weights <- as.data.frame(as.matrix(rctd@results$weights))
    weights$cell_id <- rownames(weights)
    write.csv(weights, "output/g4x_rctd_weights.csv", row.names = FALSE)
    cat(sprintf("Weights saved: %d cells x %d cell types\n",
                nrow(weights), ncol(weights) - 1))
}

# Summary stats by sample
cat("\n--- Per-Sample Summary ---\n")
# Extract sample from cell_id (e.g., E02-19 -> E02)
out_df$sample_id <- gsub("-.*", "", out_df$cell_id)

for (sample in unique(out_df$sample_id)) {
    sample_data <- out_df[out_df$sample_id == sample, ]
    cat(sprintf("\n%s (%d cells):\n", sample, nrow(sample_data)))
    type_counts <- table(sample_data$cell_type_rctd)
    top_types <- sort(type_counts, decreasing = TRUE)[1:5]
    for (i in seq_along(top_types)) {
        pct <- 100 * top_types[i] / nrow(sample_data)
        cat(sprintf("  %s: %d (%.1f%%)\n", names(top_types)[i], top_types[i], pct))
    }
}

cat("\n=== Extraction Complete ===\n")
