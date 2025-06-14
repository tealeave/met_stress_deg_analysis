# scripts/upset.R
# Purpose: To create four UpSet plots by splitting the 8 DESeq2 comparisons
# into two meaningful groups: "Treatment Effect" and "Cell Line Effect".

# --- 1. Load Libraries ---
# If you don't have these packages, install them:
# install.packages(c("UpSetR", "here"))
library(UpSetR)
library(here) # For robust file path management

# --- 2. Configuration ---
# Define significance thresholds. Genes must meet both criteria to be included.
alpha <- 0.05 # Maximum adjusted p-value (padj)
lfc_threshold <- 1 # Minimum log2 fold change (e.g., 1 means a 2-fold change)

# Define input and output directories
input_dir <- here("results", "tables_csv")
output_dir <- here("results", "upset_plots_grouped")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. Define Comparison Groups ---
# Group 1: Focuses on the treatment effect (Hcy vs Met) within each cell line.
group1_comparisons <- c(
  "468_hcy_vs_met",
  "r8_hcy_vs_met",
  "293t_hcy_vs_met",
  "r1_hcy_vs_met"
)

# Group 2: Focuses on the cell line differences under each treatment condition.
group2_comparisons <- c(
  "r8_vs_468_met",
  "r8_vs_468_hcy",
  "r1_vs_293t_met",
  "r1_vs_293t_hcy"
)

# --- 4. Load and Process Gene Lists ---
result_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(result_files) == 0) {
  stop("No result CSV files found in '", input_dir, "'. Please run the deseq_analysis.R script first.")
}

# Initialize four empty lists to hold the gene sets for each of the four plots.
g1_up_sets <- list()
g1_down_sets <- list()
g2_up_sets <- list()
g2_down_sets <- list()

print("Processing and categorizing result files into groups...")

# Loop through each DESeq2 result file.
for (file in result_files) {
  comparison_name <- gsub("\\.csv$", "", basename(file))
  results_df <- read.csv(file, row.names = 1, stringsAsFactors = FALSE)
  
  # Filter for up- and down-regulated genes based on thresholds
  up_genes <- rownames(subset(results_df, padj < alpha & log2FoldChange > lfc_threshold))
  down_genes <- rownames(subset(results_df, padj < alpha & log2FoldChange < -lfc_threshold))
  
  # Check which group the comparison belongs to and add genes to the correct list.
  if (comparison_name %in% group1_comparisons) {
    if (length(up_genes) > 0) g1_up_sets[[comparison_name]] <- up_genes
    if (length(down_genes) > 0) g1_down_sets[[comparison_name]] <- down_genes
    print(paste(" -> Group 1:", comparison_name, "processed."))
  } else if (comparison_name %in% group2_comparisons) {
    if (length(up_genes) > 0) g2_up_sets[[comparison_name]] <- up_genes
    if (length(down_genes) > 0) g2_down_sets[[comparison_name]] <- down_genes
    print(paste(" -> Group 2:", comparison_name, "processed."))
  }
}

# --- 5. Generic Plotting Function ---
# This function avoids repeating the same upset() call four times.
generate_upset_plot <- function(gene_sets, title, y_label, x_label, file_path) {
  if (length(gene_sets) > 1) {
    print(paste("Generating plot:", title))
    png(file_path, width = 12, height = 8, units = "in", res = 300)
    print(upset(
      fromList(gene_sets),
      nsets = length(gene_sets),
      nintersects = 40,
      order.by = "freq",
      mainbar.y.label = y_label,
      sets.x.label = x_label,
      text.scale = 1.5
    ))
    dev.off()
    print(paste("Plot saved to:", file_path))
  } else {
    print(paste("Skipping plot for", title, "- not enough gene sets."))
  }
}

# --- 6. Generate All Four Plots ---

# Plot 1: Group 1 (Treatment) - Upregulated
generate_upset_plot(
  gene_sets = g1_up_sets,
  title = "Group 1 (Treatment): Upregulated Genes",
  y_label = "Intersection Size (Upregulated)",
  x_label = "Total Genes in Set",
  file_path = file.path(output_dir, "upset_plot_group1_treatment_up.png")
)

# Plot 2: Group 1 (Treatment) - Downregulated
generate_upset_plot(
  gene_sets = g1_down_sets,
  title = "Group 1 (Treatment): Downregulated Genes",
  y_label = "Intersection Size (Downregulated)",
  x_label = "Total Genes in Set",
  file_path = file.path(output_dir, "upset_plot_group1_treatment_down.png")
)

# Plot 3: Group 2 (Cell Line) - Upregulated
generate_upset_plot(
  gene_sets = g2_up_sets,
  title = "Group 2 (Cell Line): Upregulated Genes",
  y_label = "Intersection Size (Upregulated)",
  x_label = "Total Genes in Set",
  file_path = file.path(output_dir, "upset_plot_group2_cellline_up.png")
)

# Plot 4: Group 2 (Cell Line) - Downregulated
generate_upset_plot(
  gene_sets = g2_down_sets,
  title = "Group 2 (Cell Line): Downregulated Genes",
  y_label = "Intersection Size (Downregulated)",
  x_label = "Total Genes in Set",
  file_path = file.path(output_dir, "upset_plot_group2_cellline_down.png")
)

print("Script finished.")
