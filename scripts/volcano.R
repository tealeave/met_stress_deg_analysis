# scripts/volcano.R
# Purpose: Generate volcano plots for each DESeq2 comparison.

# --- 1. Load Libraries ---
library(ggplot2)

# --- 2. Define Plotting Function ---
# A function to create and save a volcano plot
create_volcano_plot <- function(res, title, output_dir) {
  
  # Convert results to a data frame and add gene names as a column
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  
  # Add a column for significance status
  df$status <- "Not Significant"
  df$status[df$padj < 0.05 & !is.na(df$padj) & df$log2FoldChange > 1] <- "Up-regulated"
  df$status[df$padj < 0.05 & !is.na(df$padj) & df$log2FoldChange < -1] <- "Down-regulated"
  df$status <- factor(df$status, levels = c("Up-regulated", "Down-regulated", "Not Significant"))
  
  # Create the plot
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = status), alpha = 0.5, size = 1.5) +
    scale_color_manual(values = c("Up-regulated" = "#E41A1C", "Down-regulated" = "#377EB8", "Not Significant" = "grey")) +
    theme_bw(base_size = 14) +
    theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(title = title,
         x = bquote(~Log[2]~ "Fold Change"),
         y = bquote(~-Log[10]~italic(P)~"-value")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")

  # Define file path and save the plot
  file_path <- file.path(output_dir, paste0("volcano_", gsub(" ", "_", title), ".png"))
  ggsave(
    filename = file_path,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  print(paste("Saved volcano plot to:", file_path))
}

# --- 3. Main Script: Generate Plots ---
# Define input and output directories
rds_dir <- "results/deseq_objects"
output_dir <- "results/volcano_plots"

# Ensure output directory exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Get list of RDS files
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

# Check if there are any files to process
if (length(rds_files) == 0) {
  warning("No RDS files found in ", rds_dir)
} else {
  # Loop through each file, read the DESeq2 result, and create a plot
  for (file in rds_files) {
    tryCatch({
      res_object <- readRDS(file)
      comp_name <- gsub("\\.rds$", "", basename(file))
      create_volcano_plot(res_object, comp_name, output_dir)
    }, error = function(e) {
      warning(paste("Could not process file:", file, "\nError:", e$message))
    })
  }
}

print("All volcano plots have been generated.")
