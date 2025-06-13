# scripts/plot_volcano.R
# Purpose: Generate volcano plots for each DESeq2 comparison.

# --- 1. Load Libraries ---
library(ggplot2)

# --- 2. Define Plotting Function ---
# A function to create and save a volcano plot
create_volcano_plot <- function(res, title) {
  
  # Convert results to a data frame and add gene names as a column
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  
  # Add a column for significance status
  df$status <- "Not Significant"
  df$status[df$padj < 0.05 & df$log2FoldChange > 1] <- "Up-regulated"
  df$status[df$padj < 0.05 & df$log2FoldChange < -1] <- "Down-regulated"
  df$status <- factor(df$status, levels = c("Up-regulated", "Down-regulated", "Not Significant"))
  
  # Create the plot
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = status), alpha = 0.5, size = 1.5) +
    scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "grey")) +
    theme_bw(base_size = 14) +
    theme(legend.title = element_blank()) +
    labs(title = title,
         x = "log2(Fold Change)",
         y = "-log10(p-value)") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed")

  # Save the plot
  ggsave(
    filename = paste0("results/volcano_", gsub(" ", "_", title), ".png"),
    plot = p,
    width = 8,
    height = 6
  )
  
  print(paste("Saved volcano plot for:", title))
}

# --- 3. Main Script: Generate Plots ---
dir.create("results", showWarnings = FALSE) # Ensure results directory exists
rds_files <- list.files("results/deseq_objects", pattern = "\\.rds$", full.names = TRUE)

for (file in rds_files) {
  res_object <- readRDS(file)
  comp_name <- gsub(".rds", "", basename(file))
  create_volcano_plot(res_object, comp_name)
}

print("All volcano plots have been generated.")
