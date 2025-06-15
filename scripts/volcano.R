# scripts/volcano.R
# Purpose: Generate static volcano plots for specific DESeq2 comparisons.

# --- 1. Load Libraries ---
# Ensure ggplot2 is installed: install.packages("ggplot2")
library(ggplot2)

# --- 2. Define Plotting Function ---
# A function to create and save a static volcano plot using ggplot2
create_volcano_plot <- function(res, title, output_dir) {
  
  # Convert DESeq results to a data frame for plotting with ggplot2
  df <- as.data.frame(res)
  df$gene <- rownames(df) # Keep gene identifiers
  
  # --- Assign Significance Status ---
  # Add a column to categorize genes for coloring
  df$status <- "Not Significant"
  # Conditions for significance: adjusted p-value < 0.05 and log2 fold change > 1 or < -1
  df$status[df$padj < 0.05 & !is.na(df$padj) & df$log2FoldChange > 1] <- "Up-regulated"
  df$status[df$padj < 0.05 & !is.na(df$padj) & df$log2FoldChange < -1] <- "Down-regulated"
  
  # Convert status to a factor to ensure consistent legend order and colors
  df$status <- factor(df$status, levels = c("Up-regulated", "Down-regulated", "Not Significant"))
  
  # --- Create the Plot ---
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    # Add points, colored by significance status
    geom_point(aes(color = status), alpha = 0.5, size = 1.5) +
    # Define custom colors for each category
    scale_color_manual(values = c("Up-regulated" = "#E41A1C", "Down-regulated" = "#377EB8", "Not Significant" = "grey")) +
    # Use a clean black and white theme
    theme_bw(base_size = 14) +
    # Customize the theme to remove the legend title and gridlines for a cleaner look
    theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    # Add informative labels and title
    labs(title = title,
         x = bquote(~Log[2]~ "Fold Change"),
         y = bquote(~-Log[10]~italic(P)~"-value")) +
    # Add dashed lines for significance thresholds
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")

  # --- Save the Plot ---
  # Sanitize the title to create a valid filename
  safe_title <- gsub("[^a-zA-Z0-9_]", "_", title)
  file_path <- file.path(output_dir, paste0("volcano_", safe_title, ".png"))
  
  ggsave(
    filename = file_path,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  message("Saved static volcano plot to: ", file_path)
}

# --- 3. Main Script: Define Comparisons and Generate Plots ---
# Define input and output directories
rds_dir <- "results/deseq_objects"
output_dir <- "results/volcano_plots"

# Ensure the output directory exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- UPDATED: Define the specific list of comparisons to generate ---
comparisons_to_plot <- c(
  # Comparisons within batch 'exp2'
  "468_hcy_vs_met",
  "r8_hcy_vs_met",
  "r8_vs_468_met",
  "r8_vs_468_hcy",
  # Comparisons within batch 'exp1'
  "293t_hcy_vs_met",
  "r1_hcy_vs_met",
  "r1_vs_293t_met",
  "r1_vs_293t_hcy"
)

message("Starting static volcano plot generation for ", length(comparisons_to_plot), " comparisons.")

# --- UPDATED: Loop through the defined list instead of all files ---
for (comp_name in comparisons_to_plot) {
  # Construct the full path to the expected RDS file
  file_path <- file.path(rds_dir, paste0(comp_name, ".rds"))
  
  # Check if the file exists before trying to process it
  if (file.exists(file_path)) {
    tryCatch({
      # Read the DESeq2 results object
      res_object <- readRDS(file_path)
      
      # Use the comparison name as the plot title
      plot_title <- comp_name 
      
      # Call the function to create and save the plot
      create_volcano_plot(res_object, plot_title, output_dir)
      
    }, error = function(e) {
      # Report any errors that occur during file processing
      warning("Could not process ", file_path, ": ", e$message)
    })
  } else {
    # Warn if a specific RDS file is not found
    warning("RDS file not found for comparison: ", comp_name, ". Skipping.")
  }
}

message("All static volcano plots have been generated.")
