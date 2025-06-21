# scripts/upset_q1_q2.R
# Purpose: To create the more complex UpSet plots for Q1 & Q2,
# which require custom colors and labels. This is a standalone script.

# --- 1. Load Libraries ---
# All necessary packages are loaded here.
# install.packages(c("here", "ComplexUpset", "ggplot2", "UpSetR"))
library(here)
library(ComplexUpset)
library(ggplot2)
library(UpSetR) # fromList is used from this package


# --- 2. Configuration ---
# Define significance thresholds. Genes must meet both criteria to be included.
alpha <- 0.05 # Maximum adjusted p-value (padj)
lfc_threshold <- 1 # Minimum log2 fold change (e.g., 1 represents a 2-fold change)

# Define input and output directories
input_dir <- here("results", "tables_csv")
output_dir <- here("results", "upset_plots_by_question")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# --- 3. Define Comparison Groups for Q1 & Q2 ---
# Group for Q1 & Q2: Examines the effect of Hcy vs. Met treatment.
q1_q2_comparisons <- c(
  "468_hcy_vs_met",
  "293t_hcy_vs_met",
  "r8_hcy_vs_met",
  "r1_hcy_vs_met"
)

# Define the custom colors for the Q1 & Q2 plots.
# The order MUST match the order in q1_q2_comparisons.
q1_q2_colors <- c(
    "#8616f0", # 468
    "#018520", # 293T
    "#05bbf7", # R8
    "#de6f00"  # R1
)


# --- 4. Load and Process Gene Lists for Q1 & Q2 ---
result_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(result_files) == 0) {
  stop("No result CSV files found in '", input_dir, "'. Please run the deseq_analysis.R script first.")
}

# Initialize empty lists to hold the gene sets for the plot.
q1_q2_up_sets <- list()
q1_q2_down_sets <- list()

print("Processing and categorizing result files for Q1/Q2...")

# Loop through each DESeq2 result file.
for (file in result_files) {
  comparison_name <- gsub("\\.csv$", "", basename(file))
  
  # Only process files that are relevant to this script's purpose
  if (comparison_name %in% q1_q2_comparisons) {
    results_df <- read.csv(file, row.names = 1, stringsAsFactors = FALSE)
    
    # Filter for up- and down-regulated genes based on thresholds
    up_genes <- rownames(subset(results_df, padj < alpha & log2FoldChange > lfc_threshold))
    down_genes <- rownames(subset(results_df, padj < alpha & log2FoldChange < -lfc_threshold))
    
    # Add genes to the correct list.
    if (length(up_genes) > 0) q1_q2_up_sets[[comparison_name]] <- up_genes
    if (length(down_genes) > 0) q1_q2_down_sets[[comparison_name]] <- down_genes
    print(paste(" -> Processed:", comparison_name))
  }
}

# Reorder the lists to match the desired color and plot order.
if (length(q1_q2_up_sets) > 0) {
    q1_q2_up_sets <- q1_q2_up_sets[q1_q2_comparisons]
}
if (length(q1_q2_down_sets) > 0) {
    q1_q2_down_sets <- q1_q2_down_sets[q1_q2_comparisons]
}


# --- 5. Plotting Function for Q1 & Q2 ---
# Specialized plotting function for Q1 & Q2 using ComplexUpset for coloring.
generate_q1_q2_complex_upset <- function(gene_sets, comparisons, colors, title, file_path) {
    if (length(gene_sets) < 2) {
        print(paste("Skipping plot for '", title, "' - fewer than 2 gene sets available.", sep=""))
        return()
    }

    print(paste("Generating plot with ComplexUpset:", title))

    # Convert list to the binary matrix format that ComplexUpset needs
    upset_data <- fromList(gene_sets)

    # Define the order of sets and colors for the plot. We reverse them to get the desired layout.
    sets_in_plot_order <- rev(comparisons)
    colors_in_plot_order <- rev(colors)

    # Create a named vector for the bar colors to ensure ggplot matches them correctly
    set_colors_named <- setNames(colors, comparisons)

    # Open a PNG device to save the plot
    png(file_path, width = 12, height = 8, units = "in", res = 300)

    # Generate the plot using ComplexUpset::upset
    print(ComplexUpset::upset(
      upset_data,
      sets_in_plot_order,     # Use the specified set order
      name = "set",  
      sort_intersections = "descending",
      width_ratio = 0.2,      # Adjust space for set size bars
      base_annotations = list(
        'Intersection size' = intersection_size(
            fill = 'black',
            text = list(size = 5)
        )
      ),
      set_sizes = (
        upset_set_size(
          geom = geom_bar(aes(y = after_stat(count), fill = group), width = 0.7),
          position = 'right'
        )
          + scale_fill_manual(values = set_colors_named, guide = 'none')
          + labs(y = 'Total Genes in Set')
          + theme(
              axis.text.x = element_text(size = 12),
              axis.title.x = element_text(size = 14, face = "bold")
            )
      ),
      # This part colors the set name labels
      themes = list(
        default = theme(
          text = element_text(size = 14),
          axis.text.y = element_text(
              color = colors_in_plot_order,
              face = "bold",
              size = 15
          ),
          axis.title.y = element_blank()
        )
      )
    ))

    # Close the PNG device
    dev.off()
    print(paste("Plot with colored labels saved to:", file_path))
}


# --- 6. Generate Q1 & Q2 Plots ---
print("--- Generating Q1 & Q2 Plots ---")

# Plots for Q1 & Q2: Using the ComplexUpset function
generate_q1_q2_complex_upset(
  gene_sets = q1_q2_up_sets,
  comparisons = q1_q2_comparisons,
  colors = q1_q2_colors,
  title = "Q1 & Q2: Upregulated Genes in Hcy vs. Met",
  file_path = file.path(output_dir, "upset_q1_q2_hcy_vs_met_UP_colored.png")
)
generate_q1_q2_complex_upset(
  gene_sets = q1_q2_down_sets,
  comparisons = q1_q2_comparisons,
  colors = q1_q2_colors,
  title = "Q1 & Q2: Downregulated Genes in Hcy vs. Met",
  file_path = file.path(output_dir, "upset_q1_q2_hcy_vs_met_DOWN_colored.png")
)

print("Script finished.")
