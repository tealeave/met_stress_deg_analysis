# results/fiona_results/scripts/upset_fiona.R

# Purpose: To create UpSet plots for differentially expressed genes (DEGs)
# from Fiona's results, focusing on the Hcy treatment effect across all cell lines (Q1 & Q2).
# This script is adapted from the original upset.R to match the fiona_results data structure.

# --- 1. Load Libraries ---
# If you don't have these packages, install them:
# install.packages(c("UpSetR", "here", "biomaRt", "ggplot2", "grid"))
library(UpSetR)
library(here) # For robust file path management
library(biomaRt) # For converting gene IDs to gene names
library(ggplot2)
library(grid)  # for unit()

# --- 2. Configuration ---
# Define significance thresholds. Genes must meet both criteria to be included.
alpha <- 0.05 # Maximum adjusted p-value (padj)
lfc_threshold <- 1 # Minimum log2 fold change (e.g., 1 represents a 2-fold change)

# Define input and output directories for Fiona's results
input_dir <- here("results", "fiona_results", "data")
output_dir <- here("results", "fiona_results", "upset_plots_by_question")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. Define Comparison Groups based on Research Questions ---

# Group for Q1 & Q2: Examines the effect of Hcy vs. Met treatment in Fiona's data.
# Names are based on the subdirectories in results/fiona_results/data/
q1_q2_comparisons <- c(
  "MB468_Hcy_vs_Met",
  "HEK293T_Hcy_vs_Met",
  "R8_Hcy_vs_Met",
  "R1_Hcy_vs_Met"
)

# Define the custom colors for the Q1 & Q2 plots, maintaining the original color scheme.
# MB468 (468): purple, HEK293T (293T): green, R8: blue, R1: orange
q1_q2_colors <- c(
    "#8616f0", # MB468
    "#018520", # HEK293T
    "#05bbf7", # R8
    "#de6f00"  # R1
)

# NOTE: Q3 and Q4 comparisons are commented out as the corresponding data directories
# do not exist in the `results/fiona_results/data` folder structure.
#
# # Group for Q3: Identifies unique gene induction in revertants under Hcy.
# q3_comparisons <- c(
#  "r8_vs_468_hcy",
#  "r1_vs_293t_hcy"
# )
#
# # Group for Q4: Compares baseline expression in revertants vs. parents under Met.
# q4_comparisons <- c(
#  "r8_vs_468_met",
#  "r1_vs_293t_met"
# )

# --- 4. Load and Process Gene Lists ---
# List all files directly within the input directory.
result_files <- list.files(input_dir, full.names = TRUE)

if (length(result_files) == 0) {
  stop("No result files found in '", input_dir, "'. Please check the file structure.")
}

# Initialize empty lists to hold gene sets.
q1_q2_up_sets <- list()
q1_q2_down_sets <- list()
# q3_up_sets <- list()
# q3_down_sets <- list()
# q4_up_sets <- list()
# q4_down_sets <- list()

print("Processing and categorizing result files from fiona_results/data...")

for (file in result_files) {
  # The comparison name is the name of the file itself.
  comparison_name <- basename(file)
  # Corrected file reading to handle space-delimited format.
  results_df <- read.table(file, header = TRUE, sep = "", row.names = 1, stringsAsFactors = FALSE)
  
  up_genes <- rownames(subset(results_df, padj < alpha & log2FoldChange > lfc_threshold))
  down_genes <- rownames(subset(results_df, padj < alpha & log2FoldChange < -lfc_threshold))
  
  if (comparison_name %in% q1_q2_comparisons) {
    if (length(up_genes) > 0) q1_q2_up_sets[[comparison_name]] <- up_genes
    if (length(down_genes) > 0) q1_q2_down_sets[[comparison_name]] <- down_genes
  }
  # The following blocks are commented out as Q3/Q4 data is not available.
  # } else if (comparison_name %in% q3_comparisons) {
  #   if (length(up_genes) > 0) q3_up_sets[[comparison_name]] <- up_genes
  #   if (length(down_genes) > 0) q3_down_sets[[comparison_name]] <- down_genes
  # } else if (comparison_name %in% q4_comparisons) {
  #   if (length(up_genes) > 0) q4_up_sets[[comparison_name]] <- up_genes
  #   if (length(down_genes) > 0) q4_down_sets[[comparison_name]] <- down_genes
  # }
}

# Ensure the lists are in the correct order for plotting
q1_q2_up_sets <- q1_q2_up_sets[q1_q2_comparisons]
q1_q2_down_sets <- q1_q2_down_sets[q1_q2_comparisons]
# Remove any NULL entries that might result if a comparison had no significant genes
q1_q2_up_sets <- q1_q2_up_sets[!sapply(q1_q2_up_sets, is.null)]
q1_q2_down_sets <- q1_q2_down_sets[!sapply(q1_q2_down_sets, is.null)]


# --- 4.5. Find Unique Intersections and Convert IDs to Names ---
print("--- Calculating Unique Intersections for Q1/Q2 ---")

# Connect to Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Generic function to get gene names and print them
get_and_print_genenames <- function(gene_ids, description) {
  print(description)
  if (length(gene_ids) > 0) {
    # The Ensembl IDs from DESeq2 might have version numbers (.1, .2, etc.)
    # We need to remove these before querying biomaRt.
    gene_ids_no_version <- gsub("\\..*$", "", gene_ids)
    
    gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                       filters = 'ensembl_gene_id',
                       values = gene_ids_no_version,
                       mart = ensembl)
    
    # Create a data frame for nice printing
    results_df <- data.frame(
      ensembl_id = gene_ids,
      hgnc_symbol = gene_info$hgnc_symbol[match(gene_ids_no_version, gene_info$ensembl_gene_id)]
    )
    
    print(results_df)
  } else {
    print("None found.")
  }
  cat("\n") # Add a newline for spacing
}

# --- For Upregulated Genes ---
# Note: Updated comparison names to match Fiona's data
r1_r8_up_intersection <- intersect(q1_q2_up_sets[['R1_Hcy_vs_Met']], q1_q2_up_sets[['R8_Hcy_vs_Met']])
other_up_genes <- union(q1_q2_up_sets[['MB468_Hcy_vs_Met']], q1_q2_up_sets[['HEK293T_Hcy_vs_Met']])
unique_up_r1_r8 <- setdiff(r1_r8_up_intersection, other_up_genes)
get_and_print_genenames(unique_up_r1_r8, "Genes uniquely UPREGULATED in both R1 and R8:")

# --- For Downregulated Genes ---
# Note: Updated comparison names to match Fiona's data
r1_r8_down_intersection <- intersect(q1_q2_down_sets[['R1_Hcy_vs_Met']], q1_q2_down_sets[['R8_Hcy_vs_Met']])
other_down_genes <- union(q1_q2_down_sets[['MB468_Hcy_vs_Met']], q1_q2_down_sets[['HEK293T_Hcy_vs_Met']])
unique_down_r1_r8 <- setdiff(r1_r8_down_intersection, other_down_genes)
get_and_print_genenames(unique_down_r1_r8, "Genes uniquely DOWNREGULATED in both R1 and R8:")


# --- 5. Plotting Functions ---

# Plotting function for Q1/Q2 with custom colors and larger text.
generate_q1_q2_upset_plot <- function(gene_sets, title, file_path, set_colors, set_order) {
  if (length(gene_sets) > 1) {
    # Increased plot dimensions to better accommodate larger fonts and add padding.
    png(file_path, width = 28, height = 17, units = "in", res = 300)

    # Use the intersection of the set order and the actual names of gene sets
    # This prevents errors if a set had no significant genes and was removed.
    final_set_order <- intersect(set_order, names(gene_sets))
    final_set_colors <- set_colors[match(final_set_order, set_order)]

    print(upset(
      fromList(gene_sets),
      sets = final_set_order,
      keep.order = TRUE,
      nsets = length(gene_sets),
      nintersects = 40,
      order.by = "freq",
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Total Genes in Set",
      # Increased dot and line sizes in the matrix.
      point.size = 7,
      line.size = 3,
      # Increased text scales for better readability (~2.2x).
      text.scale = c(intersection_size_title = 4.9, intersection_size_tick_labels = 4.9, set_size_title = 4.9, set_size_tick_labels = 4.4, set_names = 5.6, main_bar_annotation = 4.9),
      main.bar.color = "black",
      sets.bar.color = final_set_colors,
      matrix.color = "gray23",
      shade.color = "whitesmoke"
    ))

    dev.off()
    
  } else {
    print(paste("Skipping plot for '", title, "' - fewer than 2 gene sets available.", sep=""))
  }
}

# --- 6. Generate All Plots ---
print("--- Generating Plots for Fiona's Q1/Q2 Data ---")
generate_q1_q2_upset_plot(gene_sets = q1_q2_up_sets, title = "Q1 & Q2: Upregulated Genes in Hcy vs. Met", file_path = file.path(output_dir, "upset_q1_q2_hcy_vs_met_UP.png"), set_colors = q1_q2_colors, set_order = q1_q2_comparisons)
generate_q1_q2_upset_plot(gene_sets = q1_q2_down_sets, title = "Q1 & Q2: Downregulated Genes in Hcy vs. Met", file_path = file.path(output_dir, "upset_q1_q2_hcy_vs_met_DOWN.png"), set_colors = q1_q2_colors, set_order = q1_q2_comparisons)

# NOTE: Q3 and Q4 plots are not generated because the data is not available.
# print("--- Generating Plots ---")
# generate_upset_plot(gene_sets = q3_up_sets, title = "Q3: Upregulated Genes in Revertants vs. Parents (Hcy)", file_path = file.path(output_dir, "upset_q3_revertant_vs_parent_hcy_UP.png"))
# generate_upset_plot(gene_sets = q3_down_sets, title = "Q3: Downregulated Genes in Revertants vs. Parents (Hcy)", file_path = file.path(output_dir, "upset_q3_revertant_vs_parent_hcy_DOWN.png"))
# generate_upset_plot(gene_sets = q4_up_sets, title = "Q4: Upregulated Genes in Revertants vs. Parents (Met)", file_path = file.path(output_dir, "upset_q4_revertant_vs_parent_met_UP.png"))
# generate_upset_plot(gene_sets = q4_down_sets, title = "Q4: Downregulated Genes in Revertants vs. Parents (Met)", file_path = file.path(output_dir, "upset_q4_revertant_vs_parent_met_DOWN.png"))

print("Script finished.")
