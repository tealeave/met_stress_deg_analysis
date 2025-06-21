# scripts/upset.R

# Purpose: To create UpSet plots for differentially expressed genes (DEGs)
# grouped by specific research questions from the analysis plan.
# Q1 & Q2: Hcy treatment effect across all cell lines.
# Q3: Unique Hcy induction in revertant lines vs. parental lines.
# Q4: Baseline expression differences between revertant and parental lines under Met.

# --- 1. Load Libraries ---
# If you don't have these packages, install them:
# install.packages(c("UpSetR", "here", "biomaRt"))
library(UpSetR)
library(here) # For robust file path management
library(biomaRt) # For converting gene IDs to gene names
library(ggplot2)
library(grid)  # for unit()

# --- 2. Configuration ---
# Define significance thresholds. Genes must meet both criteria to be included.
alpha <- 0.05 # Maximum adjusted p-value (padj)
lfc_threshold <- 1 # Minimum log2 fold change (e.g., 1 represents a 2-fold change)

# Define input and output directories
input_dir <- here("results", "tables_csv")
output_dir <- here("results", "upset_plots_by_question")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. Define Comparison Groups based on Research Questions ---

# Group for Q1 & Q2: Examines the effect of Hcy vs. Met treatment.
q1_q2_comparisons <- c(
  "468_hcy_vs_met",
  "293t_hcy_vs_met",
  "r8_hcy_vs_met",
  "r1_hcy_vs_met"
)

# Define the custom colors for the Q1 & Q2 plots.
# 468: purple, 293T: green, R8: blue, R1: orange
q1_q2_colors <- c(
    "#8616f0", # 468
    "#018520", # 293T
    "#05bbf7", # R8
    "#de6f00"  # R1
)

# Group for Q3: Identifies unique gene induction in revertants under Hcy.
q3_comparisons <- c(
  "r8_vs_468_hcy",
  "r1_vs_293t_hcy"
)

# Group for Q4: Compares baseline expression in revertants vs. parents under Met.
q4_comparisons <- c(
  "r8_vs_468_met",
  "r1_vs_293t_met"
)

# --- 4. Load and Process Gene Lists ---
result_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(result_files) == 0) {
  stop("No result CSV files found in '", input_dir, "'. Please run the deseq_analysis.R script first.")
}

# Initialize empty lists to hold gene sets.
q1_q2_up_sets <- list()
q1_q2_down_sets <- list()
q3_up_sets <- list()
q3_down_sets <- list()
q4_up_sets <- list()
q4_down_sets <- list()

print("Processing and categorizing result files...")

for (file in result_files) {
  comparison_name <- gsub("\\.csv$", "", basename(file))
  results_df <- read.csv(file, row.names = 1, stringsAsFactors = FALSE)
  
  up_genes <- rownames(subset(results_df, padj < alpha & log2FoldChange > lfc_threshold))
  down_genes <- rownames(subset(results_df, padj < alpha & log2FoldChange < -lfc_threshold))
  
  if (comparison_name %in% q1_q2_comparisons) {
    if (length(up_genes) > 0) q1_q2_up_sets[[comparison_name]] <- up_genes
    if (length(down_genes) > 0) q1_q2_down_sets[[comparison_name]] <- down_genes
  } else if (comparison_name %in% q3_comparisons) {
    if (length(up_genes) > 0) q3_up_sets[[comparison_name]] <- up_genes
    if (length(down_genes) > 0) q3_down_sets[[comparison_name]] <- down_genes
  } else if (comparison_name %in% q4_comparisons) {
    if (length(up_genes) > 0) q4_up_sets[[comparison_name]] <- up_genes
    if (length(down_genes) > 0) q4_down_sets[[comparison_name]] <- down_genes
  }
}

# Ensure the lists are in the correct order for plotting
q1_q2_up_sets <- q1_q2_up_sets[q1_q2_comparisons]
q1_q2_down_sets <- q1_q2_down_sets[q1_q2_comparisons]

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
r1_r8_up_intersection <- intersect(q1_q2_up_sets[['r1_hcy_vs_met']], q1_q2_up_sets[['r8_hcy_vs_met']])
other_up_genes <- union(q1_q2_up_sets[['468_hcy_vs_met']], q1_q2_up_sets[['293t_hcy_vs_met']])
unique_up_r1_r8 <- setdiff(r1_r8_up_intersection, other_up_genes)
get_and_print_genenames(unique_up_r1_r8, "Genes uniquely UPREGULATED in both R1 and R8:")

# --- For Downregulated Genes ---
r1_r8_down_intersection <- intersect(q1_q2_down_sets[['r1_hcy_vs_met']], q1_q2_down_sets[['r8_hcy_vs_met']])
other_down_genes <- union(q1_q2_down_sets[['468_hcy_vs_met']], q1_q2_down_sets[['293t_hcy_vs_met']])
unique_down_r1_r8 <- setdiff(r1_r8_down_intersection, other_down_genes)
get_and_print_genenames(unique_down_r1_r8, "Genes uniquely DOWNREGULATED in both R1 and R8:")


# --- 5. Plotting Functions ---

# UPDATED: Increased all visual elements and added padding by increasing plot dimensions.
generate_upset_plot <- function(gene_sets, title, file_path) {
  if (length(gene_sets) > 1) {
    # Increased plot dimensions to better accommodate larger fonts and add padding.
    png(file_path, width = 26, height = 16, units = "in", res = 300)
    print(upset(
      fromList(gene_sets), 
      nsets = length(gene_sets), 
      nintersects = 40, 
      order.by = "freq", 
      mainbar.y.label = "Intersection Size", 
      sets.x.label = "Total Genes in Set", 
      # Increased dot and line sizes in the matrix.
      point.size = 7,
      line.size = 3,
      # Increased and specified text scales for better readability (~2.2x).
      # Set "sets.x.label" (set_size_title) to be smaller than set names.
      text.scale = c(intersection_size_title = 4.4, intersection_size_tick_labels = 3.9, set_size_title = 3.9, set_size_tick_labels = 3.9, set_names = 4.9, main_bar_annotation = 4.4),
      main.bar.color = "skyblue", 
      sets.bar.color = "lightblue", 
      matrix.color = "gray23", 
      shade.color = "whitesmoke"
    ))
    dev.off()
  } else {
    print(paste("Skipping plot for '", title, "' - fewer than 2 gene sets available.", sep=""))
  }
}

# UPDATED: Increased all visual elements and added padding by increasing plot dimensions.
# FURTHER UPDATED: Increased the left margin to provide more space for the y-axis labels.
generate_q1_q2_upset_plot <- function(gene_sets, title, file_path, set_colors, set_order) {
  if (length(gene_sets) > 1) {
    # Increased plot dimensions to better accommodate larger fonts and add padding.
    png(file_path, width = 28, height = 17, units = "in", res = 300)

    print(upset(
      fromList(gene_sets),
      sets = set_order,
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
      # Set "sets.x.label" (set_size_title) to be smaller than set names.
      text.scale = c(intersection_size_title = 4.9, intersection_size_tick_labels = 4.9, set_size_title = 4.9, set_size_tick_labels = 4.4, set_names = 5.6, main_bar_annotation = 4.9),
      main.bar.color = "black",
      sets.bar.color = set_colors,
      matrix.color = "gray23",
      shade.color = "whitesmoke"
    ))

    dev.off()
    
  } else {
    print(paste("Skipping plot for '", title, "' - fewer than 2 gene sets available.", sep=""))
  }
}

# --- 6. Generate All Six Plots ---
print("--- Generating Plots ---")
generate_q1_q2_upset_plot(gene_sets = q1_q2_up_sets, title = "Q1 & Q2: Upregulated Genes in Hcy vs. Met", file_path = file.path(output_dir, "upset_q1_q2_hcy_vs_met_UP.png"), set_colors = q1_q2_colors, set_order = q1_q2_comparisons)
generate_q1_q2_upset_plot(gene_sets = q1_q2_down_sets, title = "Q1 & Q2: Downregulated Genes in Hcy vs. Met", file_path = file.path(output_dir, "upset_q1_q2_hcy_vs_met_DOWN.png"), set_colors = q1_q2_colors, set_order = q1_q2_comparisons)
generate_upset_plot(gene_sets = q3_up_sets, title = "Q3: Upregulated Genes in Revertants vs. Parents (Hcy)", file_path = file.path(output_dir, "upset_q3_revertant_vs_parent_hcy_UP.png"))
generate_upset_plot(gene_sets = q3_down_sets, title = "Q3: Downregulated Genes in Revertants vs. Parents (Hcy)", file_path = file.path(output_dir, "upset_q3_revertant_vs_parent_hcy_DOWN.png"))
generate_upset_plot(gene_sets = q4_up_sets, title = "Q4: Upregulated Genes in Revertants vs. Parents (Met)", file_path = file.path(output_dir, "upset_q4_revertant_vs_parent_met_UP.png"))
generate_upset_plot(gene_sets = q4_down_sets, title = "Q4: Downregulated Genes in Revertants vs. Parents (Met)", file_path = file.path(output_dir, "upset_q4_revertant_vs_parent_met_DOWN.png"))
print("Script finished.")
