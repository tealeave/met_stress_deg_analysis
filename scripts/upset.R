# scripts/upset.R
# Purpose: To create UpSet plots for differentially expressed genes (DEGs)
# grouped by specific research questions from the analysis plan.
# Q1 & Q2: Hcy treatment effect across all cell lines.
# Q3: Unique Hcy induction in revertant lines vs. parental lines.
# Q4: Baseline expression differences between revertant and parental lines under Met.

# --- 1. Load Libraries ---
# If you don't have these packages, install them:
# install.packages(c("UpSetR", "here"))
library(UpSetR)
library(here) # For robust file path management

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

# Initialize empty lists to hold the gene sets for each plot.
q1_q2_up_sets <- list()
q1_q2_down_sets <- list()
q3_up_sets <- list()
q3_down_sets <- list()
q4_up_sets <- list()
q4_down_sets <- list()

print("Processing and categorizing result files into groups based on research questions...")

# Loop through each DESeq2 result file.
for (file in result_files) {
  comparison_name <- gsub("\\.csv$", "", basename(file))
  results_df <- read.csv(file, row.names = 1, stringsAsFactors = FALSE)
  
  # Filter for up- and down-regulated genes based on thresholds
  up_genes <- rownames(subset(results_df, padj < alpha & log2FoldChange > lfc_threshold))
  down_genes <- rownames(subset(results_df, padj < alpha & log2FoldChange < -lfc_threshold))
  
  # Check which group the comparison belongs to and add genes to the correct list.
  if (comparison_name %in% q1_q2_comparisons) {
    if (length(up_genes) > 0) q1_q2_up_sets[[comparison_name]] <- up_genes
    if (length(down_genes) > 0) q1_q2_down_sets[[comparison_name]] <- down_genes
    print(paste(" -> Q1/Q2 Group:", comparison_name, "processed."))
  } else if (comparison_name %in% q3_comparisons) {
    if (length(up_genes) > 0) q3_up_sets[[comparison_name]] <- up_genes
    if (length(down_genes) > 0) q3_down_sets[[comparison_name]] <- down_genes
    print(paste(" -> Q3 Group:", comparison_name, "processed."))
  } else if (comparison_name %in% q4_comparisons) {
    if (length(up_genes) > 0) q4_up_sets[[comparison_name]] <- up_genes
    if (length(down_genes) > 0) q4_down_sets[[comparison_name]] <- down_genes
    print(paste(" -> Q4 Group:", comparison_name, "processed."))
  }
}

# --- 5. Generic Plotting Function ---
# This function centralizes the plotting logic to avoid repetition.
generate_upset_plot <- function(gene_sets, title, file_path) {
  if (length(gene_sets) > 1) {
    print(paste("Generating plot:", title))
    png(file_path, width = 12, height = 8, units = "in", res = 300)
    # Using print() is essential to render the plot when running the script from the command line.
    print(upset(
      fromList(gene_sets),
      nsets = length(gene_sets),
      nintersects = 40,
      order.by = "freq",
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Total Genes in Set",
      text.scale = 1.5,
      main.bar.color = "skyblue",
      sets.bar.color = "lightblue",
      matrix.color = "gray23",
      shade.color = "whitesmoke"
    ))
    dev.off()
    print(paste("Plot saved to:", file_path))
  } else {
    print(paste("Skipping plot for '", title, "' - fewer than 2 gene sets available.", sep=""))
  }
}

# --- 6. Generate All Six Plots ---

# Plots for Q1 & Q2: Hcy vs Met Treatment Effects
generate_upset_plot(
  gene_sets = q1_q2_up_sets,
  title = "Q1 & Q2: Upregulated Genes in Hcy vs. Met",
  file_path = file.path(output_dir, "upset_q1_q2_hcy_vs_met_UP.png")
)
generate_upset_plot(
  gene_sets = q1_q2_down_sets,
  title = "Q1 & Q2: Downregulated Genes in Hcy vs. Met",
  file_path = file.path(output_dir, "upset_q1_q2_hcy_vs_met_DOWN.png")
)

# Plots for Q3: Unique Revertant Response to Hcy
generate_upset_plot(
  gene_sets = q3_up_sets,
  title = "Q3: Upregulated Genes in Revertants vs. Parents (Hcy)",
  file_path = file.path(output_dir, "upset_q3_revertant_vs_parent_hcy_UP.png")
)
generate_upset_plot(
  gene_sets = q3_down_sets,
  title = "Q3: Downregulated Genes in Revertants vs. Parents (Hcy)",
  file_path = file.path(output_dir, "upset_q3_revertant_vs_parent_hcy_DOWN.png")
)

# Plots for Q4: Baseline Differences under Met
generate_upset_plot(
  gene_sets = q4_up_sets,
  title = "Q4: Upregulated Genes in Revertants vs. Parents (Met)",
  file_path = file.path(output_dir, "upset_q4_revertant_vs_parent_met_UP.png")
)
generate_upset_plot(
  gene_sets = q4_down_sets,
  title = "Q4: Downregulated Genes in Revertants vs. Parents (Met)",
  file_path = file.path(output_dir, "upset_q4_revertant_vs_parent_met_DOWN.png")
)

print("Script finished. All plots have been generated.")

