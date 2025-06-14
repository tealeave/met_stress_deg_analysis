# scripts/analytical_venn.R
# Purpose: Generate Venn diagrams specifically designed to answer the four key
# analytical questions. This script separates analyses for up- and down-regulated
# genes and saves the resulting gene lists for downstream pathway analysis.

# --- 1. Load Libraries ---
# install.packages(c("ggvenn", "ggplot2", "here", "eulerr"))
library(ggvenn)
library(ggplot2)
library(here)

# --- 2. Configuration ---
# Define significance thresholds
alpha <- 0.05       # Maximum adjusted p-value
lfc_threshold <- 1.0 # Minimum absolute log2 fold change

# Define directories
input_dir <- here("results", "tables_csv")
output_dir <- here("results", "venn_diagrams_analytical")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. Helper Functions ---
# Reads a DESeq2 result file and returns significant genes for a given direction.
get_de_genes <- function(comparison_name, direction) {
  file_path <- file.path(input_dir, paste0(comparison_name, ".csv"))
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path)); return(character(0))
  }
  results_df <- read.csv(file_path, row.names = 1)
  
  if (direction == "up") {
    sig_genes <- subset(results_df, padj < alpha & log2FoldChange > lfc_threshold)
  } else if (direction == "down") {
    sig_genes <- subset(results_df, padj < alpha & log2FoldChange < -lfc_threshold)
  }
  return(rownames(sig_genes))
}

# Saves a list of genes to a text file.
save_gene_list <- function(gene_vector, filename) {
  write.table(gene_vector, file = here(output_dir, filename),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# --- 4. Analysis for Each Question ---

# ==============================================================================
# Question 1: Is there a common response to Hcy in parental lines?
# ==============================================================================
print("--- Answering Question 1: Common parental Hcy response ---")
for (direction in c("up", "down")) {
  # Get gene lists for parental lines
  parental_293t <- get_de_genes("293t_hcy_vs_met", direction)
  parental_468 <- get_de_genes("468_hcy_vs_met", direction)
  
  # Find the intersection (common genes)
  common_genes <- intersect(parental_293t, parental_468)
  save_gene_list(common_genes, paste0("q1_common_parental_response_", direction, ".txt"))
  print(paste("Found", length(common_genes), "common", direction, "-regulated genes in parental lines."))
  
  # Visualize
  p <- ggvenn(
    list(`293T` = parental_293t, `468` = parental_468),
    fill_color = c("#0073C2FF", "#EFC000FF"), set_name_size = 5
  ) + labs(title = paste("Q1: Common Parental Hcy Response (", toupper(direction), "regulated)", sep=""))
  
  ggsave(here(output_dir, paste0("q1_venn_parental_common_response_", direction, ".png")), p, width=7, height=6)
}

# ==============================================================================
# Question 2: Which common response pathways are maintained in revertants?
# ==============================================================================
print("--- Answering Question 2: Maintained Hcy response in revertants ---")
for (direction in c("up", "down")) {
  # Load the common parental response genes we saved in Q1
  common_parental_response <- read.table(here(output_dir, paste0("q1_common_parental_response_", direction, ".txt")))$V1
  
  # Get Hcy response for revertant lines
  r1_hcy_response <- get_de_genes("r1_hcy_vs_met", direction)
  r8_hcy_response <- get_de_genes("r8_hcy_vs_met", direction)
  
  # Define the three sets for the Venn diagram
  gene_sets <- list(
    `Common Parental` = common_parental_response,
    `R1 Revertant` = r1_hcy_response,
    `R8 Revertant` = r8_hcy_response
  )
  
  # Find and save the core maintained genes (intersection of all three)
  maintained_genes <- Reduce(intersect, gene_sets)
  save_gene_list(maintained_genes, paste0("q2_maintained_core_response_", direction, ".txt"))
  print(paste("Found", length(maintained_genes), "core response genes maintained in both revertants (", direction, ")."))
  
  # Visualize
  p <- ggvenn(
    gene_sets,
    fill_color = c("#E69F00", "#56B4E9", "#009E73"), set_name_size = 5
  ) + labs(title = paste("Q2: Maintained Hcy Response (", toupper(direction), "regulated)", sep=""))

  ggsave(here(output_dir, paste0("q2_venn_maintained_response_", direction, ".png")), p, width=8, height=7)
}

# ==============================================================================
# Question 3: Are there different pathways induced in revertants?
# ==============================================================================
print("--- Answering Question 3: Unique Hcy response in revertants ---")
for (direction in c("up", "down")) {
  # --- Pair 1: R1 vs 293T ---
  r1_genes <- get_de_genes("r1_hcy_vs_met", direction)
  parent_293t_genes <- get_de_genes("293t_hcy_vs_met", direction)
  
  # Find and save genes unique to R1
  r1_unique_genes <- setdiff(r1_genes, parent_293t_genes)
  save_gene_list(r1_unique_genes, paste0("q3_unique_to_r1_response_", direction, ".txt"))
  print(paste("Found", length(r1_unique_genes), "genes unique to R1 Hcy response (", direction, ")."))
  
  # Visualize
  p1 <- ggvenn(
    list(`R1 Response` = r1_genes, `293T Response` = parent_293t_genes),
    fill_color = c("#D55E00", "#0072B2"), set_name_size = 5
  ) + labs(title = paste("Q3: R1 vs 293T Hcy Response (", toupper(direction), "regulated)", sep=""))
  ggsave(here(output_dir, paste0("q3_venn_r1_vs_parent_", direction, ".png")), p1, width=7, height=6)
  
  # --- Pair 2: R8 vs 468 ---
  r8_genes <- get_de_genes("r8_hcy_vs_met", direction)
  parent_468_genes <- get_de_genes("468_hcy_vs_met", direction)

  # Find and save genes unique to R8
  r8_unique_genes <- setdiff(r8_genes, parent_468_genes)
  save_gene_list(r8_unique_genes, paste0("q3_unique_to_r8_response_", direction, ".txt"))
  print(paste("Found", length(r8_unique_genes), "genes unique to R8 Hcy response (", direction, ")."))
  
  # Visualize
  p2 <- ggvenn(
    list(`R8 Response` = r8_genes, `468 Response` = parent_468_genes),
    fill_color = c("#CC79A7", "#F0E442"), set_name_size = 5
  ) + labs(title = paste("Q3: R8 vs 468 Hcy Response (", toupper(direction), "regulated)", sep=""))
  ggsave(here(output_dir, paste0("q3_venn_r8_vs_parent_", direction, ".png")), p2, width=7, height=6)
}

# ==============================================================================
# Question 4: Do revertants have different baseline profiles (in Met)?
# ==============================================================================
print("--- Answering Question 4: Baseline differences in revertants ---")
for (direction in c("up", "down")) {
  # Get baseline DEGs for each revertant vs its parent
  r1_baseline <- get_de_genes("r1_vs_293t_met", direction)
  r8_baseline <- get_de_genes("r8_vs_468_met", direction)
  
  print(paste("Found", length(r1_baseline), "baseline DEGs in R1 vs 293T (", direction, ")."))
  print(paste("Found", length(r8_baseline), "baseline DEGs in R8 vs 468 (", direction, ")."))
  
  # Find and save the common baseline signature
  common_signature <- intersect(r1_baseline, r8_baseline)
  save_gene_list(common_signature, paste0("q4_common_revertant_signature_", direction, ".txt"))
  print(paste("Found", length(common_signature), "common baseline DEGs in revertants (", direction, ")."))

  # Visualize
  p <- ggvenn(
    list(`R1 vs 293T` = r1_baseline, `R8 vs 468` = r8_baseline),
    fill_color = c("#868686FF", "#CD534CFF"), set_name_size = 5
  ) + labs(title = paste("Q4: Common Revertant Baseline Signature (", toupper(direction), "regulated)", sep=""))
  
  ggsave(here(output_dir, paste0("q4_venn_revertant_baseline_", direction, ".png")), p, width=7, height=6)
}

print("--- Analytical Venn diagram script finished. ---")
