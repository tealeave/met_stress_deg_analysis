# scripts/plot_venn.R
# Purpose: Generate Venn diagrams to visualize overlaps between key gene sets.

# --- 1. Load Libraries ---
# You may need to install this package: renv::install("ggvenn")
library(ggvenn)
library(ggplot2)

# --- 2. Main Script: Generate Venn Diagrams ---
dir.create("results", showWarnings = FALSE)

# Helper function to read a gene list from a file
read_gene_list <- function(filepath) {
  read.table(filepath, header = FALSE)$V1
}

# --- Venn Diagram 1: Common response to Hcy in parental lines ---
parental_293t <- read_gene_list("results/gene_lists_txt/293t_hcy_vs_met.txt")
parental_468 <- read_gene_list("results/gene_lists_txt/468_hcy_vs_met.txt")

parental_list <- list(
  `293T` = parental_293t,
  `468` = parental_468
)

p1 <- ggvenn(
  parental_list,
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
) + labs(title = "Common Response to Hcy in Parental Lines")
ggsave("results/venn_parental_hcy_response.png", plot = p1, width = 6, height = 5)
print("Saved Venn diagram for parental Hcy response.")

# --- Venn Diagram 2: Baseline differences in revertant vs parental lines ---
revertant_r1 <- read_gene_list("results/gene_lists_txt/r1_vs_293t_met.txt")
revertant_r8 <- read_gene_list("results/gene_lists_txt/r8_vs_468_met.txt")

revertant_list <- list(
  `R1 vs 293T` = revertant_r1,
  `R8 vs 468` = revertant_r8
)

p2 <- ggvenn(
  revertant_list,
  fill_color = c("#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
) + labs(title = "Baseline DEGs in Revertant Lines (in Met Media)")
ggsave("results/venn_revertant_baseline.png", plot = p2, width = 6, height = 5)
print("Saved Venn diagram for revertant baseline differences.")

print("All Venn diagrams have been generated.")
