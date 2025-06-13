# scripts/upset.R
# Purpose: Generate an Upset plot to visualize intersections of Hcy response genes.

# --- 1. Load Libraries ---
# You may need to install this: renv::install("UpSetR")
library(UpSetR)

# --- 2. Main Script: Generate Upset Plot ---
dir.create("results/upset_plot", showWarnings = FALSE, recursive = TRUE)

# Helper function to read a gene list from a file
read_gene_list <- function(filepath) {
  read.table(filepath, header = FALSE)$V1
}

# Load the four gene lists related to Hcy response
list_input <- list(
  `293t_hcy_vs_met` = read_gene_list("results/gene_lists_txt/293t_hcy_vs_met.txt"),
  `r1_hcy_vs_met`   = read_gene_list("results/gene_lists_txt/r1_hcy_vs_met.txt"),
  `468_hcy_vs_met`  = read_gene_list("results/gene_lists_txt/468_hcy_vs_met.txt"),
  `r8_hcy_vs_met`   = read_gene_list("results/gene_lists_txt/r8_hcy_vs_met.txt")
)

# Create the Upset plot
png("results/upset_plot/upset_hcy_response.png", width = 1000, height = 600, res = 100)
upset(
  fromList(list_input),
  order.by = "freq",
  nsets = 4,
  point.size = 3.5,
  line.size = 2,
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  text.scale = c(1.3, 1.3, 1, 1, 1.5, 1.2)
)
dev.off()

print("Upset plot has been generated and saved.")

