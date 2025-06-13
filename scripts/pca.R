# scripts/pca.R
# Purpose: Generate PCA plots to visualize the overall structure of the data.

# --- 1. Load Libraries ---
library(DESeq2)
library(ggplot2)

# --- 2. Load DESeq2 Data Object ---
# We need the 'dds' object which was created in the main analysis script.
# To avoid re-running the entire DESeq2 analysis, we will create it again here.
# NOTE: This part is duplicated from 'run_deseq_analysis.R' for modularity.

# Load count data
feature_counts <- read.table("data/counts.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
count_data <- feature_counts[, c(1, 7:ncol(feature_counts))]
rownames(count_data) <- count_data$Geneid
count_data$Geneid <- NULL

# Load sample information
col_data <- read.csv("data/colData.csv", header=TRUE)
rownames(col_data) <- col_data$sample_name

# *** THIS IS THE FIX: Create the 'group' column ***
col_data$group <- factor(paste0(col_data$cell_line, "_", col_data$treatment))

# Ensure data is in the correct order
count_data <- count_data[, rownames(col_data)]


# Create the full DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ 1) # A simple design is fine as we only need the object

# --- 3. Transform Data and Generate Plots ---
# For visualization, it's best to work with transformed data.
# The Variance Stabilizing Transformation (VST) is recommended for PCA.
print("Applying Variance Stabilizing Transformation (vst)...")
vsd <- vst(dds, blind=TRUE) # blind=TRUE is recommended for QC plots

# Create a directory for PCA plots
dir.create("results/pca_plots", showWarnings = FALSE, recursive = TRUE)
print("Generating PCA plots...")

# --- Plot 1: Color by 'group' (cell line + treatment) ---
p1 <- plotPCA(vsd, intgroup="group") +
  labs(title = "PCA colored by Group") +
  geom_point(size=4) + # Make points larger
  theme_bw(base_size = 14) +
  coord_fixed() # Fix aspect ratio
ggsave("results/pca_plots/pca_by_group.png", plot=p1, width=8, height=8) # Make plot square

# --- Plot 2: Color by 'cell_line' ---
p2 <- plotPCA(vsd, intgroup="cell_line") +
  labs(title = "PCA colored by Cell Line") +
  geom_point(size=4) +
  theme_bw(base_size = 14) +
  coord_fixed() # Fix aspect ratio
ggsave("results/pca_plots/pca_by_cell_line.png", plot=p2, width=8, height=8) # Make plot square

# --- Plot 3: Color by 'treatment' ---
p3 <- plotPCA(vsd, intgroup="treatment") +
  labs(title = "PCA colored by Treatment") +
  geom_point(size=4) +
  theme_bw(base_size = 14) +
  coord_fixed() # Fix aspect ratio
ggsave("results/pca_plots/pca_by_treatment.png", plot=p3, width=8, height=8) # Make plot square

# --- Plot 4: Color by 'batch' ---
p4 <- plotPCA(vsd, intgroup="batch") +
  labs(title = "PCA colored by Batch") +
  geom_point(size=4) +
  theme_bw(base_size = 14) +
  coord_fixed() # Fix aspect ratio
ggsave("results/pca_plots/pca_by_batch.png", plot=p4, width=8, height=8) # Make plot square

print("All PCA plots have been generated in 'results/pca_plots/'.")

