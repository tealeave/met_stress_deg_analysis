# scripts/deseq_analysis.R
# Purpose: Perform the core differential expression analysis and save results
# in multiple formats for downstream scripts.

# --- 1. Load Libraries ---
library(DESeq2)

# --- 2. Data Loading and Preparation ---
# Load the full featureCounts output
feature_counts <- read.table("data/counts.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Select only the Geneid (column 1) and the actual count columns (starting from column 7)
count_data <- feature_counts[, c(1, 7:ncol(feature_counts))]
rownames(count_data) <- count_data$Geneid
count_data$Geneid <- NULL

# Load the sample information and set row names
col_data <- read.csv("data/colData.csv", header=TRUE)
rownames(col_data) <- col_data$sample_name

# Make sure the column names in count_data match the row names in col_data
count_data <- count_data[, rownames(col_data)]

# Create a 'group' variable by combining cell_line and treatment
# This single factor will contain all experimental conditions.
col_data$group <- factor(paste0(col_data$cell_line, "_", col_data$treatment))

# ---
# NOTE ON DESIGN: The 'batch' variable (exp1/exp2) is completely confounded
# with the cell lines (R1/293T are only in exp1, R8/468 are only in exp2).
# Adding '~ batch + group' would cause a model matrix error.
#
# The correct approach for this design is to use '~ group' and make
# comparisons ONLY within each batch. The batch effect is controlled for
# by the nature of the contrasts themselves.
# ---
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ group)

# Pre-filter for genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
print("Running DESeq2... this may take a few minutes.")
dds <- DESeq(dds)
print("DESeq2 analysis complete.")

# --- 3. Define, Run, and Save All Comparisons ---
# Create directories to store the different output types
dir.create("results/deseq_objects", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables_csv", showWarnings = FALSE, recursive = TRUE)
dir.create("results/gene_lists_txt", showWarnings = FALSE, recursive = TRUE)

# Define a list of comparisons to make.
# Each of these comparisons is made between groups from the SAME batch,
# which is the correct way to analyze this confounded experiment.
comparisons <- list(
  # Comparisons within batch 'exp2'
  "468_hcy_vs_met"   = c("group", "468_Hcy", "468_Met"),
  "r8_hcy_vs_met"    = c("group", "R8_Hcy", "R8_Met"),
  "r8_vs_468_met"    = c("group", "R8_Met", "468_Met"),
  "r8_vs_468_hcy"    = c("group", "R8_Hcy", "468_Hcy"),
  # Comparisons within batch 'exp1'
  "293t_hcy_vs_met"  = c("group", "293T_Hcy", "293T_Met"),
  "r1_hcy_vs_met"    = c("group", "R1_Hcy", "R1_Met"),
  "r1_vs_293t_met"   = c("group", "R1_Met", "293T_Met"),
  "r1_vs_293t_hcy"   = c("group", "R1_Hcy", "293T_Hcy")
)

alpha <- 0.05 # Set the significance threshold

# Loop through each comparison, get results, and save them
for (name in names(comparisons)) {
  
  print(paste("Processing comparison:", name))
  
  # Get the results for the current comparison
  res <- results(dds, contrast=comparisons[[name]], alpha=alpha)
  res <- na.omit(res) # Remove rows with NA values
  
  # 1. Save the full result object
  saveRDS(res, file=paste0("results/deseq_objects/", name, ".rds"))
  
  # 2. Save the full results table as a CSV
  write.csv(as.data.frame(res), file=paste0("results/tables_csv/", name, ".csv"))
  
  # 3. Filter for significant genes and save the list as a TXT file
  sig_genes <- subset(res, padj < alpha)
  write.table(
    rownames(sig_genes),
    file=paste0("results/gene_lists_txt/", name, ".txt"),
    row.names=FALSE,
    col.names=FALSE,
    quote=FALSE
  )
}

print("Script finished. All results are saved in the 'results/' directory.")
