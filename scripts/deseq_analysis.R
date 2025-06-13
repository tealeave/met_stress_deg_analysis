# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(UpSetR)
library(gplots)

# --- 1. Data Loading and Preparation ---
# Load the full featureCounts output
# Note: We use "data/counts.tsv" because we are running from the project root
feature_counts <- read.table("data/counts.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Select only the Geneid (column 1) and the actual count columns (starting from column 7)
count_data <- feature_counts[, c(1, 7:ncol(feature_counts))]

# Set the row names to be the Geneids
rownames(count_data) <- count_data$Geneid
count_data$Geneid <- NULL # Remove the now-redundant Geneid column

# Load the sample information and set row names
col_data <- read.csv("data/colData.csv", header=TRUE)
rownames(col_data) <- col_data$sample_name

# IMPORTANT: Make sure the column names in your count_data match the row names in col_data
# This line ensures they are in the same order
count_data <- count_data[, rownames(col_data)]

# Create a 'group' variable by combining cell_line and treatment
col_data$group <- factor(paste0(col_data$cell_line, "_", col_data$treatment))

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ batch + group)

# Pre-filter for genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# --- 2. Answering Your Specific Questions ---

# Set the significance threshold
alpha <- 0.05

# Question 1: Common response to Hcy in 468 and 293T
res_468_hcy_vs_met <- results(dds, contrast=c("group", "468_Hcy", "468_Met"), alpha=alpha)
res_293t_hcy_vs_met <- results(dds, contrast=c("group", "293T_Hcy", "293T_Met"), alpha=alpha)

degs_468 <- rownames(res_468_hcy_vs_met[which(res_468_hcy_vs_met$padj < alpha), ])
degs_293t <- rownames(res_293t_hcy_vs_met[which(res_293t_hcy_vs_met$padj < alpha), ])

common_response_hcy <- intersect(degs_468, degs_293t)
print(paste("Found", length(common_response_hcy), "common response genes to Hcy in 468 and 293T cells."))

# Question 2: Which common response pathways are maintained in R1 and R8
res_r1_hcy_vs_met <- results(dds, contrast=c("group", "R1_Hcy", "R1_Met"), alpha=alpha)
res_r8_hcy_vs_met <- results(dds, contrast=c("group", "R8_Hcy", "R8_Met"), alpha=alpha)

degs_r1_hcy <- rownames(res_r1_hcy_vs_met[which(res_r1_hcy_vs_met$padj < alpha), ])
degs_r8_hcy <- rownames(res_r8_hcy_vs_met[which(res_r8_hcy_vs_met$padj < alpha), ])

maintained_in_r1 <- intersect(common_response_hcy, degs_r1_hcy)
maintained_in_r8 <- intersect(common_response_hcy, degs_r8_hcy)
maintained_in_both <- intersect(maintained_in_r1, maintained_in_r8)

print(paste("Found", length(maintained_in_both), "common response genes maintained in both R1 and R8 cells."))

# Question 3: Different pathways induced in Hcy in revertants vs parents
# Genes induced in R1 but not in 293T
diff_r1 <- setdiff(degs_r1_hcy, degs_293t)
# Genes induced in R8 but not in 468
diff_r8 <- setdiff(degs_r8_hcy, degs_468)

print(paste("Found", length(diff_r1), "genes induced in Hcy in R1 but not 293T."))
print(paste("Found", length(diff_r8), "genes induced in Hcy in R8 but not 468."))

# Question 4: Different expression profile in Met in R1/R8 vs parents
res_r1_vs_293t_met <- results(dds, contrast=c("group", "R1_Met", "293T_Met"), alpha=alpha)
res_r8_vs_468_met <- results(dds, contrast=c("group", "R8_Met", "468_Met"), alpha=alpha)

degs_r1_vs_293t_met <- rownames(res_r1_vs_293t_met[which(res_r1_vs_293t_met$padj < alpha), ])
degs_r8_vs_468_met <- rownames(res_r8_vs_468_met[which(res_r8_vs_468_met$padj < alpha), ])

similar_in_revertants_met <- intersect(degs_r1_vs_293t_met, degs_r8_vs_468_met)
print(paste("Found", length(similar_in_revertants_met), "genes with similar expression changes in R1 and R8 in Met media compared to parents."))


# --- 3. Saving Results ---
# Create a directory to save the results
dir.create("deseq2_results", showWarnings = FALSE)

# Save the results of each comparison to a CSV file
write.csv(as.data.frame(res_468_hcy_vs_met), file="deseq2_results/468_hcy_vs_met.csv")
write.csv(as.data.frame(res_293t_hcy_vs_met), file="deseq2_results/293t_hcy_vs_met.csv")
write.csv(as.data.frame(res_r1_hcy_vs_met), file="deseq2_results/r1_hcy_vs_met.csv")
write.csv(as.data.frame(res_r8_hcy_vs_met), file="deseq2_results/r8_hcy_vs_met.csv")
write.csv(as.data.frame(res_r1_vs_293t_met), file="deseq2_results/r1_vs_293t_met.csv")
write.csv(as.data.frame(res_r8_vs_468_met), file="deseq2_results/r8_vs_468_met.csv")

# Save the lists of genes
write.table(common_response_hcy, "deseq2_results/common_response_hcy.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(maintained_in_both, "deseq2_results/maintained_in_revertants.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(diff_r1, "deseq2_results/diff_pathways_r1.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(diff_r8, "deseq2_results/diff_pathways_r8.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(similar_in_revertants_met, "deseq2_results/similar_in_revertants_met.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)