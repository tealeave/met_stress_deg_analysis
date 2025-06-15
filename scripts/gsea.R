# scripts/gsea_analysis.R
# Purpose: Perform Gene Set Enrichment Analysis (GSEA) on the full DESeq2
# result tables to identify pathways associated with up- or down-regulation.
# --- VERSION 3: Handles p-values of 0 to prevent infinite rank scores ---

# --- 1. Load Libraries ---
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
# install.packages(c("here", "ggplot2", "dplyr", "stringr"))

library(clusterProfiler)
library(org.Hs.eg.db) # Annotation database for Human
library(enrichplot)
library(here)
library(ggplot2)
library(dplyr)
library(stringr) # For string manipulation

# --- 2. Configuration ---
pvalue_cutoff <- 0.05
output_dir <- here("results", "gsea_analysis")
input_dir <- here("results", "tables_csv")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. Helper Function for GSEA ---
perform_gsea <- function(deseq_results_df, contrast_name) {
  print(paste("--- Analyzing:", contrast_name, "---"))

  # --- A. Prepare the Ranked Gene List ---

  # 1. Remove rows with NA p-values
  ranked_list_df <- deseq_results_df %>% filter(!is.na(pvalue))
  
  ## FIX: Handle p-values of exactly 0 ##
  # Find the smallest non-zero p-value in the dataset
  min_pvalue <- min(ranked_list_df$pvalue[ranked_list_df$pvalue > 0], na.rm = TRUE)
  # Replace any 0 p-values with a number smaller than the minimum (or a system minimum)
  # This prevents log10(0) from creating 'Inf' values.
  ranked_list_df$pvalue[ranked_list_df$pvalue == 0] <- min_pvalue * 0.1
  # A more robust alternative is to use the machine's smallest representable number:
  # ranked_list_df$pvalue[ranked_list_df$pvalue == 0] <- .Machine$double.xmin

  # 2. Clean ENSEMBL IDs (remove isoform version)
  ranked_list_df$ensembl <- gsub("\\..*$", "", ranked_list_df[, 1])
  
  # 3. Create the ranking metric
  ranked_list_df$rank_metric <- sign(ranked_list_df$log2FoldChange) * -log10(ranked_list_df$pvalue)

  # 4. Map ENSEMBL IDs to Entrez IDs
  mapped_ids <- bitr(ranked_list_df$ensembl, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = TRUE)
  
  ranked_list_df <- ranked_list_df %>%
    inner_join(mapped_ids, by = c("ensembl" = "ENSEMBL")) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  # 5. Create and sort the final named vector
  gene_list_vector <- ranked_list_df$rank_metric
  names(gene_list_vector) <- ranked_list_df$ENTREZID
  gene_list_vector <- sort(gene_list_vector, decreasing = TRUE)
  
  # --- B. Run GSEA Analysis ---
  save_gsea_results <- function(gsea_result, type) {
    if (!is.null(gsea_result) && nrow(gsea_result) > 0) {
      print(paste("Found", nrow(gsea_result), "significant", type, "terms. Saving results..."))
      gsea_df <- as.data.frame(gsea_result)
      write.csv(gsea_df, file.path(output_dir, paste0(contrast_name, "_GSEA_", type, "_results.csv")))
      
      # Calculate dynamic height based on the number of pathways to plot
      n_up <- sum(gsea_df$NES > 0)
      n_down <- sum(gsea_df$NES < 0)
      pathways_to_plot <- min(n_up, 20) + min(n_down, 20)
      dynamic_height <- max(7, min(15, 5 + pathways_to_plot * 0.28)) # Clamp height

      p1 <- dotplot(gsea_result, showCategory=20, split=".sign") + facet_grid(.~.sign) + labs(title = paste(contrast_name, type, "Enrichment"))
      ggsave(file.path(output_dir, paste0(contrast_name, "_GSEA_", type, "_dotplot.png")), p1, width=12, height=dynamic_height)
      
      p2 <- ridgeplot(gsea_result, showCategory=20) + labs(title = paste(contrast_name, type, "Ridge Plot"))
      ggsave(file.path(output_dir, paste0(contrast_name, "_GSEA_", type, "_ridgeplot.png")), p2, width=12, height=10)
      
      top_up <- gsea_df %>% filter(NES > 0) %>% top_n(3, wt = abs(NES)) %>% pull(ID)
      top_down <- gsea_df %>% filter(NES < 0) %>% top_n(3, wt = abs(NES)) %>% pull(ID)
      
      if(length(c(top_up, top_down)) > 0) {
        p3 <- gseaplot2(gsea_result, geneSetID = c(top_up, top_down), pvalue_table = TRUE)
        ggsave(file.path(output_dir, paste0(contrast_name, "_GSEA_", type, "_gseaplot.png")), p3, width=12, height=10)
      }
    } else {
      print(paste("No significant", type, "terms found."))
    }
  }

  # Run for Gene Ontology (GO)
  gsea_go <- gseGO(geneList=gene_list_vector, OrgDb=org.Hs.eg.db, ont="BP", pvalueCutoff=pvalue_cutoff, minGSSize=15, maxGSSize=500, pAdjustMethod="BH", verbose=FALSE, eps=0)
  save_gsea_results(gsea_go, "GO")

  # Run for KEGG pathways
  gsea_kegg <- gseKEGG(geneList=gene_list_vector, organism='hsa', pvalueCutoff=pvalue_cutoff, minGSSize=15, maxGSSize=500, pAdjustMethod="BH", verbose=FALSE, eps=0)
  save_gsea_results(gsea_kegg, "KEGG")
}

# --- 4. Run Analysis for All DESeq2 contrasts ---
contrast_files <- list.files(path = input_dir, pattern = "\\.csv$", full.names = TRUE)

for (file_path in contrast_files) {
  contrast_name <- gsub("\\.csv$", "", basename(file_path))
  deseq_results <- read.csv(file_path)
  perform_gsea(
    deseq_results_df = deseq_results,
    contrast_name = contrast_name
  )
}

print("--- GSEA analysis script finished. ---")
print(paste("Results are saved in:", output_dir))