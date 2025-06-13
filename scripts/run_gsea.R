# scripts/run_gsea.R

# --- 1. Load Libraries ---
library(clusterProfiler)
library(org.Hs.eg.db) # Annotation database for Human
library(AnnotationDbi)
library(ggplot2)

# --- 2. Define the GSEA Function ---
# This function takes a DESeq2 result object and a comparison name,
# runs GSEA for GO and KEGG, and saves the results and plots.

perform_gsea <- function(deseq_results, comparison_name) {

  print(paste("--- Starting GSEA for:", comparison_name, "---"))

  # A. Prepare the Ranked Gene List
  res_gsea <- na.omit(deseq_results)
  gene_list <- res_gsea$stat
  names(gene_list) <- rownames(res_gsea)
  gene_list <- sort(gene_list, decreasing = TRUE)

  # B. Run Gene Ontology (GO) GSEA
  print("Running GO GSEA...")
  gse_go <- gseGO(geneList     = gene_list,
                  OrgDb        = org.Hs.eg.db,
                  keyType      = "ENSEMBL",
                  ont          = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  verbose      = FALSE) # Set to FALSE to keep output clean

  if (nrow(gse_go) > 0) {
    write.csv(as.data.frame(gse_go), file=paste0("results/gsea_go_", comparison_name, ".csv"))
    
    # Create and save dotplot
    gse_plot <- dotplot(gse_go, showCategory=20, split=".sign") + facet_grid(.~.sign)
    ggsave(paste0("results/gsea_go_", comparison_name, ".png"), plot=gse_plot, width=10, height=7)
    print("GO GSEA complete. Results and plot saved.")
  } else {
    print("No significant GO terms found.")
  }

  # C. Run KEGG Pathway GSEA
  print("Running KEGG GSEA...")
  
  # Convert ENSEMBL to ENTREZ IDs for KEGG
  entrez_ids <- mapIds(org.Hs.eg.db, keys=names(gene_list), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
  gene_list_entrez <- gene_list[!is.na(entrez_ids)]
  names(gene_list_entrez) <- entrez_ids[!is.na(entrez_ids)]

  if (length(gene_list_entrez) > 0) {
      gse_kegg <- gseKEGG(geneList     = gene_list_entrez,
                          organism     = 'hsa',
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          verbose      = FALSE)

      if (nrow(gse_kegg) > 0) {
        write.csv(as.data.frame(gse_kegg), file=paste0("results/gsea_kegg_", comparison_name, ".csv"))
        print("KEGG GSEA complete. Results saved.")
      } else {
        print("No significant KEGG pathways found.")
      }
  } else {
      print("No genes could be mapped to ENTREZ IDs for KEGG analysis.")
  }

  print(paste("--- Finished GSEA for:", comparison_name, "---"))
  
}


# --- 3. Main Script: Loop Through All Comparisons ---

# List all the result files we need to process
rds_files <- list.files("results/deseq_objects", pattern = "\\.rds$", full.names = TRUE)

# Loop over each file, load it, and run the GSEA function
for (file_path in rds_files) {
  # Extract a clean name for the comparison from the file path
  comparison_name <- gsub(".rds", "", basename(file_path))
  
  # Load the DESeq2 result object
  deseq_results <- readRDS(file_path)
  
  # Run our function
  perform_gsea(deseq_results, comparison_name)
}

print("All GSEA analyses are complete.")