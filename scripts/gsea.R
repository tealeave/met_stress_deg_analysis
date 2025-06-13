# scripts/run_gsea.R
# Purpose: Run Gene Set Enrichment Analysis (GSEA) on all DESeq2 results.
# This script automatically finds all DESeq2 output files, runs GSEA for
# Gene Ontology (GO) and KEGG pathways, and saves the results into
# separate files for up-regulated and down-regulated pathways.

# --- 1. Load Libraries ---
library(clusterProfiler)
library(org.Hs.eg.db) # Annotation database for Human
library(AnnotationDbi)
library(ggplot2)

# --- 2. Define the GSEA Function ---
# This function takes a DESeq2 result object and a comparison name,
# runs GSEA, and saves the results and plots.

perform_gsea <- function(deseq_results, comparison_name, output_dir) {

  print(paste("--- Starting GSEA for:", comparison_name, "---"))

  # A. Prepare the Ranked Gene List
  # The list is ranked by the 'stat' column from DESeq2 results.
  # A positive stat value indicates up-regulation, negative is down-regulation.
  res_gsea <- na.omit(deseq_results)
  gene_list <- res_gsea$stat
  
  # Clean ENSEMBL IDs by removing version numbers (e.g., "ENSG0000012345.6" -> "ENSG0000012345")
  names(gene_list) <- gsub("\\..*$", "", rownames(res_gsea))
  
  # Sort the gene list in decreasing order (a GSEA requirement)
  gene_list <- sort(gene_list, decreasing = TRUE)

  # B. Run Gene Ontology (GO) GSEA for Biological Processes (BP)
  print("Running GO GSEA...")
  tryCatch({
    gse_go <- gseGO(geneList      = gene_list,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = "ENSEMBL",
                    ont           = "BP", # Biological Process
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    verbose       = FALSE)

    if (!is.null(gse_go) && nrow(gse_go@result) > 0) {
      # Split results into up-regulated and down-regulated pathways
      gse_go_results <- as.data.frame(gse_go@result)
      up_regulated <- subset(gse_go_results, enrichmentScore > 0)
      down_regulated <- subset(gse_go_results, enrichmentScore < 0)

      # Save the split results to separate CSV files
      write.csv(up_regulated, file = file.path(output_dir, paste0("gsea_go_", comparison_name, "_upregulated.csv")), row.names = FALSE)
      write.csv(down_regulated, file = file.path(output_dir, paste0("gsea_go_", comparison_name, "_downregulated.csv")), row.names = FALSE)
      
      # Create and save a dotplot faceted by sign (activated vs. suppressed)
      # Increased height to prevent label overlap
      gse_plot <- dotplot(gse_go, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign)
      ggsave(file.path(output_dir, paste0("gsea_go_", comparison_name, ".png")), plot = gse_plot, width = 12, height = 12, dpi = 300)
      print("GO GSEA complete. Results and plot saved.")
    } else {
      print("No significant GO terms found.")
    }
  }, error = function(e) {
    print(paste("An error occurred during GO GSEA:", e$message))
  })

  # C. Run KEGG Pathway GSEA
  print("Running KEGG GSEA...")
  
  # Convert ENSEMBL to ENTREZ IDs for KEGG analysis
  entrez_ids <- mapIds(org.Hs.eg.db, keys = names(gene_list), column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  gene_list_entrez <- gene_list[!is.na(entrez_ids)]
  names(gene_list_entrez) <- entrez_ids[!is.na(entrez_ids)]

  if (length(gene_list_entrez) > 0) {
    tryCatch({
      gse_kegg <- gseKEGG(geneList      = gene_list_entrez,
                          organism      = 'hsa', # Homo sapiens
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          verbose       = FALSE)

      if (!is.null(gse_kegg) && nrow(gse_kegg@result) > 0) {
        # Split results into up-regulated and down-regulated pathways
        gse_kegg_results <- as.data.frame(gse_kegg@result)
        up_regulated_kegg <- subset(gse_kegg_results, enrichmentScore > 0)
        down_regulated_kegg <- subset(gse_kegg_results, enrichmentScore < 0)

        # Save the split results to separate CSV files
        write.csv(up_regulated_kegg, file = file.path(output_dir, paste0("gsea_kegg_", comparison_name, "_upregulated.csv")), row.names = FALSE)
        write.csv(down_regulated_kegg, file = file.path(output_dir, paste0("gsea_kegg_", comparison_name, "_downregulated.csv")), row.names = FALSE)

        # Create and save KEGG dotplot
        kegg_plot <- dotplot(gse_kegg, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign)
        ggsave(file.path(output_dir, paste0("gsea_kegg_", comparison_name, ".png")), plot = kegg_plot, width = 12, height = 12, dpi = 300)
        print("KEGG GSEA complete. Results and plot saved.")
      } else {
        print("No significant KEGG pathways found.")
      }
    }, error = function(e) {
      print(paste("An error occurred during KEGG GSEA:", e$message))
    })
  } else {
    print("No genes could be mapped to ENTREZ IDs for KEGG analysis.")
  }

  print(paste("--- Finished GSEA for:", comparison_name, "---"))
}


# --- 3. Main Script: Loop Through All Comparisons ---

# Define input and output directories
rds_dir <- "results/deseq_objects"
gsea_output_dir <- "results/gsea_results"

# Stop if the input directory doesn't exist
if (!dir.exists(rds_dir)) {
  stop("Directory '", rds_dir, "' not found. Please ensure it exists and contains your DESeq2 results.")
}

# Create the output directory for GSEA results
dir.create(gsea_output_dir, showWarnings = FALSE, recursive = TRUE)

# Get a list of all .rds files from the DESeq2 analysis
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(rds_files) == 0) {
  print(paste("No .rds files found in '", rds_dir, "'. No analysis will be run."))
} else {
  # Loop over each file, load it, and run the GSEA function
  for (file_path in rds_files) {
    # Extract a clean name for the comparison from the file path
    comparison_name <- gsub("\\.rds$", "", basename(file_path))
    
    # Load the DESeq2 result object
    deseq_results <- readRDS(file_path)
    
    # Run our GSEA function
    perform_gsea(deseq_results, comparison_name, gsea_output_dir)
  }
}

print("All GSEA analyses are complete.")

