# scripts/volcano_plotly.R
# Purpose: Generate interactive volcano plots for specific DESeq2 comparisons using plotly.

# --- 1. Load Libraries ---
# Ensure these libraries are installed: install.packages(c("plotly", "htmlwidgets", "AnnotationDbi"))
# BiocManager::install(c("org.Hs.eg.db", "GO.db"))
library(plotly)
library(htmlwidgets)
library(org.Hs.eg.db)      # For SYMBOL->GO mappings
library(GO.db)             # For GOID->TERM mappings
library(AnnotationDbi)     # General AnnotationDbi functions

# --- 2. Define Plotting Function ---
# A function to create and save an interactive (HTML) volcano plot using plotly
create_interactive_volcano <- function(res, title, interactive_dir) {
  
  # --- Data Preparation ---
  # Convert results to a data frame for easier manipulation
  df <- as.data.frame(res)
  # Strip version numbers from Ensembl IDs (e.g., ENSG0000012345.6 -> ENSG0000012345)
  ensembl_ids <- gsub("\\..*$", "", rownames(df))
  df$ensembl <- rownames(df)
  
  # Map ENSEMBL -> SYMBOL
  # Retrieves official gene symbols for each Ensembl ID
  df$symbol <- mapIds(
    org.Hs.eg.db,
    keys    = ensembl_ids,
    column  = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first" # If an ID maps to multiple symbols, take the first one
  )
  
  # Prepare a list of unique, non-NA symbols for GO mapping
  gene_symbols <- unique(df$symbol[!is.na(df$symbol)])
  
  # --- Map SYMBOL -> GO IDs (Biological Process only) ---
  # Fetch Gene Ontology data for our list of gene symbols
  go_data <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = gene_symbols,
    columns = c("SYMBOL", "GO", "ONTOLOGY"),
    keytype = "SYMBOL"
  )
  
  # Filter for Biological Process (BP) ontology to keep annotations relevant
  go_bp <- subset(go_data, ONTOLOGY == "BP")
  
  # Remove duplicates: ensure one GO term per gene symbol for simplicity in the plot hover
  go_bp_unique <- go_bp[!duplicated(go_bp$SYMBOL), ]
  
  # --- Map GO ID -> TERM via GO.db ---
  # Get the human-readable names for the GO IDs
  go_terms <- AnnotationDbi::select(
    GO.db,
    keys    = unique(go_bp_unique$GO),
    columns = c("GOID", "TERM"),
    keytype = "GOID"
  )
  
  # Merge the GO term names back into our BP table
  go_bp_named <- merge(
    go_bp_unique,
    go_terms,
    by.x  = "GO",
    by.y  = "GOID",
    all.x = TRUE
  )
  
  # Build a fast lookup table (named vector) to map SYMBOL -> TERM
  go_lookup <- setNames(go_bp_named$TERM, go_bp_named$SYMBOL)
  df$go_term <- go_lookup[df$symbol] # Apply the lookup to our main data frame
  
  # Fill any remaining NAs with "N/A" for a clean hover text
  df$go_term[is.na(df$go_term)] <- "N/A"
  
  # --- Create Hover Text ---
  # Custom text to display when hovering over a point in the plot
  df$hover_text <- paste0(
    "Gene: ", ifelse(is.na(df$symbol), df$ensembl, df$symbol),
    "<br>log2FC: ", round(df$log2FoldChange, 2),
    "<br>-log10(p): ", round(-log10(df$pvalue), 2),
    "<br>Process: ", df$go_term
  )
  
  # Remove rows with NA p-values, which cannot be plotted
  df <- df[!is.na(df$pvalue), ]
  
  # --- Assign Significance Status ---
  # Categorize genes based on adjusted p-value and fold change thresholds
  df$status <- "Not Significant"
  df$status[df$padj < 0.05 & df$log2FoldChange >  1] <- "Up-regulated"
  df$status[df$padj < 0.05 & df$log2FoldChange < -1] <- "Down-regulated"
  df$status <- factor(df$status,
                      levels = c("Up-regulated", "Down-regulated", "Not Significant"))
  
  # Split data by status for separate plotting traces (enables different colors/legends)
  df_up  <- subset(df, status == "Up-regulated")
  df_dn  <- subset(df, status == "Down-regulated")
  df_ns  <- subset(df, status == "Not Significant")
  
  # Calculate plot axis ranges dynamically
  max_y <- max(-log10(df$pvalue), na.rm = TRUE) * 1.05 # Add 5% buffer
  max_abs_x <- max(abs(df$log2FoldChange), na.rm = TRUE) * 1.05
  
  # --- Build Interactive Volcano with plotly ---
  p_interactive <- plot_ly(type = 'scatter', mode = 'markers') %>%
    add_trace(
      data = df_dn,
      x    = ~log2FoldChange, y = ~-log10(pvalue),
      text = ~hover_text, hoverinfo = 'text',
      name = 'Down-regulated', marker = list(color = '#377EB8', size = 6, opacity = 0.7)
    ) %>%
    add_trace(
      data = df_up,
      x    = ~log2FoldChange, y = ~-log10(pvalue),
      text = ~hover_text, hoverinfo = 'text',
      name = 'Up-regulated', marker = list(color = '#E41A1C', size = 6, opacity = 0.7)
    ) %>%
    add_trace(
      data = df_ns,
      x    = ~log2FoldChange, y = ~-log10(pvalue),
      text = ~hover_text, hoverinfo = 'text',
      name = 'Not Significant', marker = list(color = 'grey', size = 6, opacity = 0.5)
    ) %>%
    layout(
      title = list(text = title, x = 0.5), # Center title
      xaxis = list(title = "Log<sub>2</sub> Fold Change", range = c(-max_abs_x, max_abs_x)),
      yaxis = list(title = "-Log<sub>10</sub> P-value", range = c(0, max_y)),
      # Add dashed lines for significance cutoffs
      shapes = list(
        list(type='line', x0=-1, x1=-1, y0=0, y1=max_y, line=list(dash='dash', color='grey')),
        list(type='line', x0= 1, x1= 1, y0=0, y1=max_y, line=list(dash='dash', color='grey')),
        list(type='line', x0=-max_abs_x, x1=max_abs_x, y0=-log10(0.05), y1=-log10(0.05), line=list(dash='dash', color='grey'))
      ),
      hoverlabel = list(bgcolor="white", font=list(size=12)),
      legend     = list(orientation='h', xanchor="center", x=0.5, y=-0.2) # Horizontal legend below plot
    )
  
  # --- Save Plot ---
  # Construct a clean filename and save as a self-contained HTML widget
  safe_title <- gsub("[^a-zA-Z0-9_]", "_", title) # Sanitize title for filename
  interactive_path <- file.path(
    interactive_dir,
    paste0("volcano_interactive_", safe_title, ".html")
  )
  saveWidget(p_interactive, file = interactive_path, selfcontained = TRUE)
  message("Saved interactive plot to: ", interactive_path)
}

# --- 3. Main Script: Define Comparisons and Generate Plots ---
# Define directories for input objects and output plots
rds_dir <- "results/deseq_objects"
interactive_output_dir <- "results/volcano_plots_interactive"

# Create the output directory if it doesn't exist
dir.create(interactive_output_dir, showWarnings = FALSE, recursive = TRUE)

# --- UPDATED: Define the specific list of comparisons to generate ---
comparisons_to_plot <- c(
  # Comparisons within batch 'exp2'
  "468_hcy_vs_met",
  "r8_hcy_vs_met",
  "r8_vs_468_met",
  "r8_vs_468_hcy",
  # Comparisons within batch 'exp1'
  "293t_hcy_vs_met",
  "r1_hcy_vs_met",
  "r1_vs_293t_met",
  "r1_vs_293t_hcy"
)

message("Starting volcano plot generation for ", length(comparisons_to_plot), " comparisons.")

# --- UPDATED: Loop through the defined list instead of all files ---
for (comp_name in comparisons_to_plot) {
  # Construct the full path to the expected RDS file
  file_path <- file.path(rds_dir, paste0(comp_name, ".rds"))
  
  # Check if the file exists before trying to process it
  if (file.exists(file_path)) {
    tryCatch({
      # Read the DESeq2 results object
      res_object <- readRDS(file_path)
      
      # Use the comparison name as the plot title
      plot_title <- comp_name 
      
      # Call the function to create and save the plot
      create_interactive_volcano(res_object, plot_title, interactive_output_dir)
      
    }, error = function(e) {
      # Report any errors that occur during file processing
      warning("Could not process ", file_path, ": ", e$message)
    })
  } else {
    # Warn if a specific RDS file is not found
    warning("RDS file not found for comparison: ", comp_name, ". Skipping.")
  }
}

message("All interactive volcano plots have been generated.")
