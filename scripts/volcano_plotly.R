# scripts/volcano_interactive.R
# Purpose: Generate interactive volcano plots for each DESeq2 comparison using plotly.

# --- 1. Load Libraries ---
library(plotly)
library(htmlwidgets)
library(org.Hs.eg.db)     # For SYMBOL→GO mappings
library(GO.db)            # For GOID→TERM mappings
library(AnnotationDbi)    # General AnnotationDbi functions

# --- 2. Define Plotting Function ---
# A function to create and save an interactive (HTML) volcano plot using plotly
create_interactive_volcano <- function(res, title, interactive_dir) {
  
  # --- Data Preparation ---
  df <- as.data.frame(res)
  ensembl_ids <- gsub("\\..*$", "", rownames(df))
  df$ensembl <- rownames(df)
  
  # Map ENSEMBL → SYMBOL
  df$symbol <- mapIds(
    org.Hs.eg.db,
    keys    = ensembl_ids,
    column  = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # Prepare list of non-NA symbols
  gene_symbols <- unique(df$symbol[!is.na(df$symbol)])
  
  # --- Map SYMBOL → GO IDs (Biological Process only) ---
  go_data <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = gene_symbols,
    columns = c("SYMBOL", "GO", "ONTOLOGY"),
    keytype = "SYMBOL"
  )
  
  # Filter for BP ontology
  go_bp <- subset(go_data, ONTOLOGY == "BP")
  
  # Remove duplicates: keep first GO per SYMBOL
  go_bp_unique <- go_bp[!duplicated(go_bp$SYMBOL), ]
  
  # --- Map GO ID → TERM via GO.db ---
  go_terms <- AnnotationDbi::select(
    GO.db,
    keys    = unique(go_bp_unique$GO),
    columns = c("GOID", "TERM"),
    keytype = "GOID"
  )
  
  # Merge term names into the BP table
  go_bp_named <- merge(
    go_bp_unique,
    go_terms,
    by.x   = "GO",
    by.y   = "GOID",
    all.x  = TRUE
  )
  
  # Build a lookup: SYMBOL → TERM
  go_lookup <- setNames(go_bp_named$TERM, go_bp_named$SYMBOL)
  df$go_term <- go_lookup[df$symbol]
  
  # Wherever still NA, fill with "N/A"
  df$go_term[is.na(df$go_term)] <- "N/A"
  
  # --- Create Hover Text ---
  df$hover_text <- paste0(
    "Gene: ", ifelse(is.na(df$symbol), df$ensembl, df$symbol),
    "<br>log2FC: ", round(df$log2FoldChange, 2),
    "<br>-log10(p): ", round(-log10(df$pvalue), 2),
    "<br>Process: ", df$go_term
  )
  
  # Remove rows with NA p-values
  df <- df[!is.na(df$pvalue), ]
  
  # Significance status
  df$status <- "Not Significant"
  df$status[df$padj < 0.05 & df$log2FoldChange >  1] <- "Up-regulated"
  df$status[df$padj < 0.05 & df$log2FoldChange < -1] <- "Down-regulated"
  df$status <- factor(df$status,
                      levels = c("Up-regulated", "Down-regulated", "Not Significant"))
  
  # Split for plotting
  df_up  <- subset(df, status == "Up-regulated")
  df_dn  <- subset(df, status == "Down-regulated")
  df_ns  <- subset(df, status == "Not Significant")
  
  # Plot ranges
  max_y <- max(-log10(df$pvalue), na.rm = TRUE)
  min_x <- min(df$log2FoldChange, na.rm = TRUE)
  max_x <- max(df$log2FoldChange, na.rm = TRUE)
  
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
      title = title,
      xaxis = list(title = "Log<sub>2</sub> Fold Change", range = c(min_x, max_x)),
      yaxis = list(title = "-Log<sub>10</sub> P-value", range = c(0, max_y)),
      shapes = list(
        list(type='line', x0=-1, x1=-1, y0=0,    y1=max_y, line=list(dash='dash')),
        list(type='line', x0= 1, x1= 1, y0=0,    y1=max_y, line=list(dash='dash')),
        list(type='line', x0=min_x, x1=max_x, y0=-log10(0.05), y1=-log10(0.05), line=list(dash='dash'))
      ),
      hoverlabel = list(bgcolor="white", font=list(size=12)),
      legend     = list(orientation='h', xanchor="center", x=0.5)
    )
  
  # Save as self-contained HTML widget
  interactive_path <- file.path(
    interactive_dir,
    paste0("volcano_interactive_", gsub(" ", "_", title), ".html")
  )
  saveWidget(p_interactive, file = interactive_path, selfcontained = TRUE)
  message("Saved interactive plot to: ", interactive_path)
}

# --- 3. Main Script: Generate Plots ---
rds_dir <- "results/deseq_objects"
interactive_output_dir <- "results/volcano_plots_interactive"

dir.create(interactive_output_dir, showWarnings = FALSE, recursive = TRUE)

rds_files <- list.files(rds_dir, pattern="\\.rds$", full.names=TRUE)

if (length(rds_files) == 0) {
  warning("No RDS files found in ", rds_dir)
} else {
  for (file in rds_files) {
    tryCatch({
      res_object <- readRDS(file)
      comp_name  <- sub("\\.rds$", "", basename(file))
      create_interactive_volcano(res_object, comp_name, interactive_output_dir)
    }, error = function(e) {
      warning("Could not process ", file, ": ", e$message)
    })
  }
}

message("All interactive volcano plots have been generated.")
