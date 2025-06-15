# scripts/interpret_gsea_V8.R
# Purpose: A complete script to analyze and interpret GSEA results for GO and KEGG.
# --- VERSION 8: Fixed pheatmap error when clustering a single row. ---

# --- 1. Load Libraries ---
library(dplyr)
library(ggplot2)
library(here)
library(ggrepel)
library(pheatmap)

# --- 2. Configuration ---
gsea_dir <- here("results", "gsea_analysis")
output_dir <- here("results", "final_interpretation")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
label_threshold <- 0.10

# --- 3. Helper Functions ---
create_nes_scatterplot <- function(data, x_col, y_col, p_adj_x_col, p_adj_y_col, title, filename) {
  if (nrow(data) == 0) {
    print(paste("Skipping plot", basename(filename), "- no data."))
    return()
  }
  data <- data %>%
    mutate(label = if_else(!!sym(p_adj_x_col) < label_threshold & !!sym(p_adj_y_col) < label_threshold, Description, ""))
  
  p <- ggplot(data, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(aes(color = !!sym(p_adj_y_col)), alpha = 0.7) +
    geom_text_repel(aes(label = label), size = 3.5, max.overlaps = 20, box.padding = 0.5) +
    scale_color_gradient(low = "red", high = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = title, subtitle = "Pathways in top-right (common up) and bottom-left (common down)", x = paste("NES:", x_col), y = paste("NES:", y_col), color = "Adj. P-value") +
    theme_bw(base_size = 14)
    
  ggsave(filename, p, width = 12, height = 10)
  print(paste("Saved plot:", filename))
}

create_nes_heatmap <- function(data, title, filename) {
  if (nrow(data) == 0) {
    print(paste("Skipping heatmap", basename(filename), "- no data."))
    return()
  }
  
  heatmap_matrix <- data %>% select(starts_with("NES")) %>% as.matrix()
  rownames(heatmap_matrix) <- data$Description
  
  ## FIX: Disable row clustering if there's only one row, as clustering requires n >= 2.
  should_cluster_rows <- if (nrow(heatmap_matrix) > 1) TRUE else FALSE
  
  palette_breaks <- seq(-max(abs(heatmap_matrix)), max(abs(heatmap_matrix)), length.out = 101)
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  
  pheatmap(
    heatmap_matrix,
    main = title,
    filename = filename,
    width = 8,
    height = if (should_cluster_rows) 10 else 4, # Use smaller height for single-row plot
    color = color_palette,
    breaks = palette_breaks,
    cluster_cols = FALSE,
    cluster_rows = should_cluster_rows, # Use the dynamic value here
    fontsize_row = 8
  )
  print(paste("Saved heatmap:", filename))
}


# ==============================================================================
# Main Analysis Function
# ==============================================================================
run_interpretation_for_db <- function(db_type = "GO") {
  
  print(paste("############################################################"))
  print(paste("###   STARTING ANALYSIS FOR:", db_type, "PATHWAYS   ###"))
  print(paste("############################################################"))
  
  fname <- function(base_name) { paste0(base_name, "_", db_type) }

  # --- Q1: Common Parental Hcy Response ---
  print(paste("--- Q1: Common parental Hcy response (", db_type, ") ---"))
  file1_q1 <- file.path(gsea_dir, paste0("293t_hcy_vs_met_GSEA_", db_type, "_results.csv"))
  file2_q1 <- file.path(gsea_dir, paste0("468_hcy_vs_met_GSEA_", db_type, "_results.csv"))

  if (file.exists(file1_q1) && file.exists(file2_q1)) {
    parent_293t <- read.csv(file1_q1)
    parent_468 <- read.csv(file2_q1)
    common_parental <- inner_join(parent_293t, parent_468, by = "ID", suffix = c(".293t", ".468")) %>%
      filter(sign(NES.293t) == sign(NES.468)) %>%
      select(ID, Description = Description.293t, NES.293t, p.adjust.293t, NES.468, p.adjust.468) %>%
      arrange(p.adjust.293t)
    
    write.csv(common_parental, file.path(output_dir, paste0(fname("Q1_common_parental_response"), ".csv")), row.names = FALSE)
    print(paste("Found", nrow(common_parental), "commonly regulated parental pathways."))
    
    create_nes_scatterplot(
      data = common_parental, x_col = "NES.293t", y_col = "NES.468",
      p_adj_x_col = "p.adjust.293t", p_adj_y_col = "p.adjust.468",
      title = paste("Q1: Common Parental Response to Hcy (", db_type, ")"),
      filename = file.path(output_dir, paste0(fname("Q1_common_parental_response_scatterplot"), ".png"))
    )
  } else {
    print("Skipping Q1/Q2 analysis: parental GSEA result files not found.")
    common_parental <- data.frame()
  }

  # --- Q2: Maintained Response in Each Revertant ---
  print(paste("--- Q2: Maintained response in each revertant (", db_type, ") ---"))
  file_r1 <- file.path(gsea_dir, paste0("r1_hcy_vs_met_GSEA_", db_type, "_results.csv"))
  if (nrow(common_parental) > 0 && file.exists(file_r1)) {
    maintained_in_r1 <- common_parental %>%
      inner_join(read.csv(file_r1), by = "ID") %>%
      filter(sign(NES) == sign(NES.293t)) %>%
      select(ID, Description = Description.x, NES.293t, NES.468, NES) %>% rename(NES.r1 = NES)
    write.csv(maintained_in_r1, file.path(output_dir, paste0(fname("Q2a_maintained_in_R1_response"), ".csv")), row.names = FALSE)
    print(paste("Found", nrow(maintained_in_r1), "pathways maintained in R1."))
    create_nes_heatmap(maintained_in_r1, paste("Q2a: Parental Pathways Maintained in R1 (", db_type, ")"), file.path(output_dir, paste0(fname("Q2a_maintained_in_R1_heatmap"), ".png")))
  }
  
  file_r8 <- file.path(gsea_dir, paste0("r8_hcy_vs_met_GSEA_", db_type, "_results.csv"))
  if (nrow(common_parental) > 0 && file.exists(file_r8)) {
    maintained_in_r8 <- common_parental %>%
      inner_join(read.csv(file_r8), by = "ID") %>%
      filter(sign(NES) == sign(NES.293t)) %>%
      select(ID, Description = Description.x, NES.293t, NES.468, NES) %>% rename(NES.r8 = NES)
    write.csv(maintained_in_r8, file.path(output_dir, paste0(fname("Q2b_maintained_in_R8_response"), ".csv")), row.names = FALSE)
    print(paste("Found", nrow(maintained_in_r8), "pathways maintained in R8."))
    create_nes_heatmap(maintained_in_r8, paste("Q2b: Parental Pathways Maintained in R8 (", db_type, ")"), file.path(output_dir, paste0(fname("Q2b_maintained_in_R8_heatmap"), ".png")))
  }

  # --- Q3: Common Unique Response Signature ---
  print(paste("--- Q3: Common unique response signature (", db_type, ") ---"))
  file1_q3 <- file.path(gsea_dir, paste0("r1_vs_293t_hcy_GSEA_", db_type, "_results.csv"))
  file2_q3 <- file.path(gsea_dir, paste0("r8_vs_468_hcy_GSEA_", db_type, "_results.csv"))
  if (file.exists(file1_q3) && file.exists(file2_q3)) {
    common_unique_response <- inner_join(read.csv(file1_q3), read.csv(file2_q3), by = "ID", suffix = c(".r1_vs_p1", ".r8_vs_p2")) %>%
      filter(sign(NES.r1_vs_p1) == sign(NES.r8_vs_p2)) %>%
      select(ID, Description = Description.r1_vs_p1, NES.r1_vs_p1, p.adjust.r1_vs_p1, NES.r8_vs_p2, p.adjust.r8_vs_p2) %>% arrange(p.adjust.r1_vs_p1)
    write.csv(common_unique_response, file.path(output_dir, paste0(fname("Q3_common_unique_response"), ".csv")), row.names = FALSE)
    print(paste("Found", nrow(common_unique_response), "common unique response pathways."))
    create_nes_scatterplot(
      data = common_unique_response, x_col = "NES.r1_vs_p1", y_col = "NES.r8_vs_p2",
      p_adj_x_col = "p.adjust.r1_vs_p1", p_adj_y_col = "p.adjust.r8_vs_p2",
      title = paste("Q3: Common Unique Signature in Revertants (", db_type, ")"),
      filename = file.path(output_dir, paste0(fname("Q3_common_unique_response_scatterplot"), ".png"))
    )
  } else {
    print("Skipping Q3 revised analysis: one or both input files not found.")
  }

  # --- Q4: Common Baseline Differences ---
  print(paste("--- Q4: Common baseline differences in revertants (", db_type, ") ---"))
  file1_q4 <- file.path(gsea_dir, paste0("r1_vs_293t_met_GSEA_", db_type, "_results.csv"))
  file2_q4 <- file.path(gsea_dir, paste0("r8_vs_468_met_GSEA_", db_type, "_results.csv"))
  if (file.exists(file1_q4) && file.exists(file2_q4)) {
    common_baseline <- inner_join(read.csv(file1_q4), read.csv(file2_q4), by = "ID", suffix = c(".r1", ".r8")) %>%
      filter(sign(NES.r1) == sign(NES.r8)) %>%
      select(ID, Description = Description.r1, NES.r1, p.adjust.r1, NES.r8, p.adjust.r8) %>% arrange(p.adjust.r1)
    write.csv(common_baseline, file.path(output_dir, paste0(fname("Q4_common_baseline_signature"), ".csv")), row.names = FALSE)
    print(paste("Found", nrow(common_baseline), "common baseline pathways."))
    create_nes_scatterplot(
      data = common_baseline, x_col = "NES.r1", y_col = "NES.r8",
      p_adj_x_col = "p.adjust.r1", p_adj_y_col = "p.adjust.r8",
      title = paste("Q4: Common Baseline Signature in Revertants (", db_type, ")"),
      filename = file.path(output_dir, paste0(fname("Q4_common_baseline_signature_scatterplot"), ".png"))
    )
  } else {
    print("Skipping Q4 analysis: one or both input files not found.")
  }
}

# ==============================================================================
# Run the Analysis for Both GO and KEGG
# ==============================================================================
run_interpretation_for_db(db_type = "GO")
run_interpretation_for_db(db_type = "KEGG")

print("--- ALL GSEA Interpretation Scripts Finished ---")