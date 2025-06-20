# scripts/interpret_gsea_V10.R
# Purpose: A complete script to analyze and interpret GSEA results for GO and KEGG.
# --- VERSION 10: Added diverging bar charts for top pathways in Q1-Q4. ---
# --- VERSION 9: Updated scatterplot to resize to 8x8 if only one facet is present. ---
# --- VERSION 8: Fixed pheatmap error when clustering a single row. ---

# --- 1. Load Libraries ---
library(dplyr)
library(ggplot2)
library(here)
library(ggrepel)
library(pheatmap)
library(stringr) # Added for text truncation in plots

# --- 2. Configuration ---
gsea_dir <- here("results", "gsea_analysis")
output_dir <- here("results", "final_interpretation")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. Helper Functions ---

# MODIFIED: Function updated to increase point and label font size.
# --- AND to dynamically change plot size if only one facet is present. ---
create_nes_scatterplot <- function(data, x_col, y_col, p_adj_x_col, p_adj_y_col, title, filename, top_n_labels = 20) {
  if (nrow(data) == 0) {
    print(paste("Skipping plot", basename(filename), "- no data."))
    return()
  }

  # --- CLASSIFY PATHWAYS FOR FACETING ---
  data <- data %>%
    mutate(regulation_group = case_when(
      !!sym(x_col) > 0 & !!sym(y_col) > 0 ~ "Commonly Upregulated",
      !!sym(x_col) < 0 & !!sym(y_col) < 0 ~ "Commonly Downregulated",
      TRUE ~ "Other"
    )) %>%
    filter(regulation_group != "Other")

  # If, after filtering, no data remains, exit.
  if (nrow(data) == 0) {
    print(paste("Skipping plot", basename(filename), "- no commonly regulated data."))
    return()
  }

  # --- LABELING LOGIC (Adapted for facets) ---
  top_pathways_to_label <- data %>%
    group_by(regulation_group) %>%
    mutate(rank_metric = !!sym(p_adj_x_col) * !!sym(p_adj_y_col)) %>%
    slice_min(order_by = rank_metric, n = top_n_labels / 2)

  data <- data %>%
    mutate(label = if_else(ID %in% top_pathways_to_label$ID, Description, ""))

  p <- ggplot(data, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(aes(color = !!sym(p_adj_y_col)), alpha = 0.5, size = 2.4) +
    facet_wrap(~ regulation_group, scales = "free") +
    geom_text_repel(
          aes(label = label),
          size = 4,
          box.padding = 1.0,
          point.padding = 0.2,
          min.segment.length = -1,
          force = 5,
        ) +
    scale_color_gradient(low = "red", high = "blue", name = "Adj. P-value") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = title,
      subtitle = "Pathways commonly up-regulated (right) and down-regulated (left)",
      x = paste("NES:", x_col),
      y = paste("NES:", y_col)
    ) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill="grey90", color = "grey90"),
          strip.text = element_text(face = "bold"))

  # --- Dynamically set plot dimensions ---
  num_groups <- length(unique(data$regulation_group))
  plot_width <- if (num_groups == 1) 8 else 14
  plot_height <- 8
  
  ggsave(filename, p, width = plot_width, height = plot_height)
  print(paste("Saved plot:", filename, "with dimensions", plot_width, "x", plot_height))
}

create_nes_heatmap <- function(data, title, filename) {
  if (nrow(data) == 0) {
    print(paste("Skipping heatmap", basename(filename), "- no data."))
    return()
  }
  
  heatmap_matrix <- data %>% select(starts_with("NES")) %>% as.matrix()
  rownames(heatmap_matrix) <- data$Description
  
  should_cluster_rows <- if (nrow(heatmap_matrix) > 1) TRUE else FALSE
  
  palette_breaks <- seq(-max(abs(heatmap_matrix), na.rm = TRUE), max(abs(heatmap_matrix), na.rm = TRUE), length.out = 101)
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  
  pheatmap(
    heatmap_matrix,
    main = title,
    filename = filename,
    width = 8,
    height = if (should_cluster_rows) 10 else 4,
    color = color_palette,
    breaks = palette_breaks,
    cluster_cols = FALSE,
    cluster_rows = should_cluster_rows,
    fontsize_row = 8
  )
  print(paste("Saved heatmap:", filename))
}


# --- NEW FUNCTION: Diverging Barchart ---
create_diverging_barchart <- function(data, nes_col, p_adj_col, title, filename, top_n = 20) {
  if (nrow(data) == 0) {
    print(paste("Skipping barchart", basename(filename), "- no data."))
    return()
  }

  # --- 1. Select top N pathways (N/2 up, N/2 down) ---
  data_up <- data %>%
    filter(!!sym(nes_col) > 0) %>%
    slice_min(order_by = !!sym(p_adj_col), n = top_n / 2)

  data_down <- data %>%
    filter(!!sym(nes_col) < 0) %>%
    slice_min(order_by = !!sym(p_adj_col), n = top_n / 2)

  data_to_plot <- bind_rows(data_up, data_down)

  if (nrow(data_to_plot) == 0) {
    print(paste("Skipping barchart", basename(filename), "- no significant pathways to plot."))
    return()
  }

  # --- 2. Add helper columns for plotting ---
  data_to_plot <- data_to_plot %>%
    mutate(
      regulation = ifelse(!!sym(nes_col) > 0, "Upregulated", "Downregulated"),
      # Truncate long descriptions for better plotting
      short_desc = str_trunc(Description, 70)
    ) %>%
    # Order by NES for the plot
    arrange(!!sym(nes_col)) %>%
    # Make Description a factor to preserve plot order
    mutate(short_desc = factor(short_desc, levels = .$short_desc))

  # --- 3. Create the plot ---
  p <- ggplot(data_to_plot, aes(x = !!sym(nes_col), y = short_desc, fill = regulation)) +
    geom_col(alpha = 0.9) +
    # Place text labels on the opposite side of the bars
    geom_text(
      aes(x = 0, label = short_desc),
      hjust = ifelse(data_to_plot[[nes_col]] > 0, 1, 0),
      nudge_x = ifelse(data_to_plot[[nes_col]] > 0, -0.05, 0.05), # Nudge away from center
      size = 3.5
    ) +
    scale_fill_manual(values = c("Upregulated" = "#e41a1c", "Downregulated" = "#377eb8"), name = "Regulation") +
    labs(
      title = title,
      subtitle = paste("Top", nrow(data_to_plot), "most significant pathways"),
      x = "Normalized Enrichment Score (NES)",
      y = "" # Remove y-axis title
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_blank(), # Hide y-axis text since it's now inside the plot
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray30")

  # --- 4. Save the plot ---
  # Dynamic height based on number of rows
  plot_height <- max(6, nrow(data_to_plot) * 0.35)
  ggsave(filename, p, width = 11, height = plot_height, limitsize = FALSE)
  print(paste("Saved barchart:", filename))
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

    # --- Call new barchart function for Q1 ---
    common_parental_for_bar <- common_parental %>%
        mutate(NES.avg = (NES.293t + NES.468) / 2,
               p.adjust.avg = (p.adjust.293t + p.adjust.468) / 2)

    create_diverging_barchart(
        data = common_parental_for_bar, nes_col = "NES.avg", p_adj_col = "p.adjust.avg",
        title = paste("Q1: Top Common Parental Pathways (", db_type, ")"),
        filename = file.path(output_dir, paste0(fname("Q1_common_parental_response_barchart"), ".png"))
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
      select(ID, Description = Description.x, NES.293t, NES.468, NES.r1 = NES, p.adjust.r1 = p.adjust)
    
    write.csv(maintained_in_r1, file.path(output_dir, paste0(fname("Q2a_maintained_in_R1_response"), ".csv")), row.names = FALSE)
    print(paste("Found", nrow(maintained_in_r1), "pathways maintained in R1."))
    create_nes_heatmap(maintained_in_r1, paste("Q2a: Parental Pathways Maintained in R1 (", db_type, ")"), file.path(output_dir, paste0(fname("Q2a_maintained_in_R1_heatmap"), ".png")))
    
    # --- Call new barchart function for Q2a ---
    create_diverging_barchart(
        data = maintained_in_r1, nes_col = "NES.r1", p_adj_col = "p.adjust.r1",
        title = paste("Q2a: Top Maintained Pathways in R1 (", db_type, ")"),
        filename = file.path(output_dir, paste0(fname("Q2a_maintained_in_R1_barchart"), ".png"))
    )
  }
  
  file_r8 <- file.path(gsea_dir, paste0("r8_hcy_vs_met_GSEA_", db_type, "_results.csv"))
  if (nrow(common_parental) > 0 && file.exists(file_r8)) {
    maintained_in_r8 <- common_parental %>%
      inner_join(read.csv(file_r8), by = "ID") %>%
      filter(sign(NES) == sign(NES.293t)) %>%
      select(ID, Description = Description.x, NES.293t, NES.468, NES.r8 = NES, p.adjust.r8 = p.adjust)
      
    write.csv(maintained_in_r8, file.path(output_dir, paste0(fname("Q2b_maintained_in_R8_response"), ".csv")), row.names = FALSE)
    print(paste("Found", nrow(maintained_in_r8), "pathways maintained in R8."))
    create_nes_heatmap(maintained_in_r8, paste("Q2b: Parental Pathways Maintained in R8 (", db_type, ")"), file.path(output_dir, paste0(fname("Q2b_maintained_in_R8_heatmap"), ".png")))

    # --- Call new barchart function for Q2b ---
    create_diverging_barchart(
        data = maintained_in_r8, nes_col = "NES.r8", p_adj_col = "p.adjust.r8",
        title = paste("Q2b: Top Maintained Pathways in R8 (", db_type, ")"),
        filename = file.path(output_dir, paste0(fname("Q2b_maintained_in_R8_barchart"), ".png"))
    )
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

    # --- Call new barchart function for Q3 ---
    common_unique_for_bar <- common_unique_response %>%
        mutate(NES.avg = (NES.r1_vs_p1 + NES.r8_vs_p2) / 2,
               p.adjust.avg = (p.adjust.r1_vs_p1 + p.adjust.r8_vs_p2) / 2)

    create_diverging_barchart(
        data = common_unique_for_bar, nes_col = "NES.avg", p_adj_col = "p.adjust.avg",
        title = paste("Q3: Top Common Unique Pathways in Revertants (", db_type, ")"),
        filename = file.path(output_dir, paste0(fname("Q3_common_unique_response_barchart"), ".png"))
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

    # --- Call new barchart function for Q4 ---
    common_baseline_for_bar <- common_baseline %>%
        mutate(NES.avg = (NES.r1 + NES.r8) / 2,
               p.adjust.avg = (p.adjust.r1 + p.adjust.r8) / 2)

    create_diverging_barchart(
        data = common_baseline_for_bar, nes_col = "NES.avg", p_adj_col = "p.adjust.avg",
        title = paste("Q4: Top Common Baseline Pathways in Revertants (", db_type, ")"),
        filename = file.path(output_dir, paste0(fname("Q4_common_baseline_barchart"), ".png"))
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
