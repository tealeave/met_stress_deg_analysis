# main.R
#
# Purpose: This script serves as the main pipeline driver for the entire
#          metabolic stress DEG analysis project. It executes all other
#          analysis and visualization scripts in the correct order.
#
# Usage:
#   1. Set the working directory to the project root.
#   2. Run this script from the R console: source("scripts/main.R")
#   3. Or from the command line: Rscript scripts/main.R
#
# Best Practices:
#   - Each step of the pipeline is clearly demarcated.
#   - Informative messages are printed to the console and logged to a file.
#   - The script assumes that all required packages are installed.
#   - Individual scripts are responsible for loading their own libraries.
#   - Error handling is left to the individual scripts; this driver will
#     stop if any sourced script fails.
#
# -----------------------------------------------------------------------------

# --- Preamble: Set up environment and logging ---

# Load the logging library
# If you don't have log4r, install it: install.packages("log4r")
library(log4r)

# Create a directory for logs if it doesn't exist
if (!dir.exists("logs")) {
  dir.create("logs")
}

# Define the log file path with a timestamp
log_file <- paste0("logs/pipeline_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")

# FIX: Define a custom layout function. The function `default_layout` is not
# exported by the log4r package for direct use. Creating our own layout
# is the most robust solution.
default_log_layout <- function(level, ...) {
  # Formats log messages with timestamp, level, and the message content.
  paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " [", level, "] ", paste0(..., collapse = " "))
}

# Create a logger. It will write to the console and to the log file.
logger <- logger(
  threshold = "INFO",
  appenders = list(
    console_appender(layout = default_log_layout),
    file_appender(log_file, layout = default_log_layout)
  )
)

# --- Pipeline Start ---
info(logger, "===================================================================")
info(logger, "Starting the Metabolic Stress DEG Analysis Pipeline")
info(logger, paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
info(logger, "===================================================================\n")


# --- Pipeline Execution ---

# The pipeline is executed by sourcing each script in a logical sequence.
# The order is determined by the dependencies between the scripts.

# Wrap the entire pipeline in a tryCatch block for robust error logging
tryCatch({

  # 1. Core Differential Expression Analysis
  # -----------------------------------------
  info(logger, "--> STEP 1: Running DESeq2 Differential Expression Analysis...")
  source("scripts/deseq_analysis.R")
  info(logger, "... DESeq2 analysis complete.\n")

  # 2. Principal Component Analysis (PCA) for Quality Control
  # ---------------------------------------------------------
  info(logger, "--> STEP 2: Generating PCA plots for QC...")
  source("scripts/pca.R")
  info(logger, "... PCA plots generated.\n")

  # 3. Volcano Plots for Visualizing DE Results
  # -------------------------------------------
  info(logger, "--> STEP 3: Generating static and interactive volcano plots...")
  
  info(logger, "  -> Running volcano.R for static plots...")
  source("scripts/volcano.R")
  
  info(logger, "  -> Running volcano_plotly.R for interactive plots...")
  source("scripts/volcano_plotly.R")
  
  info(logger, "... Volcano plots generated.\n")

  # 4. Set-Based Analysis (Venn and UpSet Plots)
  # --------------------------------------------
  info(logger, "--> STEP 4: Generating Venn diagrams and UpSet plots...")
  source("scripts/venn.R")
  source("scripts/upset.R")
  info(logger, "... Set-based analysis plots generated.\n")

  # 5. Over-Representation Analysis (ORA)
  # -------------------------------------
  info(logger, "--> STEP 5: Running Over-Representation Analysis (ORA)...")
  source("scripts/ora.R")
  info(logger, "... ORA complete.\n")

  # 6. Gene Set Enrichment Analysis (GSEA)
  # --------------------------------------
  info(logger, "--> STEP 6: Running Gene Set Enrichment Analysis (GSEA)...")
  source("scripts/gsea.R")
  info(logger, "... GSEA complete.\n")

  # 7. Interpretation of GSEA Results
  # ---------------------------------
  info(logger, "--> STEP 7: Interpreting and visualizing GSEA results...")
  source("scripts/interpret_gsea.R")
  info(logger, "... GSEA interpretation complete.\n")

}, error = function(e) {
  # Log the error if any script fails
  error(logger, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  error(logger, "An error occurred during pipeline execution.")
  error(logger, paste("Error message:", e$message))
  error(logger, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  stop(e) # Stop the script
})


# --- Pipeline Completion ---
info(logger, "===================================================================")
info(logger, "Pipeline Finished Successfully!")
info(logger, paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
info(logger, "Check the 'results/' and 'reports/' directories for output.")
info(logger, paste("A detailed log has been saved to:", log_file))
info(logger, "===================================================================")
