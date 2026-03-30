#' @file entropy_correlation_analysis.R
#' @title Analyze Entropy-Indel Correlations by Region
#' @description This script performs a correlation analysis between entropy metrics (Mean, Median, Max, StDev)
#' across different regions (SL1, SL2, SL3, TL) and Indel scores.
#' It reads an input TSV file, cleans the data, calculates Spearman correlations for
#' each region and metric, and exports the results to separate TSV files.
#' @author [Author Name Here]
#' @date 2023-10-27
#' @keywords correlation, entropy, indel, Spearman, data analysis, R
#'

#' @section Libraries:
#' Loads necessary R packages for data manipulation and analysis.
library(dplyr)
library(tidyr)

#' @section 1. Setup Data:
#' Initializes input file path, loads the dataset, and defines regions and metrics
#' for subsequent analysis. Includes error handling for missing input file.
input_tsv <- "entropy_advanced_strict.tsv"

if (!file.exists(input_tsv)) {
  stop(sprintf("Error: %s not found in the current directory.", input_tsv))
}

df <- read.delim(input_tsv, sep = "\t", check.names = FALSE, na.strings = c("NA", "NaN", ""))

regions <- c("SL1", "SL2", "SL3", "TL")
metrics <- c("Mean", "Median", "Max", "StDev")

#' @section 2. Clean Data:
#' Processes the raw dataframe by converting relevant columns (Indel and all entropy features)
#' to numeric type and removing rows with missing values (NA) in these critical columns.
df$Indel <- as.numeric(as.character(df$Indel))
# Calculate all features for processing
all_features <- c()
for (r in regions) {
  for (m in metrics) {
    all_features <- c(all_features, paste(r, m, sep = "_"))
    df[[paste(r, m, sep = "_")]] <- as.numeric(as.character(df[[paste(r, m, sep = "_")]]))
  }
}
df_clean <- drop_na(df, all_of(c(all_features, "Indel")))

#' @section 3. Calculate and Export by Region:
#' Iterates through each defined region (SL1, SL2, SL3, TL) and for each region,
#' calculates the Spearman correlation between each entropy metric (Mean, Median, Max, StDev)
#' and the 'Indel' score. Results include Rho values and a descriptive 'Mechanism_Role'
#' based on the metric type. Outputs are combined, ranked by absolute Rho, and exported
#' to region-specific TSV files.
for (r in regions) {
  region_results <- list()

  for (m in metrics) {
    feat_name <- paste(r, m, sep = "_")

    # Calculate correlation
    stat_test <- suppressWarnings(cor.test(df_clean[[feat_name]], df_clean$Indel, method = "spearman"))

    # Store with a descriptive explanation for your discussion
    region_results[[m]] <- data.frame(
      Metric = m,
      Rho = round(stat_test$estimate, 4),
      Mechanism_Role = case_when(
        m == "Max"    ~ "Single-Point Weak Link/Hinge",
        m == "Median" ~ "Baseline Structural Bulk",
        m == "StDev"  ~ "Rigid-to-Flexible Contrast",
        m == "Mean"   ~ "Average Global Looseness"
      ),
      stringsAsFactors = FALSE
    )
  }

  # Combine results for this specific region
  results_df <- do.call(rbind, region_results)

  # Rank by absolute Rho strength
  results_df <- results_df[order(-abs(results_df$Rho)), ]

  # Generate filename and export
  output_filename <- paste0("summary_", r, "_correlation.txt")
  write.table(results_df, file = output_filename, sep = "\t", row.names = FALSE, quote = FALSE)

  cat(sprintf("Generated Summary for %s: %s\n", r, output_filename))
}
