#' @file entropy_cumulative_scorecard.R
#' @title Analyze Cumulative Entropy Metrics and Indel Efficiency
#' @description This script performs a comprehensive analysis of cumulative entropy metrics
#'   against Indel efficiency. It calculates Spearman correlations, generates a leaderboard
#'   of features by correlation strength, and produces a multi-panel scatter plot
#'   visualizing these relationships.
#'
#' @details The script processes an input TSV file containing various sequence-level
#'   entropy statistics and 'Indel' efficiency data. Key outputs include a summary
#'   table of correlation results and a PDF plot. This analysis is designed to
#'   characterize the mechanical properties of proteins based on overall sequence entropy.
#'
#' @author [Your Name/Organization Here]
#'
#' @section Dependencies:
#'   Requires the following R packages: `ggplot2`, `dplyr`, and `tidyr`.
#'
#' @section Input Data:
#'   `entropy_cumulative_stats.tsv`: A tab-separated file expected to contain columns
#'   for various entropy metrics (e.g., Ent_Sum, Ent_Median) and 'Indel' efficiency.
#'   This file is assumed to be an output from a preceding Python merge step.
#'
#' @section Output Files:
#'   `summary_cumulative_advanced.txt`: A tab-separated file summarizing the Spearman
#'     correlation results, including Rho and P-value for each entropy feature,
#'     along with its proposed mechanical role. The results are ordered by the
#'     absolute value of Rho.
#'   `entropy_cumulative_scorecard_plot.pdf`: A PDF document containing a multi-panel
#'     scatter plot. Each panel displays the relationship between a single entropy
#'     feature and 'Indel' efficiency, including a linear regression line.
#'
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Setup Data
#' @section 1. Setup Data:
#'   This section defines the paths for input and output files. It performs a check
#'   to ensure the input data file exists before proceeding. The primary dataset is
#'   loaded, and initial handling of NA values (including explicit "NA", "NaN", "")
#'   is performed during the read operation.
input_tsv <- "entropy_cumulative_stats.tsv"
output_pdf <- "entropy_cumulative_scorecard_plot.pdf"

if (!file.exists(input_tsv)) {
  stop(sprintf("Error: %s not found in the current directory.", input_tsv))
}

# Load and handle NAs
df <- read.delim(input_tsv, sep = "\t", check.names = FALSE, na.strings = c("NA", "NaN", ""))

#' @description `features`: A character vector explicitly listing the 7 target entropy
#'   features that will be analyzed. These are cumulative metrics calculated across sequences.
features <- c("Ent_Sum", "Ent_Median", "Ent_Max", "Ent_Min", "Ent_StDev", "Ent_IQR", "Ent_Range")

# 2. Force Data Types & Clean
#' @section 2. Force Data Types & Clean:
#'   This section ensures that the 'Indel' column and all specified entropy `features`
#'   are coerced into numeric data types. It then cleans the dataset by dropping
#'   any rows that have missing values in the 'Indel' column, preparing the data
#'   for correlation analysis.
df$Indel <- as.numeric(as.character(df$Indel))
for (feat in features) {
  if (feat %in% colnames(df)) {
    df[[feat]] <- as.numeric(as.character(df[[feat]]))
  }
}

# Drop rows missing Indel data
df_clean <- drop_na(df, all_of(c("Indel")))

# 3. Calculate Correlations
#' @section 3. Calculate Correlations:
#'   This section performs Spearman correlation tests between each entropy feature
#'   and 'Indel' efficiency. A conceptual mechanical role is assigned to each feature.
#'   The correlation coefficients (Rho) and p-values are calculated, compiled into
#'   a data frame, and then exported to a summary file.
#'
#' @description `roles`: A named character vector that maps each entropy feature
#'   to a descriptive mechanical role, providing contextual understanding of its
#'   potential significance.
roles <- c(Ent_Max    = "Single-Point Weak Link/Hinge",
           Ent_Median = "Baseline Structural Bulk",
           Ent_StDev  = "Rigid-to-Flexible Contrast",
           Ent_Sum    = "Total Sequence Looseness",
           Ent_Min    = "Rigidity Anchor Point",
           Ent_Range  = "Mechanical Stretch/Gradient",
           Ent_IQR    = "Consistency of Bulk Flexibility")

results_list <- list()

for (feat in features) {
  if (feat %in% colnames(df_clean)) {
    # Calculate Spearman correlation
    stat_test <- suppressWarnings(cor.test(df_clean[[feat]], df_clean$Indel, method = "spearman"))

    res_row <- data.frame(
      Feature = feat,
      Rho = round(stat_test$estimate, 4),
      P_value = signif(stat_test$p.value, 4),
      Mechanism_Role = roles[feat],
      stringsAsFactors = FALSE
    )

    results_list[[feat]] <- res_row
  }
}

# Export Global Summary File
domain_df <- do.call(rbind, results_list) %>%
  arrange(desc(abs(Rho)))

write.table(domain_df, file = "summary_cumulative_advanced.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Generated: summary_cumulative_advanced.txt\n")

# 4. Global Leaderboard output to console
#' @section 4. Global Leaderboard Output to Console:
#'   This section prints a neatly formatted leaderboard of the correlation results
#'   to the R console. It displays the 'Feature', 'Rho' (Spearman's correlation
#'   coefficient), and 'P_value' for each analyzed entropy metric, ordered by
#'   the absolute value of Rho, providing a quick summary of the most influential features.
cat("\n======================================================\n")
cat(sprintf("   CUMULATIVE SCORECARD LEADERBOARD (N = %d)\n", nrow(df_clean)))
cat("======================================================\n")
print(domain_df %>% select(Feature, Rho, P_value), row.names = FALSE)
cat("======================================================\n\n")

# 5. Generate the Scatter Plot (using facet_wrap for a 1D list of features)
#' @section 5. Generate Scatter Plot:
#'   This final section prepares the cleaned data (`df_clean`) by transforming it
#'   into a "long" format, which is ideal for creating multi-panel plots using
#'   `ggplot2::facet_wrap`. It then generates a scatter plot grid, where each panel
#'   visualizes the relationship between one entropy feature and 'Indel' efficiency.
#'   A linear regression line is overlaid on each scatter plot to highlight trends.
#'   The factor levels for `Feature` are explicitly set to ensure a consistent
#'   plotting order. The resulting plot is saved as a high-resolution PDF.
df_long <- pivot_longer(df_clean, cols = any_of(features), names_to = "Feature", values_to = "Value")

# Lock the factor levels so the plot outputs in the exact order of the 'features' vector
df_long$Feature <- factor(df_long$Feature, levels = features)

p <- ggplot(df_long, aes(x = Value, y = Indel)) +
  geom_point(alpha = 0.6, color = "black", size = 1.5) +
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = FALSE) +
  facet_wrap(~ Feature, scales = "free_x", ncol = 4) +
  theme_bw() +
  labs(title = "Cumulative Mechanical Scorecard: Global Entropy vs Indel Efficiency",
       subtitle = "Whole-sequence structural distributions",
       x = "Positional Entropy Metric",
       y = "Indel Frequency (%)") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14))

# Saving with optimized dimensions for a 7-panel grid (4 columns x 2 rows)
ggsave(output_pdf, plot = p, width = 12, height = 8, dpi = 300)
cat(sprintf("Saved plot grid to: %s\n", output_pdf))
