#' @title Full Mechanical Scorecard Analysis
#' @description This script performs an advanced structural analysis relating positional entropy (B-factor) features from different RNA regions to indel efficiency.
#'   It calculates Spearman correlations between 28 structural features (from 4 RNA regions and 7 metrics) and indel frequency,
#'   generates region-specific summary files, and visualizes the relationships using a 4x7 scatter plot grid.
#' @param input_tsv A character string specifying the path to the input TSV file containing entropy data and Indel frequencies.
#'   Defaults to "entropy_advanced_strict.tsv".
#' @param output_pdf A character string specifying the path for the output PDF plot.
#'   Defaults to "entropy_advanced_scorecard_plot.pdf".
#' @return This script does not return a value but generates several output files:
#'   \itemize{
#'     \item Tab-separated text files (e.g., \code{summary_TL_advanced.txt}) for each region, summarizing Spearman correlations.
#'     \item A PDF file (\code{entropy_advanced_scorecard_plot.pdf}) containing a 4x7 grid of scatter plots.
#'   }
#' @details
#'   The script processes structural features derived from positional entropy (B-factor) for four RNA regions (TL, SL1, SL2, SL3)
#'   and seven statistical metrics (Mean, Median, Max, Min, StDev, Range, IQR).
#'   Each metric is associated with a specific 'Mechanism_Role':
#'   \itemize{
#'     \item Max: "Single-Point Weak Link/Hinge"
#'     \item Median: "Baseline Structural Bulk"
#'     \item StDev: "Rigid-to-Flexible Contrast"
#'     \item Mean: "Average Global Looseness"
#'     \item Min: "Rigidity Anchor Point"
#'     \item Range: "Mechanical Stretch/Gradient"
#'     \item IQR: "Consistency of Bulk Flexibility"
#'   }
#'   Correlations are calculated using Spearman's rank correlation coefficient.
#' @importFrom ggplot2 ggplot geom_point geom_smooth facet_grid theme_bw labs theme ggsave element_text
#' @importFrom dplyr %>% select arrange all_of
#' @importFrom tidyr pivot_longer separate drop_na
#' @importFrom stats cor.test
#' @importFrom utils read.delim write.table head
#' @examples
#' \dontrun{
#' # Assuming 'entropy_advanced_strict.tsv' is in the current directory:
#' # source("your_script_name.R")
#' }
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[stats]{cor.test}}, \code{\link[tidyr]{pivot_longer}}
#' @keywords entropy structural_biology RNA Indel correlation plot
#' @export
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Setup Data
input_tsv <- "entropy_advanced_strict.tsv"
output_pdf <- "entropy_advanced_scorecard_plot.pdf"

if (!file.exists(input_tsv)) {
  stop(sprintf("Error: %s not found in the current directory.", input_tsv))
}

# Load and handle NAs
df <- read.delim(input_tsv, sep = "\t", check.names = FALSE, na.strings = c("NA", "NaN", ""))

regions <- c("TL", "SL1", "SL2", "SL3")
# Expanded Metric Set
metrics <- c("Mean", "Median", "Max", "Min", "StDev", "Range", "IQR")

# Define all 28 features (4 regions * 7 metrics)
features <- c()
for (r in regions) {
  for (m in metrics) {
    features <- c(features, paste(r, m, sep = "_"))
  }
}

# 2. Force Data Types & Clean
df$Indel <- as.numeric(as.character(df$Indel))
for (feat in features) {
  if (feat %in% colnames(df)) {
    df[[feat]] <- as.numeric(as.character(df[[feat]]))
  }
}

df_clean <- drop_na(df, all_of(c("Indel"))) # Clean based on Indel, then handle feature NAs individually

# 3. Calculate Correlations and Export Domain Files
all_results <- list()
roles <- c(Max    = "Single-Point Weak Link/Hinge", 
           Median = "Baseline Structural Bulk", 
           StDev  = "Rigid-to-Flexible Contrast", 
           Mean   = "Average Global Looseness",
           Min    = "Rigidity Anchor Point",
           Range  = "Mechanical Stretch/Gradient",
           IQR    = "Consistency of Bulk Flexibility")

for (r in regions) {
  region_list <- list()
  
  for (m in metrics) {
    feat_name <- paste(r, m, sep = "_")
    
    # Check if column exists (safety for partial extraction)
    if (feat_name %in% colnames(df_clean)) {
      stat_test <- suppressWarnings(cor.test(df_clean[[feat_name]], df_clean$Indel, method = "spearman"))
      
      res_row <- data.frame(
        Feature = feat_name,
        Metric = m,
        Rho = round(stat_test$estimate, 4),
        P_value = signif(stat_test$p.value, 4),
        Mechanism_Role = roles[m],
        stringsAsFactors = FALSE
      )
      
      region_list[[m]] <- res_row
      all_results[[feat_name]] <- res_row
    }
  }
  
  # Export Domain-Specific Summary
  domain_df <- do.call(rbind, region_list) %>%
    select(Metric, Rho, P_value, Mechanism_Role) %>%
    arrange(desc(abs(Rho)))
  
  write.table(domain_df, file = paste0("summary_", r, "_advanced.txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("Generated: summary_%s_advanced.txt\n", r))
}

# 4. Global Leaderboard (Top 10 strongest features)
leaderboard <- do.call(rbind, all_results) %>%
  arrange(desc(abs(Rho)))

cat("\n======================================================\n")
cat(sprintf("   STRUCTURAL SCORECARD LEADERBOARD (N = %d)\n", nrow(df_clean)))
cat("======================================================\n")
print(head(leaderboard %>% select(Feature, Rho, P_value), 10), row.names = FALSE)
cat("======================================================\n\n")

# 5. Generate the 4x7 Scatter Plot Grid
df_long <- pivot_longer(df_clean, cols = any_of(features), names_to = "Feature", values_to = "Value") %>%
  separate(Feature, into = c("Region", "Metric"), sep = "_")

df_long$Region <- factor(df_long$Region, levels = regions)
df_long$Metric <- factor(df_long$Metric, levels = metrics)

p <- ggplot(df_long, aes(x = Value, y = Indel)) +
  geom_point(alpha = 0.6, color = "black", size = 1.5) +
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = FALSE) +
  facet_grid(Region ~ Metric, scales = "free_x") +
  theme_bw() +
  labs(title = "Full Mechanical Scorecard: 28 Structural Features vs Indel Efficiency",
       subtitle = "Rows: RNA Regions | Columns: Distributional Metrics",
       x = "Positional Entropy (B-factor)",
       y = "Indel Frequency (%)") +
  theme(strip.text = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14))

# Saving with expanded dimensions for the 7 columns
ggsave(output_pdf, plot = p, width = 16, height = 10, dpi = 300)
cat(sprintf("Saved 4x7 scorecard grid to: %s\n", output_pdf))
