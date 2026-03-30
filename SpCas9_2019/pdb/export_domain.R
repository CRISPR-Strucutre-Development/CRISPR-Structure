# Load required library
library(dplyr)
library(tidyr)

# 1. Setup Data
input_tsv <- "entropy_advanced_strict.tsv"

if (!file.exists(input_tsv)) {
  stop(sprintf("Error: %s not found in the current directory.", input_tsv))
}

df <- read.delim(input_tsv, sep = "\t", check.names = FALSE, na.strings = c("NA", "NaN", ""))

regions <- c("SL1", "SL2", "SL3", "TL")
metrics <- c("Mean", "Median", "Max", "StDev")

# 2. Clean Data
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

# 3. Calculate and Export by Region
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
