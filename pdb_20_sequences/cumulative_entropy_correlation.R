# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Setup Data
# Ensure this matches the output file from the Python merge step
input_tsv <- "entropy_cumulative_stats.tsv" 
output_pdf <- "entropy_cumulative_scorecard_plot.pdf"

if (!file.exists(input_tsv)) {
  stop(sprintf("Error: %s not found in the current directory.", input_tsv))
}

# Load and handle NAs
df <- read.delim(input_tsv, sep = "\t", check.names = FALSE, na.strings = c("NA", "NaN", ""))

# Explicitly define the 7 target features (excluding Residue_Count and Indel)
features <- c("Ent_Sum", "Ent_Median", "Ent_Max", "Ent_Min", "Ent_StDev", "Ent_IQR", "Ent_Range")

# 2. Force Data Types & Clean
df$Indel <- as.numeric(as.character(df$Indel))
for (feat in features) {
  if (feat %in% colnames(df)) {
    df[[feat]] <- as.numeric(as.character(df[[feat]]))
  }
}

# Drop rows missing Indel data
df_clean <- drop_na(df, all_of(c("Indel"))) 

# 3. Calculate Correlations
# Mapping mechanical roles directly to the cumulative metrics
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
cat("\n======================================================\n")
cat(sprintf("   CUMULATIVE SCORECARD LEADERBOARD (N = %d)\n", nrow(df_clean)))
cat("======================================================\n")
print(domain_df %>% select(Feature, Rho, P_value), row.names = FALSE)
cat("======================================================\n\n")

# 5. Generate the Scatter Plot (using facet_wrap for a 1D list of features)
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
