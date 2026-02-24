# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Setup Data
input_tsv <- "entropy_cumulative_stats.tsv" 

if (!file.exists(input_tsv)) {
  stop(sprintf("Error: %s not found.", input_tsv))
}

df <- read.delim(input_tsv, sep = "\t", check.names = FALSE, na.strings = c("NA", "NaN", ""))
df$Indel <- as.numeric(as.character(df$Indel))
df$Ent_Max <- as.numeric(as.character(df$Ent_Max))

df_clean <- drop_na(df, all_of(c("Indel", "Ent_Max")))

# 2. Define Thresholds
biological_thresh <- 30
structural_thresh <- 2.940

# Assign Quadrant Categories for coloring
df_clean <- df_clean %>%
  mutate(Classification = case_when(
    Indel >= biological_thresh & Ent_Max < structural_thresh ~ "True Positive (TP)",
    Indel < biological_thresh & Ent_Max >= structural_thresh ~ "True Negative (TN)",
    Indel < biological_thresh & Ent_Max < structural_thresh ~ "False Positive (FP)",
    Indel >= biological_thresh & Ent_Max >= structural_thresh ~ "False Negative (FN)"
  ))

# 3. Generate the Quadrant Plot
p_quad <- ggplot(df_clean, aes(x = Ent_Max, y = Indel, fill = Classification)) +
  geom_point(size = 5, shape = 21, color = "black", alpha = 0.8) +
  
  # The Crosshairs
  geom_vline(xintercept = structural_thresh, linetype = "dashed", color = "darkred", linewidth = 1.2) +
  geom_hline(yintercept = biological_thresh, linetype = "dashed", color = "darkblue", linewidth = 1.2) +
  
  # Define Custom Colors for the Quadrants
  scale_fill_manual(values = c(
    "True Positive (TP)" = "#377EB8",  # Blue (Success)
    "True Negative (TN)" = "#4DAF4A",  # Green (Success)
    "False Positive (FP)" = "#FF7F00", # Orange (Error)
    "False Negative (FN)" = "#E41A1C"  # Red (Error)
  )) +
  
  theme_bw() +
  labs(title = "In Silico Gatekeeper vs. Biological Ground Truth",
       subtitle = sprintf("Crosshairs: Ent_Max Threshold (%.3f) | Indel Threshold (%.0f%%)", 
                          structural_thresh, biological_thresh),
       x = "Maximum Positional Entropy (Ent_Max)",
       y = "Experimental Indel Efficiency (%)") +
  
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

# Save the plot
ggsave("EntMax_Quadrant_Analysis.pdf", plot = p_quad, width = 10, height = 7, dpi = 300)
cat("Saved Quadrant Plot to: EntMax_Quadrant_Analysis.pdf\n")
