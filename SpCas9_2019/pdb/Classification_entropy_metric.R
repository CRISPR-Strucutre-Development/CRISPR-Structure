# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ROCR)
library(RColorBrewer)

# 1. Setup Data
# Ensure the file 'entropy_cumulative_stats.tsv' contains the DG_UNFOLD column
input_tsv <- "entropy_cumulative_stats.tsv" 

if (!file.exists(input_tsv)) {
  stop(sprintf("Error: %s not found. Please verify the file exists.", input_tsv))
}

df <- read.delim(input_tsv, sep = "\t", check.names = FALSE, na.strings = c("NA", "NaN", ""))

# FEATURE LIST: Now including DG_UNFOLD and Indel
features <- c("Indel", "DG_UNFOLD", "Ent_Sum", "Ent_Median", "Ent_Max", 
              "Ent_Min", "Ent_StDev", "Ent_IQR", "Ent_Range")

# 2. Clean and Binarize the Data
for (feat in features) {
  if (feat %in% colnames(df)) {
    df[[feat]] <- as.numeric(as.character(df[[feat]]))
  }
}

# Remove rows where Indel is missing (Ground Truth is essential)
df_clean <- drop_na(df, all_of(c("Indel")))

# Binarize: Active (>30%), Inactive (<30%)
df_clean <- df_clean %>%
  mutate(
    Activity_State = ifelse(Indel > 30, 1, 0),
    Class_Label = factor(ifelse(Indel > 30, "Active (>30%)", "Inactive (<30%)"), 
                         levels = c("Inactive (<30%)", "Active (>30%)"))
  )

# 3. Calculate ROC-AUC & Mann-Whitney P-Values
auc_results <- list()
pred_objects <- list() 
facet_labels <- c() 

for (feat in features) {
  if (feat %in% colnames(df_clean)) {
    
    # Check if data exists for this feature
    if (all(is.na(df_clean[[feat]]))) next
    
    # ROCR Prediction and AUC
    pred_raw <- prediction(df_clean[[feat]], df_clean$Activity_State)
    raw_auc <- performance(pred_raw, "auc")@y.values[[1]]
    
    # Directionality adjustment: 
    # For DG_UNFOLD and Ent_Max, lower values usually predict better activity
    if (raw_auc < 0.5) {
      adj_auc <- 1 - raw_auc
      pred_plot <- prediction(-df_clean[[feat]], df_clean$Activity_State) 
    } else {
      adj_auc <- raw_auc
      pred_plot <- pred_raw
    }
    
    # Mann-Whitney U Test (Wilcoxon Rank Sum)
    # Tests if the distributions differ significantly between Active and Inactive
    pval <- wilcox.test(df_clean[[feat]] ~ df_clean$Activity_State, exact = FALSE)$p.value
    
    auc_results[[feat]] <- data.frame(
      Feature = feat, 
      AUC = round(adj_auc, 4), 
      MW_P_Value = pval, 
      stringsAsFactors = FALSE
    )
    
    pred_objects[[feat]] <- pred_plot
    
    # Labels for Faceted Boxplots
    if(feat == "Indel") {
       facet_labels[feat] <- sprintf("%s (Ground Truth)\n(AUC=%.2f | p=%.1e)", feat, adj_auc, pval)
    } else {
       facet_labels[feat] <- sprintf("%s\n(AUC=%.2f | p=%.1e)", feat, adj_auc, pval)
    }
  }
}

leaderboard <- do.call(rbind, auc_results) %>% arrange(desc(AUC))

cat("\n=======================================================\n")
cat("   CLASSIFICATION LEADERBOARD (Sequence vs. Structure)\n")
cat("=======================================================\n")
print(leaderboard, row.names = FALSE)
cat("=======================================================\n\n")

# 4. Generate the Comparative ROC Plot
pdf("entropy_ROC_curves.pdf", width = 8, height = 8)
# Colors: Updated to handle 9 potential features
colors <- brewer.pal(n = min(9, length(leaderboard$Feature)), name = "Set1") 
sorted_features <- leaderboard$Feature

plot(performance(pred_objects[[sorted_features[1]]], "tpr", "fpr"), 
     col = colors[1], lwd = 4, main = "ROC Curves: Thermodynamics & Entropy vs. Indel",
     xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray60", lwd = 2) 

for (i in 2:length(sorted_features)) {
  plot(performance(pred_objects[[sorted_features[i]]], "tpr", "fpr"), 
        col = colors[i], lwd = 2, add = TRUE)
}

legend_labels <- sapply(1:length(sorted_features), function(i) {
  sprintf("%s (AUC = %.3f)", sorted_features[i], leaderboard$AUC[i])
})
legend("bottomright", legend = legend_labels, col = colors, lwd = c(4, rep(2, length(sorted_features)-1)), bty = "n", cex = 0.8)
dev.off()

# 5. Generate the AUC Bar Plot
leaderboard$Fill_Color <- ifelse(leaderboard$Feature == "Indel", "Baseline", "Predictor")

p_bar <- ggplot(leaderboard, aes(x = reorder(Feature, AUC), y = AUC, fill = Fill_Color)) +
  geom_col(color = "black", alpha = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("Baseline" = "#E41A1C", "Predictor" = "#377EB8")) +
  geom_text(aes(label = sprintf("%.3f", AUC)), hjust = -0.2, fontface = "bold", size = 4) +
  theme_bw() +
  labs(title = "AUC Comparison: DG_UNFOLD vs. Entropy Metrics",
       x = "Structural Feature",
       y = "Adjusted AUC") +
  ylim(0, 1.2) + 
  theme(legend.position = "none",
        axis.text.y = element_text(face = "bold", size = 11))

ggsave("entropy_AUC_barplot.pdf", plot = p_bar, width = 8, height = 5)

# 6. Generate the Boxplots Grid
df_long <- pivot_longer(df_clean, cols = any_of(features), names_to = "Feature", values_to = "Value")
df_long$Facet_Name <- facet_labels[df_long$Feature]

# Set order based on AUC leaderboard
ordered_facets <- facet_labels[sorted_features]
df_long$Facet_Name <- factor(df_long$Facet_Name, levels = ordered_facets)

p_box <- ggplot(df_long, aes(x = Class_Label, y = Value, fill = Class_Label)) +
  geom_boxplot(width = 0.5, color = "black", outlier.shape = 21, outlier.fill = "white", alpha = 0.8) +
  facet_wrap(~ Facet_Name, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Inactive (<30%)" = "#E41A1C", "Active (>30%)" = "#377EB8")) +
  theme_bw() +
  labs(title = "Distributions: Thermodynamic Stability (DG) vs. Entropy Predictors",
       x = "sgRNA Efficacy Category",
       y = "Metric Value") +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(face = "bold", size = 9))

ggsave("entropy_boxplots_all.pdf", plot = p_box, width = 12, height = 10)