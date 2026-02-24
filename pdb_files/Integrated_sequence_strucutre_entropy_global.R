# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ROCR)
library(RColorBrewer)

# 1. Load and Prepare Data
# Ensure your file matches the name: entropy_cumulative_stats_update_final.tsv
df <- read.delim("entropy_cumulative_stats_update_final.tsv", sep="\t", check.names=FALSE)

# Binarize: Active (>30%), Inactive (<30%) 
df <- df %>%
  mutate(
    Activity_State = ifelse(Indel > 30, 1, 0),
    Class_Label = factor(ifelse(Indel > 30, "Active (>30%)", "Inactive (<30%)"),
                         levels = c("Inactive (<30%)", "Active (>30%)"))
  )

# 2. Performance Metrics (Structural Threshold = 2.94)
structural_thresh <- 2.940
df <- df %>%
  mutate(
    Predicted_State = ifelse(Ent_Max < structural_thresh, 1, 0),
    Classification = case_when(
      Activity_State == 1 & Predicted_State == 1 ~ "True Positive (TP)",
      Activity_State == 0 & Predicted_State == 0 ~ "True Negative (TN)",
      Activity_State == 0 & Predicted_State == 1 ~ "False Positive (FP)",
      Activity_State == 1 & Predicted_State == 0 ~ "False Negative (FN)"
    )
  )

# Calculate NPV, PPV, Sensitivity, Specificity
conf_matrix <- table(Predicted = df$Predicted_State, Actual = df$Activity_State)
TN <- conf_matrix["0", "0"]; FN <- conf_matrix["0", "1"]
TP <- conf_matrix["1", "1"]; FP <- conf_matrix["1", "0"]

NPV <- (TN / (TN + FN)) * 100
PPV <- (TP / (TP + FP)) * 100
Sens <- (TP / (TP + FN)) * 100
Spec <- (TN / (TN + FP)) * 100

cat("\n=======================================================\n")
cat("   DIAGNOSTIC PERFORMANCE (Ent_Max < 2.94)\n")
cat("=======================================================\n")
cat(sprintf("NPV (Filter Reliability): %.1f%%\n", NPV))
cat(sprintf("PPV (Selection Reliability): %.1f%%\n", PPV))
cat(sprintf("Sensitivity: %.1f%%\n", Sens))
cat(sprintf("Specificity: %.1f%%\n", Spec))
cat("=======================================================\n")

# 3. AUC Leaderboard
features_to_test <- c("Indel", "SSC", "Ent_Max")
auc_results <- list()
pred_objects <- list()

for (feat in features_to_test) {
  if (feat == "Ent_Max") {
    pred_raw <- prediction(-df[[feat]], df$Activity_State)
  } else {
    pred_raw <- prediction(df[[feat]], df$Activity_State)
  }
  
  raw_auc <- performance(pred_raw, "auc")@y.values[[1]]
  auc_results[[feat]] <- data.frame(Feature = feat, AUC = round(raw_auc, 4))
  pred_objects[[feat]] <- pred_raw
}

leaderboard <- do.call(rbind, auc_results) %>% arrange(desc(AUC))
print(leaderboard)

# 4. Generate Plots
# 4.1 AUC Bar Plot
p_bar <- ggplot(leaderboard, aes(x=reorder(Feature, AUC), y=AUC, fill=Feature == "Indel")) +
  geom_col(color="black", alpha=0.8) +
  coord_flip() +
  scale_fill_manual(values=c("TRUE"="#E41A1C", "FALSE"="#377EB8")) +
  geom_text(aes(label=sprintf("%.3f", AUC)), hjust=-0.2, fontface="bold") +
  theme_bw() +
  labs(title="AUC Comparison: Sequence vs Structure", x="Metric", y="AUC") +
  theme(legend.position="none")
ggsave("AUC_Comparison_Barplot.pdf", plot=p_bar, width=8, height=4)

# 4.2 Annotated Quadrant Plot
p_quad <- ggplot(df, aes(x=Ent_Max, y=Indel, fill=Classification)) +
  geom_point(size=5, shape=21, color="black", alpha=0.8) +
  geom_vline(xintercept=structural_thresh, linetype="dashed", color="darkred") +
  geom_hline(yintercept=30, linetype="dashed", color="darkblue") +
  scale_fill_manual(values=c("True Positive (TP)"="#377EB8", "True Negative (TN)"="#4DAF4A",
                             "False Positive (FP)"="#FF7F00", "False Negative (FN)"="#E41A1C")) +
  theme_bw() +
  labs(title="In Silico Gatekeeper Analysis",
       subtitle=sprintf("NPV: %.1f%% | PPV: %.1f%%", NPV, PPV),
       x="Maximum Positional Entropy (Ent_Max)", y="Indel Efficiency (%)")
ggsave("Quadrant_Analysis_Final.pdf", plot=p_quad, width=8, height=6)

# 4.3 ADDED: The Structural Veto Plot
# Categorize based on SSC median and Structural Threshold
ssc_med <- median(df$SSC, na.rm=TRUE)
df_veto <- df %>%
  mutate(
    SSC_Grp = factor(ifelse(SSC > ssc_med, "High SSC", "Low SSC"), levels=c("High SSC", "Low SSC")),
    Ent_Grp = factor(ifelse(Ent_Max < structural_thresh, "Stable (Low Ent)", "Unstable (High Ent)"),
                     levels=c("Stable (Low Ent)", "Unstable (High Ent)"))
  ) %>%
  group_by(SSC_Grp, Ent_Grp) %>%
  summarise(Mean_Indel = mean(Indel, na.rm=TRUE), N = n(), .groups = 'drop')

p_veto <- ggplot(df_veto, aes(x=SSC_Grp, y=Mean_Indel, fill=Ent_Grp)) +
  geom_col(position=position_dodge(width=0.8), color="black", width=0.7) +
  geom_text(aes(label=sprintf("%.1f%%\n(N=%d)", Mean_Indel, N)), 
            position=position_dodge(width=0.8), vjust=-0.5, fontface="bold", size=3.5) +
  scale_fill_manual(values=c("Stable (Low Ent)"="#377EB8", "Unstable (High Ent)"="#E41A1C")) +
  theme_bw() +
  labs(title="Sequence vs. Structure: The 'Structural Veto' Effect",
       x="Sequence Score Category (SSC)", y="Mean Indel Efficiency (%)", fill="Structural State") +
  ylim(0, 100) +
  theme(plot.title=element_text(face="bold"), legend.position="top")

ggsave("SSC_vs_Entropy_Veto_Plot.pdf", plot=p_veto, width=8, height=6)

# 4.4 SSC Residual Analysis
fit_ssc <- lm(Indel ~ SSC, data=df)
df$SSC_Residuals <- resid(fit_ssc)

p_resid <- ggplot(df, aes(x=Ent_Max, y=SSC_Residuals)) +
  geom_point(aes(color=Class_Label), size=4) +
  geom_smooth(method="lm", color="black", linetype="dashed", se=FALSE) +
  theme_bw() +
  labs(title="SSC Residuals vs Structural Entropy",
       subtitle="Negative residuals = performed WORSE than sequence score predicted",
       x="Ent_Max", y="Residuals (Observed - SSC Predicted)")
ggsave("SSC_Residual_Analysis.pdf", plot=p_resid, width=8, height=6)

# 5. Export Summary 2x2 Table
write.table(df_veto, "Summary_Veto_Stats.tsv", sep="\t", row.names=FALSE)
cat("\nFinal Summary Table saved to 'Summary_Veto_Stats.tsv'\n")
