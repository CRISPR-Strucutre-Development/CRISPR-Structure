# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ROCR)
library(RColorBrewer)
library(lmtest) # Required for the Likelihood Ratio Test

# 1. Load and Prepare Data
df <- read.delim("entropy_cumulative_stats_update_final.tsv", sep="\t", check.names=FALSE)

# Binarize: Active (>30%), Inactive (<30%) 
df <- df %>%
  mutate(
    Activity_State = ifelse(Indel > 30, 1, 0),
    Class_Label = factor(ifelse(Indel > 30, "Active (>30%)", "Inactive (<30%)"),
                         levels = c("Inactive (<30%)", "Active (>30%)"))
  )

# ---------------------------------------------------------
# 2. UNBIASED THRESHOLD CALCULATION (Youden's J Statistic)
# ---------------------------------------------------------
# Since lower Ent_Max predicts Activity (1), we use negative Ent_Max for ROCR
pred_obj <- prediction(-df$Ent_Max, df$Activity_State)
perf_obj <- performance(pred_obj, "tpr", "fpr")

# Extract True Positive Rate (Sensitivity) and False Positive Rate (1 - Specificity)
tpr <- perf_obj@y.values[[1]]
fpr <- perf_obj@x.values[[1]]

# Calculate Youden's J for all possible thresholds
youdens_j <- tpr - fpr

# Find the index of the maximum J value
optimal_idx <- which.max(youdens_j)

# Extract the optimal threshold (reversing the negative sign)
structural_thresh <- -perf_obj@alpha.values[[1]][optimal_idx]

cat("\n=======================================================\n")
cat("   UNBIASED THRESHOLD CALCULATION (Youden's Index)\n")
cat("=======================================================\n")
cat(sprintf("Optimal Ent_Max Threshold: %.3f\n", structural_thresh))
cat(sprintf("Maximized Sensitivity:     %.1f%%\n", tpr[optimal_idx] * 100))
cat(sprintf("Maximized Specificity:     %.1f%%\n", (1 - fpr[optimal_idx]) * 100))
cat("=======================================================\n")

# 3. Diagnostic Metrics using the Objective Threshold
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

conf_matrix <- table(Predicted = df$Predicted_State, Actual = df$Activity_State)
TN <- conf_matrix["0", "0"]; FN <- conf_matrix["0", "1"]
TP <- conf_matrix["1", "1"]; FP <- conf_matrix["1", "0"]

NPV <- (TN / (TN + FN)) * 100
PPV <- (TP / (TP + FP)) * 100

cat("\n=======================================================\n")
cat(sprintf("   DIAGNOSTIC PERFORMANCE (Ent_Max < %.3f)\n", structural_thresh))
cat("=======================================================\n")
cat(sprintf("NPV (Filter Reliability):    %.1f%%\n", NPV))
cat(sprintf("PPV (Selection Reliability): %.1f%%\n", PPV))
cat("=======================================================\n")

# 4. AUC Leaderboard
features_to_test <- c("Indel", "SSC", "Ent_Max", "DG_UNFOLD")
auc_results <- list()
pred_objects <- list()

for (feat in features_to_test) {
  if (!(feat %in% colnames(df))) next
  
  pred_raw <- prediction(df[[feat]], df$Activity_State)
  raw_auc <- performance(pred_raw, "auc")@y.values[[1]]
  
  # Directionality correction: if AUC < 0.5, lower values are better (invert it)
  if (raw_auc < 0.5) {
    adj_auc <- 1 - raw_auc
    pred_plot <- prediction(-df[[feat]], df$Activity_State)
  } else {
    adj_auc <- raw_auc
    pred_plot <- pred_raw
  }
  
  auc_results[[feat]] <- data.frame(Feature = feat, AUC = round(adj_auc, 4))
  pred_objects[[feat]] <- pred_plot
}

leaderboard <- do.call(rbind, auc_results) %>% arrange(desc(AUC))

# ---------------------------------------------------------
# 5. NESTED LIKELIHOOD RATIO TESTS (Benchmarked to SSC)
# ---------------------------------------------------------
library(lmtest)

# Build the models
model_seq    <- glm(Activity_State ~ SSC, data = df, family = binomial)
model_local  <- glm(Activity_State ~ SSC + Ent_Max, data = df, family = binomial)
model_global <- glm(Activity_State ~ SSC + DG_UNFOLD, data = df, family = binomial)
model_all    <- glm(Activity_State ~ SSC + Ent_Max + DG_UNFOLD, data = df, family = binomial)

# Test all structural additions strictly against M1 (Sequence Only)
lrt_local  <- lrtest(model_seq, model_local)
lrt_global <- lrtest(model_seq, model_global)
lrt_all    <- lrtest(model_seq, model_all)

cat("\n=======================================================\n")
cat("   LIKELIHOOD RATIO TEST: COMPARED TO SEQUENCE BASELINE\n")
cat("=======================================================\n")
cat(sprintf("M1: INDEL ~ SSC                           LL = %.3f\n", logLik(model_seq)))
cat(sprintf("M2: INDEL ~ SSC + Ent_Max                 LL = %.3f | p-value vs M1 = %.5f\n", logLik(model_local), lrt_local$`Pr(>Chisq)`[2]))
cat(sprintf("M3: INDEL ~ SSC + DG_UNFOLD               LL = %.3f | p-value vs M1 = %.5f\n", logLik(model_global), lrt_global$`Pr(>Chisq)`[2]))
cat(sprintf("M4: INDEL ~ SSC + Ent_Max + DG_UNFOLD     LL = %.3f | p-value vs M1 = %.5f\n", logLik(model_all), lrt_all$`Pr(>Chisq)`[2]))
cat("=======================================================\n")
# ---------------------------------------------------------
# 6. Generate Plots
# ---------------------------------------------------------

# 6.1 AUC Bar Plot (Revamped & Categorized)
leaderboard <- leaderboard %>%
  mutate(
    Category = case_when(
      Feature == "Indel" ~ "Ground Truth (Perfect)",
      Feature == "SSC" ~ "1D Sequence Predictor",
      Feature %in% c("Ent_Max", "DG_UNFOLD") ~ "3D Biophysical Predictor"
    ),
    Category = factor(Category, levels = c("Ground Truth (Perfect)", "3D Biophysical Predictor", "1D Sequence Predictor"))
  )

p_bar <- ggplot(leaderboard, aes(x=reorder(Feature, AUC), y=AUC, fill=Category)) +
  geom_col(color="black", width=0.7, alpha=0.9) +
  geom_text(aes(label=sprintf("%.3f", AUC)), hjust=-0.2, size=4.5, fontface="bold") +
  coord_flip() +
  scale_fill_manual(values=c("Ground Truth (Perfect)"="#E41A1C", 
                             "3D Biophysical Predictor"="#377EB8", 
                             "1D Sequence Predictor"="#4DAF4A")) +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.2)) + # Expanded limit prevents text cutoff
  theme_minimal() +
  theme(
    axis.text.y = element_text(face="bold", size=12, color="black"),
    axis.text.x = element_text(size=11, color="black"),
    plot.title = element_text(face="bold", size=14),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  labs(title="Predictive Power: Sequence vs. Biophysics",
       subtitle="Area Under the ROC Curve (AUC) for individual features",
       x="Predictive Feature", y="AUC Score")

ggsave("AUC_Comparison_Extended.pdf", plot=p_bar, width=8, height=5)

# 6.2 Annotated Quadrant Plot (Structural Veto)
p_quad <- ggplot(df, aes(x=Ent_Max, y=Indel, fill=Classification)) +
  geom_point(size=5, shape=21, color="black", alpha=0.8) +
  geom_vline(xintercept=structural_thresh, linetype="dashed", color="darkred", size=1) +
  theme_bw() + labs(title="Structural Gatekeeper Analysis", 
                    subtitle=sprintf("Optimal Cutoff Threshold: %.3f", structural_thresh), 
                    x="Ent_Max", y="Indel %")
ggsave("Quadrant_Analysis_Final.pdf", plot=p_quad, width=8, height=6)

# 6.3 The Structural Veto Plot (SSC vs Ent_Max)
ssc_med <- median(df$SSC, na.rm=TRUE)
df_veto <- df %>%
  mutate(
    SSC_Grp = factor(ifelse(SSC > ssc_med, "High SSC", "Low SSC"), levels=c("High SSC", "Low SSC")),
    Ent_Grp = factor(ifelse(Ent_Max < structural_thresh, "Stable (Low Ent)", "Unstable (High Ent)"))
  ) %>%
  group_by(SSC_Grp, Ent_Grp) %>%
  summarise(Mean_Indel = mean(Indel, na.rm=TRUE), N = n(), .groups = 'drop')

p_veto <- ggplot(df_veto, aes(x=SSC_Grp, y=Mean_Indel, fill=Ent_Grp)) +
  geom_col(position=position_dodge(width=0.8), color="black") +
  scale_fill_manual(values=c("Stable (Low Ent)"="#377EB8", "Unstable (High Ent)"="#E41A1C")) +
  theme_bw() + labs(title="The Structural Veto Effect", y="Mean Indel %")
ggsave("SSC_vs_Entropy_Veto_Plot.pdf", plot=p_veto, width=8, height=6)

# 6.4 Thermodynamic Influence Plot (DG_UNFOLD vs Indel)
p_dg <- ggplot(df, aes(x=DG_UNFOLD, y=Indel, color=Class_Label)) +
  geom_point(size=4, alpha=0.8) +
  geom_smooth(method="lm", color="black", linetype="dashed", se=FALSE) +
  scale_color_manual(values = c("Inactive (<30%)" = "#E41A1C", "Active (>30%)" = "#377EB8")) +
  theme_bw() + labs(title="Thermodynamic Influence: Global Energy vs Editing", 
                     x="Global Free Energy: DG_UNFOLD (kcal/mol)", y="Indel %")
ggsave("DG_Unfold_vs_Indel.pdf", plot=p_dg, width=8, height=6)

# Export Final Table
write.table(df_veto, "Summary_Veto_Stats.tsv", sep="\t", row.names=FALSE)