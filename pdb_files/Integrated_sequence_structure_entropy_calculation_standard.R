#' @title Analyze Biophysical Predictors of Protein Activity
#' @description This script performs a comprehensive analysis of various biophysical features
#'   (e.g., Ent_Max, DG_UNFOLD) and sequence features (e.g., SSC) to predict protein activity,
#'   which is defined by an Indel percentage threshold.
#' @details The analysis pipeline encompasses several key steps:
#'   \itemize{
#'     \item Data loading and initial binarization of activity states based on Indel percentage.
#'     \item Calculation of an unbiased optimal threshold for 'Ent_Max' using Youden's J statistic
#'           to maximize diagnostic performance.
#'     \item Evaluation of diagnostic metrics such as Negative Predictive Value (NPV) and
#'           Positive Predictive Value (PPV) based on the derived threshold.
#'     \item Benchmarking of individual feature predictive power through Area Under the ROC Curve (AUC) analysis.
#'     \item Application of Nested Likelihood Ratio Tests (LRTs) to quantitatively assess the
#'           incremental predictive value of biophysical features when added to a sequence-based model.
#'     \item Generation of a suite of publication-quality plots to visually represent the findings,
#'           including AUC comparisons, a structural gatekeeper quadrant plot, an illustration of
#'           the "structural veto" effect, and the influence of global free energy on activity.
#'     \item Export of summary statistics tables derived from the analysis.
#'   }
#'   This script is designed for a specific dataset (`entropy_cumulative_stats_update_final.tsv`)
#'   and assumes its presence in the working directory. It does not contain user-defined functions
#'   intended for external export, but rather executes a complete analytical workflow.
#' @author (Your Name Here) <your.email@example.com>
#' @date 2023-10-27
#' @keywords data-analysis biophysics protein-activity machine-learning ROC AUC LRT
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import ROCR
#' @import RColorBrewer
#' @import lmtest
# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ROCR)
library(RColorBrewer)
library(lmtest) # Required for the Likelihood Ratio Test

#' @section 1. Load and Prepare Data
#' @description This section loads the raw experimental data from a TSV file and performs
#'   initial data transformations.
#' @details The 'Indel' column is used to binarize protein activity into 'Active (>30%)'
#'   and 'Inactive (<30%)' states, creating 'Activity_State' (0/1) and 'Class_Label' factors.
df <- read.delim("entropy_cumulative_stats_update_final.tsv", sep="\t", check.names=FALSE)

# Binarize: Active (>30%), Inactive (<30%)
df <- df %>%
  mutate(
    Activity_State = ifelse(Indel > 30, 1, 0),
    Class_Label = factor(ifelse(Indel > 30, "Active (>30%)", "Inactive (<30%)"),
                         levels = c("Inactive (<30%)", "Active (>30%)"))
  )

# ---------------------------------------------------------
#' @section 2. Unbiased Threshold Calculation (Youden's J Statistic)
#' @description This section determines an optimal, unbiased threshold for the 'Ent_Max'
#'   feature to differentiate between active and inactive protein states.
#' @details Youden's J statistic (Sensitivity + Specificity - 1) is used to find the
#'   threshold that maximizes diagnostic accuracy. The 'Ent_Max' values are negated
#'   for compatibility with the ROCR package, as lower 'Ent_Max' values are indicative of activity.
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

#' @section 3. Diagnostic Metrics using the Objective Threshold
#' @description This section applies the previously calculated 'structural_thresh'
#'   to classify proteins and then computes key diagnostic performance metrics.
#' @details A confusion matrix is generated, and from it, Negative Predictive Value (NPV)
#'   and Positive Predictive Value (PPV) are calculated to assess the reliability of
#'   the filter and selection processes, respectively.
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

#' @section 4. AUC Leaderboard
#' @description This section evaluates the individual predictive power of several features
#'   by calculating their Area Under the ROC Curve (AUC).
#' @details Each feature's AUC is computed, and if an AUC is less than 0.5, it is
#'   inverted (1-AUC) to ensure that higher AUC always indicates better predictive performance,
#'   regardless of the inherent directionality of the feature. The results are compiled
#'   into a leaderboard sorted by AUC.
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
#' @section 5. Nested Likelihood Ratio Tests (Benchmarked to SSC)
#' @description This section performs nested Likelihood Ratio Tests (LRTs) to statistically
#'   assess the added predictive value of biophysical features ('Ent_Max', 'DG_UNFOLD')
#'   when introduced to a baseline logistic regression model that only includes 'SSC'.
#' @details Logistic regression models are built with increasing complexity. LRTs are then
#'   used to compare the more complex models (M2, M3, M4) against the simpler baseline model (M1),
#'   providing p-values that indicate significant improvements in model fit.
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
#' @section 6. Generate Plots
#' @description This section creates several analytical plots to visualize the results
#'   of the predictive modeling and feature analysis.
#' @details Four distinct plots are generated and saved as PDF files:
#'   \itemize{
#'     \item \strong{AUC Bar Plot}: A categorized bar chart comparing the AUC scores of different features.
#'     \item \strong{Annotated Quadrant Plot}: A scatter plot showing 'Ent_Max' vs 'Indel',
#'           with quadrants defined by the optimal 'Ent_Max' threshold.
#'     \item \strong{Structural Veto Plot}: A bar chart illustrating the effect of 'Ent_Max'
#'           (structural stability) on 'Indel' across different 'SSC' groups.
#'     \item \strong{Thermodynamic Influence Plot}: A scatter plot of 'DG_UNFOLD' vs 'Indel',
#'           highlighting the relationship between global free energy and protein activity.
#'   }
# ---------------------------------------------------------

# 6.1 AUC Bar Plot (Revamped & Categorized)
#' @subsection 6.1 AUC Bar Plot (Revamped & Categorized)
#' @description Generates a bar plot visualizing the AUC scores of various features,
#'   categorized by their type (Ground Truth, 3D Biophysical, 1D Sequence Predictor).
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
#' @subsection 6.2 Annotated Quadrant Plot (Structural Veto)
#' @description Generates a scatter plot of 'Ent_Max' vs 'Indel', illustrating the
#'   classification of proteins into quadrants based on activity state and the
#'   calculated optimal 'Ent_Max' threshold.
p_quad <- ggplot(df, aes(x=Ent_Max, y=Indel, fill=Classification)) +
  geom_point(size=5, shape=21, color="black", alpha=0.8) +
  geom_vline(xintercept=structural_thresh, linetype="dashed", color="darkred", size=1) +
  theme_bw() + labs(title="Structural Gatekeeper Analysis", 
                    subtitle=sprintf("Optimal Cutoff Threshold: %.3f", structural_thresh), 
                    x="Ent_Max", y="Indel %")
ggsave("Quadrant_Analysis_Final.pdf", plot=p_quad, width=8, height=6)

# 6.3 The Structural Veto Plot (SSC vs Ent_Max)
#' @subsection 6.3 The Structural Veto Plot (SSC vs Ent_Max)
#' @description Creates a bar plot demonstrating the "structural veto" effect,
#'   showing how structural stability (Ent_Max) influences protein activity (Indel)
#'   within different sequence complexity (SSC) groups.
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
#' @subsection 6.4 Thermodynamic Influence Plot (DG_UNFOLD vs Indel)
#' @description Generates a scatter plot exploring the relationship between global free
#'   energy of unfolding (DG_UNFOLD) and protein activity (Indel %), colored by activity class.
p_dg <- ggplot(df, aes(x=DG_UNFOLD, y=Indel, color=Class_Label)) +
  geom_point(size=4, alpha=0.8) +
  geom_smooth(method="lm", color="black", linetype="dashed", se=FALSE) +
  scale_color_manual(values = c("Inactive (<30%)" = "#E41A1C", "Active (>30%)" = "#377EB8")) +
  theme_bw() + labs(title="Thermodynamic Influence: Global Energy vs Editing", 
                     x="Global Free Energy: DG_UNFOLD (kcal/mol)", y="Indel %")
ggsave("DG_Unfold_vs_Indel.pdf", plot=p_dg, width=8, height=6)

# Export Final Table
#' @section 7. Export Summary Table
#' @description Exports the 'df_veto' dataframe, which contains summarized statistics
#'   related to the "structural veto" analysis, into a TSV file.
write.table(df_veto, "Summary_Veto_Stats.tsv", sep="\t", row.names=FALSE)
