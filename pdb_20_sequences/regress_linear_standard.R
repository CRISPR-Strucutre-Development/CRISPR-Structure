#' @title Analyze the relationship between structural metrics and Indel frequency, controlling for SSC.
#' @description This script performs partial regression analysis to assess the independent contribution of various structural metrics to Indel frequency, while strictly regressing out the influence of SSC (Short Self-Complementarity).
#' @details
#' The script first loads data from 'entropy_averages_strict_final.tsv'.
#' For each structural metric, it performs a multiple linear regression where Indel is predicted by the structural metric and SSC.
#' It then uses the Frisch-Waugh-Lovell theorem to visually represent the partial regression, plotting residuals of Indel (after regressing out SSC) against residuals of the structural metric (after regressing out SSC).
#' The coefficients and p-values from the full multiple regression are displayed on the partial regression plots.
#' Finally, it generates a section on sample size justification, comparing the current sample size to a proposed n=150 to reduce standard error.
#' All plots and a summary table of coefficients are saved to a PDF report and a CSV file, respectively.
#' @author [Your Name/Organization Here]
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[gridExtra]{grid.arrange}}, \code{\link[dplyr]{dplyr}}, \code{\link[stats]{lm}}, \code{\link[base]{scale}}
#' @keywords analysis regression visualization report
#' @examples
#' # This is a script, so no direct function calls.
#' # To run, ensure 'entropy_averages_strict_final.tsv' is in the working directory.
#' # source("your_script_name.R")
NULL
# Load required libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(dplyr))

# 1. Load the data
#' @name df
#' @title Input Data: Entropy Averages and Experimental Metrics
#' @description A data frame loaded from 'entropy_averages_strict_final.tsv' containing various structural metrics, Indel frequencies, and SSC values.
#' @format A data frame with N rows and P columns, where N is the number of observations and P is the number of variables, including:
#' \describe{
#'   \item{ID}{Unique identifier.}
#'   \item{sgRNA}{sgRNA sequence.}
#'   \item{TargetSequence(RNAversion)}{RNA target sequence.}
#'   \item{Indel}{Indel frequency (dependent variable).}
#'   \item{SSC}{Short Self-Complementarity (covariate).}
#'   \item{Spacer.Sequence}{Spacer sequence.}
#'   \item{Scaffold.Sequence}{Scaffold sequence.}
#'   \item{Indel_Binary}{Binary indel outcome.}
#'   \item{[Structural Metrics]}{Various columns representing structural metrics (independent variables), e.g., entropy values, minimum free energy, etc.}
#' }
#' @source entropy_averages_strict_final.tsv
NULL
df <- read.delim('entropy_averages_strict_final.tsv', sep='\t', check.names=FALSE)
n_current <- nrow(df)

# Define columns to ignore
ignore_cols <- c('ID', 'sgRNA', 'TargetSequence(RNAversion)', 'Indel', 'SSC', 
                 'Spacer.Sequence', 'Scaffold.Sequence', 'Indel_Binary')

# Identify the structural metric columns
structure_cols <- setdiff(colnames(df), ignore_cols)

# Create a dataframe to store summary statistics
#' @name results_df
#' @title Summary Statistics of Partial Regression Coefficients
#' @description A data frame storing the estimated coefficients and p-values for each structural metric from the partial regression analysis.
#' @format A data frame with M rows and 3 columns, where M is the number of structural metrics analyzed:
#' \describe{
#'   \item{Metric}{The name of the structural metric analyzed.}
#'   \item{Coefficient}{The estimated regression coefficient for the metric (Indel_scaled ~ Metric_scaled, with SSC_scaled regressed out).}
#'   \item{P_Value}{The p-value associated with the coefficient.}
#' }
NULL
results_df <- data.frame(Metric=character(), Coefficient=numeric(), P_Value=numeric(), stringsAsFactors=FALSE)

cat(sprintf("Generating Partial Regression Plots for %d metrics...\n", length(structure_cols)))

# Initialize the PDF device
pdf('linear_ssc_regressed_out_report_R.pdf', width=12, height=6)

# 2. Iterate through each structural metric
for (col in structure_cols) {
  
  # Isolate data and drop NAs for this specific metric
  sub_df <- df[!is.na(df[[col]]), c("Indel", "SSC", col)]
  
  # Standardize (scale) the variables to get comparable coefficients
  sub_df$Indel_scaled <- scale(sub_df$Indel)
  sub_df$SSC_scaled <- scale(sub_df$SSC)
  sub_df$Metric_scaled <- scale(sub_df[[col]])
  
  # --- Full Multiple Regression (for precise coefficients/p-values) ---
  formula_full <- as.formula(paste("Indel_scaled ~ SSC_scaled + Metric_scaled"))
  model_full <- lm(formula_full, data=sub_df)
  
  coef_val <- coef(summary(model_full))["Metric_scaled", "Estimate"]
  pval_val <- coef(summary(model_full))["Metric_scaled", "Pr(>|t|)"]
  
  results_df <- rbind(results_df, data.frame(Metric=col, Coefficient=coef_val, P_Value=pval_val))
  
  # --- Partial Regression (Frisch-Waugh-Lovell Theorem) ---
  mod_y_ssc <- lm(Indel_scaled ~ SSC_scaled, data=sub_df)
  resid_y <- residuals(mod_y_ssc)
  
  mod_x_ssc <- lm(Metric_scaled ~ SSC_scaled, data=sub_df)
  resid_x <- residuals(mod_x_ssc)
  
  # Step 3: Plot the residuals against each other
  plot_data <- data.frame(Resid_X = resid_x, Resid_Y = resid_y)
  
  p1 <- ggplot(plot_data, aes(x=Resid_X, y=Resid_Y)) +
    geom_point(color="blue", alpha=0.6, size=2) +
    geom_smooth(method="lm", color="red", se=FALSE) +
    geom_hline(yintercept=0, linetype="dashed", color="black", alpha=0.3) +
    geom_vline(xintercept=0, linetype="dashed", color="black", alpha=0.3) +
    theme_bw() +
    labs(title=sprintf("Partial Regression: Indel vs %s\n(SSC strictly regressed out)", col),
         subtitle=sprintf("Coefficient: %.3f | P-Value: %.4f", coef_val, pval_val),
         x=sprintf("%s Residuals (Isolated from SSC)", col),
         y="Indel Residuals (Isolated from SSC)") +
    theme(plot.title = element_text(face="bold"))
  
  grid.arrange(p1, ncol=1)
}

# 3. Sort results and execute outputs
results_df <- results_df[order(results_df$P_Value), ]

# --- CSV EXPORT ADDED HERE ---
write.csv(results_df, "linear_coefficients_results.csv", row.names=FALSE)
cat("\nSaved numerical results to 'linear_coefficients_results.csv'\n")

best_metric <- results_df$Metric[1]

# --- Section 2: Sample Size Justification (The N=150 Argument) ---
p2 <- ggplot(df, aes_string(x=paste0("`", best_metric, "`"), y="Indel")) +
  geom_point(color="blue", alpha=0.6, size=2) +
  geom_smooth(method="lm", color="red", se=TRUE, fill="grey70") +
  theme_bw() +
  labs(title=sprintf("Systemic Uncertainty at n=%d", n_current),
       subtitle="95% CI bands obscure the signal",
       x=sprintf("Raw %s", best_metric),
       y="Indel %") +
  theme(plot.title = element_text(face="bold"))

n_range <- seq(10, 200, by=5)
se_decay <- 1 / sqrt(n_range)
se_data <- data.frame(n = n_range, SE = se_decay)

se_current <- 1 / sqrt(n_current)
se_150 <- 1 / sqrt(150)
reduction <- ((se_current - se_150) / se_current) * 100

p3 <- ggplot(se_data, aes(x=n, y=SE)) +
  geom_line(size=1) +
  geom_vline(xintercept=n_current, color="red", linetype="dashed", size=1) +
  geom_vline(xintercept=150, color="green4", linetype="dashed", size=1) +
  annotate("text", x=n_current+5, y=max(se_decay)*0.9, label=sprintf("Current: n=%d", n_current), color="red", hjust=0) +
  annotate("text", x=155, y=max(se_decay)*0.8, label="Proposed: n=150", color="green4", hjust=0) +
  theme_bw() +
  labs(title=sprintf("Standard Error Decay (~%.0f%% Drop at n=150)", reduction),
       x="Sample Size (n)",
       y="Relative Standard Error Multiplier") +
  theme(plot.title = element_text(face="bold"))

grid.arrange(p2, p3, ncol=2)
dev.off()

cat("Report complete. Saved as 'linear_ssc_regressed_out_report_R.pdf'\n")
