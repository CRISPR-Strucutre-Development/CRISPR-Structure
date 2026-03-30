# Load required libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(dplyr))

# 1. Load the data
df <- read.delim('entropy_averages_strict_final.tsv', sep='\t', check.names=FALSE)
n_current <- nrow(df)

# Define columns to ignore
ignore_cols <- c('ID', 'sgRNA', 'TargetSequence(RNAversion)', 'Indel', 'SSC', 
                 'Spacer.Sequence', 'Scaffold.Sequence', 'Indel_Binary')

# Identify the structural metric columns
structure_cols <- setdiff(colnames(df), ignore_cols)

# Create a dataframe to store summary statistics
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
