#' @file
#' @title Power Analysis for CRISPR-Cas9 Entropy Metrics
#' @description This script performs a power analysis to determine the required sample size for various CRISPR-Cas9 entropy metrics, based on observed effect sizes from the `entropy_averages_strict_final.tsv` dataset. It calculates the f2 effect size for each metric in a multiple regression model (predicting 'Indel' from 'SSC' and the metric) and determines the sample size needed for 80% power. Finally, it generates a CSV file with the results and a power curve visualization for the top 5 most promising metrics.

#' @section Load Libraries:
#' Loads necessary R packages for statistical analysis, plotting, and data manipulation.
library(pwr)     # For power analysis calculations (e.g., pwr.f2.test)
library(ggplot2) # For creating visualizations (e.g., power curves)
library(dplyr)   # For data manipulation and piping operations

#' @section 1. Load Data:
#' Loads the primary dataset and identifies columns representing structural metrics for analysis.
#' The dataset 'entropy_averages_strict_final.tsv' contains experimental data, including
#' Indel frequencies, SSC scores, and various entropy/structure-related metrics.
df <- read.delim('entropy_averages_strict_final.tsv', sep='\t', check.names=FALSE)
#' Defines columns to be ignored or treated as outcome/covariates,
#' allowing identification of columns representing structural metrics.
ignore_cols <- c('ID', 'sgRNA', 'TargetSequence(RNAversion)', 'Indel', 'SSC', 
                 'Spacer.Sequence', 'Scaffold.Sequence', 'Indel_Binary')
#' Identifies columns that represent structural metrics by excluding `ignore_cols` from all columns.
structure_cols <- setdiff(colnames(df), ignore_cols)

#' @section 2. Calculate Observed Effect Sizes and Required N for ALL Metrics:
#' Iterates through each structural metric, performs a multiple regression to calculate
#' the f2 effect size, and determines the sample size needed for 80% statistical power.
power_data <- data.frame()
#' Loop through each identified structural column to perform power analysis.
for (col in structure_cols) {
  #' Subset the data to exclude rows with NA values for the current metric.
  sub_df <- df[!is.na(df[[col]]), ]
  
  #' @subsection Standardized Multiple Regression:
  #' Fits a linear model with standardized variables to calculate the effect size (R-squared).
  #' The model predicts 'Indel' frequency based on 'SSC' and the current structural metric.
  fit <- lm(scale(Indel) ~ scale(SSC) + scale(sub_df[[col]]), data=sub_df)
  r2 <- summary(fit)$r.squared
  f2 <- r2 / (1 - r2) 
  
  #' @subsection Calculate Required Sample Size (N) for 80% Power:
  #' Determines the minimum sample size needed to detect the observed effect size with 80% power
  #' and a significance level of 0.05. Uses `tryCatch` to handle cases where `f2` is zero or very small,
  #' which might lead to errors in `pwr.f2.test`.
  required_n <- tryCatch({
    res <- pwr.f2.test(u=2, f2=f2, sig.level=0.05, power=0.80)
    ceiling(res$v + res$u + 1) # Calculate total N from v (degrees of freedom for error) and u (degrees of freedom for numerator)
  }, error = function(e) NA) # Returns NA if an error occurs (e.g., f2 is too small)
  
  #' Appends the calculated metric, its effect size, and required N to the `power_data` dataframe.
  power_data <- rbind(power_data, data.frame(
    Metric = col, 
    Effect_Size_f2 = f2, 
    Required_N_for_80_Power = required_n
  ))
}

#' @section 3. Save Full Power Analysis Results:
#' Saves the comprehensive table of all metrics, their observed effect sizes, and the
#' required sample sizes for 80% power to a CSV file. The table is sorted by the
#' required sample size, indicating the metrics with the strongest signal first.
#' Sorts the `power_data` dataframe by `Required_N_for_80_Power` in ascending order.
power_data <- power_data %>% arrange(Required_N_for_80_Power)
#' Writes the sorted power analysis results to a CSV file.
write.csv(power_data, "power_analysis_full_table.csv", row.names=FALSE)
#' Prints a success message to the console.
cat("SUCCESS: Full results saved to 'power_analysis_full_table.csv'\n")

#' @section 4. Filter for Top 5 Promising Metrics to Plot:
#' Selects the top 5 metrics from the sorted `power_data` table for visualization.
#' These are the metrics that require the smallest sample size for 80% power.
top_metrics <- head(power_data, 5)

#' @section 5. Generate Data for Power Curves:
#' Calculates statistical power across a range of sample sizes for each of the top 5
#' most promising metrics. This data will be used to plot power curves.
#' Defines the range of sample sizes (N) to evaluate for the power curves.
n_range <- seq(10, 250, by=5)
plot_df <- data.frame()

#' Loop through each of the top 5 metrics to calculate power across `n_range`.
for (i in 1:nrow(top_metrics)) {
  m_f2 <- top_metrics$Effect_Size_f2[i] # Observed f2 effect size for the current metric
  m_name <- top_metrics$Metric[i]     # Name of the current metric
  #' Calculates power for each `n` in `n_range` using `pwr.f2.test`.
  powers <- sapply(n_range, function(n) {
    pwr.f2.test(u=2, v=n-3, f2=m_f2, sig.level=0.05)$power
  })
  #' Appends the calculated N, Power, and Metric name to `plot_df`.
  plot_df <- rbind(plot_df, data.frame(N=n_range, Power=powers, Metric=m_name))
}

#' @section 6. Create Power Curve Visualization:
#' Generates a PDF plot visualizing the power curves for the top 5 metrics.
#' The plot shows how statistical power changes with sample size (N) for each metric,
#' highlighting current and target sample sizes and the 80% power threshold.
#' Opens a PDF device for saving the plot.
pdf('power_analysis_metrics.pdf', width=10, height=7)
#' Creates the ggplot visualization:
#' - `geom_line`: Draws power curves for each metric.
#' - `geom_vline`: Adds vertical lines for current (N=19) and target (N=150) sample sizes.
#' - `geom_hline`: Adds a horizontal line for the 80% power threshold.
#' - `labs`: Sets the title, subtitle, and axis labels.
#' - `theme_minimal`: Applies a clean, minimal theme to the plot.
ggplot(plot_df, aes(x=N, y=Power, color=Metric)) +
  geom_line(size=1.2) +
  geom_vline(xintercept=19, linetype="dashed", color="red") +
  geom_vline(xintercept=150, linetype="dashed", color="darkgreen") +
  geom_hline(yintercept=0.8, linetype="dotted", color="black") +
  labs(title="Power Analysis: Probability of Success",
       subtitle="Dashed lines indicate current sample (19) vs target sample (150)",
       y="Statistical Power (Probability)", x="Sample Size (N)") +
  theme_minimal()
#' Closes the PDF device, saving the plot to the specified file.
dev.off()
