# Load libraries
library(pwr)
library(ggplot2)
library(dplyr)

# 1. Load Data
df <- read.delim('entropy_averages_strict_final.tsv', sep='\t', check.names=FALSE)
ignore_cols <- c('ID', 'sgRNA', 'TargetSequence(RNAversion)', 'Indel', 'SSC', 
                 'Spacer.Sequence', 'Scaffold.Sequence', 'Indel_Binary')
structure_cols <- setdiff(colnames(df), ignore_cols)

# 2. Calculate observed effect sizes and Required N for ALL metrics
power_data <- data.frame()
for (col in structure_cols) {
  sub_df <- df[!is.na(df[[col]]), ]
  
  # Standardized Multiple Regression
  fit <- lm(scale(Indel) ~ scale(SSC) + scale(sub_df[[col]]), data=sub_df)
  r2 <- summary(fit)$r.squared
  f2 <- r2 / (1 - r2) 
  
  # Calculate exactly how many samples (N) are needed for 80% Power
  # we use tryCatch to handle cases where f2 is effectively 0
  required_n <- tryCatch({
    res <- pwr.f2.test(u=2, f2=f2, sig.level=0.05, power=0.80)
    ceiling(res$v + res$u + 1)
  }, error = function(e) NA)
  
  power_data <- rbind(power_data, data.frame(
    Metric = col, 
    Effect_Size_f2 = f2, 
    Required_N_for_80_Power = required_n
  ))
}

# 3. SAVE THE FULL CSV (This is what was missing)
# Sort by best signal (lowest Required N)
power_data <- power_data %>% arrange(Required_N_for_80_Power)
write.csv(power_data, "power_analysis_full_table.csv", row.names=FALSE)
cat("SUCCESS: Full results saved to 'power_analysis_full_table.csv'\n")

# 4. Filter for Top 5 Promising Metrics to Plot
top_metrics <- head(power_data, 5)

# 5. Generate Power Curves (Visualizing the top ones)
n_range <- seq(10, 250, by=5)
plot_df <- data.frame()

for (i in 1:nrow(top_metrics)) {
  m_f2 <- top_metrics$Effect_Size_f2[i]
  m_name <- top_metrics$Metric[i]
  powers <- sapply(n_range, function(n) {
    pwr.f2.test(u=2, v=n-3, f2=m_f2, sig.level=0.05)$power
  })
  plot_df <- rbind(plot_df, data.frame(N=n_range, Power=powers, Metric=m_name))
}

# 6. Create the Visualization
pdf('power_analysis_metrics.pdf', width=10, height=7)
ggplot(plot_df, aes(x=N, y=Power, color=Metric)) +
  geom_line(size=1.2) +
  geom_vline(xintercept=19, linetype="dashed", color="red") +
  geom_vline(xintercept=150, linetype="dashed", color="darkgreen") +
  geom_hline(yintercept=0.8, linetype="dotted", color="black") +
  labs(title="Power Analysis: Probability of Success",
       subtitle="Dashed lines indicate current sample (19) vs target sample (150)",
       y="Statistical Power (Probability)", x="Sample Size (N)") +
  theme_minimal()
dev.off()