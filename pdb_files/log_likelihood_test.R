# 1. Setup Data
df <- read.delim("entropy_cumulative_stats_update_final.tsv", sep="\t")

# Binarize the target (The "All-or-Nothing" outcome)
df$Activity_State <- ifelse(df$Indel > 30, 1, 0)

# 2. Build Logistic Models (family = binomial)
# Model 1: Probability of Activity based ONLY on Sequence
model_seq <- glm(Activity_State ~ SSC, data = df, family = binomial)

# Model 2: Probability of Activity based on Sequence + Structure
model_combined <- glm(Activity_State ~ SSC + Ent_Max, data = df, family = binomial)

# 3. Likelihood Ratio Test (LRT)
# This calculates the exact p-value for adding structure to the sequence model
library(lmtest)
lrt_result <- lrtest(model_seq, model_combined)

# 4. Results
cat("\n=======================================================\n")
cat("   LOGISTIC LIKELIHOOD RATIO TEST (Bimodal Data)\n")
cat("=======================================================\n")
cat(sprintf("Log-Likelihood (SSC Only):    %.3f\n", logLik(model_seq)))
cat(sprintf("Log-Likelihood (SSC + Ent):   %.3f\n", logLik(model_combined)))
cat(sprintf("Exact P-Value for Structure:  %.5f\n", lrt_result$`Pr(>Chisq)`[2]))
cat("=======================================================\n")
