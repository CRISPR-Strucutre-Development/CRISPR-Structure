library(dplyr)
library(lmtest)

# 1. Load Data
df <- read.delim("entropy_cumulative_stats_update_final.tsv", sep="\t", check.names=FALSE)

# Binarize Activity
df$Activity_State <- ifelse(df$Indel > 30, 1, 0)

# 2. Build the Nested Model Hierarchy
# Model 1: Sequence Only (The baseline)
model_1 <- glm(Activity_State ~ SSC, data = df, family = binomial)

# Model 2: Sequence + Local Structure (Local flexibility veto)
model_2 <- glm(Activity_State ~ SSC + Ent_Max, data = df, family = binomial)

# Model 3: Sequence + Local Structure + Global Thermodynamics (Full Model)
model_3 <- glm(Activity_State ~ SSC + Ent_Max + DG_UNFOLD, data = df, family = binomial)

# Model 4: Sequence + Global Thermodynamics (Direct comparison)
model_4 <- glm(Activity_State ~ SSC + DG_UNFOLD, data = df, family = binomial)

# 3. Perform Likelihood Ratio Tests
# Test A: Does Local Entropy (Ent_Max) improve the Sequence model?
lrt_ent <- lrtest(model_1, model_2)

# Test B: Does Global Stability (DG_UNFOLD) improve the Sequence + Entropy model?
lrt_dg_added <- lrtest(model_2, model_3)

# Test C: Does Global Stability (DG_UNFOLD) improve the Sequence model alone?
lrt_dg_only <- lrtest(model_1, model_4)

# 4. Print Results
cat("\n=======================================================\n")
cat("   NESTED LOGISTIC LIKELIHOOD RATIO ANALYSIS\n")
cat("=======================================================\n")
cat(sprintf("1. SSC Baseline LL:          %.3f\n", logLik(model_1)))
cat(sprintf("2. SSC + Ent_Max LL:         %.3f (p = %.5f)\n", logLik(model_2), lrt_ent$`Pr(>Chisq)`[2]))
cat(sprintf("3. SSC + Ent + DG_UNFOLD LL: %.3f (p = %.5f)\n", logLik(model_3), lrt_dg_added$`Pr(>Chisq)`[2]))
cat(sprintf("4. SSC + DG_UNFOLD LL:       %.3f (p = %.5f)\n", logLik(model_4), lrt_dg_only$`Pr(>Chisq)`[2]))
cat("=======================================================\n")
