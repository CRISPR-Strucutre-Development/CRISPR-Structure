library(ggplot2)
library(dplyr)

# 1. Load Data
df <- read.delim("entropy_cumulative_stats_update_final.tsv", sep="\t", check.names=FALSE)

# Binarize for the Logistic Model
df$Activity_State <- ifelse(df$Indel > 30, 1, 0)
df$Class_Label <- factor(ifelse(df$Indel > 30, "Active (>30%)", "Inactive (<30%)"),
                         levels = c("Inactive (<30%)", "Active (>30%)"))

# 2. Fit the Logistic "Switch" Model
model_logit <- glm(Activity_State ~ Ent_Max, data = df, family = binomial)

# Generate a smooth curve for the S-shape
range_ent <- seq(min(df$Ent_Max, na.rm=TRUE), max(df$Ent_Max, na.rm=TRUE), length.out = 100)
predict_df <- data.frame(Ent_Max = range_ent)
predict_df$Prob <- predict(model_logit, newdata = predict_df, type = "response")

# 3. Create the Plot
p_switch <- ggplot(df, aes(x = Ent_Max, y = Activity_State)) +
  # Actual experimental points (0 and 1)
  # jitter adds a tiny bit of vertical space so points don't overlap
  geom_jitter(aes(color = Class_Label), width = 0, height = 0.04, size = 4, alpha = 0.8) +
  # The Logistic "Switch" Curve
  geom_line(data = predict_df, aes(x = Ent_Max, y = Prob), size = 1.5, color = "black") +
  # The 2.94 Veto Threshold
  geom_vline(xintercept = 2.94, linetype = "dashed", color = "darkred", size = 1) +
  # Labeling and Theme
  scale_color_manual(values = c("Inactive (<30%)" = "#E41A1C", "Active (>30%)" = "#377EB8")) +
  theme_bw() +
  labs(title = "The Structural Veto Switch",
       subtitle = "Probability of guide activity crashes as structural entropy exceeds 2.94",
       x = "Structural Entropy (Ent_Max)",
       y = "Probability of Editing Success (0 to 1)") +
  annotate("text", x = 3.1, y = 0.5, label = "VETO ZONE", color = "darkred", fontface = "bold", angle = 90) +
  theme(legend.position = "top", plot.title = element_text(face = "bold", size = 14))

ggsave("Logistic_Switch_Plot.pdf", plot = p_switch, width = 8, height = 6)
