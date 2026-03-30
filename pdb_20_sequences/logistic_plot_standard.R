#' @title Analyze Structural Entropy and Editing Success
#' @description This script loads entropy data, binarizes an activity state,
#'   fits a logistic regression model to predict editing success based on structural entropy,
#'   and generates a ggplot visualization of the "Structural Veto Switch".
#' @details The script processes a TSV file containing entropy and indel statistics.
#'   It defines an 'Activity_State' (binary: 0 or 1) and 'Class_Label' (factor)
#'   based on a 30% indel threshold. A logistic regression model (`model_logit`)
#'   is then fitted to predict 'Activity_State' using 'Ent_Max'.
#'   Finally, a `ggplot` object (`p_switch`) is created to visualize this relationship,
#'   highlighting a "Veto Threshold" at Ent_Max = 2.94. The plot is saved as a PDF.
#' @return This script implicitly saves a PDF file named "Logistic_Switch_Plot.pdf"
#'   to the working directory, containing the generated visualization.
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[stats]{glm}}, \code{\link[base]{read.delim}}
#' @import ggplot2
#' @import dplyr
library(ggplot2)
library(dplyr)

# 1. Load Data
df <- read.delim("entropy_cumulative_stats_update_final.tsv", sep="\t", check.names=FALSE)

# Binarize for the Logistic Model
df$Activity_State <- ifelse(df$Indel > 30, 1, 0)
df$Class_Label <- factor(ifelse(df$Indel > 30, "Active (>30%)", "Inactive (<30%)"),
                         levels = c("Inactive (<30%)", "Active (>30%)"))

# 2. Fit the Logistic "Switch" Model
#' @title Logistic Regression Model for Structural Veto Switch
#' @name model_logit
#' @description A generalized linear model (GLM) fitted to predict 'Activity_State'
#'   (binary: 0 or 1) based on 'Ent_Max' (Structural Entropy), using a binomial family.
#' @details This model is crucial for understanding the "structural veto switch"
#'   mechanism, as it quantifies the probability of a guide being 'Active'
#'   (Indel > 30%) as a function of its maximum structural entropy. The fitted
#'   model is subsequently used to generate the smooth S-shaped curve in the visualization.
#' @examples
#' # After running the script, inspect the model summary:
#' # summary(model_logit)
#' # Extract model coefficients:
#' # coef(model_logit)
model_logit <- glm(Activity_State ~ Ent_Max, data = df, family = binomial)

# Generate a smooth curve for the S-shape
range_ent <- seq(min(df$Ent_Max, na.rm=TRUE), max(df$Ent_Max, na.rm=TRUE), length.out = 100)
predict_df <- data.frame(Ent_Max = range_ent)
predict_df$Prob <- predict(model_logit, newdata = predict_df, type = "response")

# 3. Create the Plot
#' @title Structural Veto Switch Visualization
#' @name p_switch
#' @description A `ggplot2` object visualizing the relationship between
#'   Structural Entropy (Ent_Max) and the Probability of Editing Success,
#'   highlighting a "Structural Veto Switch".
#' @details This plot displays several key components:
#'   \itemize{
#'     \item Jittered experimental points representing the observed 'Activity_State'
#'           (0 or 1) for each `Ent_Max` value, colored by 'Class_Label'.
#'     \item A smooth logistic "switch" curve, derived from `model_logit`,
#'           showing the predicted probability of editing success.
#'     \item A vertical dashed line at `Ent_Max` = 2.94, marking the "Veto Threshold".
#'     \item Annotations and labels clarifying the "Veto Zone" and overall theme.
#'   }
#' @return A `ggplot` object that can be printed, further modified, or saved.
#' @seealso \code{\link{model_logit}}, \code{\link[ggplot2]{geom_jitter}}, \code{\link[ggplot2]{geom_line}}, \code{\link[ggplot2]{geom_vline}}
#' @examples
#' # After generation by the script, the plot can be displayed:
#' # print(p_switch)
#' # The script also saves this plot to "Logistic_Switch_Plot.pdf".
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
