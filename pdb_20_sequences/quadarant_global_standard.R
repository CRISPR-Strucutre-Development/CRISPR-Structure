#' @title Analyze and Plot Indel Efficiency vs. Positional Entropy
#' @description
#' This script loads experimental indel efficiency and maximum positional entropy data,
#' cleans it, applies biological and structural thresholds to classify data points
#' into True Positive, True Negative, False Positive, and False Negative categories,
#' and then generates a quadrant plot to visualize these classifications.
#'
#' The plot provides an in silico gatekeeper vs. biological ground truth comparison.
#'
#' Assumes the input TSV file 'entropy_cumulative_stats.tsv' is in the same directory.
#' If the file is not found, the script will terminate with an error.
#'
#' The script performs the following steps:
#' 1.  **Data Setup**: Loads data from a TSV, converts relevant columns to numeric,
#'     and removes rows with missing 'Indel' or 'Ent_Max' values.
#' 2.  **Define Thresholds**: Sets `biological_thresh` and `structural_thresh`
#'     for classification.
#' 3.  **Classification**: Assigns each data point to one of four quadrants
#'     (True Positive, True Negative, False Positive, False Negative) based on
#'     its 'Indel' and 'Ent_Max' values relative to the defined thresholds.
#' 4.  **Generate Quadrant Plot**: Creates a scatter plot using ggplot2,
#'     with points colored according to their classification.
#'     Crosshairs indicate the defined thresholds. The plot is saved as a PDF.
#' @importFrom readr read_tsv cols
#' @importFrom dplyr mutate filter case_when
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline labs scale_fill_manual theme_bw
#' @importFrom ggplot2 theme element_text element_rect element_blank ggsave
library(readr)
library(dplyr)
library(ggplot2)

# 1. Setup Data
input_tsv <- "entropy_cumulative_stats.tsv"

if (!file.exists(input_tsv)) {
  stop(paste0("Error: ", input_tsv, " not found."))
}

# Using readr::read_tsv for robust TSV reading, handling NA values
df <- read_tsv(input_tsv, na = c("NA", "NaN", ""), col_types = cols())

# Convert relevant columns to numeric, coercing errors to NaN (handled by as.numeric implicitly)
df <- df %>%
  mutate(
    Indel = as.numeric(Indel),
    Ent_Max = as.numeric(Ent_Max)
  )

# Drop rows where 'Indel' or 'Ent_Max' are NaN (NA in R)
df_clean <- df %>%
  filter(!is.na(Indel) & !is.na(Ent_Max))

# 2. Define Thresholds
biological_thresh <- 30
structural_thresh <- 2.940

# Assign Quadrant Categories for coloring using dplyr::case_when
df_clean <- df_clean %>%
  mutate(
    Classification = case_when(
      (Indel >= biological_thresh) & (Ent_Max < structural_thresh) ~ "True Positive (TP)",
      (Indel < biological_thresh) & (Ent_Max >= structural_thresh) ~ "True Negative (TN)",
      (Indel < biological_thresh) & (Ent_Max < structural_thresh) ~ "False Positive (FP)",
      (Indel >= biological_thresh) & (Ent_Max >= structural_thresh) ~ "False Negative (FN)",
      TRUE ~ "Unknown" # Default for any unhandled cases
    )
  )

# Ensure Classification is a factor for consistent plotting order
df_clean$Classification <- factor(df_clean$Classification, levels = c(
  "True Positive (TP)",
  "True Negative (TN)",
  "False Positive (FP)",
  "False Negative (FN)"
))

# 3. Generate the Quadrant Plot

# Define Custom Colors for the Quadrants
colors <- c(
  "True Positive (TP)" = "#377EB8",  # Blue (Success)
  "True Negative (TN)" = "#4DAF4A",  # Green (Success)
  "False Positive (FP)" = "#FF7F00", # Orange (Error)
  "False Negative (FN)" = "#E41A1C"  # Red (Error)
)

# Create the plot using ggplot2
p <- ggplot(df_clean, aes(x = Ent_Max, y = Indel, fill = Classification)) + # Use fill for Classification
  geom_point(shape = 21, color = "black", size = 5, stroke = 1, alpha = 0.8, # shape=21 for filled circle with border
             show.legend = TRUE) + # Ensure points contribute to legend
  # The Crosshairs
  geom_vline(xintercept = structural_thresh, linetype = "dashed", color = "darkred", size = 1.2) +
  geom_hline(yintercept = biological_thresh, linetype = "dashed", color = "darkblue", size = 1.2) +
  # Set custom fill colors
  scale_fill_manual(values = colors) +
  # Set plot title and labels
  labs(
    title = "In Silico Gatekeeper vs. Biological Ground Truth",
    subtitle = paste0("Crosshairs: Ent_Max Threshold (", sprintf("%.3f", structural_thresh),
                      ") | Indel Threshold (", sprintf("%.0f", biological_thresh), "%)"),
    x = "Maximum Positional Entropy (Ent_Max)",
    y = "Experimental Indel Efficiency (%)",
    fill = "" # Hide legend title for fill aesthetic
  ) +
  # Set theme for better aesthetics, mimicking the Python example
  theme_bw() + # Starting with a black and white theme
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5), # Centered title
    plot.subtitle = element_text(size = 11, hjust = 0.5), # Centered subtitle
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_blank(), # Remove legend title explicitly
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(linetype = "dotted", color = "grey", size = 0.5), # Lighter grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1) # Add a border around the plot area
  )

# Adjust plot layout (done automatically by ggplot normally, but can be controlled with ggsave params)
# plt.tight_layout(rect=[0, 0.03, 1, 0.9]) in Python is largely handled by ggplot2's defaults and theme settings.

# Save the plot
output_filename <- "EntMax_Quadrant_Analysis.pdf"
ggsave(output_filename, plot = p, width = 10, height = 7, dpi = 300)
message(paste0("Saved Quadrant Plot to: ", output_filename))
