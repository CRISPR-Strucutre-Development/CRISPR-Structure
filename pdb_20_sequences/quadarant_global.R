"""
This script loads experimental indel efficiency and maximum positional entropy data,
cleans it, applies biological and structural thresholds to classify data points
into True Positive, True Negative, False Positive, and False Negative categories,
and then generates a quadrant plot to visualize these classifications.

The plot provides an in silico gatekeeper vs. biological ground truth comparison.

Assumes the input TSV file 'entropy_cumulative_stats.tsv' is in the same directory.
If the file is not found, the script will terminate with an error.

The script performs the following steps:
1.  **Data Setup**: Loads data from a TSV, converts relevant columns to numeric,
    and removes rows with missing 'Indel' or 'Ent_Max' values.
2.  **Define Thresholds**: Sets `biological_thresh` and `structural_thresh`
    for classification.
3.  **Classification**: Assigns each data point to one of four quadrants
    (True Positive, True Negative, False Positive, False Negative) based on
    its 'Indel' and 'Ent_Max' values relative to the defined thresholds.
4.  **Generate Quadrant Plot**: Creates a scatter plot using Matplotlib
    and Seaborn, with points colored according to their classification.
    Crosshairs indicate the defined thresholds. The plot is saved as a PDF.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 1. Setup Data
input_tsv = "entropy_cumulative_stats.tsv" 

if not os.path.exists(input_tsv):
  raise FileNotFoundError(f"Error: {input_tsv} not found.")

df = pd.read_csv(input_tsv, sep="\t", na_values=["NA", "NaN", ""])

# Convert relevant columns to numeric, coercing errors to NaN
df['Indel'] = pd.to_numeric(df['Indel'], errors='coerce')
df['Ent_Max'] = pd.to_numeric(df['Ent_Max'], errors='coerce')

# Drop rows where 'Indel' or 'Ent_Max' are NaN
df_clean = df.dropna(subset=['Indel', 'Ent_Max']).copy() # Using .copy() to avoid SettingWithCopyWarning

# 2. Define Thresholds
biological_thresh = 30
structural_thresh = 2.940

# Assign Quadrant Categories for coloring using numpy.select (similar to R's case_when)
conditions = [
    (df_clean['Indel'] >= biological_thresh) & (df_clean['Ent_Max'] < structural_thresh),
    (df_clean['Indel'] < biological_thresh) & (df_clean['Ent_Max'] >= structural_thresh),
    (df_clean['Indel'] < biological_thresh) & (df_clean['Ent_Max'] < structural_thresh),
    (df_clean['Indel'] >= biological_thresh) & (df_clean['Ent_Max'] >= structural_thresh)
]
choices = [
    "True Positive (TP)",
    "True Negative (TN)",
    "False Positive (FP)",
    "False Negative (FN)"
]
df_clean['Classification'] = np.select(conditions, choices, default='Unknown')

# 3. Generate the Quadrant Plot

# Define Custom Colors for the Quadrants
colors = {
    "True Positive (TP)": "#377EB8",  # Blue (Success)
    "True Negative (TN)": "#4DAF4A",  # Green (Success)
    "False Positive (FP)": "#FF7F00", # Orange (Error)
    "False Negative (FN)": "#E41A1C"  # Red (Error)
}

# Set Seaborn style for better aesthetics, mimicking ggplot's theme_bw
sns.set_style("whitegrid")

plt.figure(figsize=(10, 7))
ax = sns.scatterplot(
    data=df_clean,
    x="Ent_Max",
    y="Indel",
    hue="Classification",
    palette=colors,
    s=200, # Size of points (s is area, roughly matches R's size=5 for prominence)
    edgecolor="black",
    alpha=0.8,
    linewidth=1 # Border width for filled circles
)

# The Crosshairs
ax.axvline(x=structural_thresh, linestyle="dashed", color="darkred", linewidth=1.2)
ax.axhline(y=biological_thresh, linestyle="dashed", color="darkblue", linewidth=1.2)

# Set plot title and labels
ax.set_title("In Silico Gatekeeper vs. Biological Ground Truth", fontsize=15, fontweight='bold')
ax.set_xlabel("Maximum Positional Entropy (Ent_Max)", fontsize=12, fontweight='bold')
ax.set_ylabel("Experimental Indel Efficiency (%)", fontsize=12, fontweight='bold')

# Set subtitle using plt.suptitle for placement control
plt.suptitle(f"Crosshairs: Ent_Max Threshold ({structural_thresh:.3f}) | Indel Threshold ({biological_thresh:.0f}%)",
             y=0.92, x=0.5, ha='center', fontsize=11)

# Configure legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels, title="", loc='right', fontsize=11)

# Adjust plot layout to make space for titles
plt.tight_layout(rect=[0, 0.03, 1, 0.9])

# Save the plot
output_filename = "EntMax_Quadrant_Analysis.pdf"
plt.savefig(output_filename, dpi=300)
print(f"Saved Quadrant Plot to: {output_filename}")
