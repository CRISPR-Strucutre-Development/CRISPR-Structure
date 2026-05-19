import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') # Headless mode for HPC
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

def load_and_process(file_path, label):
    df = pd.read_csv(file_path)
    # Identify value columns by their prefix
    val_cols = [c for c in df.columns if any(p in c for p in ['TL_', 'SL1_', 'SL2_', 'SL3_'])]
    
    df_melted = df.melt(id_vars=['sgRNA', 'Rank'], 
                        value_vars=val_cols, 
                        var_name='Region', 
                        value_name='Binary_Value')
    
    # Standardize names to TL, SL1, SL2, SL3
    df_melted['Region'] = df_melted['Region'].str.split('_').str[0]
    df_melted['Dataset'] = label
    return df_melted

def get_asterisks(p):
    if p <= 0.001: return "***"
    if p <= 0.01: return "**"
    if p <= 0.05: return "*"
    return "ns"

if len(sys.argv) < 3:
    print("Usage: python script.py <file1.csv> <file2.csv>")
    sys.exit(1)

# 1. Load datasets with requested labels
label_a = "spacer:scaffold [2020:2019]"
label_b = "spacer:scaffold[2019:2019]"

df1 = load_and_process(sys.argv[1], label_a)
df2 = load_and_process(sys.argv[2], label_b)
df_all = pd.concat([df1, df2])

# Set order for X-axis
regions = ['TL', 'SL1', 'SL2', 'SL3']
df_all['Region'] = pd.Categorical(df_all['Region'], categories=regions, ordered=True)

# 2. Setup Plot (Vertical Orientation)
plt.figure(figsize=(10, 8))

# x is now 'Region', y is 'Binary_Value' (the proportion)
ax = sns.barplot(data=df_all, 
                 x='Region', 
                 y='Binary_Value', 
                 hue='Dataset',
                 errorbar=('ci', 95), 
                 capsize=.1, 
                 palette=['grey', 'royalblue'])

# 3. Statistical Testing and Annotation
stats_results = []
for i, region in enumerate(regions):
    sub_a = df1[df1['Region'] == region]['Binary_Value']
    sub_b = df2[df2['Region'] == region]['Binary_Value']
    
    ones_a, zeros_a = int(sub_a.sum()), int(len(sub_a) - sub_a.sum())
    ones_b, zeros_b = int(sub_b.sum()), int(len(sub_b) - sub_b.sum())
    
    _, p_val = fisher_exact([[ones_a, zeros_a], [ones_b, zeros_b]])
    stars = get_asterisks(p_val)
    
    stats_results.append({'Region': region, 'P_Value': p_val, 'Significance': stars})
    
    # Calculate Y-position for the annotation (above the highest bar/error bar)
    # Using a 95% CI approximation for height placement
    mean_a = sub_a.mean()
    mean_b = sub_b.mean()
    # Adding a small buffer above the tallest bar
    y_max = max(mean_a, mean_b) + 0.12 
    
    # Place stars above the group
    plt.text(i, y_max, stars, ha='center', va='bottom', 
             fontweight='bold', fontsize=14)

# 4. Save and Finalize
pd.DataFrame(stats_results).to_csv('fisher_test_results.csv', index=False)

plt.ylabel('Proportion of 1s (Mean Density)', fontsize=12)
plt.xlabel('Structural Region', fontsize=12)
plt.title('sgRNA Feature Density Comparison', fontsize=14, pad=25)
plt.ylim(0, 1.2) # Extended for annotations
plt.legend(title='Scaffold Type', loc='upper right')
plt.grid(axis='y', linestyle='--', alpha=0.3)

plt.savefig('sgrna_vertical_comparison_stacking.png', dpi=300, bbox_inches='tight')
plt.savefig('sgrna_vertical_comparison_stacking.pdf', bbox_inches='tight')

print("Final vertical comparison plots and CSV generated.")
