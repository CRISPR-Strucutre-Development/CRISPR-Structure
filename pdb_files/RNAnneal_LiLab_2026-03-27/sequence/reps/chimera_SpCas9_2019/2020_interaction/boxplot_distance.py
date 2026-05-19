import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') # Headless mode for HPC
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import shapiro, ttest_ind, mannwhitneyu

def load_and_process(file_path, label):
    df = pd.read_csv(file_path, na_values=['NA', 'nan', 'NaN'])
    val_cols = [c for c in df.columns if any(p in c for p in ['TL_', 'SL1_', 'SL2_', 'SL3_'])]
    df_melted = df.melt(id_vars=['sgRNA', 'Rank'], 
                        value_vars=val_cols, 
                        var_name='Region', 
                        value_name='Distance')
    df_melted['Region'] = df_melted['Region'].str.split('_').str[0]
    df_melted['Dataset'] = label
    return df_melted

def get_asterisks(p):
    if pd.isna(p): return "ns"
    if p <= 0.001: return "***"
    if p <= 0.01: return "**"
    if p <= 0.05: return "*"
    return "ns"

if len(sys.argv) < 3:
    print("Usage: python boxplot_distance_stats.py <file1.csv> <file2.csv>")
    sys.exit(1)

label_a = "spacer:scaffold [2020:2019]"
label_b = "spacer:scaffold[2019:2019]"

df1 = load_and_process(sys.argv[1], label_a)
df2 = load_and_process(sys.argv[2], label_b)
df_all = pd.concat([df1, df2], ignore_index=True)

regions = ['TL', 'SL1', 'SL2', 'SL3']
df_all['Region'] = pd.Categorical(df_all['Region'], categories=regions, ordered=True)

# Plotting
plt.figure(figsize=(12, 8))
ax = sns.boxplot(data=df_all, x='Region', y='Distance', hue='Dataset', 
                 palette=['grey', 'royalblue'], showfliers=False)
sns.stripplot(data=df_all, x='Region', y='Distance', hue='Dataset', 
              dodge=True, alpha=0.3, palette=['black', 'black'], ax=ax, legend=False)

# Statistics Logic
stats_results = []
print(f"{'Region':<8} | {'Norm_A_p':<10} | {'Norm_B_p':<10} | {'Test_Used':<15} | {'Final_p':<10} | {'Sig'}")
print("-" * 75)

for i, region in enumerate(regions):
    sub_a = df1[df1['Region'] == region]['Distance'].dropna()
    sub_b = df2[df2['Region'] == region]['Distance'].dropna()
    
    p_norm_a, p_norm_b, p_val, test_used = np.nan, np.nan, np.nan, "N/A"
    
    if len(sub_a) >= 3 and len(sub_b) >= 3:
        # 1. Normality Check
        _, p_norm_a = shapiro(sub_a)
        _, p_norm_b = shapiro(sub_b)
        
        # 2. Choice of Test (Blunt Truth Logic)
        if p_norm_a > 0.05 and p_norm_b > 0.05:
            _, p_val = ttest_ind(sub_a, sub_b)
            test_used = "t-test"
        else:
            _, p_val = mannwhitneyu(sub_a, sub_b, alternative='two-sided')
            test_used = "Mann-Whitney"
    
    stars = get_asterisks(p_val)
    
    # Store for CSV
    stats_results.append({
        'Region': region,
        'Shapiro_A_p': p_norm_a,
        'Shapiro_B_p': p_norm_b,
        'Test_Used': test_used,
        'Final_P_Value': p_val,
        'Significance': stars
    })
    
    # Print to Command Line
    print(f"{region:<8} | {p_norm_a:<10.4f} | {p_norm_b:<10.4f} | {test_used:<15} | {p_val:<10.4f} | {stars}")
    
    # Annotation on Plot
    y_max = df_all[df_all['Region'] == region]['Distance'].max()
    if not pd.isna(y_max):
        plt.text(i, y_max * 1.05, stars, ha='center', va='bottom', fontweight='bold', fontsize=14)

# Save Stats to CSV
pd.DataFrame(stats_results).to_csv('statistical_analysis_report.csv', index=False)

plt.title('sgRNA Distance Comparison (Shapiro-Wilk + t-test/Mann-Whitney)')
plt.ylabel('Distance (Å)')
plt.ylim(0, df_all['Distance'].max() * 1.25)
plt.legend(title='Scaffold Type', loc='upper right')
plt.savefig('distance_comparison_final.png', dpi=300, bbox_inches='tight')

print("\nSuccess: 'statistical_analysis_report.csv' and 'distance_comparison_final.png' saved.")
