import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, ks_2samp

def run_percentile_tail_squeeze(merged_csv_path):
    df = pd.read_csv(merged_csv_path)
    df.columns = df.columns.str.strip()
    
    # Sort strictly by Indel Rate to establish flawless rankings
    df_sorted = df.sort_values(by='Indel_Rate', ascending=False).reset_index(drop=True)
    
    # Extract the exact top 10% (15 guides) and bottom 10% (15 guides)
    N = 15
    top_10_percent = df_sorted.head(N).copy()
    bottom_10_percent = df_sorted.tail(N).copy()
    
    top_10_percent['Tail_Group'] = 'Top 10% Elite Elite'
    bottom_10_percent['Tail_Group'] = 'Bottom 10% Near Dead'
    
    squeeze_df = pd.concat([top_10_percent, bottom_10_percent])
    
    g1 = top_10_percent['Max_True_OffTarget_Score'].values
    g2 = bottom_10_percent['Max_True_OffTarget_Score'].values
    
    print("="*80)
    print(f"      STRICT PERCENTILE SQUEEZE ANALYSIS (N = {N} vs N = {N})")
    print("="*80)
    print(f"Top 10% Indel Range   : {top_10_percent['Indel_Rate'].min():.2f}% to {top_10_percent['Indel_Rate'].max():.2f}%")
    print(f"Bottom 10% Indel Range: {bottom_10_percent['Indel_Rate'].min():.2f}% to {bottom_10_percent['Indel_Rate'].max():.2f}%\n")
    
    # Mann-Whitney U Test
    _, mw_p = mannwhitneyu(g1, g2, alternative='two-sided')
    print(f"[1] MANN-WHITNEY U TEST P-VALUE: {mw_p:.6f}")
    
    # KS Test
    _, ks_p = ks_2samp(g1, g2)
    print(f"[2] KOLMOGOROV-SMIRNOV P-VALUE : {ks_p:.6f}")
    print("="*80)

if __name__ == "__main__":
    run_percentile_tail_squeeze('merged_structural_offtarget_analysis.csv')