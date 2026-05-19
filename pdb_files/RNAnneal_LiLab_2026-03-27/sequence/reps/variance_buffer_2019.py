import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

def prove_structural_buffer(distance_file, segregated_indel_file):
    # 1. Load data
    df_dist = pd.read_csv(distance_file)
    df_indel = pd.read_csv(segregated_indel_file)
    
    # Clean ranks
    df_dist['Rank_Num'] = df_dist['Rank'].astype(str).str.extract(r'(\d+)').astype(int)
    regions = ['TL_32_37', 'SL1_54_60', 'SL2_71_76', 'SL3_88_90']
    dist_cols = [c for c in df_dist.columns if any(r in c for r in regions)]
    
    # 2. Compute Fluctuation (Standard Deviation across the 5 Ranks)
    # This measures the "Wiggle Room". High variance = flexible/buffered lock.
    ensemble_fluctuation = df_dist.groupby('sgRNA')[dist_cols].std().add_suffix('_Fluctuation')
    
    # 3. Merge with your segregated groups
    df_indel['sgRNA'] = df_indel['sgRNA'].astype(str).str.strip()
    ensemble_fluctuation = ensemble_fluctuation.reset_index()
    ensemble_fluctuation['sgRNA'] = ensemble_fluctuation['sgRNA'].astype(str).str.strip()
    
    master_df = df_indel.merge(ensemble_fluctuation, on='sgRNA')
    
    # Identify your high/low groups from the file
    group_col = [c for c in df_indel.columns if 'ind' not in c.lower() and 'sgrna' not in c.lower()][0]
    groups = master_df[group_col].unique()
    high_label = [g for g in groups if 'high' in str(g).lower() or 'good' in str(g).lower()][0]
    low_label = [g for g in groups if 'low' in str(g).lower() or 'bad' in str(g).lower()][0]
    
    high_group = master_df[master_df[group_col] == high_label]
    low_group = master_df[master_df[group_col] == low_label]
    
    print("="*115)
    print("                      PROVING THE STRUCTURAL BUFFER: ENSEMBLE FLUCTUATION (VARIANCE)")
    print("="*115)
    print(f"{'Structural Region':<30} | {f'Fluctuation {high_label}':<22} | {f'Fluctuation {low_label}':<22} | {'P-Value':<12}")
    print("-"*115)
    
    for col in [c for c in master_df.columns if '_Fluctuation' in c]:
        high_vals = high_group[col].dropna()
        low_vals = low_group[col].dropna()
        
        # Mann-Whitney to see if the flexibility profile is significantly different
        stat, p_val = mannwhitneyu(high_vals, low_vals, alternative='two-sided')
        sig_flag = " *SIGNIFICANT BUFFER*" if p_val < 0.05 else ""
        
        print(f"{col:<30} | {high_vals.mean():<22.2f} | {low_vals.mean():<22.2f} | {p_val:<12.4e}{sig_flag}")
    print("="*115)

if __name__ == "__main__":
    prove_structural_buffer('tail-heavy-atom.csv', 'filtered_indels_clean.csv')
