import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

def run_2019_exact_descriptive(pairing_file, stacking_file, distance_file, segregated_indel_file, distance_cutoff=3.2):
    # Load files
    df_pair = pd.read_csv(pairing_file)
    df_stack = pd.read_csv(stacking_file)
    df_dist = pd.read_csv(distance_file)
    df_indel = pd.read_csv(segregated_indel_file)
    
    # Standardize ranks
    for df in [df_pair, df_stack, df_dist]:
        df['Rank_Num'] = df['Rank'].astype(str).str.extract(r'(\d+)').astype(int)
        
    regions = ['TL_32_37', 'SL1_54_60', 'SL2_71_76', 'SL3_88_90']
    
    # 1. Compute Pairing & Stacking Persistence (0 to 5)
    pair_features = df_pair.groupby('sgRNA')[regions].sum().add_suffix('_Pair_Persist')
    stack_features = df_stack.groupby('sgRNA')[regions].sum().add_suffix('_Stack_Persist')
    
    # 2. Compute Distance Persistence using 'tail-heavy-atom.csv'
    # Finds columns ending with '_Dist' or containing the region names
    dist_cols = [c for c in df_dist.columns if any(r in c for r in regions)]
    df_dist_binary = df_dist[['sgRNA', 'Rank_Num']].copy()
    
    for col in dist_cols:
        # 1 if the interaction is tightly closed, 0 if open
        df_dist_binary[f"{col}_Is_Closed"] = (df_dist[col] <= distance_cutoff).astype(int)
        
    closed_cols = [c for c in df_dist_binary.columns if '_Is_Closed' in c]
    dist_features = df_dist_binary.groupby('sgRNA')[closed_cols].sum()
    
    # Standardize names to match the table format
    dist_features.columns = [c.split('_')[0] + '_' + c.split('_')[1] + '_' + c.split('_')[2] + '_Dist_Persist' for c in dist_features.columns]

    # 3. Master Merge
    structural_matrix = pair_features.merge(stack_features, on='sgRNA').merge(dist_features, on='sgRNA').reset_index()
    df_indel['sgRNA'] = df_indel['sgRNA'].astype(str).str.strip()
    structural_matrix['sgRNA'] = structural_matrix['sgRNA'].astype(str).str.strip()
    
    master_df = df_indel.merge(structural_matrix, on='sgRNA')
    
    # 4. Separate based on your pre-existing categories in the file
    # Dynamically identifies the group column (handles 'Group', 'Phenotype', 'Class', etc.)
    group_col = [c for c in df_indel.columns if 'ind' not in c.lower() and 'sgrna' not in c.lower()][0]
    
    groups = master_df[group_col].unique()
    high_label = [g for g in groups if 'high' in str(g).lower() or 'good' in str(g).lower()][0]
    low_label = [g for g in groups if 'low' in str(g).lower() or 'bad' in str(g).lower()][0]
    
    high_group = master_df[master_df[group_col] == high_label]
    low_group = master_df[master_df[group_col] == low_label]
    
    all_features = [c for c in master_df.columns if '_Persist' in c]
    
    print("="*115)
    print("                      CORRECTED DESCRIPTIVE STATISTICS & ENSEMBLE SIGNALS (2019 SCAFFOLD)")
    print("="*115)
    print(f"{'Structural Feature':<30} | {f'Mean {high_label}':<22} | {f'Mean {low_label}':<22} | {'U-Statistic':<12} | {'P-Value':<12}")
    print("-"*115)
    
    for feat in all_features:
        high_vals = high_group[feat].dropna()
        low_vals = low_group[feat].dropna()
        
        stat, p_val = mannwhitneyu(high_vals, low_vals, alternative='two-sided')
        sig_flag = " *SIGNIFICANT*" if p_val < 0.05 else ""
        
        print(f"{feat:<30} | {high_vals.mean():<22.2f} | {low_vals.mean():<22.2f} | {stat:<12.1f} | {p_val:<12.4e}{sig_flag}")
    print("="*115)

if __name__ == "__main__":
    run_2019_exact_descriptive(
        pairing_file='tail-base-pairing.csv',
        stacking_file='tail-base-stacking.csv',
        distance_file='tail-heavy-atom.csv',
        segregated_indel_file='filtered_indels_clean.csv',
        distance_cutoff=3.2
    )