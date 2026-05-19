import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

def analyze_per_region(dist_file, pair_file, stack_file, indel_file):
    # 1. Load and Standardize
    df_dist = pd.read_csv(dist_file, na_values=['NA'])
    df_pair = pd.read_csv(pair_file)
    df_stack = pd.read_csv(stack_file)
    df_indel = pd.read_csv(indel_file)
    df_indel['Indel_Rate'] = pd.to_numeric(df_indel['Indel_Rate'], errors='coerce')
    df_indel = df_indel.dropna(subset=['Indel_Rate'])

    # 2. Define High and Low Cutoffs
    high_ids = df_indel[df_indel['Indel_Rate'] > 80]['sgRNA'].tolist()
    low_ids = df_indel[df_indel['Indel_Rate'] < 15]['sgRNA'].tolist()

    # 3. Filter structural data for these specific IDs
    regions = ['TL_32_37', 'SL1_54_60', 'SL2_71_76', 'SL3_88_90']
    
    # Process each region separately
    for reg in regions:
        print(f"\n{'='*10} ANALYZING REGION: {reg} {'='*10}")
        
        # Aggregations
        dist_reg = df_dist.groupby('sgRNA')[reg].min().reset_index()
        pair_reg = df_pair.groupby('sgRNA')[reg].max().reset_index()
        stack_reg = df_stack.groupby('sgRNA')[reg].max().reset_index()

        # --- PART A: DISTANCE (Placeholder Logic per user request) ---
        dist_reg['Group'] = 'Middle'
        dist_reg.loc[dist_reg['sgRNA'].isin(high_ids), 'Group'] = 'High (>80%)'
        dist_reg.loc[dist_reg['sgRNA'].isin(low_ids), 'Group'] = 'Low (<15%)'
        
        reg_dist_data = dist_reg[dist_reg['Group'] != 'Middle'].merge(df_indel, on='sgRNA')
        
        high_dist = reg_dist_data[reg_dist_data['Group'].str.contains('High')][reg]
        low_dist = reg_dist_data[reg_dist_data['Group'].str.contains('Low')][reg]
        
        p_dist = 1.0 # Default
        if not high_dist.empty and not low_dist.empty:
            _, p_dist = mannwhitneyu(high_dist, low_dist)
            print(f"Distance p-value (High vs Low): {p_dist:.4e}")

        plt.figure(figsize=(6, 5))
        sns.boxplot(data=reg_dist_data, x='Group', y=reg, hue='Group', palette='Set2', legend=False)
        plt.title(f'Heavy-Atom Distance: {reg}\np = {p_dist:.4f}')
        plt.ylabel('Min Distance (Å)')
        plt.savefig(f'Tail_Dist_{reg}.png', dpi=300)
        plt.close()

        # --- PART B: PAIRING INTERACTION STATS ---
        pair_reg_merged = pair_reg[pair_reg['sgRNA'].isin(high_ids + low_ids)].merge(df_indel, on='sgRNA')
        
        # Calculate p-value comparing Indel Rates of 'Locked' (1) vs 'Open' (0)
        p_pair = 1.0
        grp0 = pair_reg_merged[pair_reg_merged[reg] == 0]['Indel_Rate']
        grp1 = pair_reg_merged[pair_reg_merged[reg] == 1]['Indel_Rate']
        
        if len(grp0) > 0 and len(grp1) > 0:
            _, p_pair = mannwhitneyu(grp0, grp1)
            print(f"Pairing p-value for {reg}: {p_pair:.4e}")
        else:
            print(f"Skipping Pairing p-value for {reg}: One group is empty.")

        plt.figure(figsize=(6, 5))
        sns.boxplot(data=pair_reg_merged, x=reg, y='Indel_Rate', hue=reg, palette='Reds', legend=False)
        plt.title(f'Pairing Lock Impact: {reg}\np = {p_pair:.4f}')
        plt.xlabel(f'Pairing in {reg} (0=No, 1=Yes)')
        plt.ylabel('Experimental Indel %')
        plt.savefig(f'Tail_Pairing_{reg}.png', dpi=300)
        plt.close()

        # --- PART C: STACKING INTERACTION STATS ---
        stack_reg_merged = stack_reg[stack_reg['sgRNA'].isin(high_ids + low_ids)].merge(df_indel, on='sgRNA')
        
        p_stack = 1.0
        s_grp0 = stack_reg_merged[stack_reg_merged[reg] == 0]['Indel_Rate']
        s_grp1 = stack_reg_merged[stack_reg_merged[reg] == 1]['Indel_Rate']
        
        if len(s_grp0) > 0 and len(s_grp1) > 0:
            _, p_stack = mannwhitneyu(s_grp0, s_grp1)
            print(f"Stacking p-value for {reg}: {p_stack:.4e}")
        else:
            print(f"Skipping Stacking p-value for {reg}: One group is empty.")

        plt.figure(figsize=(6, 5))
        sns.boxplot(data=stack_reg_merged, x=reg, y='Indel_Rate', hue=reg, palette='Oranges', legend=False)
        plt.title(f'Stacking Lock Impact: {reg}\np = {p_stack:.4f}')
        plt.xlabel(f'Stacking in {reg} (0=No, 1=Yes)')
        plt.ylabel('Experimental Indel %')
        plt.savefig(f'Tail_Stacking_{reg}.png', dpi=300)
        plt.close()

    print("\nAll regional tail plots (Dist, Pairing, Stacking) saved.")

if __name__ == "__main__":
    analyze_per_region('sgrna_interactions_distance.csv', 
                       'sgrna_interactions_pair.csv', 
                       'sgrna_interactions_stack.csv', 
                       'experimental_indels.csv')