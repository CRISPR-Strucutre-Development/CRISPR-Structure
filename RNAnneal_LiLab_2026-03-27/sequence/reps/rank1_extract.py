import pandas as pd

def extract_rank1_truth_table(tail_pairing_file, tail_stacking_file, filtered_indels_file):
    # 1. Load Data
    df_p = pd.read_csv(tail_pairing_file)
    df_s = pd.read_csv(tail_stacking_file)
    df_i = pd.read_csv(filtered_indels_file)
    
    # 2. Filter for Rank 1 ONLY
    # Assuming your Rank column contains strings like 'rank01' or integers
    # Adjust the filter condition if your rank naming is different
    df_p_r1 = df_p[df_p['Rank'].astype(str).str.contains('1')].copy()
    df_s_r1 = df_s[df_s['Rank'].astype(str).str.contains('1')].copy()

    regions = ['TL_32_37', 'SL1_54_60', 'SL2_71_76', 'SL3_88_90']
    
    # 3. Rename columns to indicate Rank 1
    df_p_r1 = df_p_r1[['sgRNA'] + regions].rename(columns={c: f"{c}_R1_Pair" for c in regions})
    df_s_r1 = df_s_r1[['sgRNA'] + regions].rename(columns={c: f"{c}_R1_Stack" for c in regions})

    # 4. Merge with Indel Phenotypes
    master_r1 = df_i.merge(df_p_r1, on='sgRNA').merge(df_s_r1, on='sgRNA')
    
    master_r1.to_csv('sgRNA_Rank1_Truth_Table.csv', index=False)
    print("Rank 1 Truth Table generated: sgRNA_Rank1_Truth_Table.csv")

if __name__ == "__main__":
    extract_rank1_truth_table('tail-base-pairing.csv', 'tail-base-stacking.csv', 'filtered_indels_clean.csv')
