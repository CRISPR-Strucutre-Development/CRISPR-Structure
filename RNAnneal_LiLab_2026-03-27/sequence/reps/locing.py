import pandas as pd

def generate_master_file(tail_pairing_file, tail_stacking_file, filtered_indels_file):
    try:
        # 1. Load Data
        df_p = pd.read_csv(tail_pairing_file)
        df_s = pd.read_csv(tail_stacking_file)
        df_i = pd.read_csv(filtered_indels_file)
        
        # Clean IDs
        for df in [df_p, df_s, df_i]:
            df['sgRNA'] = df['sgRNA'].astype(str).str.strip()
        
        regions = ['TL_32_37', 'SL1_54_60', 'SL2_71_76', 'SL3_88_90']
        
        def get_binary_votes(df, suffix):
            # Sum interactions across 5 ranks
            votes = df.groupby('sgRNA')[regions].sum()
            # Apply Rule: 1 if interaction exists in >= 3 ranks, else 0
            binary = (votes >= 3).astype(int)
            return binary.add_suffix(suffix)

        # 2. Extract Majority Binary Flags
        pair_binary = get_binary_votes(df_p, '_Pair')
        stack_binary = get_binary_votes(df_s, '_Stack')

        # 3. Merge into the Master Table
        # Starting with the filtered Indel data as the base
        master_df = df_i.merge(pair_binary, on='sgRNA', how='left')
        master_df = master_df.merge(stack_binary, on='sgRNA', how='left')

        # 4. Final Formatting
        # Organizing columns: sgRNA, Group, Indel_Rate, then all structural flags
        cols = ['sgRNA', 'Group', 'Indel_Rate'] + \
               [c for c in master_df.columns if '_Pair' in c or '_Stack' in c]
        
        master_df = master_df[cols]

        # 5. Save the file
        output_name = 'sgRNA_Structural_Truth_Table.csv'
        master_df.to_csv(output_name, index=False)
        
        print(f"--- Master File Generated ---")
        print(f"Rows: {len(master_df)}")
        print(f"Columns: {list(master_df.columns)}")
        print(f"Saved as: {output_name}")

    except Exception as e:
        print(f"Error generating master file: {e}")

if __name__ == "__main__":
    generate_master_file(
        'tail-base-pairing.csv', 
        'tail-base-stacking.csv', 
        'filtered_indels_clean.csv'
    )