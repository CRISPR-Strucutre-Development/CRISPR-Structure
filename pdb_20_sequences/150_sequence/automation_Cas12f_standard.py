"""Generates an expanded sgRNA design file based on target and scaffold truncation series.

This script reads two input files:
    - 'Cas12f1_tracr_spacer_truncation.xlsx': Contains unique target definitions (Target, Target (RNA), MS2).
    - 'truncation_series.csv': Contains various scaffold sequences along with their truncation labels and step numbers.

It then iterates through each scaffold in the truncation series and combines it with every unique target to generate a comprehensive list of sgRNA designs. For each combination, it creates a unique ID and includes all relevant target and scaffold information.

The final expanded data is saved to a CSV file named 'Cas12f1_expanded_final.csv'.

Input Files:
    - Cas12f1_tracr_spacer_truncation.xlsx (Excel file)
    - truncation_series.csv (CSV file)

Output File:
    - Cas12f1_expanded_final.csv (CSV file)
"""
import pandas as pd

# 1. Load the input files
# Note: reading .xlsx requires the 'openpyxl' library (pip install openpyxl)
main_df = pd.read_excel('Cas12f1_tracr_spacer_truncation.xlsx')
series_df = pd.read_csv('truncation_series.csv')

# 2. Extract the unique targets (t1, t2, t3)
# These are the three target types that need to be repeated for every scaffold
unique_targets = main_df[['Target', 'Target (RNA)', 'MS2']].drop_duplicates().reset_index(drop=True)

expanded_rows = []

# 3. Iterate through the truncation series
for _, s_row in series_df.iterrows():
    step = s_row['Step']
    trunc_label = s_row['Tuncation'] # Column name in truncation_series.csv
    scaffold_seq = s_row['Scaffold']
    
    # Generate the sc ID (Step 0 = sc1, Step 1 = sc2, etc.)
    sc_idx = step + 1
    
    # 4. Create a row for each target (t1, t2, t3) using the current scaffold
    for t_idx, t_row in unique_targets.iterrows():
        t_num = t_idx + 1
        expanded_rows.append({
            'ID': f"sgRNA_t{t_num}_sc{sc_idx}",
            'Target': t_row['Target'],
            'Scaffold (Direct Repeat + tracRNA)': scaffold_seq,
            'Target (RNA)': t_row['Target (RNA)'],
            'MS2': t_row['MS2'],
            'Trunction': trunc_label
        })

# 5. Create final DataFrame and save as .csv
final_df = pd.DataFrame(expanded_rows)
final_df.to_csv('Cas12f1_expanded_final.csv', index=False)

print("File 'Cas12f1_expanded_final.csv' has been created.")
