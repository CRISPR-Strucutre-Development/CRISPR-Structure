import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

def calculate_structural_enrichment(truth_table_file):
    df = pd.read_csv(truth_table_file)
    
    # Define categorical outcomes
    df['is_failure'] = df['Indel_Rate'] < 10
    
    regions = [c for c in df.columns if '_R1_' in c]
    results = []

    for col in regions:
        # Create 2x2 Contingency Table
        # Rows: Locked (1) vs Open (0)
        # Cols: Failure (True) vs Success (False)
        contingency_table = pd.crosstab(df[col], df['is_failure'])
        
        # If the table isn't 2x2 (missing a category), skip
        if contingency_table.shape == (2, 2):
            odds_ratio, p_value = fisher_exact(contingency_table)
            
            # Calculate % of failures in each group
            fail_rate_locked = (contingency_table.iloc[1, 1] / contingency_table.iloc[1].sum()) * 100
            fail_rate_open = (contingency_table.iloc[0, 1] / contingency_table.iloc[0].sum()) * 100
            
            results.append({
                'Feature': col,
                'Odds_Ratio': round(odds_ratio, 2),
                'P_Value_Fisher': round(p_value, 4),
                'Fail_Rate_Locked_%': round(fail_rate_locked, 1),
                'Fail_Rate_Open_%': round(fail_rate_open, 1),
                'Enrichment_Factor': round(fail_rate_locked / fail_rate_open, 2) if fail_rate_open > 0 else 0
            })

    results_df = pd.DataFrame(results).sort_values('P_Value_Fisher')
    print("--- Structural Enrichment (Fisher's Exact Test) ---")
    print(results_df.to_string(index=False))
    results_df.to_csv('Structural_Enrichment_Results.csv', index=False)

if __name__ == "__main__":
    calculate_structural_enrichment('sgRNA_Rank1_Truth_Table.csv')
