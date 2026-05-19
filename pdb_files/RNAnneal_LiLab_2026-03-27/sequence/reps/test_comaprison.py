import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind, ks_2samp

def run_triple_stats(truth_table_file):
    # Load the master file we generated
    df = pd.read_csv(truth_table_file)
    
    regions = ['TL_32_37', 'SL1_54_60', 'SL2_71_76', 'SL3_88_90']
    types = ['Pair', 'Stack']
    
    results = []

    for i_type in types:
        for reg in regions:
            col_name = f"{reg}_{i_type}"
            
            # Split Indel rates based on majority lock (0 vs 1)
            group0 = df[df[col_name] == 0]['Indel_Rate']
            group1 = df[df[col_name] == 1]['Indel_Rate']
            
            if len(group0) > 2 and len(group1) > 2:
                # 1. Mann-Whitney U (Non-parametric)
                _, p_mwu = mannwhitneyu(group0, group1, alternative='two-sided')
                
                # 2. T-test (Parametric - comparing means)
                _, p_ttest = ttest_ind(group0, group1, equal_var=False)
                
                # 3. KS-test (Comparing distributions)
                _, p_ks = ks_2samp(group0, group1)
                
                results.append({
                    'Interaction': col_name,
                    'N_Open': len(group0),
                    'N_Locked': len(group1),
                    'P_MannWhitney': round(p_mwu, 5),
                    'P_Ttest': round(p_ttest, 5),
                    'P_KSTest': round(p_ks, 5)
                })

    # Convert to DataFrame for a clean summary
    results_df = pd.DataFrame(results)
    results_df.to_csv('Structural_Stats_Comparison.csv', index=False)
    
    print("--- Statistical Blunt Truth Table ---")
    print(results_df.to_string(index=False))

if __name__ == "__main__":
    run_triple_stats('sgRNA_Structural_Truth_Table.csv')
