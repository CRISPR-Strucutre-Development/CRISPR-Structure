import pandas as pd
from scipy.stats import spearmanr

def run_whole_library_correlation(struct_file, indel_file):
    print("Loading structural profiles and continuous Indel data...")
    
    try:
        df_struct = pd.read_csv(struct_file)
        df_indel = pd.read_csv(indel_file)
    except FileNotFoundError as e:
        print(f"CRITICAL ERROR: {e}")
        return

    # Clean string identifiers
    df_struct['sgRNA'] = df_struct['sgRNA'].astype(str).str.strip()
    df_indel['sgRNA'] = df_indel['sgRNA'].astype(str).str.strip()
    
    # Merge the physics with the biology
    master_df = df_indel.merge(df_struct, on='sgRNA')
    
    if 'Indel_Rate' not in master_df.columns:
        print("CRITICAL ERROR: 'Indel_Rate' column not found in your CSV.")
        print(f"Available columns: {master_df.columns.tolist()}")
        return
        
    # Drop any rows missing either an Indel rate or an RMSD value
    clean_df = master_df.dropna(subset=['Indel_Rate', 'Mean_Scaffold_RMSD'])
    total_n = len(clean_df)
    
    if total_n < 20:
        print(f"WARNING: Low sample size detected (n={total_n}). Ensure your Indel file contains all guides.")

    metrics = ['Mean_Scaffold_RMSD', 'Max_Scaffold_RMSD', 'Ensemble_RMSD_Fluctuation']
    
    print("\n" + "="*110)
    print(f"        WHOLE-LIBRARY CORRELATION: SCAFFOLD DISTORTION VS CONTINUOUS INDEL RATE (n={total_n})")
    print("="*110)
    print(f"{'Structural Distortion Metric':<35} | {'Spearman Rho (ρ)':<20} | {'P-Value':<12}")
    print("-"*110)
    
    for metric in metrics:
        if metric not in clean_df.columns:
            continue
            
        # Calculate continuous rank correlation
        rho, p_val = spearmanr(clean_df[metric], clean_df['Indel_Rate'])
        
        sig_flag = " <-- *STATISTICALLY SIGNIFICANT DRIVER*" if p_val < 0.05 else ""
        print(f"{metric:<35} | {rho:<20.4f} | {p_val:<12.4e}{sig_flag}")
        
    print("="*110)
    
    # Save the merged master file for graphing in R/Python later
    clean_df.to_csv('Master_150_Guide_Correlation_Matrix.csv', index=False)
    print("Exported fully merged dataset to 'Master_150_Guide_Correlation_Matrix.csv' for plotting.")

if __name__ == "__main__":
    run_whole_library_correlation(
        struct_file='Scaffold_Ensemble_RMSD_Profiles.csv',
        indel_file='experimental_indels.csv'
    )
