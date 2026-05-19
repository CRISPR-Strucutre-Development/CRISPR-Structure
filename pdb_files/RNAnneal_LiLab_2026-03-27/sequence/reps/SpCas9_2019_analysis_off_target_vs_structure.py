import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

def analyze_scaffold_crisis(file1_path, file2_path, output_csv='merged_structural_offtarget_analysis.csv'):
    print("Loading datasets...")
    # Read the two files
    df1 = pd.read_csv(file1_path)
    df2 = pd.read_csv(file2_path)
    
    # Standardize join columns
    df1.columns = df1.columns.str.strip()
    df2.columns = df2.columns.str.strip()
    
    # Rename join columns to match perfectly
    df1 = df1.rename(columns={'sgRNA': 'sgRNA_ID'})
    
    # Merge datasets on guide ID
    merged_df = pd.merge(df1, df2, on='sgRNA_ID', how='inner')
    print(f"Successfully merged {len(merged_df)} guides with complete structural and off-target metrics.")
    
    # Save the consolidated matrix
    merged_df.to_csv(output_csv, index=False)
    
    # -------------------------------------------------------------
    # STATISTICAL ANALYSIS: THE BLUNT TRUTH CORRELATIONS
    # -------------------------------------------------------------
    print("\n" + "="*80)
    print("      MATHEMATICAL CORRELATION MATRIX (STRUCTURAL VS. OFF-TARGET ANALYSIS)")
    print("="*80)
    
    structural_metrics = ['Mean_Scaffold_RMSD', 'Max_Scaffold_RMSD', 'Ensemble_RMSD_Fluctuation']
    offtarget_metrics = ['Max_True_OffTarget_Score', 'True_OffTarget_Matches', 'High_Risk_Sites_Count', 'Indel_Rate']
    
    results = []
    for struct in structural_metrics:
        for oft in offtarget_metrics:
            if struct in merged_df.columns and oft in merged_df.columns:
                # Drop NaNs just in case
                valid_data = merged_df[[struct, oft]].dropna()
                if len(valid_data) > 2:
                    p_coeff, p_val = pearsonr(valid_data[struct], valid_data[oft])
                    s_coeff, s_val = spearmanr(valid_data[struct], valid_data[oft])
                    results.append({
                        'Structural_Metric': struct,
                        'OffTarget_Metric': oft,
                        'Pearson_r': p_coeff,
                        'Pearson_p': p_val,
                        'Spearman_rho': s_coeff,
                        'Spearman_p': s_val
                    })
                    
    corr_df = pd.DataFrame(results)
    print(corr_df.to_string(index=False))
    print("="*80)
    
    # -------------------------------------------------------------
    # PLOT 1: THE SMOKING GUN SCATTER (RMSD vs Max True Off-Target Score)
    # -------------------------------------------------------------
    sns.set_theme(style='ticks')
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter plot with Indel_Rate mapping to color
    sc = ax.scatter(
        merged_df['Mean_Scaffold_RMSD'], 
        merged_df['Max_True_OffTarget_Score'], 
        c=merged_df['Indel_Rate'], 
        cmap='YlOrRd', 
        s=100, 
        edgecolors='black', 
        alpha=0.85
    )
    
    # Add trend line
    if len(merged_df) > 1:
        m, b = np.polyfit(merged_df['Mean_Scaffold_RMSD'], merged_df['Max_True_OffTarget_Score'], 1)
        ax.plot(merged_df['Mean_Scaffold_RMSD'], m*merged_df['Mean_Scaffold_RMSD'] + b, color='red', linestyle='--', alpha=0.7, label='Linear Trend')
    
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('On-Target Indel Rate (%)', fontsize=12, fontweight='bold', labelpad=10)
    
    ax.set_xlabel('Mean Scaffold Structural Distortion (RMSD $\AA$)', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_ylabel('Max True Genomic Off-Target Score (CFD)', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_title('The SpCas9 2019 Scaffold Specificity Crisis:\nStructural Deformation Decoupled From Off-Target Regulation', fontsize=13, fontweight='bold', pad=15)
    
    # Highlight quadrants
    ax.axhline(y=0.20, color='gray', linestyle=':', alpha=0.5)
    ax.text(ax.get_xlim()[0] + 0.5, 0.22, 'High-Risk Threshold (CFD ≥ 0.20)', color='gray', fontsize=10, style='italic')
    
    plt.tight_layout()
    plt.savefig('Scaffold_RMSD_vs_Max_OffTarget_Risk.png', dpi=300)
    plt.close()
    print("Saved 'Scaffold_RMSD_vs_Max_OffTarget_Risk.png'")
    
    # -------------------------------------------------------------
    # PLOT 2: THE PROMISCUITY MATRIX (Fluctuation vs Total Matches)
    # -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Log scale for True_OffTarget_Matches if variation is huge
    y_data = merged_df['True_OffTarget_Matches']
    y_label = 'Total Genomic Off-Target Matches (Count)'
    if y_data.max() / (y_data.min() + 1) > 10:
        y_data = np.log10(y_data + 1)
        y_label = 'Total Genomic Off-Target Matches ($\log_{10}$ Count)'
        
    sc2 = ax.scatter(
        merged_df['Ensemble_RMSD_Fluctuation'], 
        y_data, 
        c=merged_df['Max_True_OffTarget_Score'], 
        cmap='coolwarm', 
        s=120, 
        edgecolors='black', 
        alpha=0.9
    )
    
    cbar2 = plt.colorbar(sc2, ax=ax)
    cbar2.set_label('Max Individual Off-Target Score', fontsize=12, fontweight='bold', labelpad=10)
    
    ax.set_xlabel('Ensemble Scaffold Backbone Fluctuation ($\AA$)', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_ylabel(y_label, fontsize=12, fontweight='bold', labelpad=10)
    ax.set_title('Genomic Escape Footprint vs Scaffold Elasticity', fontsize=13, fontweight='bold', pad=15)
    
    plt.tight_layout()
    plt.savefig('Scaffold_Fluctuation_vs_Genomic_Footprint.png', dpi=300)
    plt.close()
    print("Saved 'Scaffold_Fluctuation_vs_Genomic_Footprint.png'")
    
    print("\nAnalysis complete. Consolidated metrics saved to CSV.")

if __name__ == "__main__":
    # Update filenames as needed for your pipeline execution
    analyze_scaffold_crisis(
        file1_path='Master_150_Guide_Correlation_Matrix.csv', 
        file2_path='Corrected_Guide_Level_Elevation_Cumulative_Risk.csv'
    )
