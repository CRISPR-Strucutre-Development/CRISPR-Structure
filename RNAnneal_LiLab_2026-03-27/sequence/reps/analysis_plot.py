import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, ks_2samp

def generate_rank1_boxplots(truth_table_file):
    # 1. Load the Rank 1 Truth Table
    df = pd.read_csv(truth_table_file)
    
    # 2. Define the new Rank 1 specific column names
    regions = ['TL_32_37', 'SL1_54_60', 'SL2_71_76', 'SL3_88_90']
    types = ['Pair', 'Stack']
    
    fig, axes = plt.subplots(2, 4, figsize=(22, 12), sharey=True)
    
    for i, i_type in enumerate(types):
        for j, reg in enumerate(regions):
            # Target the R1 columns specifically
            col_name = f"{reg}_R1_{i_type}"
            ax = axes[i, j]
            
            # Check if column exists to avoid crash
            if col_name not in df.columns:
                ax.set_title(f"Missing: {col_name}")
                continue

            # 3. Statistical Comparison
            group0 = df[df[col_name] == 0]['Indel_Rate']
            group1 = df[df[col_name] == 1]['Indel_Rate']
            
            p_mwu, p_ks = "N/A", "N/A"
            if len(group0) > 1 and len(group1) > 1:
                _, p_mwu = mannwhitneyu(group0, group1)
                _, p_ks = ks_2samp(group0, group1)
                p_mwu = f"{p_mwu:.4f}"
                p_ks = f"{p_ks:.4f}"
            
            # 4. Plotting
            sns.boxplot(ax=ax, data=df, x=col_name, y='Indel_Rate', 
                        palette='vlag', hue=col_name, legend=False)
            sns.stripplot(ax=ax, data=df, x=col_name, y='Indel_Rate', 
                          color='black', alpha=0.4, size=5)
            
            # 5. Descriptive Stats for the Title
            ax.set_title(f"{col_name}\nMWU-p: {p_mwu} | KS-p: {p_ks}", fontsize=11, fontweight='bold')
            ax.set_xlabel("Rank 1 Lock (0=No, 1=Yes)")
            if j == 0:
                ax.set_ylabel("Indel Rate (%)")
            else:
                ax.set_ylabel("")

    plt.suptitle("CRITICAL BIOPHYSICS: Impact of Rank 1 (Ground State) Locks on Indel Efficiency", fontsize=16, y=0.96)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    output_name = 'Rank1_Structural_Boxplots.png'
    plt.savefig(output_name, dpi=300)
    print(f"Extraction successful. Plot saved as {output_name}")

if __name__ == "__main__":
    generate_rank1_boxplots('sgRNA_Rank1_Truth_Table.csv')