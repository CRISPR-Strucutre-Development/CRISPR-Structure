import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_interaction_heatmaps(pair_file, stack_file):
    # 1. Load the interaction data
    df_pair = pd.read_csv(pair_file)
    df_stack = pd.read_csv(stack_file)

    # 2. Aggregate: Sum the 1s across the 5 ranks per sgRNA
    # This gives us a 'Frequency Score' from 0 to 5
    pair_freq = df_pair.groupby('sgRNA')[['TL_32_37', 'SL1_54_60', 'SL2_71_76', 'SL3_88_90']].sum()
    stack_freq = df_stack.groupby('sgRNA')[['TL_32_37', 'SL1_54_60', 'SL2_71_76', 'SL3_88_90']].sum()

    # 3. Create the Visualization
    # We use a very tall figure (30 inches) so 150 rows are actually readable
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 30))

    # Plot Pairing Heatmap
    sns.heatmap(pair_freq, ax=ax1, cmap="Reds", vmin=0, vmax=5, 
                cbar_kws={'label': 'Frequency (Number of Ranks)'},
                linewidths=0.05, linecolor='gray')
    ax1.set_title('Base-Pairing Frequency (Ensemble 1-5)', fontsize=16, pad=20)
    ax1.set_ylabel('sgRNA ID', fontsize=12)

    # Plot Stacking Heatmap
    sns.heatmap(stack_freq, ax=ax2, cmap="Oranges", vmin=0, vmax=5, 
                cbar_kws={'label': 'Frequency (Number of Ranks)'},
                linewidths=0.05, linecolor='gray')
    ax2.set_title('Base-Stacking Frequency (Ensemble 1-5)', fontsize=16, pad=20)
    ax2.set_ylabel('') # Remove label on second plot for cleanliness

    plt.tight_layout()
    output_name = 'Interaction_Ensemble_Heatmap.png'
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    print(f"Heatmap successfully saved as {output_name}")

if __name__ == "__main__":
    plot_interaction_heatmaps('sgrna_interactions_pair.csv', 'sgrna_interactions_stack.csv')
