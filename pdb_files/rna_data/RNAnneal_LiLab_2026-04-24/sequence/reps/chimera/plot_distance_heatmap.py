import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Suppress GUI display for headless HPC nodes
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_distance_heatmap_hpc.py <input_csv_file>")
        sys.exit(1)

    input_file = sys.argv[1]

    try:
        df = pd.read_csv(input_file, na_values=['NA'])
        
        # Combine sgRNA and rank for row labels
        df['sgRNA_Rank'] = df['sgRNA'].astype(str) + '_' + df['Rank'].astype(str)
        df.set_index('sgRNA_Rank', inplace=True)
        df.drop(columns=['sgRNA', 'Rank'], inplace=True)

        # Convert data to float
        df = df.astype(float)

        # Plot heatmap
        plt.figure(figsize=(10, 14))
        sns.heatmap(
            df,
            cmap='viridis_r',
            annot=True,
            fmt=".1f",
            mask=df.isnull(),
            linewidths=0.5,
            linecolor='lightgray',
            cbar_kws={'label': 'Minimum Heavy-Atom Distance (Å)'}
        )

        plt.title('Heavy-Atom Minimum Distance Matrix across Ranks', fontsize=14, pad=15)
        plt.ylabel('sgRNA, Rank', fontsize=12)
        plt.xlabel('Interaction Regions', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()

        # Save output in both PNG and PDF formats
        plt.savefig('distance_heatmap.png', dpi=300)
        plt.savefig('distance_heatmap.pdf', dpi=300)
        
        print("Success! Plots saved as 'distance_heatmap.png' and 'distance_heatmap.pdf'")

    except Exception as e:
        print(f"Error processing the file: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
