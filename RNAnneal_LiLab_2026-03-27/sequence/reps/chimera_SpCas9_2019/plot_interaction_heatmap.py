import sys
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')  # Suppress GUI display for headless HPC nodes
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import re

def get_sgrna_num(s):
    """Extracts the integer number from sgRNA identifier for natural sorting."""
    match = re.search(r'\d+', s)
    return int(match.group()) if match else 0

def process_heatmap(csv_file_path, output_suffix):
    """Processes, plots, and saves interaction heatmaps."""
    try:
        df = pd.read_csv(csv_file_path)
    except FileNotFoundError:
        print(f"Error: The file '{csv_file_path}' was not found.")
        sys.exit(1)
        
    interaction_cols = ['TL_32_37', 'SL1_54_60', 'SL2_71_76', 'SL3_88_90']
    
    # Check if all columns exist
    missing_cols = [col for col in interaction_cols if col not in df.columns]
    if missing_cols:
        print(f"Missing columns in {csv_file_path}: {missing_cols}")
        sys.exit(1)
        
    # Extract unique sgRNAs and sort them naturally (e.g., 1, 2, ..., 137)
    unique_sgrnas = sorted(df['sgRNA'].unique(), key=get_sgrna_num)
    
    # Convert 'sgRNA' to categorical with the correct numerical sequence
    df['sgRNA'] = pd.Categorical(df['sgRNA'], categories=unique_sgrnas, ordered=True)
    
    # Sort values by sgRNA and Rank
    df = df.sort_values(by=['sgRNA', 'Rank'])
    
    heatmap_data = df.set_index(['sgRNA', 'Rank'])
    heatmap_data = heatmap_data[interaction_cols]
    
    # Define a custom colormap: 0 = blue, 1 = red
    cmap = ListedColormap(['blue', 'red'])
    
    fig, ax = plt.subplots(figsize=(10, 12))
    sns.heatmap(
        heatmap_data, 
        cmap=cmap, 
        linewidths=0.5, 
        linecolor='lightgray', 
        cbar_kws={'label': 'Interaction Present (1 = Red, 0 = Blue)', 'ticks': [0, 1]},
        ax=ax
    )
    
    # Formatting the plot
    plt.title(f'sgRNA Interaction ({output_suffix.upper()}) Heatmap across Ranks', fontsize=14, pad=15)
    plt.ylabel('sgRNA, Rank', fontsize=12)
    plt.xlabel('Interaction Regions', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    
    plt.tight_layout()
    
    # Save files to both formats
    png_filename = f'interaction_heatmap_{output_suffix}.png'
    pdf_filename = f'interaction_heatmap_{output_suffix}.pdf'
    
    plt.savefig(png_filename, dpi=300)
    plt.savefig(pdf_filename, dpi=300)
    
    print(f"Success! Heatmap saved as {png_filename} and {pdf_filename}")
    plt.close('all')

def main():
    if len(sys.argv) < 3:
        print("Usage: python plot_interaction_heatmap.py <path_to_pair_csv> <path_to_stack_csv>")
        sys.exit(1)
        
    pair_file_path = sys.argv[1]
    stack_file_path = sys.argv[2]
    
    # Process base pairing
    print("Processing Base Pairing interactions...")
    process_heatmap(pair_file_path, output_suffix="pair")
    
    # Process base stacking
    print("\nProcessing Base Stacking interactions...")
    process_heatmap(stack_file_path, output_suffix="stack")

if __name__ == "__main__":
    main()
