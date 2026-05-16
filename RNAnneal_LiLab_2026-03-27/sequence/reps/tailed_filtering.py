import pandas as pd
import sys

def extract_tail_structures(filtered_indel_file, pair_file, stack_file, dist_file):
    try:
        # 1. Load the Filtered IDs (The 1-10% and >80% cohorts)
        df_tails = pd.read_csv(filtered_indel_file)
        tail_ids = df_tails['sgRNA'].unique().tolist()
        
        print(f"Loaded {len(tail_ids)} target sgRNAs from filtered indel file.")

        # 2. Define File Mapping for the extraction
        source_files = {
            'tail-base-pairing.csv': pair_file,
            'tail-base-stacking.csv': stack_file,
            'tail-heavy-atom.csv': dist_file
        }

        # 3. Perform the Extraction
        for output_name, input_path in source_files.items():
            print(f"Processing {input_path}...")
            
            # Load the large interaction file
            df_full = pd.read_csv(input_path)
            
            # Ensure ID column is standardized to 'sgRNA' for the filter
            # (Adjust 'ID' to 'sgRNA' if your raw files use 'ID')
            if 'sgRNA' not in df_full.columns and 'ID' in df_full.columns:
                df_full = df_full.rename(columns={'ID': 'sgRNA'})

            # Filter for only the tail IDs
            df_filtered = df_full[df_full['sgRNA'].isin(tail_ids)].copy()
            
            # Save the new tail-specific file
            df_filtered.to_csv(output_name, index=False)
            print(f"Saved {len(df_filtered)} interaction rows to {output_name}")

        print("\n--- Structural Extraction Complete ---")
        print("Files generated: tail-base-pairing.csv, tail-base-stacking.csv, tail-heavy-atom.csv")

    except Exception as e:
        print(f"Error during extraction: {e}")

if __name__ == "__main__":
    # Ensure all file paths are correct for your ihc-grid-1-1-1 environment
    extract_tail_structures(
        'filtered_indels_clean.csv',        # The file from v3 script
        'sgrna_interactions_pair.csv',      # Raw pairing file
        'sgrna_interactions_stack.csv',     # Raw stacking file
        'sgrna_interactions_distance.csv'   # Raw distance file
    )
