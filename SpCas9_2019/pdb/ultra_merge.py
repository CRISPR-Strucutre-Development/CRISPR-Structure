import pandas as pd
import sys

def main():
    if len(sys.argv) != 3:
        print("Usage: python merge_by_id.py <file1.txt> <file2.txt>")
        sys.exit(1)

    file1, file2 = sys.argv[1], sys.argv[2]
    output_name = "Final_Merged_CRISPR_Data.txt"

    try:
        # Load files as tab-separated
        df1 = pd.read_csv(file1, sep='\t')
        df2 = pd.read_csv(file2, sep='\t')

        # Clean the ID column: remove potential whitespace to ensure match
        df1['ID'] = df1['ID'].astype(str).str.strip()
        df2['ID'] = df2['ID'].astype(str).str.strip()

        # Merge strictly on ID
        # how='inner' keeps only IDs present in BOTH files.
        # Use how='left' if you want to keep all IDs from File 1 regardless of File 2.
        merged_df = pd.merge(df1, df2, on='ID', how='inner', suffixes=('', '_drop'))

        # Define the exact columns you want for the final answer
        # Checking for both SL2 and SL1 variations based on your previous snippet
        target_cols = ['ID', 'SSC_score', 'indel', 
                       'SL2_mean', 'SL2_median', 'SL2_Mean', 'SL2_Median']
        
        # Filter only for columns that actually exist in the merged data
        existing_cols = [c for c in target_cols if c in merged_df.columns]
        final_df = merged_df[existing_cols]

        # Export to tab-separated format
        final_df.to_csv(output_name, sep='\t', index=False)
        
        print(f"--- Merge Summary ---")
        print(f"Rows matched on ID: {len(final_df)}")
        print(f"Columns extracted: {existing_cols}")
        print(f"Output saved to: {output_name}")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
