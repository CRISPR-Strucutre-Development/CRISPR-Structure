import pandas as pd
import sys

def main():
    if len(sys.argv) != 3:
        print("Usage: python merge_crispr.py <file1.txt> <file2.txt>")
        sys.exit(1)

    file1, file2 = sys.argv[1], sys.argv[2]
    output_name = "SSC_sequence_indel.txt"

    try:
        # Load files
        df1 = pd.read_csv(file1, sep=',')
        df2 = pd.read_csv(file2, sep='\t')

        # Clean sequence columns: Remove spaces and standardize to DNA (U->T) Uppercase
        # We create a 'merge_key' column so we don't destroy the original sequence formatting
        for df in [df1, df2]:
            if 'sequence' in df.columns:
                df['merge_key'] = df['sequence'].str.strip().str.upper().str.replace('U', 'T')
            else:
                print(f"Error: Column 'sequence' not found in one of the files.")
                sys.exit(1)

        # Perform Outer Join on the normalized merge_key
        # suffixes helps distinguish columns like 'ID' if they exist in both files
        merged_df = pd.merge(df1, df2, on='merge_key', how='outer', suffixes=('_f1', '_f2'))

        # Final cleanup: 
        # 1. Prioritize the original sequence from the first file, fill with second if missing
        if 'sequence_f1' in merged_df.columns:
            merged_df['sequence'] = merged_df['sequence_f1'].fillna(merged_df['sequence_f2'])
            merged_df.drop(columns=['sequence_f1', 'sequence_f2'], inplace=True)
        
        # 2. Remove the helper key
        merged_df.drop(columns=['merge_key'], inplace=True)

        # Save to tab-separated file
        merged_df.to_csv(output_name, sep='\t', index=False)
        print(f"Successfully merged. Results saved to {output_name}")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
