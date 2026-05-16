import pandas as pd
import sys

def direct_extract(input_tsv, output_csv='experimental_indels.csv'):
    try:
        # Load the file with auto-detection
        df = pd.read_csv(input_tsv, sep=None, engine='python')
        
        # Strip hidden spaces from headers
        df.columns = [str(col).strip() for col in df.columns]
        
        # Direct extraction of the two columns as they are
        if 'ID' in df.columns and 'indel' in df.columns:
            clean_df = df[['ID', 'indel']].copy()
            
            # Rename for the next script's compatibility
            clean_df.columns = ['sgRNA', 'Indel_Rate']
            
            # Save the raw strings
            clean_df.to_csv(output_csv, index=False)
            
            print(f"Success! Extracted {len(clean_df)} rows.")
            print("\n--- RAW PREVIEW ---")
            print(clean_df.head(10))
        else:
            print(f"FAILED: ID or indel not found. Columns are: {df.columns.tolist()}")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        direct_extract(sys.argv[1])
    else:
        print("Usage: python script_name.py <input_file>")
