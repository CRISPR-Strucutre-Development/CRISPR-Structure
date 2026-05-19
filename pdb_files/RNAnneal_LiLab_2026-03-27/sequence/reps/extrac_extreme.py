import pandas as pd
import sys

def extract_cleaned_tails(input_file, output_file='filtered_indels_clean.csv'):
    try:
        # Load file with auto-separator detection
        df = pd.read_csv(input_file, sep=None, engine='python')
        
        # Clean headers
        df.columns = [str(col).strip() for col in df.columns]
        
        # Verify required columns
        if 'sgRNA' not in df.columns or 'Indel_Rate' not in df.columns:
            print(f"Error: Required columns 'sgRNA' and 'Indel_Rate' not found.")
            return

        # Ensure numeric conversion and handle non-numeric entries
        df['Indel_Rate'] = pd.to_numeric(df['Indel_Rate'], errors='coerce')
        df = df.dropna(subset=['Indel_Rate'])

        # Apply the Refined Filters:
        # Success: > 80%
        # Failure: > 1% and < 10% (avoiding negative/zero noise)
        high_perf = df[df['Indel_Rate'] > 78].copy()
        low_perf = df[(df['Indel_Rate'] < 10) & (df['Indel_Rate'] > 1)].copy()
        
        # Label groups for analysis
        high_perf['Group'] = 'High'
        low_perf['Group'] = 'Low'
        
        # Merge and finalize
        final_df = pd.concat([high_perf, low_perf]).reset_index(drop=True)
        final_df.to_csv(output_file, index=False)
        
        print(f"--- Cleaned Tail Extraction Successful ---")
        print(f"High Efficiency (>75%):  {len(high_perf)}")
        print(f"Low Efficiency (1-10%): {len(low_perf)}")
        print(f"Total sequences saved:  {len(final_df)}")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        extract_cleaned_tails(sys.argv[1])
    else:
        print("Usage: python extract_tails_v3.py <input_file>")
