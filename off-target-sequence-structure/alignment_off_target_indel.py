import pandas as pd
import sys
import os

def to_rna(dna_seq):
    """Converts DNA sequence string to RNA by replacing T with U."""
    if pd.isna(dna_seq):
        return ""
    return str(dna_seq).upper().replace('T', 'U')

def find_mismatches(ref, query):
    """
    Performs a position-wise comparison.
    Returns count and 0-indexed positions.
    """
    ref = str(ref).upper()
    query = str(query).upper()
    
    length = min(len(ref), len(query))
    mismatch_indices = []
    
    for i in range(length):
        if ref[i] != query[i]:
            mismatch_indices.append(str(i))
            
    count = len(mismatch_indices)
    count += abs(len(ref) - len(query)) # Handle length differences
    
    indices_str = ";".join(mismatch_indices) if mismatch_indices else "None"
    return count, indices_str

def main():
    if len(sys.argv) != 2:
        print("Usage: python align_targets.py <input_file.csv>")
        sys.exit(1)

    input_file = sys.argv[1]
    
    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)

    try:
        # Load the CSV
        df = pd.read_csv(input_file)
        
        # Ensure the column exists
        target_col = "Target (3' -> 5')"
        if target_col not in df.columns:
            print(f"Error: Column '{target_col}' not found in CSV.")
            sys.exit(1)

        # 1. Generate the RNA Version Column
        df['Target sequence (RNA version)'] = df[target_col].apply(to_rna)

        # 2. Identify On-Target Reference for alignment
        on_target_data = df[df['ID'] == 'on-target-1']
        if on_target_data.empty:
            print("Critical Error: 'on-target-1' not found in the ID column.")
            sys.exit(1)
        
        ref_seq = on_target_data[target_col].values[0]

        # 3. Perform Alignment and Mismatch Calculation
        alignment_results = df[target_col].apply(lambda x: find_mismatches(ref_seq, x))
        
        df['total_mismatches'] = [r[0] for r in alignment_results]
        df['mismatch_positions_0_idx'] = [r[1] for r in alignment_results]

        # 4. Save and Output
        output_filename = "processed_alignment_with_rna.csv"
        df.to_csv(output_filename, index=False)
        
        print("Processing successful.")
        print(f"Results saved to: {output_filename}")
        print(df[['ID', 'Target sequence (RNA version)', 'total_mismatches']].head())

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
