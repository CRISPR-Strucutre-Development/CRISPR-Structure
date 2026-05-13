import sys
import RNA
import pandas as pd
import csv

def get_pairings(dot_bracket):
    """Parse dot-bracket string to extract paired residues (1-based index)."""
    stack = []
    pairs = []
    for i, char in enumerate(dot_bracket):
        if char == '(':
            stack.append(i + 1)
        elif char == ')':
            left = stack.pop()
            pairs.append((left, i + 1))
    return pairs

def check_region_interaction(pairs, region_a, region_b):
    """Check if any base pairing exists between two interaction regions."""
    for p in pairs:
        r1, r2 = p
        if (r1 in region_a and r2 in region_b) or (r1 in region_b and r2 in region_a):
            return 1
    return 0

def main():
    if len(sys.argv) < 2:
        print("Usage: python rnafold_analysis.py <input_csv_file>")
        sys.exit(1)
        
    input_file = sys.argv[1]
    output_filename = 'rnafold_interaction_matrix.csv'
    
    # Corrected dictionary keys to match the accessed keys in the loop
    regions = {
        'TL_32_37': range(32, 38),
        'SL1_54_60': range(54, 61),
        'SL2_71_76': range(71, 77),
        'SL3_88_90': range(88, 91)
    }
    
    spacer_region = range(1, 21)
    rank_limit = 5
    
    try:
        df_input = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
        sys.exit(1)
        
    # Open output file for streaming data directly to disk
    with open(output_filename, mode='w', newline='') as f_out:
        writer = csv.writer(f_out)
        header = ['sgRNA', 'Rank', 'TL_32_37', 'SL1_54_60', 'SL2_71_76', 'SL3_88_90']
        writer.writerow(header)
        
        for index, row in df_input.iterrows():
            sgrna_id = row['ID']
            sequence = row['sgRNA']
            
            print(f"Processing {sgrna_id}...")
            
            fc = RNA.fold_compound(sequence)
            mfe_struct, mfe_en = fc.mfe()
            
            try:
                # Restrict to 2 kcal/mol to avoid exponential memory growth
                subopts = fc.subopt(1000)
            except Exception as e:
                print(f"Error processing {sgrna_id}: {e}", file=sys.stderr)
                continue
                
            rank_idx = 1
            for s in subopts:
                if rank_idx > rank_limit:
                    break
                    
                pairs = get_pairings(s.structure)
                
                out_row = [
                    sgrna_id,
                    f'rank{rank_idx:02d}',
                    check_region_interaction(pairs, regions['TL_32_37'], spacer_region),
                    check_region_interaction(pairs, regions['SL1_54_60'], spacer_region),
                    check_region_interaction(pairs, regions['SL2_71_76'], spacer_region),
                    check_region_interaction(pairs, regions['SL3_88_90'], spacer_region)
                ]
                
                writer.writerow(out_row)
                rank_idx += 1
                
    print(f"Success! Matrix saved to {output_filename}")

if __name__ == '__main__':
    main()
