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
    
    regions = {
        'TL_37_42': range(37, 43),
        'SL1_64_70': range(64, 71),
        'SL2_82_87': range(82, 88),
        'SL3_98_100': range(98, 101)
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
        header = ['sgRNA', 'Rank', 'TL_37_42', 'SL1_64_70', 'SL2_82_87', 'SL3_98_100']
        writer.writerow(header)
        
        for index, row in df_input.iterrows():
            sgrna_id = row['sgRNA']
            sequence = row['Sequence']
            
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
                    check_region_interaction(pairs, regions['TL_37_42'], spacer_region),
                    check_region_interaction(pairs, regions['SL1_64_70'], spacer_region),
                    check_region_interaction(pairs, regions['SL2_82_87'], spacer_region),
                    check_region_interaction(pairs, regions['SL3_98_100'], spacer_region)
                ]
                
                writer.writerow(out_row)
                rank_idx += 1
                
    print(f"Success! Matrix saved to {output_filename}")

if __name__ == '__main__':
    main()
