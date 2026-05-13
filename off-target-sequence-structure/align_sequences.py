import sys
import csv

def to_rna(dna_seq):
    """Replaces T with U to convert DNA to RNA."""
    return dna_seq.upper().replace('T', 'U')

def calculate_mismatches(ref_dna, query_dna):
    """
    Compares the off-target against the on-target.
    Returns the total count and 1-based indices of mismatches.
    """
    mismatches = 0
    indices = []
    
    # Compare position by position
    length = min(len(ref_dna), len(query_dna))
    for i in range(length):
        if ref_dna[i] != query_dna[i]:
            mismatches += 1
            indices.append(str(i + 1)) # 1-based index
            
    # Account for length differences as mismatches
    diff = abs(len(ref_dna) - len(query_dna))
    mismatches += diff
    
    idx_string = ";".join(indices) if indices else "None"
    return mismatches, idx_string

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <off_target_file.txt> <on_target_dna>")
        sys.exit(1)

    off_target_file = sys.argv[1]
    on_target_dna = sys.argv[2].upper()
    
    # DEFINITIVE LOGIC: 
    # Guide RNA is ONLY the first 20nt of the On-Target DNA, converted to RNA.
    guide_rna_fixed = to_rna(on_target_dna[:20])
    
    results = []

    # 1. Process the On-Target row
    results.append({
        "DNA_Sequence": on_target_dna,
        "RNA_Sequence": to_rna(on_target_dna),
        "Guide_RNA": guide_rna_fixed,
        "Mismatch_Count": 0,
        "Mismatch_Indices": "NA"
    })

    # 2. Process Off-Target file
    try:
        with open(off_target_file, 'r') as f:
            for line in f:
                off_dna = line.strip().upper()
                if not off_dna:
                    continue
                
                count, idx_str = calculate_mismatches(on_target_dna, off_dna)
                
                results.append({
                    "DNA_Sequence": off_dna,
                    "RNA_Sequence": to_rna(off_dna),
                    "Guide_RNA": guide_rna_fixed, # Constant from On-Target
                    "Mismatch_Count": count,
                    "Mismatch_Indices": idx_str
                })
    except FileNotFoundError:
        sys.exit(f"Error: File '{off_target_file}' not found.")

    # 3. Save Results
    output_name = "alignment_results.csv"
    fields = ["DNA_Sequence", "RNA_Sequence", "Guide_RNA", "Mismatch_Count", "Mismatch_Indices"]
    
    with open(output_name, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)

    print(f"File generated: {output_name}")
    print(f"Guide RNA used for all rows: {guide_rna_fixed}")

if __name__ == "__main__":
    main()
