import sys
import csv

def to_rna(dna_seq):
    """Replaces T with U to convert DNA to RNA."""
    return dna_seq.upper().replace('T', 'U')

def reverse_complement_rna(rna_seq):
    """Generates the reverse complement of an RNA sequence."""
    complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    # Complement then reverse
    comp_seq = "".join([complement.get(base, base) for base in rna_seq.upper()])
    return comp_seq[::-1]

def calculate_mismatches(ref_dna, query_dna):
    """
    Compares the off-target against the on-target.
    Returns the total count and 1-based indices of mismatches.
    """
    mismatches = 0
    indices = []
    
    length = min(len(ref_dna), len(query_dna))
    for i in range(length):
        if ref_dna[i] != query_dna[i]:
            mismatches += 1
            indices.append(str(i + 1))
            
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
    
    guide_rna_fixed = to_rna(on_target_dna[:20])
    results = []

    # Sequence processing logic
    def process_sequence(dna):
        rna = to_rna(dna)
        rev_comp_rna = reverse_complement_rna(rna)
        return dna, rna, rev_comp_rna

    # 1. Process On-Target
    dna, rna, rev_comp = process_sequence(on_target_dna)
    results.append({
        "Target Sequence (DNA)": dna,
        "Target Sequence (RNA)": rna,
        "Target Sequence: RNA complementary": rev_comp,
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
                dna_off, rna_off, rev_comp_off = process_sequence(off_dna)
                
                results.append({
                    "Target Sequence (DNA)": dna_off,
                    "Target Sequence (RNA)": rna_off,
                    "Target Sequence: RNA complementary": rev_comp_off,
                    "Guide_RNA": guide_rna_fixed,
                    "Mismatch_Count": count,
                    "Mismatch_Indices": idx_str
                })
    except FileNotFoundError:
        sys.exit(f"Error: File '{off_target_file}' not found.")

    # 3. Save Results
    output_name = "off_target_analysis.csv"
    fields = [
        "Target Sequence (DNA)", 
        "Target Sequence (RNA)", 
        "Target Sequence: RNA complementary", 
        "Guide_RNA", 
        "Mismatch_Count", 
        "Mismatch_Indices"
    ]
    
    with open(output_name, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)

    print(f"File generated: {output_name}")
    print(f"Guide RNA: {guide_rna_fixed}")

if __name__ == "__main__":
    main()
