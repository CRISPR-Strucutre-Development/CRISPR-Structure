import csv

def generate_3prime_truncation_csv(baseline, target, filename="truncation_series.csv"):
    start_idx = baseline.find(target)
    end_idx = start_idx + len(target)
    
    if start_idx == -1:
        print("Error: Target segment not found in baseline.")
        return

    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Step", "Omitted_Nucleotides", "Edited_Sequence"])
        
        # Loop to truncate from the RIGHT (3' end) of the target segment
        for i in range(len(target) + 1):
            # Identifying the string being removed (e.g., A, AA, CAA...)
            omitted = target[len(target)-i:] if i > 0 else "Baseline"
            
            # Sequence = Prefix + (Target minus the last i bases) + Suffix
            edited_seq = baseline[:start_idx] + target[:len(target)-i] + baseline[end_idx:]
            
            writer.writerow([i, omitted, edited_seq])
    
    print(f"File '{filename}' generated successfully.")

# Your Data
base_seq = "ACCGCUUCACCAAAAGCUGUCCCUUAGGGGAUUAGAACUUGAGUGAAGGUGGGCUGCUUGCAUCAGCCUAAUGUCGAGAAGUGCUUUCUUCGGAAAGUAACCCUCGAAACAAAUUCAGUGCUCCUCUCCAAUUCUGCACAAGAAAGUUGCAGAACCCGAAUAGAGCAAUGAAGGAAUGCAAC"
tracr_leader = "CAAAUUCAGUGCUCCUCUCCAAUUCUGCACAA"

# Execute
generate_3prime_truncation_csv(base_seq, tracr_leader)
