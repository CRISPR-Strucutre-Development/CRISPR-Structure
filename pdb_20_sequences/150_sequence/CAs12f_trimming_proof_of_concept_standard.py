import csv

def generate_3prime_truncation_csv(baseline, target, filename="truncation_series.csv"):
    """Generates a CSV file with 3' truncations of a target segment within a baseline sequence.

    This function finds the first occurrence of the `target` segment within the `baseline`
    sequence. It then generates a series of new sequences where the `target` segment
    is progressively truncated from its 3' (right) end, one nucleotide at a time,
    until the entire segment is removed. These truncated sequences, along with
    the number of nucleotides omitted and the omitted string, are written to a CSV file.

    Args:
        baseline (str): The full sequence string in which the `target` segment is searched and modified.
        target (str): The specific subsequence to be found within `baseline` and then progressively
                      truncated from its 3' end.
        filename (str, optional): The name of the CSV file to be created.
                                  Defaults to "truncation_series.csv".

    Returns:
        None: The function prints status messages to the console and generates a CSV file
              as a side effect. If the target segment is not found, an error message
              is printed, and no file is generated.
    """
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
