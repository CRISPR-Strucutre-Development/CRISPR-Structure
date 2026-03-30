import csv

# Define the constants
baseline_scaffold = "ACCGCUUCACCAAAAGCUGUCCCUUAGGGGAUUAGAACUUGAGUGAAGGUGGGCUGCUUGCAUCAGCCUAAUGUCGAGAAGUGCUUUCUUCGGAAAGUAACCCUCGAAACAAAUUCAGUGCUCCUCUCCAAUUCUGCACAAGAAAGUUGCAGAACCCGAAUAGAGCAAUGAAGGAAUGCAAC"
spacers = [
    "UUUGCACACACACAGUGGGCUACC", 
    "UUUGCAUCCCCAGGACACACACAC",
    "UUUAAGAACACAUACCCCUGGGCC"
]
ms2 = "UUUUAUUUU"

output_data = []

# 1. Process all baselines first
for spacer in spacers:
    full_baseline_seq = baseline_scaffold + spacer + ms2
    output_data.append([baseline_scaffold, spacer, ms2, "baseline", full_baseline_seq])
    
# 2. Process deletions step-by-step, cycling through all spacers at each step
for deletion_length in range(12, 25): 
    truncated_scaffold = baseline_scaffold[deletion_length:]
    label = f"{deletion_length}nt-deleted"
    
    for spacer in spacers:
        full_seq = truncated_scaffold + spacer + ms2
        output_data.append([truncated_scaffold, spacer, ms2, label, full_seq])

# Write the collected data to a CSV file
filename = "scaffold_sequences_interleaved.csv"
with open(filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Scaffold", "Spacer", "MS2", "Label", "Full_Sequence"])
    writer.writerows(output_data)

print(f"Data generation complete. File saved as: {filename}")
