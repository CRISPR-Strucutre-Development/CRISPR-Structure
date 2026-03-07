import pandas as pd

# 1. Load the data
input_file = 'SpCas9_2020_sequence.csv'
df = pd.read_csv(input_file)

# 2. Function to analyze the PAM
def analyze_pam(row):
    """Analyzes the Protospacer Adjacent Motif (PAM) sequence from a given row.

    This function extracts a target sequence, converts it to DNA, and validates
    if the critical 2-nucleotide PAM region (indices 25-26) is "GG". It also
    checks for sequence length validity based on Azimuth indexing requirements.

    Args:
        row (pd.Series): A row from a pandas DataFrame, expected to contain
                         a 'TargetSequence(RNAversion)' column.

    Returns:
        tuple: A tuple containing:
            - str: The validation status (e.g., "Valid (GG)", "Invalid (CA)", "Sequence Too Short").
            - str: The 2-nucleotide sequence extracted from the PAM region, or "N/A" if the sequence is too short.
    """
    # Convert RNA to DNA string for consistent indexing
    seq = str(row['TargetSequence(RNAversion)']).upper().replace('U', 'T')
    
    # Validation: Sequence must be 30nt for Azimuth indexing
    if len(seq) < 27:
        return "Sequence Too Short", "N/A"
    
    # Extract the critical 2 nucleotides that Azimuth checks (indices 25 and 26)
    critical_pam = seq[25:27]
    
    if critical_pam == "GG":
        return "Valid (GG)", critical_pam
    else:
        return "Invalid (" + critical_pam + ")", critical_pam

# 3. Apply the logic
results = df.apply(analyze_pam, axis=1)
df['PAM_Status'], df['Extracted_GG_Region'] = zip(*results)

# 4. Save to a separate diagnostic file
output_file = 'PAM_Validation_Results.tsv'
df[['TargetSequence(RNAversion)', 'Extracted_GG_Region', 'PAM_Status']].to_csv(output_file, sep='\t', index=False)

# 5. Summary Report for the user
print "Analysis Complete."
print "Valid (GG) count:   ", len(df[df['PAM_Status'] == "Valid (GG)"])
print "Invalid count:      ", len(df[df['PAM_Status'] != "Valid (GG)"])
print "Results saved to:    " + output_file
