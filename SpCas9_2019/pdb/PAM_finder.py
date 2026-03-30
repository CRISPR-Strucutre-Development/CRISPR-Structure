import pandas as pd

# 1. Load the data
input_file = 'SpCas9_2020_sequence.csv'
df = pd.read_csv(input_file)

# 2. Function to analyze the PAM
def analyze_pam(row):
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
