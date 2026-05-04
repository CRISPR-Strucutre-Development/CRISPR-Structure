import pandas as pd

# Load the dataset
# Assuming your file is named 'sgrna_data.csv' and is tab-separated
# If it is space-separated, use sep='\s+'
df = pd.read_csv('sgrna_data.csv', sep='\s+')

# Extract the first 20 nucleotides
df['sgRNA_20nt'] = df['sgRNA'].str[:20]

# Select only the ID and the new 20nt sequence column
output_df = df[['sgRNA_20nt']]

# Save to a new CSV file
output_df.to_csv('extracted_sgRNAs.csv', index=False)

print("Extraction complete. Results saved to 'extracted_sgRNAs.csv'.")
