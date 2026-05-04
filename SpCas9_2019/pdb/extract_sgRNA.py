import pandas as pd

# Load the source CSV
df = pd.read_csv('SpCas9_2019_150_sequence.csv')

# 1. Translate U to T in the sequence column
df['sequence'] = df['sequence'].str.replace('U', 'T')

# 2. Filter to keep ONLY the sequence column
final_df = df[['sequence']]

# 3. Save to the final CSV
final_df.to_csv('extracted_sgRNA.csv', index=False)
