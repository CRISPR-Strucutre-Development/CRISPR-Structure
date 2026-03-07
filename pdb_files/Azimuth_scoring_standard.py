"""This script processes a CRISPR guide RNA sequence dataset to generate Azimuth off-target
prediction scores.

It performs the following steps:
1. Loads a CSV file containing guide RNA sequences.
2. Transforms the RNA sequences to DNA format (uppercase, 'U' to 'T') for
   compatibility with the Azimuth model.
3. Generates Azimuth off-target prediction scores using the transformed sequences.
4. Constructs an output DataFrame containing the original RNA sequences and their
   corresponding Azimuth scores.
5. Exports the results to a TSV file named 'Azimuth.tsv'.
"""

import pandas as pd
import numpy as np
import azimuth.model_comparison

# 1. Load the dataset
input_file = 'SpCas9_2020_sequence.csv'
df = pd.read_csv(input_file)

# 2. Transform sequences: Uppercase and convert RNA (U) to DNA (T)
# This ensures compatibility with the Azimuth model requirements
df['TransformedSequence'] = df['TargetSequence(RNAversion)'].str.upper().str.replace('U', 'T')

# 3. Generate Predictions
# Passing None for amino_acid_cut_positions and percent_peptides as requested
sequences = df['TransformedSequence'].values
scores = azimuth.model_comparison.predict(
    sequences, 
    None, 
    None
)

# 4. Construct output DataFrame
# Retains original RNA sequence and maps the new score
output_df = pd.DataFrame({
    'TargetSequence(RNAversion)': df['TargetSequence(RNAversion)'],
    'Score': scores
})

# 5. Export to TSV format
output_df.to_csv('Azimuth.tsv', sep='\t', index=False)

print("Process complete. Output saved to Azimuth.tsv.")
