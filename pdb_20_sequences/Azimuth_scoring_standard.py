"""Predicts CRISPR/Cas9 cutting efficiency scores using the Azimuth model.

This script reads target RNA sequences from a CSV dataset, converts them into 
DNA sequences by substituting 'U' with 'T', and evaluates them using the 
Azimuth predictive model. The original RNA sequences and their corresponding 
predicted scores are then exported to a TSV file.

Attributes:
    input_file (str): Path to the input CSV file containing sequence data.
    df (pandas.DataFrame): DataFrame containing the loaded and transformed sequences.
    sequences (numpy.ndarray): Array of formatted DNA sequences for prediction.
    scores (numpy.ndarray): Array of Azimuth model prediction scores.
    output_df (pandas.DataFrame): DataFrame mapping original sequences to predicted scores.
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
