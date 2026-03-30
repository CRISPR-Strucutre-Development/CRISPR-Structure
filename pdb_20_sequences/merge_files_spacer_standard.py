import pandas as pd

def update_tsv_with_spacer(tsv_filepath, csv_filepath, output_filepath):
    """Updates a TSV file by merging it with spacer sequences from a CSV file.

    Reads a primary TSV file and a secondary CSV file, extracts the 'sgRNA' 
    and 'TargetSequence(RNAversion)' columns from the CSV, formats the RNA 
    sequence into a DNA format (uppercase, replaces 'U' with 'T'), and 
    performs a left join on the 'sgRNA' column. The resulting DataFrame is 
    exported as a new TSV file.

    Args:
        tsv_filepath (str): The file path to the input TSV file.
        csv_filepath (str): The file path to the input CSV file containing 
            the 'sgRNA' and 'TargetSequence(RNAversion)' columns.
        output_filepath (str): The file path where the updated TSV file 
            will be saved.

    Returns:
        None
    """
    # 1. Load the data
    tsv_df = pd.read_csv(tsv_filepath, sep='\t')
    csv_df = pd.read_csv(csv_filepath, sep=',')

    # 2. Isolate target variable and the new merge key (sgRNA)
    # Using .copy() prevents SettingWithCopyWarning when modifying the subset
    csv_subset = csv_df[['sgRNA', 'TargetSequence(RNAversion)']].copy()

    # 3. Transform Spacer.Sequence (Caps and U -> T)
    # Vectorized string operations for high performance
    csv_subset['TargetSequence(RNAversion)'] = csv_subset['TargetSequence(RNAversion)'].str.upper().str.replace('U', 'T')

    # 4. Execute Left Join on 'sgRNA'
    updated_tsv_df = tsv_df.merge(csv_subset, on='sgRNA', how='left')

    # 5. Write to disk
    updated_tsv_df.to_csv(output_filepath, sep='\t', index=False)

# Execution
update_tsv_with_spacer(
    tsv_filepath='entropy_cumulative_stats.tsv', 
    csv_filepath='SpCas9_2020_sequence.csv', 
    output_filepath='entropy_cumulative_stats_updated.tsv'
)
