import pandas as pd

import pandas as pd

def update_tsv_with_spacer(tsv_filepath, csv_filepath, output_filepath):
    # 1. Load the data
    tsv_df = pd.read_csv(tsv_filepath, sep='\t')
    csv_df = pd.read_csv(csv_filepath, sep=',')

    # 2. Isolate relevant columns and prep for merge
    # Transform CSV sgRNA to DNA format (U -> T) BEFORE merging
    csv_subset = csv_df[['ID', 'sgRNA', 'indel']].copy()
    csv_subset['sgRNA'] = csv_subset['sgRNA'].str.upper().str.replace('U', 'T')

    # 3. Ensure TSV sgRNA is also cleaned/standardized if necessary
    tsv_df['sgRNA'] = tsv_df['sgRNA'].str.upper().str.replace('U', 'T')

    # 4. Execute Left Join
    # This keeps all rows from TSV and adds the 'indel' column where ID and sgRNA match
    updated_tsv_df = tsv_df.merge(csv_subset, on=['ID', 'sgRNA'], how='left')

    # 5. Write the full DataFrame to disk
    updated_tsv_df.to_csv(output_filepath, sep='\t', index=False)

# Execution
if __name__ == "__main__":
    update_tsv_with_spacer(
        tsv_filepath='entropy_advanced_strict.tsv', 
        csv_filepath='SpCas9_2019_150_sequence.csv', 
        output_filepath='entropy_advanced_strict_updated.tsv'
    )