import pandas as pd

def update_tsv_with_spacer(tsv_filepath, csv_filepath, output_filepath):
    # 1. Load the data
    tsv_df = pd.read_csv(tsv_filepath, sep='\t')
    csv_df = pd.read_csv(csv_filepath, sep=',')

    # 2. Clean TSV ID: Remove everything from the last underscore onwards (e.g., _rank01)
    # This regex r'_(?!.*_).*' matches the last underscore and everything following it
    # Or more simply, if the format is consistent: .str.split('_rank').str[0]
    tsv_df['ID'] = tsv_df['ID'].str.replace(r'_rank\d+', '', regex=True)

    # 3. Standardize sgRNA sequences (U -> T)
    csv_subset = csv_df[['ID', 'sgRNA', 'indel']].copy()
    csv_subset['sgRNA'] = csv_subset['sgRNA'].str.upper().str.replace('U', 'T')
    tsv_df['sgRNA'] = tsv_df['sgRNA'].str.upper().str.replace('U', 'T')

    # 4. Execute Left Join
    # Now ID and sgRNA formats match between both files
    updated_tsv_df = tsv_df.merge(csv_subset, on=['ID', 'sgRNA'], how='left')

    # 5. Diagnostic Check
    matched_count = updated_tsv_df['indel'].notna().sum()
    print(f"Successfully matched {matched_count} rows out of {len(tsv_df)}.")

    # 6. Write to disk
    updated_tsv_df.to_csv(output_filepath, sep='\t', index=False)

if __name__ == "__main__":
    update_tsv_with_spacer(
        tsv_filepath='entropy_advanced_strict.tsv', 
        csv_filepath='SpCas9_2019_150_sequence.csv', 
        output_filepath='entropy_advanced_strict_updated.tsv'
    )
