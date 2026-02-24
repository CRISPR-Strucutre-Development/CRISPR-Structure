import pandas as pd

def append_ssc_to_tsv(updated_tsv_filepath, ssc_filepath, output_filepath):
    # 1. Load the existing TSV and the SSC dataset
    tsv_df = pd.read_csv(updated_tsv_filepath, sep='\t')
    
    # .out files frequently use variable whitespace rather than strict tabs.
    # sep='\s+' handles any combination of spaces or tabs automatically.
    ssc_df = pd.read_csv(ssc_filepath, sep='\s+') 

    # 2. Isolate required columns from the SSC file
    ssc_subset = ssc_df[['TargetSequence(RNAversion)', 'SSC']].copy()

    # 3. Standardize the merge key in both dataframes (Defensive Programming)
    # Applying this to both ensures a 100% match rate, even if the upstream 
    # TSV wasn't fully standardized in a previous step.
    tsv_df['TargetSequence(RNAversion)'] = tsv_df['TargetSequence(RNAversion)'].str.upper().str.replace('U', 'T')
    ssc_subset['TargetSequence(RNAversion)'] = ssc_subset['TargetSequence(RNAversion)'].str.upper().str.replace('U', 'T')

    # 4. Execute Left Join
    final_df = tsv_df.merge(ssc_subset, on='TargetSequence(RNAversion)', how='left')

    # 5. Write to disk
    final_df.to_csv(output_filepath, sep='\t', index=False)

# Execution
append_ssc_to_tsv(
    updated_tsv_filepath='entropy_cumulative_stats_updated.tsv', 
    ssc_filepath='SSC_sequence.out',
    output_filepath='entropy_cumulative_stats_update_final.tsv' # Or overwrite the updated one
)
