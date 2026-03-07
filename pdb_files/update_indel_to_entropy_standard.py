import os
import pandas as pd

def append_indel_data(ref_file, tsv_files):
    """Appends Indel data from a reference file to multiple target TSV files.

    This function reads a reference file (CSV or TSV) containing 'sgRNA' and 'Indel'
    columns. It then iterates through a list of target TSV files, loads each,
    merges the Indel data based on matching 'sequence(sgRNA)' values, and
    overwrites the original target TSV files with the updated data.

    Args:
        ref_file (str): The path to the reference file (CSV or TSV) containing
            the 'sgRNA' and 'Indel' columns.
        tsv_files (list[str]): A list of paths to the target TSV files to be
            updated with Indel data.
    """
    # 1. Load the reference file containing the Indel data
    if not os.path.exists(ref_file):
        print(f"Error: Reference file '{ref_file}' not found.")
        return
        
    print(f"Loading reference data from: {ref_file}...")
    
    # Detect separator based on file extension (handles both .csv and .tsv)
    sep = '\t' if ref_file.endswith('.tsv') else ','
    
    try:
        # Load only the required columns to save memory
        df_ref = pd.read_csv(ref_file, sep=sep, usecols=['sgRNA', 'Indel'])
    except ValueError as e:
        print(f"Error: Missing expected columns in reference file. {e}")
        return

    # Drop potential duplicates in the reference file to prevent row multiplication during merge
    df_ref = df_ref.drop_duplicates(subset=['sgRNA'])
    
    # Rename 'sgRNA' to 'sequence(sgRNA)' so the columns match perfectly for the merge
    df_ref = df_ref.rename(columns={'sgRNA': 'sequence(sgRNA)'})

    # 2. Iterate through your three generated TSV files and update them
    for tsv in tsv_files:
        if not os.path.exists(tsv):
            print(f"Warning: {tsv} not found. Skipping.")
            continue
            
        # Load the PDB extraction TSV
        df_target = pd.read_csv(tsv, sep='\t')
        
        # If the script was already run previously, drop the old Indel column to avoid duplicates
        if 'Indel' in df_target.columns:
            df_target = df_target.drop(columns=['Indel'])
        
        # Perform a Left Join: Keep all PDB data, attach Indel where sequences match
        df_merged = pd.merge(df_target, df_ref, on='sequence(sgRNA)', how='left')
        
        # Count how many sequences from the PDB files did NOT have a match in your reference file
        missing_count = df_merged['Indel'].isna().sum()
        
        # Fill missing values with "NA" for a clean text output
        df_merged['Indel'] = df_merged['Indel'].fillna("NA")
        
        # Overwrite the original TSV with the new merged data
        df_merged.to_csv(tsv, sep='\t', index=False)
        
        print(f"Updated: {tsv} | (Sequences lacking Indel matches: {missing_count})")

if __name__ == "__main__":
    # ---> CHANGE THIS to the actual name of your file containing the Indel column <---
    REFERENCE_FILE = "SpCas9_2020_sequence.csv" 
    
    TSV_FILES = [
        "entropy_averages_super_relaxed.tsv",
        "entropy_averages_relaxed.tsv",
        "entropy_averages_strict.tsv"
    ]
    
    append_indel_data(REFERENCE_FILE, TSV_FILES)
