def merge_sequence_files(file1_path, file2_path, output_path):
    # Load datasets
    df1 = pd.read_csv(file1_path)
    df2 = pd.read_csv(file2_path)

    # Logic: Normalize RNA (U) to DNA (T) for comparison
    # This ensures 'AUGC' becomes 'ATGC' to match File 2
    df1['sequence'] = df1['sequence'].str.replace('U', 'T', case=True)

    # Merge based on the 'sequence' column
    # 'inner' join ensures only matching sequences are kept
    merged_df = pd.merge(df1, df2, on='sequence', how='inner')

    # Export final result
    merged_df.to_csv(output_path, index=False)
    print(f"Success: Merged file saved to {output_path}")

# Usage
# merge_sequence_files('file1.csv', 'file2.csv', 'merged_output.csv')
