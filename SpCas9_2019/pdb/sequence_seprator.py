import pandas as pd

def extract_sequence_to_seq(input_filepath, output_filepath):
    # 1. Load the updated dataset
    df = pd.read_csv(input_filepath, sep='\t')

    # 2. Isolate the final column
    # .iloc stands for integer-location based indexing. 
    # [:, -1] translates to: [all rows, strictly the last column]
    final_column = df.iloc[:, -1]

    # 3. Write to .seq file
    # header=False prevents writing the column name at the top of the file
    # index=False prevents writing row numbers
    final_column.to_csv(output_filepath, index=False, header=False)

# Execution
extract_sequence_to_seq(
    input_filepath='entropy_averages_strict_updated.tsv',
    output_filepath='SSC_sequence.seq'
)
