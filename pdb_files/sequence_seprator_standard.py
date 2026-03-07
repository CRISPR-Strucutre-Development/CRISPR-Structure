import pandas as pd

def extract_sequence_to_seq(input_filepath, output_filepath):
    """Extracts the last column from a TSV file and saves it as a sequence file.

    Loads a tab-separated values (TSV) file, isolates its final column, and writes
    the content of this column to a new file. The output file will not contain
    any headers or row indices, making it suitable for sequence data where only
    the raw column values are needed.

    Args:
        input_filepath (str): The path to the input TSV file.
        output_filepath (str): The path where the extracted sequence will be saved
            (e.g., a `.seq` file).

    Returns:
        None: The function writes the data directly to the `output_filepath` as a
            side effect.
    """
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
