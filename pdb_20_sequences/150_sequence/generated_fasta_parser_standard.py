import pandas as pd

def main():
    """Processes CRISPR-Cas9 indel data, performs stratified sampling, and preprocesses sequences.

    Loads indel data from 'SpCas9_Indel_2020.csv', filters sequences based on length
    and PAM motif (NGG), and performs a value-based stratified sampling to select
    150 sequences. The selected sequences are then processed by extracting guide
    sequences and converting DNA (T) to RNA (U) bases. Finally, the processed
    sequences, indel values, and guide sequences are saved to
    'SpCas9_2020_150_sequence.csv'.

    The script exits prematurely if fewer than 150 sequences remain after initial
    PAM filtering.
    """
    # Load data
    df = pd.read_csv('SpCas9_Indel_2020.csv')

    # PAM filtering
    df_filtered = df[(df['sequence'].str.len() >= 27) & (df['sequence'].str.slice(25, 27).str.upper() == 'GG')].copy()

    # Sequence counting and check
    num_sequences = len(df_filtered)
    print(f'Number of sequences present after NGG filtering: {num_sequences}')
    if num_sequences < 150:
        exit()

    # Value-based sampling and top-off
    df_filtered['bin'] = pd.cut(df_filtered['indel'], bins=4, include_lowest=True)
    targets = [45, 22, 38, 45]
    sampled_df = pd.DataFrame()
    for target, cat in zip(targets, df_filtered['bin'].cat.categories):
        bin_df = df_filtered[df_filtered['bin'] == cat]
        if len(bin_df) == 0:
            continue
        try:
            sampled_df = pd.concat([sampled_df, bin_df.sample(n=target, replace=False)])
        except ValueError:
            sampled_df = pd.concat([sampled_df, bin_df.sample(n=target, replace=True)])

    # Adjust sample size if needed
    if len(sampled_df) < 150:
        deficit = 150 - len(sampled_df)
        sampled_df = pd.concat([sampled_df, df_filtered.sample(n=deficit, replace=False)])

    # Sequence processing
    sampled_df['guide'] = sampled_df['sequence'].str.slice(4, 24)
    sampled_df['sequence'] = sampled_df['sequence'].str.replace('T', 'U', case=False)
    sampled_df['guide'] = sampled_df['guide'].str.replace('T', 'U', case=False)

    # Output
    sampled_df[['sequence', 'indel', 'guide']].to_csv('SpCas9_2020_150_sequence.csv', index=False)

if __name__ == "__main__":
    main()
