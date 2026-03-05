import pandas as pd

def main():
    # Load data
    df = pd.read_csv('SpCas9_Indel_2020.csv')

    # PAM Filtering
    df_filtered = df[(df['Sequence'].str.len() >= 27) & (df['Sequence'].str.slice(25, 27).str.upper() == 'GG')].copy()

    # Sequence Counting
    print(f"Number of sequences present after NGG filtering: {len(df_filtered)}")
    if len(df_filtered) < 150:
        exit()

    # Value-Based Sampling & Top-Off
    df_filtered['bin'] = pd.cut(df_filtered['indel'], bins=4, include_lowest=True)
    targets = [45, 22, 38, 45] 
    sampled_df = pd.DataFrame()
    for target, cat in zip(targets, df_filtered['bin'].cat.categories):
        bin_df = df_filtered[df_filtered['bin'] == cat]
        if len(bin_df) == 0: continue
        try:
            sampled_df = pd.concat([sampled_df, bin_df.sample(n=target, replace=False)])
        except ValueError:
            sampled_df = pd.concat([sampled_df, bin_df.sample(n=target, replace=True)])
    if len(sampled_df) < 150:
        deficit = 150 - len(sampled_df)
        sampled_df = pd.concat([sampled_df, df_filtered.sample(n=deficit, replace=False)])

    # Sequence Processing
    sampled_df['guide'] = sampled_df['Sequence'].str.slice(4, 24)
    sampled_df['Sequence'] = sampled_df['Sequence'].str.replace('T', 'U', case=False)
    sampled_df['guide'] = sampled_df['guide'].str.replace('T', 'U', case=False)

    # Output Formatting
    sampled_df[['Sequence', 'indel', 'guide']].to_csv('SpCas9_2020_150_sequence.csv', index=False)

if __name__ == "__main__":
    main()
