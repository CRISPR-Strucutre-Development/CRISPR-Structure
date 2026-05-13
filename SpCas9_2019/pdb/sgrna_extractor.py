import pandas as pd

# Load the input CSV file
# Replace 'your_input_file.csv' with the actual path to your file
df = pd.read_csv('SpCas9_2019_150_sequence.csv')

# Generate the list of IDs you want to extract: sgRNA_1 to sgRNA_5 and sgRNA_137 to sgRNA_141
target_ids = [f'sgRNA_{i}' for i in range(1, 6)] + [f'sgRNA_{i}' for i in range(137, 142)]

# Filter the DataFrame for these specific IDs
filtered_df = df[df['ID'].isin(target_ids)]

# Select only the required columns
output_df = filtered_df[['sgRNA', 'ID']]

# Save the filtered data to a new CSV file
output_df.to_csv('output_file.csv', index=False)

print("Extraction complete. The new file 'output_file.csv' has been saved.")
