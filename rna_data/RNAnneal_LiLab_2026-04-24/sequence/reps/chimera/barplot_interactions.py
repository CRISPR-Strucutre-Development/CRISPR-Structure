import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 1. Load your data
df = pd.read_csv("sgrna_interactions_pair.csv")

# 2. Reshape data from wide to long format
val_cols = ['TL_37_42', 'SL1_64_70', 'SL2_82_87', 'SL3_98_100']
df_melted = df.melt(id_vars=['sgRNA', 'Rank'], 
                    value_vars=val_cols, 
                    var_name='Region', 
                    value_name='Binary_Value')

# 3. Clean Region names (remove the underscore and numbers)
df_melted['Region'] = df_melted['Region'].str.split('_').str[0]

# 4. Generate the Plot
plt.figure(figsize=(10, 6))
sns.barplot(data=df_melted, 
            y='Region', 
            x='Binary_Value', 
            errorbar=('ci', 95), 
            capsize=.1, 
            palette='magma')

plt.xlabel('Proportion of 1s (Mean Density)')
plt.ylabel('Structural Region')
plt.title('Comparison of sgRNA Feature Density across Regions')
plt.xlim(0, 1)
plt.grid(axis='x', linestyle='--', alpha=0.6)

# SAVE the file instead of plt.show()
plt.savefig('base_pair_comparison_plot.png', dpi=300, bbox_inches='tight')
print("Plot saved as sgrna_comparison_plot.png")
