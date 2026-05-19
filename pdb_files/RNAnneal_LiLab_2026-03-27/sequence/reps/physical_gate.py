import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

def calculate_effective_entropy(entropy_csv, pair_csv, stack_csv, indel_csv):
    # 1. Load Data
    df_s = pd.read_csv(entropy_csv)
    df_p = pd.read_csv(pair_csv)
    df_st = pd.read_csv(stack_csv)
    df_i = pd.read_csv(indel_csv)
    
    # 2. Pivot Entropy to have regions as columns per sgRNA/Rank
    df_s_pivot = df_s.pivot(index=['sgRNA', 'Rank'], columns='Region', values='Entropy').reset_index()
    
    # 3. Merge with Interactions (Ir)
    # Ir = 1 if pairing OR stacking exists
    data = df_s_pivot.merge(df_p, on=['sgRNA', 'Rank'], suffixes=('', '_p'))
    data = data.merge(df_st, on=['sgRNA', 'Rank'], suffixes=('', '_st'))

    # 4. Apply the Equation: (S_spacer + S_hinge) * (1 - Ir)
    target_hinges = {'TL': 'TL_32_37', 'SL1': 'SL1_54_60', 'SL2': 'SL2_71_76', 'SL3': 'SL3_88_90'}
    
    for short_name, col_name in target_hinges.items():
        # Binary Gate: Ir = 1 if either pairing or stacking is detected
        ir = data[[f'{col_name}', f'{col_name}_st']].max(axis=1)
        
        # Calculate Term: (S_spacer + S_hinge) * (1 - Ir)
        data[f'term_{short_name}'] = (data['spacer'] + data[short_name]) * (1 - ir)

    # 5. Average across Ranks (1/R * Sum)
    ensemble_results = data.groupby('sgRNA')[[f'term_{h}' for h in target_hinges]].mean().reset_index()
    ensemble_results.columns = ['sgRNA'] + [f'S_eff_{h}' for h in target_hinges]

    # 6. Final Comparison with Indel Tails
    final_df = ensemble_results.merge(df_i, on='sgRNA')
    high = final_df[final_df['Indel_Rate'] > 80]
    low = final_df[final_df['Indel_Rate'] < 15]

    for h in target_hinges:
        col = f'S_eff_{h}'
        _, p = mannwhitneyu(high[col], low[col])
        print(f"Region {h} | Effective Entropy p-value: {p:.4e}")

        plt.figure(figsize=(6, 5))
        sns.boxplot(data=pd.concat([high, low]), x='Indel_Rate', y=col, palette='vlag')
        plt.title(f'Effective Entropy Analysis: {h}\np = {p:.4f}')
        plt.savefig(f'Final_Seff_{h}.png')

if __name__ == "__main__":
    calculate_effective_entropy('pdb_entropy_extracted.csv', 
                                'sgrna_interactions_pair.csv', 
                                'sgrna_interactions_stack.csv', 
                                'experimental_indels.csv')
