import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

# 1. Load your generated files
dist_df = pd.read_csv('sgrna_interactions_distance.csv', na_values=['NA'])
pair_df = pd.read_csv('sgrna_interactions_pair.csv')
stack_df = pd.read_csv('sgrna_interactions_stack.csv')
indel_df = pd.read_csv('experimental_indels.csv') # Ensure columns: 'sgRNA', 'Indel_Rate'

# 2. Aggregation: Find the Structural Bottleneck
# Take MIN distance (the tightest spot) across Ranks 01-05
dist_summary = dist_df.groupby('sgRNA')[['SL1_54_60', 'SL2_71_76']].min().reset_index()
dist_summary['Min_Hinge_Dist'] = dist_summary[['SL1_54_60', 'SL2_71_76']].min(axis=1)

# Take MAX Interaction (1 if a 'Lock' appears in any of the top 5 ranks)
pair_summary = pair_df.groupby('sgRNA')[['SL1_54_60', 'SL2_71_76']].max().reset_index()
pair_summary['Any_Pairing'] = pair_summary[['SL1_54_60', 'SL2_71_76']].max(axis=1)

stack_summary = stack_df.groupby('sgRNA')[['SL1_54_60', 'SL2_71_76']].max().reset_index()
stack_summary['Any_Stacking'] = stack_summary[['SL1_54_60', 'SL2_71_76']].max(axis=1)

# 3. Merge with Experimental Data
data = indel_df.merge(dist_summary, on='sgRNA')\
               .merge(pair_summary, on='sgRNA')\
               .merge(stack_summary, on='sgRNA')
data = data.dropna(subset=['Min_Hinge_Dist', 'Indel_Rate'])

# 4. Calculate Spearman Rho and p-values
print("--- Statistical Results ---")
for feat in ['Min_Hinge_Dist', 'Any_Pairing', 'Any_Stacking']:
    rho, p = spearmanr(data[feat], data['Indel_Rate'])
    print(f"{feat:15} | Rho: {rho:6.3f} | p-value: {p:.2e}")

# 5. Visualization for PI Presentation
plt.figure(figsize=(10, 6))
sns.regplot(data=data, x='Min_Hinge_Dist', y='Indel_Rate', 
            scatter_kws={'alpha':0.5, 'color':'teal'}, line_kws={'color':'red'})
plt.title('Hinge Accessibility vs. Editing Efficiency (n=150)', fontsize=14)
plt.xlabel('Min Heavy-Atom Distance (Å) at SL1/SL2', fontsize=12)
plt.ylabel('Experimental Indel %', fontsize=12)
plt.savefig('Correlation_Distance_Indels.png', dpi=300)

plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
sns.boxplot(data=data, x='Any_Pairing', y='Indel_Rate', palette='vlag')
plt.title('Impact of Spacer-Scaffold Pairing')
plt.subplot(1, 2, 2)
sns.boxplot(data=data, x='Any_Stacking', y='Indel_Rate', palette='vlag')
plt.title('Impact of Spacer-Scaffold Stacking')
plt.savefig('Interaction_Impact_Boxplots.png', dpi=300)
