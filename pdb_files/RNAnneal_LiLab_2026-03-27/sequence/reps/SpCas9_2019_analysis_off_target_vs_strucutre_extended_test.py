import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, ks_2samp, fisher_exact

def run_extreme_tail_analysis(merged_csv_path):
    # Load consolidated dataframe
    df = pd.read_csv(merged_csv_path)
    df.columns = df.columns.str.strip()
    
    # Define strict extreme boundaries (excluding the middle 10% to 78% noise)
    high_condition = df['Indel_Rate'] > 78.0
    low_condition = df['Indel_Rate'] < 10.0
    
    # Filter dataset to isolate the tails
    df_tails = df[high_condition | low_condition].copy()
    
    # Assign categorical efficiency bins based on the new limits
    df_tails['Efficiency_Group'] = np.where(df_tails['Indel_Rate'] > 78.0, 'High_Efficiency', 'Low_Efficiency')
    
    # Maintain standard empirical off-target risk cutoff (CFD >= 0.20)
    OFFTARGET_CUTOFF = 0.20
    df_tails['Risk_Group'] = np.where(df_tails['Max_True_OffTarget_Score'] >= OFFTARGET_CUTOFF, 'High_Risk', 'Low_Risk')
    
    high_eff_scores = df_tails[df_tails['Efficiency_Group'] == 'High_Efficiency']['Max_True_OffTarget_Score'].dropna()
    low_eff_scores = df_tails[df_tails['Efficiency_Group'] == 'Low_Efficiency']['Max_True_OffTarget_Score'].dropna()
    
    print("="*80)
    print("   EXTREME TAILS STATISTICAL PROFILE FOR THE 2019 SCAFFOLD CRISIS")
    print("="*80)
    print(f"Sample Sizes: High-Efficiency (>78.0%): {len(high_eff_scores)} guides")
    print(f"              Low-Efficiency  (<10.0%): {len(low_eff_scores)} guides\n")
    
    if len(high_eff_scores) == 0 or len(low_eff_scores) == 0:
        print("CRITICAL ERROR: One of the extreme groups has 0 samples with these cutoffs.")
        print("Verify your data ranges before proceeding.")
        return

    # 1. MANN-WHITNEY U TEST (Medians)
    u_stat, mw_p = mannwhitneyu(high_eff_scores, low_eff_scores, alternative='two-sided')
    print(f"[1] MANN-WHITNEY U TEST (Comparing Off-Target Medians between Extremes)")
    print(f"    U-Statistic: {u_stat:.2f}")
    print(f"    P-Value    : {mw_p:.6f}")
    if mw_p > 0.05:
        print("    STATUS     : NON-SIGNIFICANT. Even when comparing hyper-active guides to near-dead")
        print("                 guides, the median off-target risk remains completely indistinguishable.")
    else:
        print("    STATUS     : SIGNIFICANT. A true localized median shift exists between the tails.")
    print("-"*80)

    # 2. KOLMOGOROV-SMIRNOV (KS) TEST (Distribution Shapes)
    ks_stat, ks_p = ks_2samp(high_eff_scores, low_eff_scores)
    print(f"[2] KOLMOGOROV-SMIRNOV (KS) TEST (Comparing Extreme Risk Distribution Shapes)")
    print(f"    KS-Statistic: {ks_stat:.4f}")
    print(f"    P-Value     : {ks_p:.6f}")
    if ks_p > 0.05:
        print("    STATUS      : NON-SIGNIFICANT. The cumulative distribution shape of off-target risk")
        print("                  is statistically identical between both extreme populations.")
    else:
        print("    STATUS      : SIGNIFICANT. The global probability distributions diverge.")
    print("-"*80)

    # 3. FISHER'S EXACT TEST (Categorical Association)
    contingency_matrix = pd.crosstab(df_tails['Efficiency_Group'], df_tails['Risk_Group'])
    
    # Explicitly re-index to force a full 2x2 layout if any cell is zero
    for group in ['High_Efficiency', 'Low_Efficiency']:
        if group not in contingency_matrix.index:
            contingency_matrix.loc[group] = [0, 0]
    for risk in ['High_Risk', 'Low_Risk']:
        if risk not in contingency_matrix.columns:
            contingency_matrix[risk] = 0
    contingency_matrix = contingency_matrix.loc[['High_Efficiency', 'Low_Efficiency'], ['High_Risk', 'Low_Risk']]
    
    print("[3] FISHER'S EXACT TEST CONTINGENCY TABLE (EXTREMES ONLY):")
    print(contingency_matrix)
    print("")
    
    odds_ratio, fisher_p = fisher_exact(contingency_matrix.values)
    print(f"    Odds Ratio : {odds_ratio:.4f}")
    print(f"    P-Value    : {fisher_p:.6f}")
    if fisher_p > 0.05:
        print("    STATUS     : NON-SIGNIFICANT. Hyper-efficiency vs. complete failure does not associate")
        print("                 with off-target risk. Off-target propensity is an inherent structural flaw.")
    else:
        print("    STATUS     : SIGNIFICANT. Extreme efficiency boundaries alter risk proportions.")
    print("="*80)

if __name__ == "__main__":
    run_extreme_tail_analysis('merged_structural_offtarget_analysis.csv')