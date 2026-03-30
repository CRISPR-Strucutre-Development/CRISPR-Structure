import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import chi2

def load_and_preprocess_data(filepath):
    """Loads data from a TSV file and preprocesses it.

    Args:
        filepath: The path to the TSV file.

    Returns:
        A pandas DataFrame with an added 'Activity_State' column.
    """
    df = pd.read_csv(filepath, sep='\t')
    # Binarize Activity
    df['Activity_State'] = (df['Indel'] > 30).astype(int)
    return df

def build_nested_models(df):
    """Builds and fits a hierarchy of nested logistic regression models.

    Args:
        df: The input DataFrame containing 'Activity_State', 'SSC', 'Ent_Max', and 'DG_UNFOLD'.

    Returns:
        A tuple containing the four fitted GLM models: (model_1, model_2, model_3, model_4).
    """
    # Model 1: Sequence Only (The baseline)
    model_1 = smf.glm("Activity_State ~ SSC", data=df, family=sm.families.Binomial()).fit()

    # Model 2: Sequence + Local Structure (Local flexibility veto)
    model_2 = smf.glm("Activity_State ~ SSC + Ent_Max", data=df, family=sm.families.Binomial()).fit()

    # Model 3: Sequence + Local Structure + Global Thermodynamics (Full Model)
    model_3 = smf.glm("Activity_State ~ SSC + Ent_Max + DG_UNFOLD", data=df, family=sm.families.Binomial()).fit()

    # Model 4: Sequence + Global Thermodynamics (Direct comparison)
    model_4 = smf.glm("Activity_State ~ SSC + DG_UNFOLD", data=df, family=sm.families.Binomial()).fit()

    return model_1, model_2, model_3, model_4

def perform_likelihood_ratio_tests(model_1, model_2, model_3, model_4):
    """Performs likelihood ratio tests between nested models.

    Args:
        model_1: The baseline GLM model (Sequence Only).
        model_2: GLM model with Sequence + Local Structure.
        model_3: GLM model with Sequence + Local Structure + Global Thermodynamics.
        model_4: GLM model with Sequence + Global Thermodynamics.

    Returns:
        A dictionary containing the p-values for each test.
        Keys are 'lrt_ent', 'lrt_dg_added', 'lrt_dg_only'.
    """
    results = {}

    # Test A: Does Local Entropy (Ent_Max) improve the Sequence model?
    # H0: model_1 (SSC) vs H1: model_2 (SSC + Ent_Max)
    ll_1 = model_1.llf
    ll_2 = model_2.llf
    # df_diff is the number of additional parameters in the unrestricted model (model_2)
    df_diff_a = len(model_2.params) - len(model_1.params)
    lrt_stat_a = -2 * (ll_1 - ll_2)
    results['lrt_ent'] = chi2.sf(lrt_stat_a, df_diff_a)

    # Test B: Does Global Stability (DG_UNFOLD) improve the Sequence + Entropy model?
    # H0: model_2 (SSC + Ent_Max) vs H1: model_3 (SSC + Ent_Max + DG_UNFOLD)
    ll_3 = model_3.llf
    # df_diff is the number of additional parameters in the unrestricted model (model_3)
    df_diff_b = len(model_3.params) - len(model_2.params)
    lrt_stat_b = -2 * (ll_2 - ll_3)
    results['lrt_dg_added'] = chi2.sf(lrt_stat_b, df_diff_b)

    # Test C: Does Global Stability (DG_UNFOLD) improve the Sequence model alone?
    # H0: model_1 (SSC) vs H1: model_4 (SSC + DG_UNFOLD)
    ll_4 = model_4.llf
    # df_diff is the number of additional parameters in the unrestricted model (model_4)
    df_diff_c = len(model_4.params) - len(model_1.params)
    lrt_stat_c = -2 * (ll_1 - ll_4)
    results['lrt_dg_only'] = chi2.sf(lrt_stat_c, df_diff_c)

    return results

def print_analysis_results(model_1, model_2, model_3, model_4, lrt_p_values):
    """Prints the log-likelihoods of the models and the p-values from LRTs.

    Args:
        model_1: The baseline GLM model (Sequence Only).
        model_2: GLM model with Sequence + Local Structure.
        model_3: GLM model with Sequence + Local Structure + Global Thermodynamics.
        model_4: GLM model with Sequence + Global Thermodynamics.
        lrt_p_values: A dictionary of p-values from the likelihood ratio tests.
                      Expected keys: 'lrt_ent', 'lrt_dg_added', 'lrt_dg_only'.
    """
    print("\n=======================================================")
    print("   NESTED LOGISTIC LIKELIHOOD RATIO ANALYSIS")
    print("=======================================================")
    print(f"1. SSC Baseline LL:          {model_1.llf:.3f}")
    print(f"2. SSC + Ent_Max LL:         {model_2.llf:.3f} (p = {lrt_p_values['lrt_ent']:.5f})")
    print(f"3. SSC + Ent + DG_UNFOLD LL: {model_3.llf:.3f} (p = {lrt_p_values['lrt_dg_added']:.5f})")
    print(f"4. SSC + DG_UNFOLD LL:       {model_4.llf:.3f} (p = {lrt_p_values['lrt_dg_only']:.5f})")
    print("=======================================================")

# Main execution
# Define the path to the data file
data_filepath = "entropy_cumulative_stats_update_final.tsv"

# 1. Load Data
df = load_and_preprocess_data(data_filepath)

# 2. Build the Nested Model Hierarchy
model_1, model_2, model_3, model_4 = build_nested_models(df)

# 3. Perform Likelihood Ratio Tests
lrt_p_values = perform_likelihood_ratio_tests(model_1, model_2, model_3, model_4)

# 4. Print Results
print_analysis_results(model_1, model_2, model_3, model_4, lrt_p_values)
