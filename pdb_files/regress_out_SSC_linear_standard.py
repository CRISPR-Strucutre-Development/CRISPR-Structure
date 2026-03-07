import pandas as pd
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler

def regress_out_sequence(filepath):
    """Performs ordinary least squares (OLS) regression to evaluate the contribution of RNA structure metrics to 'Indel' rates, while controlling for 'SSC'.

    The function loads data from a specified TSV file, standardizes relevant features,
    and then fits a separate OLS model for each RNA structure metric. Each model
    predicts 'Indel' based on 'SSC' and one specific RNA structure metric.
    The coefficients and p-values for the structure metrics are extracted and
    printed in a sorted table.

    Args:
        filepath (str): The path to the tab-separated (TSV) data file containing
                        'Indel', 'SSC', and various RNA structure metrics.

    Returns:
        None: The function prints a formatted pandas DataFrame to the console
              showing the coefficients and p-values for each RNA structure metric,
              sorted by p-value.
    """
    # 1. Load data
    df = pd.read_csv(filepath, sep='\t')
    
    # 2. Isolate the variables
    y = df['Indel']
    ssc = df['SSC']
    
    # Define columns to ignore (Identifiers and the Sequence/Outcome columns)
    ignore_cols = ['ID', 'sgRNA', 'TargetSequence(RNAversion)', 'Indel', 'SSC', 
                   'Spacer.Sequence', 'Scaffold.Sequence']
    
    # Get all the RNA structure columns (TL_Mean, SL1_Max, etc.)
    structure_cols = [col for col in df.columns if col not in ignore_cols]
    
    # 3. Setup lists to store our results
    results_data = []
    
    # 4. Standardize the data to compare coefficients apples-to-apples
    scaler = StandardScaler()
    
    # 5. Iterate through each structure metric one by one
    for col in structure_cols:
        # Create a mini-dataframe with just SSC and the current structure metric
        X = df[['SSC', col]].dropna() 
        y_subset = y.loc[X.index] # Ensure y matches X after dropping any NaNs
        
        # Standardize the features
        X_scaled = pd.DataFrame(scaler.fit_transform(X), columns=['SSC', col], index=X.index)
        X_scaled = sm.add_constant(X_scaled) # Add the beta_0 intercept
        
        # Fit the model: Indel ~ SSC + Structure_Metric
        model = sm.OLS(y_subset, X_scaled).fit()
        
        # Extract the metrics specifically for the Structure variable
        struct_coef = model.params[col]
        struct_pval = model.pvalues[col]
        
        results_data.append({
            'Structure_Metric': col,
            'Structure_Coefficient': struct_coef,
            'P-Value': struct_pval
        })
        
    # 6. Format and sort the results
    results_df = pd.DataFrame(results_data)
    # Sort by the most statistically significant (lowest p-value)
    results_df = results_df.sort_values(by='P-Value').reset_index(drop=True)
    
    print("--- RNA Structure Contribution (Controlling for SSC) ---")
    print(results_df.to_string(index=False))

# Execution
regress_out_sequence('entropy_averages_strict_final.tsv')
