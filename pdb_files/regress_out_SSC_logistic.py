import pandas as pd
import statsmodels.api as sm
import numpy as np
from sklearn.preprocessing import StandardScaler

def run_logistic_regression(filepath):
    df = pd.read_csv(filepath, sep='\t')
    
    # 1. Binarize the Dependent Variable (Indel)
    # 1 = Active (Indel > 10%), 0 = Inactive (Indel <= 10%)
    df['Indel_Binary'] = (df['Indel'] > 10).astype(int)
    
    y = df['Indel_Binary']
    
    # 2. Define features to ignore
    ignore_cols = ['ID', 'sgRNA', 'TargetSequence(RNAversion)', 'Indel', 'Indel_Binary', 'SSC', 
                   'Spacer.Sequence', 'Scaffold.Sequence']
    structure_cols = [col for col in df.columns if col not in ignore_cols]
    
    results_data = []
    scaler = StandardScaler()
    
    # 3. Iterate and run Logistic Regression (Logit)
    for col in structure_cols:
        # Isolate SSC and the target structure metric
        X = df[['SSC', col]].dropna()
        y_subset = y.loc[X.index]
        
        # Standardize features for comparable coefficients
        X_scaled = pd.DataFrame(scaler.fit_transform(X), columns=['SSC', col], index=X.index)
        X_scaled = sm.add_constant(X_scaled)
        
        try:
            # Fit the Logistic Regression model
            # disp=0 suppresses the iterative solver output text
            model = sm.Logit(y_subset, X_scaled).fit(disp=0) 
            
            struct_coef = model.params[col]
            struct_pval = model.pvalues[col]
            
            # In logistic regression, exponentiating the coefficient gives the Odds Ratio
            odds_ratio = np.exp(struct_coef)
            
            results_data.append({
                'Structure_Metric': col,
                'Coefficient (Log-Odds)': struct_coef,
                'Odds Ratio': odds_ratio,
                'P-Value': struct_pval
            })
            
        except np.linalg.LinAlgError:
            # Catches perfect separation errors common in very small sample sizes
            pass 

    # 4. Format and display
    results_df = pd.DataFrame(results_data)
    results_df = results_df.sort_values(by='P-Value').reset_index(drop=True)
    
    print("--- Logistic Regression: Structure Metrics (Controlling for SSC) ---")
    print(results_df.to_string(index=False))

# Execution
run_logistic_regression('entropy_averages_strict_final.tsv')
