import pandas as pd
import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler

def analyze_and_plot(filepath):
    """Analyzes CRISPR-Cas9 indel data and generates diagnostic plots for model fits.

    This function loads experimental data from a TSV file, preprocesses it to
    create a binary indel outcome, and then fits ordinary least squares (OLS)
    and logistic regression models for various structural metrics. It identifies
    the 'best' metric for both linear and logistic models based on the lowest
    p-value and generates a comprehensive 2x2 diagnostic figure. The figure
    includes a linear fit plot, a residuals plot, a logistic fit plot, and a
    volcano plot summarizing all linear model results. The generated plots are
    saved as high-resolution PNG and PDF files.

    Args:
        filepath (str): The path to the tab-separated (TSV) CSV file containing
            experimental data. This file is expected to include 'Indel' percentages,
            'SSC' scores, and various structural metrics.

    Returns:
        None: This function saves diagnostic plots to disk and prints save messages,
            but does not return any values.

    Raises:
        FileNotFoundError: If the `filepath` does not point to an existing file.
        pd.errors.EmptyDataError: If the `filepath` points to an empty file.
        KeyError: If essential columns like 'Indel' or 'SSC' are missing from the dataframe.
    """
    # 1. Load Data
    df = pd.read_csv(filepath, sep='\t')
    df['Indel_Binary'] = (df['Indel'] > 10).astype(int)
    
    ignore_cols = ['ID', 'sgRNA', 'TargetSequence(RNAversion)', 'Indel', 'Indel_Binary', 'SSC', 
                   'Spacer.Sequence', 'Scaffold.Sequence']
    structure_cols = [col for col in df.columns if col not in ignore_cols]
    
    scaler = StandardScaler()
    
    lin_results, log_results = [], []
    best_lin_model, best_log_model = None, None
    best_lin_metric, best_log_metric = None, None
    lowest_lin_pval, lowest_log_pval = float('inf'), float('inf')
    best_lin_data, best_log_data = None, None

    # 2. Fit Models and Find the Best Metric
    for col in structure_cols:
        X = df[['SSC', col]].dropna()
        y_lin = df['Indel'].loc[X.index]
        y_log = df['Indel_Binary'].loc[X.index]
        
        X_scaled = pd.DataFrame(scaler.fit_transform(X), columns=['SSC', col], index=X.index)
        X_scaled = sm.add_constant(X_scaled)
        
        # --- Linear ---
        lin_model = sm.OLS(y_lin, X_scaled).fit()
        pval = lin_model.pvalues[col]
        lin_results.append({'Metric': col, 'Coef': lin_model.params[col], 'P-Value': pval})
        
        if pval < lowest_lin_pval:
            lowest_lin_pval = pval
            best_lin_metric = col
            best_lin_model = lin_model
            best_lin_data = (X_scaled, y_lin, X[col]) # Store raw feature for plotting

        # --- Logistic ---
        try:
            log_model = sm.Logit(y_log, X_scaled).fit(disp=0)
            log_pval = log_model.pvalues[col]
            log_results.append({'Metric': col, 'Coef': log_model.params[col], 'P-Value': log_pval})
            
            if log_pval < lowest_log_pval:
                lowest_log_pval = log_pval
                best_log_metric = col
                best_log_model = log_model
                best_log_data = (X_scaled, y_log, X[col])
        except np.linalg.LinAlgError:
            pass

    # 3. Generate Diagnostics Figure for the PI
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    sns.set_theme(style="whitegrid")

    # --- Plot 1: Linear Fit (Best Metric) ---
    ax = axes[0, 0]
    sns.regplot(x=best_lin_data[2], y=best_lin_data[1], ax=ax, scatter_kws={'alpha':0.5}, line_kws={'color':'red'})
    ax.set_title(f'Linear Fit: Indel vs {best_lin_metric}\n(OLS line in red)', fontsize=12, fontweight='bold')
    ax.set_xlabel(f'Raw {best_lin_metric}')
    ax.set_ylabel('Indel %')

    # --- Plot 2: Linear Residuals (The Issue Checker) ---
    ax = axes[0, 1]
    fitted_vals = best_lin_model.fittedvalues
    residuals = best_lin_model.resid
    sns.scatterplot(x=fitted_vals, y=residuals, ax=ax, alpha=0.5)
    ax.axhline(0, color='red', linestyle='--')
    ax.set_title(f'Residuals vs Fitted (OLS)\nCheck for non-random patterns/funneling', fontsize=12, fontweight='bold')
    ax.set_xlabel('Fitted Indel Values')
    ax.set_ylabel('Residuals')

    # --- Plot 3: Logistic Fit (Best Metric) ---
    ax = axes[1, 0]
    # Generate smooth curve for predicted probabilities
    X_log_scaled, y_log, raw_log_feat = best_log_data
    x_range = np.linspace(X_log_scaled[best_log_metric].min(), X_log_scaled[best_log_metric].max(), 100)
    # Hold SSC at its mean (0 since standardized) to plot the isolated effect of the metric
    X_pred = sm.add_constant(pd.DataFrame({'SSC': np.zeros(100), best_log_metric: x_range}))
    y_pred_prob = best_log_model.predict(X_pred)
    
    # Map scaled x_range back to raw for intuitive x-axis
    raw_x_range = np.linspace(raw_log_feat.min(), raw_log_feat.max(), 100)
    
    sns.scatterplot(x=raw_log_feat, y=y_log, ax=ax, alpha=0.5, label='Actual (0 or 1)')
    ax.plot(raw_x_range, y_pred_prob, color='red', label='Predicted Probability')
    ax.set_title(f'Logistic Fit: Probability of Indel > 10% vs {best_log_metric}', fontsize=12, fontweight='bold')
    ax.set_xlabel(f'Raw {best_log_metric}')
    ax.set_ylabel('Probability / Actual Class')
    ax.legend()

    # --- Plot 4: Summary Volcano Plot of All Metrics ---
    ax = axes[1, 1]
    lin_df = pd.DataFrame(lin_results)
    lin_df['-log10(P-Value)'] = -np.log10(lin_df['P-Value'])
    sns.scatterplot(data=lin_df, x='Coef', y='-log10(P-Value)', ax=ax, color='blue', alpha=0.7)
    
    # Highlight top 3 most significant
    top_3 = lin_df.nsmallest(3, 'P-Value')
    for _, row in top_3.iterrows():
        ax.text(row['Coef'], row['-log10(P-Value)'], row['Metric'], fontsize=9, verticalalignment='bottom')
        
    ax.axhline(-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
    ax.set_title('Linear Model Summary: Effect Size vs Significance', fontsize=12, fontweight='bold')
    ax.set_xlabel('Coefficient (Standardized)')
    ax.set_ylabel('-log10(P-Value)')
    ax.legend()

    plt.tight_layout()
    
    # --- HPC Save Execution ---
    # Save as high-resolution PNG
    plt.savefig('model_diagnostics.png', dpi=300, bbox_inches='tight')
    print("Saved model_diagnostics.png")
    
    # Save as vector PDF
    plt.savefig('model_diagnostics.pdf', format='pdf', bbox_inches='tight')
    print("Saved model_diagnostics.pdf")
    
    # Close the figure to free up memory on the compute node
    plt.close(fig)

# Execution
analyze_and_plot('entropy_averages_strict_final.tsv')
