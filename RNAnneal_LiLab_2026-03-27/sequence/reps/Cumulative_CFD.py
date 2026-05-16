import pandas as pd
import numpy as np

def recalculate_true_offtarget_risk(scored_csv, output_csv):
    print(f"Loading data from {scored_csv}...")
    df = pd.read_csv(scored_csv)
    
    # CRITICAL STEP: Filter out the on-target sites (mismatches must be > 0)
    # This removes the 1.0 CFD scores belonging to the intended targets
    df_off_only = df[df['Mismatches'] > 0]
    
    print(f"Filtered out on-target sites. Processing {len(df_off_only)} true off-target alignments...")
    
    grouped = df_off_only.groupby('sgRNA_ID')
    summary_records = []
    
    for sgRNA_id, group in grouped:
        cfd_scores = group['CFD_Score'].values
        
        # Calculate Elevation probability of cutting AT LEAST ONE TRUE off-target site
        prob_no_cut = 1.0 - cfd_scores
        prob_miss_all = np.prod(prob_no_cut)
        cumulative_risk = 1.0 - prob_miss_all
        
        total_off_sites = len(group)
        high_risk_sites = np.sum(cfd_scores >= 0.2)
        max_true_offtarget = np.max(cfd_scores) if len(cfd_scores) > 0 else 0
        
        summary_records.append({
            'sgRNA_ID': sgRNA_id,
            'True_OffTarget_Matches': total_off_sites,
            'High_Risk_Sites_Count': high_risk_sites,
            'Max_True_OffTarget_Score': max_true_offtarget,
            'Cumulative_OffTarget_Risk': cumulative_risk
        })
        
    summary_df = pd.DataFrame(summary_records)
    
    # Sort by cumulative risk to isolate the most promiscuous guides
    summary_df = summary_df.sort_values(by='Cumulative_OffTarget_Risk', ascending=False)
    
    summary_df.to_csv(output_csv, index=False)
    print(f"SUCCESS: Corrected summary saved to '{output_csv}'")
    
    print("\n" + "="*80)
    print("      CORRECTED ELEVATION CUMULATIVE RISK RANKING (ON-TARGETS EXCLUDED)")
    print("="*80)
    print(summary_df.head(15).to_string(index=False))
    print("="*80)

if __name__ == "__main__":
    recalculate_true_offtarget_risk(
        scored_csv='Master_CFD_Scored_OffTargets.csv',
        output_csv='Corrected_Guide_Level_Elevation_Cumulative_Risk.csv'
    )
