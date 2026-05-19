import os
import pickle
import argparse
import pandas as pd
import numpy as np

def get_parser():
    parser = argparse.ArgumentParser(description='Unified true CFD and Cumulative Elevation Scorer')
    parser.add_argument('--input', type=str, default='cas_offinder_output.txt',
                        help='Raw Cas-OFFinder 7-column output file')
    parser.add_argument('--output_all', type=str, default='Master_CFD_Scored_OffTargets.csv',
                        help='Output path for all scored rows')
    parser.add_argument('--output_summary', type=str, default='Corrected_Guide_Level_Elevation_Cumulative_Risk.csv',
                        help='Output path for cumulative whole-guide metrics')
    return parser

def load_cfd_matrices():
    """Loads and verifies empirical matrix picklings."""
    try:
        with open('mismatch_score.pkl', 'rb') as f:
            mm_scores = pickle.load(f, encoding='latin1')
        with open('pam_scores.pkl', 'rb') as f:
            pam_scores = pickle.load(f, encoding='latin1')
        return mm_scores, pam_scores
    except FileNotFoundError as e:
        raise FileNotFoundError(
            "CRITICAL ERROR: 'mismatch_score.pkl' or 'pam_scores.pkl' missing from directory."
        ) from e

def calculate_single_cfd(wt_spacer, off_spacer, off_pam, mm_scores, pam_scores):
    """Calculates the exact CFD score by mapping to target-strand DNA complements."""
    score = 1.0
    
    # Standardize to RNA representation (Uracil) for matrix indexing consistency
    wt_rna = wt_spacer.upper().replace('T', 'U')
    off_rna = off_spacer.upper().replace('T', 'U')
    off_pam = off_pam.upper()
    
    # Map protospacer strand bases to their target complementary strand equivalents
    complement_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
    
    # 1. Process Spacer Mismatches
    for i in range(len(wt_rna)):
        if wt_rna[i] != off_rna[i]:
            # Translate protospacer base to target strand base
            dna_target_base = complement_map.get(off_rna[i], off_rna[i])
            
            # Format expected by official Broad matrix: 'rA:dG,5'
            key = f"r{wt_rna[i]}:d{dna_target_base},{i+1}"
            
            if key in mm_scores:
                score *= float(mm_scores[key])
                
    # 2. Process PAM Suffix
    pam_core = off_pam[-2:]
    if pam_core in pam_scores:
        score *= float(pam_scores[pam_core])
    else:
        score *= 0.01 # Missing/non-canonical PAM penalty
        
    return score

def execute_pipeline():
    args = get_parser().parse_args()
    
    print("Initializing Python 3 matrix variables...")
    mm_scores, pam_scores = load_cfd_matrices()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found.")
        return
        
    print(f"Parsing Cas-OFFinder results from {args.input}...")
    raw_data = []
    with open(args.input, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                raw_data.append(parts[:7])
                
    df = pd.DataFrame(raw_data, columns=[
        'Query_Seq', 'Chromosome', 'Position', 'Full_OffTarget_Seq', 'Direction', 'Mismatches', 'sgRNA_ID'
    ])
    df['Mismatches'] = df['Mismatches'].astype(int)
    
    print(f"Loaded {len(df)} lines. Computing individual CFD scores...")
    
    cfd_scores = []
    for _, row in df.iterrows():
        full_target = row['Full_OffTarget_Seq'].upper()
        off_spacer = full_target[:20]
        off_pam = full_target[20:]
        wt_spacer = row['Query_Seq'].upper()[:20]
        
        score = calculate_single_cfd(wt_spacer, off_spacer, off_pam, mm_scores, pam_scores)
        cfd_scores.append(score)
        
    df['CFD_Score'] = cfd_scores
    df.to_csv(args.output_all, index=False)
    print(f"Saved individual row details to '{args.output_all}'")
    
    # -------------------------------------------------------------
    # CUMULATIVE SUMMARY GENERATION (ON-TARGETS EXCLUDED)
    # -------------------------------------------------------------
    print("Isolating true off-target events (Mismatches > 0)...")
    df_off_only = df[df['Mismatches'] > 0]
    
    grouped = df_off_only.groupby('sgRNA_ID')
    summary_records = []
    
    print("Computing cumulative whole-guide Elevation risks...")
    for sgRNA_id, group in grouped:
        scores = group['CFD_Score'].values
        
        # Elevation Cumulative Formula: 1 - Product(1 - CFD_i)
        prob_no_cut = 1.0 - scores
        prob_miss_all = np.prod(prob_no_cut)
        cumulative_risk = 1.0 - prob_miss_all
        
        summary_records.append({
            'sgRNA_ID': sgRNA_id,
            'True_OffTarget_Matches': len(group),
            'High_Risk_Sites_Count': np.sum(scores >= 0.2),
            'Max_True_OffTarget_Score': np.max(scores) if len(scores) > 0 else 0.0,
            'Cumulative_OffTarget_Risk': cumulative_risk
        })
        
    summary_df = pd.DataFrame(summary_records)
    summary_df = summary_df.sort_values(by='Cumulative_OffTarget_Risk', ascending=False)
    summary_df.to_csv(args.output_summary, index=False)
    
    print("\n" + "="*80)
    print("      CORRECTED ELEVATION CUMULATIVE RISK RANKING (ON-TARGETS EXCLUDED)")
    print("="*80)
    print(summary_df.head(15).to_string(index=False))
    print("="*80)

if __name__ == '__main__':
    execute_pipeline()
