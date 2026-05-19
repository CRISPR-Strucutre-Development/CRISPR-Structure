import os
import pandas as pd
import glob

def extract_pdb_entropy(pdb_dir, output_csv='pdb_entropy_extracted.csv'):
    # Define residues of interest
    regions = {
        'spacer': list(range(1, 21)),
        'TL': list(range(32, 38)),
        'SL1': list(range(54, 61)),
        'SL2': list(range(71, 77)),
        'SL3': list(range(88, 91))
    }
    
    results = []
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    
    print(f"Processing {len(pdb_files)} PDB files...")

    for pdb_path in pdb_files:
        # Expected filename format: sgRNA_ID_rank01.pdb
        filename = os.path.basename(pdb_path)
        parts = filename.replace(".pdb", "").split("_")
        sgrna_id = "_".join(parts[:-1])
        rank = parts[-1]

        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    res_seq = int(line[22:26].strip())
                    entropy_val = float(line[60:66].strip()) # 10th Column (B-factor)
                    
                    for region_name, res_list in regions.items():
                        if res_seq in res_list:
                            results.append({
                                'sgRNA': sgrna_id,
                                'Rank': rank,
                                'Region': region_name,
                                'Residue': res_seq,
                                'Entropy': entropy_val
                            })
                            
    df = pd.DataFrame(results)
    # Average per region per rank per sgRNA
    df_summary = df.groupby(['sgRNA', 'Rank', 'Region'])['Entropy'].mean().reset_index()
    df_summary.to_csv(output_csv, index=False)
    print(f"Extraction complete. Data saved to {output_csv}")

if __name__ == "__main__":
    extract_pdb_entropy('/local/projects-t3/lilab/vmenon/SpCas9-Strucutre/pdb/RNAnneal_LiLab_2026-03-27/sequence/reps/')
