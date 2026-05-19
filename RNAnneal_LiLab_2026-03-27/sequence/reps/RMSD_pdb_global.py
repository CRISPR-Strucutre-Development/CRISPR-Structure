import os
import re
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, Superimposer

def build_library_rmsd_profiles(pdb_root_dir):
    print(f"Scanning flat PDB repository: {pdb_root_dir}")
    if not os.path.exists(pdb_root_dir):
        print(f"CRITICAL ERROR: Path '{pdb_root_dir}' does not exist.")
        return

    parser = PDBParser(QUIET=True)
    superimposer = Superimposer()
    
    # Match flat folder schema exactly: sgRNA_99_rank05.pdb (case-insensitive)
    file_pattern = re.compile(r'^(sgrna)_(\d+)_rank(\d+)\.pdb$', re.IGNORECASE)
    all_files = os.listdir(pdb_root_dir)
    
    valid_matches = []
    for f in all_files:
        match = file_pattern.match(f)
        if match:
            sgrna_id = f"sgRNA_{match.group(2)}"
            valid_matches.append((f, sgrna_id, int(match.group(3)), int(match.group(2))))
                
    df_files = pd.DataFrame(valid_matches, columns=['Filename', 'sgRNA_ID', 'RankNum', 'NumericID'])
    df_files = df_files[df_files['RankNum'] <= 5]
    
    if df_files.empty:
        print("CRITICAL ERROR: Found 0 files matching 'sgRNA_{ID}_rank{num}.pdb'.")
        return
        
    unique_ids = sorted(df_files['sgRNA_ID'].unique(), key=lambda x: int(x.split('_')[1]))
    print(f"Processing {len(df_files)} rank files across {len(unique_ids)} unique sgRNAs...")
    
    # Structural boundaries matching your sequential PDB layout (Scaffold = 21-100)
    scaffold_start, scaffold_end = 21, 100
    
    # Set sgRNA_1 Rank 1 as the absolute geometric baseline for the library alignment
    ref_rows = df_files[(df_files['NumericID'] == 1) & (df_files['RankNum'] == 1)]
    ref_file_name = ref_rows['Filename'].values[0] if not ref_rows.empty else df_files.sort_values(['NumericID', 'RankNum']).iloc[0]['Filename']
        
    ref_path = os.path.join(pdb_root_dir, ref_file_name)
    ref_structure = parser.get_structure('REF', ref_path)
    ref_model = ref_structure[0]
    ref_atoms = [a for r in ref_model.get_residues() if scaffold_start <= r.get_id()[1] <= scaffold_end for a in r.get_atoms() if a.get_name() == "C4'"]
    
    records = []
    
    for sgrna in unique_ids:
        guide_files = df_files[df_files['sgRNA_ID'] == sgrna]
        rmsd_values = []
        
        for _, row in guide_files.iterrows():
            pdb_file = os.path.join(pdb_root_dir, row['Filename'])
            try:
                structure = parser.get_structure(sgrna, pdb_file)
                model = structure[0]
                target_atoms = [a for r in model.get_residues() if scaffold_start <= r.get_id()[1] <= scaffold_end for a in r.get_atoms() if a.get_name() == "C4'"]
                
                min_len = min(len(ref_atoms), len(target_atoms))
                if min_len == 0:
                    continue
                    
                superimposer.set_atoms(ref_atoms[:min_len], target_atoms[:min_len])
                rmsd_values.append(superimposer.rms)
            except Exception:
                continue
                
        if rmsd_values:
            records.append({
                'sgRNA': sgrna,
                'Mean_Scaffold_RMSD': np.mean(rmsd_values),
                'Max_Scaffold_RMSD': np.max(rmsd_values),
                'Ensemble_RMSD_Fluctuation': np.std(rmsd_values)
            })
            
    master_df = pd.DataFrame(records)
    master_df.to_csv('Scaffold_Ensemble_RMSD_Profiles.csv', index=False)
    print("\n" + "="*65)
    print(" STEP 1 COMPLETE: Master structural profile file generated:")
    print(" --> 'Scaffold_Ensemble_RMSD_Profiles.csv'")
    print("="*65)

if __name__ == "__main__":
    build_library_rmsd_profiles(
        pdb_root_dir='/local/projects-t3/lilab/vmenon/SpCas9-Strucutre/pdb/RNAnneal_LiLab_2026-03-27/sequence/reps/pdb_files/'
    )
