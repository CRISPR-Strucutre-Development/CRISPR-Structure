import os
import re
import pandas as pd
import numpy as np
import barnaba as bb
from Bio.PDB import PDBParser

def get_base_atoms_coords(residue):
    """
    Extracts 3D coordinates of heavy atoms strictly inside the nitrogenous base ring.
    Excludes the sugar-phosphate backbone to isolate the base-pairing/stacking geometry.
    """
    backbone_atoms = ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]
    coords = [a.get_coord() for a in residue.get_atoms() if a.get_name() not in backbone_atoms and a.element != 'H']
    if not coords:
        coords = [a.get_coord() for a in residue.get_atoms() if a.element != 'H']
    return coords

def extract_residue_number(res_string):
    match = re.search(r'\d+', res_string)
    if match:
        return int(match.group())
    raise ValueError(f"Could not parse numeric ID from Barnaba: '{res_string}'")

def build_structural_master_files(pdb_root_dir):
    print(f"Scanning flat PDB repository: {pdb_root_dir}")
    if not os.path.exists(pdb_root_dir):
        print(f"CRITICAL ERROR: The path '{pdb_root_dir}' does not exist.")
        return
        
    parser = PDBParser(QUIET=True)
    
    # 1. Match files exactly like 'sgRNA_99_rank05.pdb' case-insensitively
    file_pattern = re.compile(r'^(sgrna)_(\d+)_rank(\d+)\.pdb$', re.IGNORECASE)
    
    all_files = os.listdir(pdb_root_dir)
    valid_matches = []
    
    for f in all_files:
        match = file_pattern.match(f)
        if match:
            # Group by sgRNA ID to aggregate ranks later
            sgrna_id = f"sgRNA_{match.group(2)}"
            valid_matches.append((f, sgrna_id, int(match.group(3))))
            
    if not valid_matches:
        print("CRITICAL ERROR: Found 0 files matching 'sgRNA_{ID}_rank{num}.pdb' pattern.")
        print("Double-check filenames using: ls " + pdb_root_dir)
        return
        
    df_files = pd.DataFrame(valid_matches, columns=['Filename', 'sgRNA_ID', 'RankNum'])
    # Restrict to Ranks 1-5
    df_files = df_files[df_files['RankNum'] <= 5]
    unique_ids = df_files['sgRNA_ID'].unique()
    print(f"Identified {len(df_files)} total rank files covering {len(unique_ids)} unique sgRNAs.")

    scaffold_start, scaffold_end = 21, 100
    
    pairing_records = []
    stacking_records = []
    distance_records = []
    
    # 2. Iterate through each unique guide sequence pool
    for sgrna in sorted(unique_ids):
        guide_files = df_files[df_files['sgRNA_ID'] == sgrna]
        
        pair_tracker = {}
        stack_tracker = {}
        dist_tracker = {}
        
        valid_files_processed = 0
        
        for _, row in guide_files.iterrows():
            pdb_file = os.path.join(pdb_root_dir, row['Filename'])
            
            try:
                # Top-tier geometric parsing via Barnaba
                stackings, pairings, res_list = bb.annotate(pdb_file)
                
                # Direct distance tracking via Bio.PDB
                structure = parser.get_structure(sgrna, pdb_file)
                model = structure[0]
                pdb_residues = {r.get_id()[1]: r for r in model.get_residues() if scaffold_start <= r.get_id()[1] <= scaffold_end}
                
                all_interactions = []
                if len(pairings) > 0 and pairings[0] is not None:
                    all_interactions.extend([(p[0], p[1], 'pair') for p in pairings[0][0]])
                if len(stackings) > 0 and stackings[0] is not None:
                    all_interactions.extend([(s[0], s[1], 'stack') for s in stackings[0][0]])
                
                if all_interactions:
                    valid_files_processed += 1
                    
                for idx1, idx2, int_type in all_interactions:
                    res_num1 = extract_residue_number(res_list[idx1])
                    res_num2 = extract_residue_number(res_list[idx2])
                    
                    # Enforce the scaffold footprint bounds (21-100)
                    if scaffold_start <= res_num1 <= scaffold_end and scaffold_start <= res_num2 <= scaffold_end:
                        id1, id2 = sorted([res_num1, res_num2])
                        pair_key = f"Res_{id1}_{id2}"
                        
                        if id1 in pdb_residues and id2 in pdb_residues:
                            coords1 = get_base_atoms_coords(pdb_residues[id1])
                            coords2 = get_base_atoms_coords(pdb_residues[id2])
                            min_euclidean_d = min([np.linalg.norm(c1 - c2) for c1 in coords1 for c2 in coords2])
                        else:
                            min_euclidean_d = 99.0
                            
                        if int_type == 'pair':
                            pair_tracker[pair_key] = pair_tracker.get(pair_key, 0) + 1
                        if int_type == 'stack':
                            stack_tracker[pair_key] = stack_tracker.get(pair_key, 0) + 1
                            
                        if min_euclidean_d <= 3.5:
                            dist_tracker[pair_key] = dist_tracker.get(pair_key, 0) + 1
                            
            except Exception as e:
                continue
                
        if valid_files_processed > 0:
            p_row = {'sgRNA': sgrna}; p_row.update(pair_tracker)
            s_row = {'sgRNA': sgrna}; s_row.update(stack_tracker)
            d_row = {'sgRNA': sgrna}; d_row.update(dist_tracker)
            
            pairing_records.append(p_row)
            stacking_records.append(s_row)
            distance_records.append(d_row)

    if not pairing_records:
        print("CRITICAL ERROR: Files parsed but zero cross-residue interactions matched the Barnaba parameters.")
        return

    # 3. Save Wide Master Matrices
    df_pair_master = pd.DataFrame(pairing_records).fillna(0)
    df_stack_master = pd.DataFrame(stacking_records).fillna(0)
    df_dist_master = pd.DataFrame(distance_records).fillna(0)
    
    df_pair_master.to_csv('scaffold_base_pairing.csv', index=False)
    df_stack_master.to_csv('scaffold_base_stacking.csv', index=False)
    df_dist_master.to_csv('scaffold_heavy_distances.csv', index=False)
    
    print("\n" + "="*75)
    print(f"Success! Master structural matrices generated for {len(pairing_records)} active IDs:")
    print(" 1. scaffold_base_pairing.csv")
    print(" 2. scaffold_base_stacking.csv")
    print(" 3. scaffold_heavy_distances.csv")
    print("="*75)

if __name__ == "__main__":
    build_structural_master_files(
        pdb_root_dir='/local/projects-t3/lilab/vmenon/SpCas9-Strucutre/pdb/RNAnneal_LiLab_2026-03-27/sequence/reps/pdb_files/'
    )