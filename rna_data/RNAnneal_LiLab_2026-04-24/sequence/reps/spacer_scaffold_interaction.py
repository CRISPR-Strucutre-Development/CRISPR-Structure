import os
import csv
import re
import barnaba as bb
from Bio.PDB import PDBParser

def get_res_num(res_str):
    """Extracts the integer residue number from the residue identifier string."""
    match = re.search(r'\d+', res_str)
    return int(match.group()) if match else None

def check_region_interaction(interactions_data, res, region_a, region_b):
    """Checks if any residue in region_a interacts with any residue in region_b."""
    if interactions_data and len(interactions_data) > 0 and len(interactions_data[0][0]) > 0:
        for p in range(len(interactions_data[0][0])):
            res1_idx = interactions_data[0][0][p][0]
            res2_idx = interactions_data[0][0][p][1]
            
            r1 = get_res_num(res[res1_idx])
            r2 = get_res_num(res[res2_idx])
            
            if (r1 in region_a and r2 in region_b) or (r1 in region_b and r2 in region_a):
                return 1
    return 0

def get_min_heavy_atom_distance(pdb_path, region_a, region_b):
    """Calculates the minimum heavy-atom distance between two regions in Ångströms."""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('rna', pdb_path)
        atoms_a = []
        atoms_b = []
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    res_id = residue.id[1]
                    if res_id in region_a:
                        for atom in residue:
                            if atom.element != 'H':
                                atoms_a.append(atom)
                    elif res_id in region_b:
                        for atom in residue:
                            if atom.element != 'H':
                                atoms_b.append(atom)
                                
        min_dist = float('inf')
        for a in atoms_a:
            for b in atoms_b:
                dist = a - b
                if dist < min_dist:
                    min_dist = dist
                    
        return round(min_dist, 3) if min_dist != float('inf') else None
    except Exception:
        return None

def process_sgrna_structures(folder_path=".", output_prefix="sgrna_"):
    """
    Scans the directory for PDB files and writes three output CSV files:
    Pairings, Stackings, and Heavy-Atom Minimum Distances.
    """
    spacer_region = range(1, 21)
    regions = {
        "TL_37_42": range(37, 43),
        "SL1_64_70": range(64, 71),
        "SL2_82_87": range(82, 88),
        "SL3_98_100": range(98, 101)
    }
    
    files = os.listdir(folder_path)
    sgrna_prefixes = set()
    
    for f in files:
        if f.endswith('.pdb') and '_rank' in f:
            prefix = f.split('_rank')[0]
            sgrna_prefixes.add(prefix)
            
    # Initialize three file handles
    pair_file_path = os.path.join(folder_path, f"{output_prefix}interactions_pair.csv")
    stack_file_path = os.path.join(folder_path, f"{output_prefix}interactions_stack.csv")
    dist_file_path = os.path.join(folder_path, f"{output_prefix}interactions_distance.csv")
    
    print(f"Target Directory: {os.path.abspath(folder_path)}")
    print(f"Writing metrics to {pair_file_path}, {stack_file_path}, and {dist_file_path}")
    
    with open(pair_file_path, mode='w', newline='') as f_pair, \
         open(stack_file_path, mode='w', newline='') as f_stack, \
         open(dist_file_path, mode='w', newline='') as f_dist:
        
        pair_writer = csv.writer(f_pair)
        stack_writer = csv.writer(f_stack)
        dist_writer = csv.writer(f_dist)
        
        header = ['sgRNA', 'Rank', 'TL_37_42', 'SL1_64_70', 'SL2_82_87', 'SL3_98_100']
        
        pair_writer.writerow(header)
        stack_writer.writerow(header)
        dist_writer.writerow(header)
        
        for prefix in sorted(list(sgrna_prefixes)):
            for i in range(1, 6):
                rank_str = f"rank{i:02d}"
                pdb_name = f"{prefix}_{rank_str}.pdb"
                pdb_path = os.path.join(folder_path, pdb_name)
                
                if os.path.exists(pdb_path):
                    try:
                        stackings, pairings, res = bb.annotate(pdb_path)
                        
                        # Extract 1/0 for pairings
                        p_tl = check_region_interaction(pairings, res, spacer_region, regions["TL_37_42"])
                        p_sl1 = check_region_interaction(pairings, res, spacer_region, regions["SL1_64_70"])
                        p_sl2 = check_region_interaction(pairings, res, spacer_region, regions["SL2_82_87"])
                        p_sl3 = check_region_interaction(pairings, res, spacer_region, regions["SL3_98_100"])
                        
                        # Extract 1/0 for stackings
                        s_tl = check_region_interaction(stackings, res, spacer_region, regions["TL_37_42"])
                        s_sl1 = check_region_interaction(stackings, res, spacer_region, regions["SL1_64_70"])
                        s_sl2 = check_region_interaction(stackings, res, spacer_region, regions["SL2_82_87"])
                        s_sl3 = check_region_interaction(stackings, res, spacer_region, regions["SL3_98_100"])
                        
                        # Calculate continuous heavy atom minimum distance
                        d_tl = get_min_heavy_atom_distance(pdb_path, spacer_region, regions["TL_37_42"])
                        d_sl1 = get_min_heavy_atom_distance(pdb_path, spacer_region, regions["SL1_64_70"])
                        d_sl2 = get_min_heavy_atom_distance(pdb_path, spacer_region, regions["SL2_82_87"])
                        d_sl3 = get_min_heavy_atom_distance(pdb_path, spacer_region, regions["SL3_98_100"])
                        
                        pair_writer.writerow([prefix, rank_str, p_tl, p_sl1, p_sl2, p_sl3])
                        stack_writer.writerow([prefix, rank_str, s_tl, s_sl1, s_sl2, s_sl3])
                        dist_writer.writerow([prefix, rank_str, d_tl, d_sl1, d_sl2, d_sl3])
                        
                    except Exception as e:
                        print(f"Error processing structure {pdb_name}: {e}")
                        pair_writer.writerow([prefix, rank_str, 0, 0, 0, 0])
                        stack_writer.writerow([prefix, rank_str, 0, 0, 0, 0])
                        dist_writer.writerow([prefix, rank_str, 'NA', 'NA', 'NA', 'NA'])
                else:
                    pair_writer.writerow([prefix, rank_str, 0, 0, 0, 0])
                    stack_writer.writerow([prefix, rank_str, 0, 0, 0, 0])
                    dist_writer.writerow([prefix, rank_str, 'NA', 'NA', 'NA', 'NA'])
                    
    print("Processing complete.")

if __name__ == "__main__":
    target_directory = "/local/projects-t3/lilab/vmenon/SpCas9-Strucutre/pdb/rna_data/RNAnneal_LiLab_2026-04-24/sequence/reps/chimera/"
    process_sgrna_structures(folder_path=target_directory)
