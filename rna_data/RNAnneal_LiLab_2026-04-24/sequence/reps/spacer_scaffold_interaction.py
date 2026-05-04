import os
import csv
import re
import barnaba as bb

def get_res_num(res_str):
    """Extracts the integer residue number from the residue identifier string."""
    match = re.search(r'\d+', res_str)
    return int(match.group()) if match else None

def check_region_interaction(pairings, res, region_a, region_b):
    """Checks if any residue in region_a interacts with any residue in region_b."""
    if pairings and len(pairings) > 0 and len(pairings[0][0]) > 0:
        for p in range(len(pairings[0][0])):
            res1_idx = pairings[0][0][p][0]
            res2_idx = pairings[0][0][p][1]
            
            r1 = get_res_num(res[res1_idx])
            r2 = get_res_num(res[res2_idx])
            
            if (r1 in region_a and r2 in region_b) or (r1 in region_b and r2 in region_a):
                return 1
    return 0

def process_sgrna_structures(folder_path=".", output_csv="sgrna_interactions.csv"):
    """
    Scans the directory for PDB files, calculates base-pair interactions
    for ranks 1 through 5 only, and writes the output to a CSV file.
    """
    # Define the regions
    spacer_region = range(1, 21)
    regions = {
        "TL_37_42": range(37, 43),
        "SL1_64_70": range(64, 71),
        "SL2_82_87": range(82, 88),
        "SL3_98_100": range(98, 101)
    }
    
    # Discover unique sgRNA files by grouping prefixes
    files = os.listdir(folder_path)
    sgrna_prefixes = set()
    
    for f in files:
        if f.endswith('.pdb') and '_rank' in f:
            prefix = f.split('_rank')[0]
            sgrna_prefixes.add(prefix)
            
    output_file_path = os.path.join(folder_path, output_csv)
    print(f"Target Directory: {os.path.abspath(folder_path)}")
    print(f"Writing CSV file to: {output_file_path}")
    
    with open(output_file_path, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        # Headers now accurately reflect the defined dictionary keys
        header = ['sgRNA', 'TL_37_42', 'SL1_64_70', 'SL2_82_87', 'SL3_98_100', 'Rank']
        writer.writerow(header)
        
        for prefix in sorted(list(sgrna_prefixes)):
            # Process up to rank 05 for all sgRNAs
            for i in range(1, 6):
                rank_str = f"rank{i:02d}"
                pdb_name = f"{prefix}_{rank_str}.pdb"
                pdb_path = os.path.join(folder_path, pdb_name)
                
                # If the rank file exists, read and process it
                if os.path.exists(pdb_path):
                    try:
                        stackings, pairings, res = bb.annotate(pdb_path)
                        
                        tl_val = check_region_interaction(pairings, res, spacer_region, regions["TL_37_42"])
                        sl1_val = check_region_interaction(pairings, res, spacer_region, regions["SL1_64_70"])
                        sl2_val = check_region_interaction(pairings, res, spacer_region, regions["SL2_82_87"])
                        sl3_val = check_region_interaction(pairings, res, spacer_region, regions["SL3_98_100"])
                        
                        writer.writerow([prefix, tl_val, sl1_val, sl2_val, sl3_val, rank_str])
                    except Exception as e:
                        print(f"Error reading structure {pdb_name}: {e}")
                        writer.writerow([prefix, 0, 0, 0, 0, rank_str])
                else:
                    # Missing rank file; populate with 0
                    writer.writerow([prefix, 0, 0, 0, 0, rank_str])
                    
    print("Processing complete.")

if __name__ == "__main__":
    # Specify the directory containing your PDB files
    target_directory = "/local/projects-t3/lilab/vmenon/SpCas9-Strucutre/rna_data/RNAnneal_LiLab_2026-04-24/sequence/reps/chimera/"
    process_sgrna_structures(folder_path=target_directory)
