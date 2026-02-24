import os
import numpy as np
from Bio import PDB

def process_cumulative_entropy_data(pdb_dir, out_file):
    parser = PDB.PDBParser(QUIET=True)
    nuc_map = {'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T', 'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U'}
    
    with open(out_file, 'w') as out_f:
        # Write the new header with all requested statistical metrics
        header = "ID\tsequence(sgRNA)\tResidue_Count\tEnt_Sum\tEnt_Median\tEnt_Max\tEnt_Min\tEnt_StDev\tEnt_IQR\tEnt_Range\n"
        out_f.write(header)
        
        for filename in os.listdir(pdb_dir):
            if not filename.endswith('.pdb'):
                continue
                
            filepath = os.path.join(pdb_dir, filename)
            clean_id = filename.replace('.pdb', '')
            
            try:
                structure = parser.get_structure(clean_id, filepath)
            except Exception as e:
                print(f"Error parsing {filename}: {e}")
                continue
            
            sequence = []
            entropy_values = []
            
            # Extract sequence and B-factors from the first model
            for model in structure:
                for chain in model:
                    for residue in chain:
                        res_name = residue.get_resname().strip()
                        
                        if res_name in nuc_map:
                            sequence.append(nuc_map[res_name])
                            atoms = list(residue.get_atoms())
                            if atoms:
                                # Grab the B-factor of the first atom (typically the backbone/C4' or P)
                                entropy_values.append(atoms[0].get_bfactor())
                break  # Only process the first model
            
            seq_string = "".join(sequence)
            
            # If no valid entropy values were extracted, fill with NA
            if not entropy_values:
                print(f"ALERT [{clean_id}]: No valid nucleotide coordinates found.")
                out_f.write(f"{clean_id}\t{seq_string}\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
                continue

            # Convert to numpy array for efficient statistical operations
            arr = np.array(entropy_values)
            res_count = len(arr)
            
            # Calculate metrics
            ent_sum = round(np.sum(arr), 4)
            ent_median = round(np.median(arr), 4)
            ent_max = round(np.max(arr), 4)
            ent_min = round(np.min(arr), 4)
            ent_stdev = round(np.std(arr), 4)
            
            q75, q25 = np.percentile(arr, [75, 25])
            ent_iqr = round(q75 - q25, 4)
            ent_range = round(ent_max - ent_min, 4)
            
            # Write to file
            out_f.write(f"{clean_id}\t{seq_string}\t{res_count}\t{ent_sum}\t{ent_median}\t{ent_max}\t{ent_min}\t{ent_stdev}\t{ent_iqr}\t{ent_range}\n")

if __name__ == "__main__":
    PDB_DIRECTORY = "." 
    OUTPUT_FILE = "entropy_cumulative_stats.tsv"
    
    if os.path.exists(PDB_DIRECTORY):
        print("Starting global entropy extraction and statistical analysis...\n")
        process_cumulative_entropy_data(PDB_DIRECTORY, OUTPUT_FILE)
        print(f"\nExtraction complete.")
        print(f"- Data saved to: {OUTPUT_FILE}")
    else:
        print(f"Directory not found: {PDB_DIRECTORY}")
