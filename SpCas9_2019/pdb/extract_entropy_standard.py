import os
from Bio import PDB

def process_tri_entropy_data(pdb_dir, out_super_rel_txt, out_rel_txt, out_str_txt):
    """
    Process PDB files to calculate and save average entropy data for structural regions.

    This function reads all PDB files in a given directory, extracts the nucleotide
    sequence and the corresponding B-factor (entropy) for each residue, and calculates
    the average entropy for four regions (TL, SL1, SL2, SL3) across three different
    stringency thresholds (Super Relaxed, Relaxed, Strict). The results are saved into
    three corresponding TSV files.

    Args:
        pdb_dir (str): Path to the directory containing input PDB files.
        out_super_rel_txt (str): File path for the output 'super relaxed' TSV file.
        out_rel_txt (str): File path for the output 'relaxed' TSV file.
        out_str_txt (str): File path for the output 'strict' TSV file.

    Returns:
        None
    """
    parser = PDB.PDBParser(QUIET=True)
    nuc_map = {'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T', 'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U'}
    
    # Open all three output files simultaneously
    with open(out_super_rel_txt, 'w') as out_srel, open(out_rel_txt, 'w') as out_rel, open(out_str_txt, 'w') as out_str:
        
        header = "ID\tsequence(sgRNA)\tTL\tSL1\tSL2\tSL3\n"
        out_srel.write(header)
        out_rel.write(header)
        out_str.write(header)
        
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
            entropy_dict = {}
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        res_name = residue.get_resname().strip()
                        res_id = residue.get_id()[1] 
                        
                        if res_name in nuc_map:
                            sequence.append(nuc_map[res_name])
                            atoms = list(residue.get_atoms())
                            if atoms:
                                entropy_dict[res_id] = atoms[0].get_bfactor()
                break 

            seq_string = "".join(sequence)
            
            def get_avg(start, end, region_name, dataset_type):
                """
                Calculate the average entropy for a specific range of residues.

                Args:
                    start (int): The starting residue ID (inclusive).
                    end (int): The ending residue ID (inclusive).
                    region_name (str): The name of the region being processed (e.g., 'TL(32-47)').
                    dataset_type (str): The stringency level or dataset type for logging purposes
                        (e.g., 'Super Relaxed').

                Returns:
                    float or str: The average entropy rounded to 4 decimal places, or "NA"
                        if no residues are found in the specified range.
                """
                missing_res = []
                vals = []
                
                for i in range(start, end + 1):
                    if i in entropy_dict:
                        vals.append(entropy_dict[i])
                    else:
                        missing_res.append(i)
                
                if missing_res:
                    print(f"ALERT [{clean_id} | {dataset_type}]: Missing {region_name} coordinates -> {missing_res}")
                    
                if not vals:
                    return "NA"
                return round(sum(vals) / len(vals), 4)

            # 1. Calculate Super Relaxed Coordinates
            tl_srel = get_avg(32, 47, "TL(32-47)", "Super Relaxed")
            sl1_srel = get_avg(63, 72, "SL1(63-72)", "Super Relaxed")
            sl2_srel = get_avg(79, 90, "SL2(79-90)", "Super Relaxed")
            sl3_srel = get_avg(92, 106, "SL3(92-106)", "Super Relaxed")
            
            out_srel.write(f"{clean_id}\t{seq_string}\t{tl_srel}\t{sl1_srel}\t{sl2_srel}\t{sl3_srel}\n")

            # 2. Calculate Relaxed Coordinates
            tl_rel = get_avg(36, 43, "TL(36-43)", "Relaxed")
            sl1_rel = get_avg(63, 71, "SL1(63-71)", "Relaxed")
            sl2_rel = get_avg(80, 90, "SL2(80-90)", "Relaxed")
            sl3_rel = get_avg(92, 101, "SL3(92-101)", "Relaxed")
            
            out_rel.write(f"{clean_id}\t{seq_string}\t{tl_rel}\t{sl1_rel}\t{sl2_rel}\t{sl3_rel}\n")

            # 3. Calculate Strict Coordinates
            tl_str = get_avg(38, 41, "TL(38-41)", "Strict")
            sl1_str = get_avg(65, 69, "SL1(65-69)", "Strict")
            sl2_str = get_avg(83, 86, "SL2(83-86)", "Strict")
            sl3_str = get_avg(98, 100, "SL3(98-100)", "Strict")
            
            out_str.write(f"{clean_id}\t{seq_string}\t{tl_str}\t{sl1_str}\t{sl2_str}\t{sl3_str}\n")

if __name__ == "__main__":
    PDB_DIRECTORY = "." 
    SUPER_RELAXED_FILE = "entropy_averages_super_relaxed.tsv"
    RELAXED_FILE = "entropy_averages_relaxed.tsv"
    STRICT_FILE = "entropy_averages_strict.tsv"
    
    if os.path.exists(PDB_DIRECTORY):
        print("Starting triple entropy extraction (Super Relaxed, Relaxed, & Strict parameters)...\n")
        process_tri_entropy_data(PDB_DIRECTORY, SUPER_RELAXED_FILE, RELAXED_FILE, STRICT_FILE)
        print(f"\nExtraction complete.")
        print(f"- Super Relaxed data saved to: {SUPER_RELAXED_FILE}")
        print(f"- Relaxed data saved to:       {RELAXED_FILE}")
        print(f"- Strict data saved to:        {STRICT_FILE}")
    else:
        print(f"Directory not found: {PDB_DIRECTORY}")
