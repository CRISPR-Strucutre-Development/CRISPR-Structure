import os
import statistics
import numpy as np
from Bio import PDB

def process_advanced_metrics(pdb_dir, output_txt):
    parser = PDB.PDBParser(QUIET=True)
    nuc_map = {'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T', 'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U'}
    
    with open(output_txt, 'w') as out:
        # 1. Expand regions and metrics
        regions = ["TL", "SL1", "SL2", "SL3"]
        metrics = ["Mean", "Median", "Max", "Min", "StDev", "Range", "IQR"]
        
        header_cols = ["ID", "sequence(sgRNA)"]
        for r in regions:
            for m in metrics:
                header_cols.append(f"{r}_{m}")
                
        out.write("\t".join(header_cols) + "\n")
        
        for filename in os.listdir(pdb_dir):
            if not filename.endswith('.pdb'):
                continue
                
            filepath = os.path.join(pdb_dir, filename)
            clean_id = filename.replace('.pdb', '')
            
            try:
                structure = parser.get_structure(clean_id, filepath)
            except Exception:
                continue
            
            sequence = []
            entropy_dict = {}
            
            # Extract first model data
            for model in structure:
                for chain in model:
                    for residue in chain:
                        res_name = residue.get_resname().strip()
                        res_id = residue.get_id()[1] 
                        if res_name in nuc_map:
                            sequence.append(nuc_map[res_name])
                            atoms = list(residue.get_atoms())
                            if atoms:
                                # Use B-factor as the entropy metric
                                entropy_dict[res_id] = atoms[0].get_bfactor()
                break 

            seq_string = "".join(sequence)
            
            def get_stats(start, end):
                vals = [entropy_dict[i] for i in range(start, end + 1) if i in entropy_dict]
                if not vals:
                    return ["NA"] * 7
                
                # Minimum of 2 points required for StDev and IQR to make sense
                if len(vals) == 1:
                    v = round(vals[0], 4)
                    return [v, v, v, v, 0.0, 0.0, 0.0]
                
                # Arithmetic and Central Tendency
                mean_v = statistics.mean(vals)
                med_v = statistics.median(vals)
                max_v = max(vals)
                min_v = min(vals)
                
                # Dispersion Metrics
                stdev_v = statistics.stdev(vals)
                range_v = max_v - min_v
                
                # Interquartile Range (IQR) - using numpy for accuracy on small sets
                q75, q25 = np.percentile(vals, [75, 25])
                iqr_v = q75 - q25
                
                return [round(x, 4) for x in [mean_v, med_v, max_v, min_v, stdev_v, range_v, iqr_v]]

            # Use your validated "Strict" coordinates
            stats_tl  = get_stats(38, 41)
            stats_sl1 = get_stats(65, 69)
            stats_sl2 = get_stats(83, 86)
            stats_sl3 = get_stats(98, 100)
            
            row_data = [clean_id, seq_string] + stats_tl + stats_sl1 + stats_sl2 + stats_sl3
            out.write("\t".join(map(str, row_data)) + "\n")

if __name__ == "__main__":
    PDB_DIRECTORY = "." 
    OUTPUT_FILE = "entropy_advanced_strict.tsv"
    
    if os.path.exists(PDB_DIRECTORY):
        process_advanced_metrics(PDB_DIRECTORY, OUTPUT_FILE)
        print(f"Extraction complete. Results saved to {OUTPUT_FILE} (30 columns total).")
    else:
        print("Directory not found.")