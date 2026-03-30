"""Extracts nucleotide sequences from PDB files and saves them to a text file.

This script iterates through a specified directory of PDB files, parses the
first model of each structure, and extracts the standard DNA and RNA
nucleotide sequences for each chain. The results are written to an output file.
"""
import os
from Bio import PDB

def extract_pdb_nucleotides(input_directory, output_txt_file):
    """Extracts nucleotide sequences from all PDB files in a directory.

    Iterates over all `.pdb` files in the specified input directory, reads the 
    first structural model using Bio.PDB, and translates standard nucleotide 
    residue names into their single-letter sequence codes. The extracted 
    sequences are then written to the specified output text file.

    Args:
        input_directory (str): Path to the directory containing `.pdb` files.
        output_txt_file (str): Path to the output text file where the extracted
            sequences will be written.
    """
    # Suppress PDB construction warnings for cleaner execution
    parser = PDB.PDBParser(QUIET=True)
    
    # Standard PDB nomenclature for nucleotides
    nucleotide_map = {
        'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T',  # DNA
        'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U'       # RNA
    }

    with open(output_txt_file, 'w') as out_file:
        for filename in os.listdir(input_directory):
            if not filename.endswith('.pdb'):
                continue
                
            filepath = os.path.join(input_directory, filename)
            
            try:
                # Load the structure. 'filename' acts as the structure ID.
                structure = parser.get_structure(filename, filepath)
            except Exception as e:
                out_file.write(f">{filename} - ERROR PARSING: {e}\n\n")
                continue
            
            out_file.write(f">{filename}\n")
            has_nucleotides = False
            
            # PDB files can have multiple models (e.g., NMR structures). 
            # We only extract from the first model (index 0) to avoid duplicates.
            model = structure[0]
            
            for chain in model:
                sequence = []
                for residue in chain:
                    # Strip whitespace from residue names (e.g., ' DA ' -> 'DA')
                    res_name = residue.get_resname().strip()
                    
                    if res_name in nucleotide_map:
                        sequence.append(nucleotide_map[res_name])
                
                # If nucleotides were found in this chain, write them out
                if sequence:
                    joined_sequence = "".join(sequence)
                    out_file.write(f"Chain_{chain.id}: {joined_sequence}\n")
                    has_nucleotides = True
            
            if not has_nucleotides:
                out_file.write("No standard nucleotide chains detected.\n")
            
            out_file.write("\n")

if __name__ == "__main__":
    # Define your paths here
    PDB_DIRECTORY = "/local/projects-t3/lilab/vmenon/SpCas9-Strucutre/pdb/pdb_files/" 
    OUTPUT_FILE = "extracted_nucleotides.txt"
    
    # Ensure the directory exists before running
    if os.path.exists(PDB_DIRECTORY):
        extract_pdb_nucleotides(PDB_DIRECTORY, OUTPUT_FILE)
        print(f"Extraction complete. Results saved to {OUTPUT_FILE}")
    else:
        print(f"Error: Directory '{PDB_DIRECTORY}' not found.")
