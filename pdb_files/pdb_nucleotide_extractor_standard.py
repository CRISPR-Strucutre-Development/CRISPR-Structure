import os
from Bio import PDB

def extract_pdb_nucleotides(input_directory, output_txt_file):
    """Extracts nucleotide sequences from PDB files in a given directory.

    This function iterates through all `.pdb` files in the specified input directory,
    parses each structure, and extracts sequences of standard DNA (DA, DC, DG, DT)
    and RNA (A, C, G, U) nucleotides. The extracted sequences are written to an
    output text file, with each PDB file's sequences grouped together.

    For each PDB file:
    - It attempts to parse the structure using Biopython's PDBParser.
    - Only the first model (index 0) of the structure is processed.
    - It iterates through chains and residues, mapping recognized nucleotide
      residue names to their single-letter codes.
    - If a chain contains nucleotides, its sequence is written to the output file.
    - If no standard nucleotide chains are detected in a PDB file, a corresponding
      message is written.
    - If a PDB file fails to parse, an error message is logged in the output file.

    Args:
        input_directory (str): The path to the directory containing the PDB files.
        output_txt_file (str): The path to the output text file where the
                                extracted nucleotide sequences will be written.

    Returns:
        None: This function does not return any value; it writes its output
              directly to the specified `output_txt_file`.
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
