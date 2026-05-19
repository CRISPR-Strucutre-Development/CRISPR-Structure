import os

def parse_fasta_and_format(fasta_file, output_txt, genome_path, mismatch_limit=4):
    print(f"Reading full guide sequences from: {fasta_file}")
    
    if not os.path.exists(fasta_file):
        print(f"CRITICAL ERROR: '{fasta_file}' not found.")
        return

    # 1. Parse the FASTA file manually (no external dependencies needed)
    sequences = {}
    with open(fasta_file, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Save the previous sequence before starting a new one
                if header:
                    sequences[header] = "".join(seq_lines)
                # Clean the header to use as the ID (take the first word to avoid spaces)
                header = line[1:].split()[0]
                seq_lines = []
            elif line:
                seq_lines.append(line)
        # Catch the final sequence in the file
        if header:
            sequences[header] = "".join(seq_lines)

    # 2. Extract the spacer and write the Cas-OFFinder input file
    print(f"Extracting 20nt spacers and generating {output_txt}...")
    
    valid_count = 0
    with open(output_txt, 'w') as out_f:
        # Line 1: Genome FASTA directory path
        out_f.write(f"{genome_path}\n")
        
        # Line 2: SpCas9 Template (20 Ns + NRG for canonical/non-canonical PAMs)
        out_f.write("NNNNNNNNNNNNNNNNNNNNNRG\n")
        
        # Line 3+: Process each guide
        for guide_id, full_sequence in sequences.items():
            full_sequence = full_sequence.upper()
            
            # Slice the first 20 nucleotides (indices 0 through 19)
            spacer = full_sequence[0:20]
            
            # Strict quality control
            if len(spacer) != 20:
                print(f"WARNING: Skipping '{guide_id}'. Extracted sequence is {len(spacer)}nt, expected 20nt.")
                continue
                
            # Write exactly: [20nt Spacer] [Mismatches] [ID]
            out_f.write(f"{spacer} {mismatch_limit} {guide_id}\n")
            valid_count += 1

    print("\n" + "="*70)
    print(f" SUCCESS: Processed {valid_count} sequences from FASTA.")
    print(f" Scrapped scaffold data. Kept only the 20nt 5'-end spacer.")
    print(f" Ready for GPU run: ./cas-offinder {output_txt} G0 cas_output.txt")
    print("="*70)

if __name__ == "__main__":
    # --- UPDATE THESE THREE VARIABLES BEFORE RUNNING ---
    
    FASTA_INPUT = "input.fa"  # Your raw FASTA file with full sequences
    TXT_OUTPUT = "150_cas_off_finder_input.txt"         # The file to feed into Cas-OFFinder
    GENOME_DIR = "/local/projects-t3/lilab/vmenon/CRISPRi/Genome/chromosomes/" # Absolute path to genome
    
    parse_fasta_and_format(
        fasta_file=FASTA_INPUT,
        output_txt=TXT_OUTPUT,
        genome_path=GENOME_DIR,
        mismatch_limit=4
    )
