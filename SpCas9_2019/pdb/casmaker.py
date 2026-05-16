import io
import pandas as pd

def generate_23nt_input_from_rna(csv_file_path, output_txt_path, genome_dir_path):
    print(f"Reading raw transcript data from: {csv_file_path}")
    
    try:
        df = pd.read_csv(csv_file_path)
    except FileNotFoundError:
        print(f"CRITICAL ERROR: File '{csv_file_path}' not found.")
        return

    # Clean column names to eliminate accidental whitespace
    df.columns = df.columns.str.strip()
    
    # Track which columns we are using
    id_col = 'ID' if 'ID' in df.columns else 'sgRNA'
    guide_col = 'guide'
    
    if guide_col not in df.columns:
        print(f"CRITICAL ERROR: Column '{guide_col}' missing from your CSV.")
        print(f"Available columns: {list(df.columns)}")
        return

    print("Extracting spacers, converting to DNA, and appending NGG templates...")
    
    valid_count = 0
    with open(output_txt_path, 'w') as f:
        # Line 1: Absolute path to your hg38 chromosome files
        f.write(f"{genome_dir_path.strip()}\n")
        
        # Line 2: The global template (20 Ns + NRG to capture NAG/NGG variants)
        f.write("NNNNNNNNNNNNNNNNNNNNNRG\n")
        
        # Line 3+: Process individual sequences
        for index, row in df.iterrows():
            guide_id = str(row[id_col]).strip()
            rna_spacer = str(row[guide_col]).strip().upper()
            
            # Extract exactly 20nt just in case there is trailing noise
            rna_spacer = rna_spacer[:20]
            
            if len(rna_spacer) != 20:
                print(f"Skipping {guide_id}: Extracted spacer length is {len(rna_spacer)}nt instead of 20nt.")
                continue
                
            # Convert RNA to DNA string (Uracil -> Thymine)
            dna_spacer = rna_spacer.replace('U', 'T')
            
            # Append canonical 'NGG' to construct the full 23nt target sequence
            # This fills the 23-char slot without breaking Cas-OFFinder syntax
            target_23mer = f"{dna_spacer}NGG"
            
            # Format: [23-char Target] [Mismatch Limit] [sgRNA_ID]
            f.write(f"{target_23mer} 4 {guide_id}\n")
            valid_count += 1
            
    print("\n" + "="*70)
    print(f" SUCCESS: Generated perfect 23nt targets for {valid_count} guides.")
    print(f" File saved to: '{output_txt_path}'")
    print("="*70)

# --- Test execution block using your sample data stream ---
if __name__ == "__main__":
    # Simulated file path variables
    INPUT_CSV = "SpCas9_2019_150_sequence.csv" 
    OUTPUT_TXT = "150_cas_off_finder_input.txt"
    GENOME_DIR = "/local/projects-t3/lilab/vmenon/CRISPRi/Genome/chromosomes/" # Must be your absolute path
    
    # Execute
    generate_23nt_input_from_rna(
        csv_file_path=INPUT_CSV,
        output_txt_path=OUTPUT_TXT,
        genome_dir_path=GENOME_DIR
    )
