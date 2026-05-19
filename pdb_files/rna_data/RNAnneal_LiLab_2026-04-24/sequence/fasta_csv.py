import sys
import csv

def fasta_to_csv(fasta_path, csv_path):
    with open(fasta_path, 'r') as f_in, open(csv_path, 'w', newline='') as f_out:
        writer = csv.writer(f_out)
        # Write CSV header
        writer.writerow(['sgRNA', 'Sequence'])
        
        current_id = None
        current_seq = []
        
        for line in f_in:
            line = line.strip()
            if line.startswith('>'):
                # Process the previous entry if it exists
                if current_id is not None:
                    writer.writerow([current_id, ''.join(current_seq)])
                    current_seq = []
                # Clean up the identifier by stripping the '>' character
                current_id = line[1:]
            else:
                current_seq.append(line)
                
        # Flush the last sequence in the file
        if current_id is not None:
            writer.writerow([current_id, ''.join(current_seq)])
            
    print(f"Successfully converted {fasta_path} to {csv_path}")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python fasta_to_csv.py <input.fasta> <output.csv>")
        sys.exit(1)
        
    fasta_file = sys.argv[1]
    csv_file = sys.argv[2]
    
    fasta_to_csv(fasta_file, csv_file)
