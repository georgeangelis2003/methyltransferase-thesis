import pandas as pd
import glob
import os

ID_COLUMN = 'subject_id'
SEQUENCE_COLUMN = 'subject_sequence'

# Paths
INPUT_DIR = '/home/angelis/thesis/reference_proteomes/new_filtered_diamond_tsv'
QUERY_DIR = '/home/angelis/thesis/fasta_files'
OUTPUT_DIR = os.path.join(INPUT_DIR, 'fasta_for_mafft')

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

tsv_files = glob.glob(os.path.join(INPUT_DIR, '*_filtered.tsv'))

if not tsv_files:
    print(f"No files found in {INPUT_DIR}")
else:
    for tsv_file in tsv_files:
        try:
            # 1. Identify the base ID (e.g., A0A0K1YW34)
            base_name = os.path.basename(tsv_file).split('_')[0]
            
            # 2. Look for the reference query file
            query_fasta_path = os.path.join(QUERY_DIR, f"{base_name}.fasta")
            
            # Load the TSV data
            df = pd.read_csv(tsv_file, sep='\t')
            fasta_filename = os.path.join(OUTPUT_DIR, f"{base_name}.fasta")

            with open(fasta_filename, 'w') as f_out:
                # 3. First, copy the content of the reference query FASTA
                if os.path.exists(query_fasta_path):
                    with open(query_fasta_path, 'r') as f_ref:
                        ref_content = f_ref.read()
                        # Ensure it ends with a newline
                        f_out.write(ref_content.strip() + "\n")
                else:
                    print(f"Warning: Reference file {query_fasta_path} not found.")

                # 4. Then, append the subjects from the TSV
                for index, row in df.iterrows():
                    seq_id = row[ID_COLUMN]
                    sequence = row[SEQUENCE_COLUMN]
                    if pd.notna(sequence):
                        f_out.write(f">{seq_id}\n")
                        f_out.write(f"{sequence}\n")
            
            print(f"Converted: {base_name}.fasta (Reference sequence included)")

        except Exception as e:
            print(f"Error in {tsv_file}: {e}")

    print("Done.")