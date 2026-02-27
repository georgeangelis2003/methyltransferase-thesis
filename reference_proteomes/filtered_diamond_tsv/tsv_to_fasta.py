import pandas as pd
import glob
import os

ID_COLUMN = 'subject_id'

SEQUENCE_COLUMN = 'trimmed_sequence'

OUTPUT_DIR = 'fasta_for_mafft'

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

tsv_files = glob.glob('*_dmnd_e02_filtered.tsv')

if not tsv_files:
    print("Δεν βρέθηκαν αρχεία .tsv με το μοτίβο *_dmnd_e02_filtered.tsv.")
else:
    print(f"Βρέθηκαν {len(tsv_files)} αρχεία για μετατροπή")
    
    for tsv_file in tsv_files:
        try:
            df = pd.read_csv(tsv_file, sep='\t')
            
            base_name = tsv_file.split('_')[0]
            fasta_filename = os.path.join(OUTPUT_DIR, f"{base_name}.fasta")

            with open(fasta_filename, 'w') as f_out:
                for index, row in df.iterrows():
  
                    seq_id = row[ID_COLUMN]
                    sequence = row[SEQUENCE_COLUMN]

                    f_out.write(f">{seq_id}\n")
                    f_out.write(f"{sequence}\n")
            
            print(f"Μετατροπή επιτυχής: {fasta_filename}")

        except Exception as e:
            print(f"Σφάλμα κατά την επεξεργασία του αρχείου {tsv_file}: {e}")

    print("Όλα τα αρχεία μετατράπηκαν.")