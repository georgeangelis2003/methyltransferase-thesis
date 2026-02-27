import pandas as pd
import os
from Bio import SeqIO
from pathlib import Path
import glob


FASTA_DIR = '/home/angelis/thesis/fasta_files'
INPUT_DIR = '/home/angelis/thesis/reference_proteomes/enriched_tsv_files'
OUTPUT_DIR = '/home/angelis/thesis/reference_proteomes/new_filtered_diamond_tsv'
LENGTH_THRESHOLD = 30

column_names = [
    'query_id', 'subject_id', 'subject_description', 'pident', 'evalue',
    'bitscore', 'subject_length', 'query_coverage', 'subject_coverage',
    'qstart', 'qend', 'sstart', 'send', 'has_pdb_structure', 'subject_sequence'
]

def get_query_length(fasta_path):
    """Διαβάζει την αλληλουχία από το FASTA και επιστρέφει το μήκος της."""
    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
        return len(record.seq)
    except:
        return None

def load_blast_results(input_file, col_names):
    """Διαβάζει το αρχείο TSV των αποτελεσμάτων BLAST σε ένα DataFrame."""
    try:
        df = pd.read_csv(
            input_file, 
            sep='\t', 
            header=None, 
            names=col_names
        )
        
       
        numeric_cols = ['pident', 'evalue', 'bitscore', 'subject_length', 
                       'query_coverage', 'subject_coverage', 'qstart', 
                       'qend', 'sstart', 'send']
        
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
    
        critical_cols = ['sstart', 'send', 'subject_length', 'subject_sequence']
        df_clean = df.dropna(subset=critical_cols)
        
        if len(df_clean) < len(df):
            print(f'  Προειδοποίηση: Αφαιρέθηκαν {len(df) - len(df_clean)} γραμμές με missing values.')
        
        print(f'  Loaded {len(df_clean)} valid lines.')
        return df_clean
    except FileNotFoundError:
        print(f'  ERROR: File not found')
        return None
    except Exception as e:
        print(f'  ERROR: {e}')
        return None

def filter_by_length(df, query_length, threshold=30):
    """Κρατάει μόνο τις αλληλουχίες που είναι query_length +/- threshold."""
    if query_length is None:
        return df
    
    min_len = query_length - threshold
    max_len = query_length + threshold
    
    initial_count = len(df)
    df_filtered = df[(df['subject_length'] >= min_len) & (df['subject_length'] <= max_len)].copy()
    
    print(f"  Φιλτράρισμα μήκους (Target: {query_length}): Κρατήθηκαν {len(df_filtered)} από {initial_count} hits.")
    return df_filtered

def remove_duplicate_sequences(df):
    print('  Removing duplicate sequences.')
    
    df_sorted = df.sort_values(by='evalue', ascending=True)
    
    df_unique = df_sorted.drop_duplicates(
        subset=['subject_sequence'], 
        keep='first'
    ).reset_index(drop=True)
    
    print(f"  Αρχικά hits: {len(df)}. Μοναδικές αλληλουχίες: {len(df_unique)}.")
    return df_unique

def process_single_file(input_tsv, output_tsv, col_names, threshold, protein_id):
    
    print(f"Processing: {protein_id}")
    
    fasta_path = os.path.join(FASTA_DIR, f"{protein_id}.fasta")
    query_length = get_query_length(fasta_path)

    df = load_blast_results(input_tsv, col_names)
    if df is None or len(df) == 0:
        print(f"  SKIPPED: No valid data in {input_tsv}")
        return False

    df_filtered_len = filter_by_length(df, query_length, threshold)
    
    df_filtered_len['Organism'] = df_filtered_len['subject_description'].str.extract(r'OS=(.*?) OX=')

    df_final_filtered = remove_duplicate_sequences(df_filtered_len)

    final_columns = [
        'query_id', 
        'subject_id', 
        'Organism', 
        'pident', 
        'evalue', 
        'bitscore',
        'subject_length',
        'query_coverage',
        'subject_coverage',
        'qstart',
        'qend',
        'sstart',
        'send',
        'has_pdb_structure',
        'subject_sequence' 
    ]
    
    df_final = df_final_filtered[final_columns]
    df_final.to_csv(output_tsv, sep='\t', index=False)
    
    print(f"Saved: {output_tsv}")
    print(f"  Final sequences (unique): {len(df_final)}")
    return True

def main():
    """
    Κύρια συνάρτηση που επεξεργάζεται όλα τα enriched.tsv files.
    """
  
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}")
    
   
    
    all_enriched = glob.glob(os.path.join(INPUT_DIR, '*_enriched.tsv'))
    

    input_files = all_enriched 
    
    
    print(f"\nFound {len(input_files)} enriched.tsv files to process.")
    
    print(f"(Excluded {len(all_enriched) - len(input_files)} _dmnd_e02_enriched.tsv files)")
    
    if len(input_files) == 0:
        print("No files found! Check the INPUT_DIR path.")
        return
    
    success_count = 0
    failed_count = 0
    
    for input_file in sorted(input_files):
  
        basename = os.path.basename(input_file)
        protein_id = basename.replace('_enriched.tsv', '')
        

        output_filename = f"{protein_id}_filtered.tsv"
        output_file = os.path.join(OUTPUT_DIR, output_filename)
        

        try:
            if process_single_file(input_file, output_file, column_names, LENGTH_THRESHOLD, protein_id):
                success_count += 1
            else:
                failed_count += 1
        except Exception as e:
            print(f"  ERROR processing {protein_id}: {e}")
            failed_count += 1
    

    print(f"BATCH PROCESSING COMPLETE")
    print(f" Successfully processed: {success_count}")
    print(f" Failed/Skipped: {failed_count}")
    print(f"Total files: {len(input_files)}")
    print(f"\nFiltered files saved in: {OUTPUT_DIR}")

if __name__ == '__main__':
    main()