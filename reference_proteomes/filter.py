import pandas as pd
import os
from Bio import SeqIO
from pathlib import Path
import glob

# Ρυθμίσεις μονοπατιών
FASTA_DIR = '/home/angelis/thesis/fasta_files'
INPUT_DIR = '/home/angelis/thesis/reference_proteomes'
OUTPUT_DIR = '/home/angelis/thesis/reference_proteomes/filtered_diamond_tsv'
LEN_THRESHOLD = 30  # +/- αμινοξέα

column_names = [
    'query_id', 'subject_id', 'subject_description', 'pident', 'evalue',
    'bitscore', 'subject_length', 'query_coverage', 'subject_coverage',
    'qstart', 'qend', 'sstart', 'send', 'has_pdb_structure', 'subject_sequence'
]

def get_query_length(protein_id):
    """Βρίσκει το FASTA αρχείο και επιστρέφει το μήκος της query αλληλουχίας."""
    fasta_path = os.path.join(FASTA_DIR, f"{protein_id}.fasta")
    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
        return len(record.seq)
    except Exception as e:
        print(f"  ΣΦΑΛΜΑ: Δεν βρέθηκε το FASTA για το {protein_id} στο {fasta_path}")
        return None

def load_blast_results(input_file, col_names):
    """Διαβάζει το αρχείο TSV των αποτελεσμάτων BLAST."""
    try:
        df = pd.read_csv(input_file, sep='\t', header=None, names=col_names)
        
        numeric_cols = ['pident', 'evalue', 'bitscore', 'subject_length', 
                       'query_coverage', 'subject_coverage', 'qstart', 
                       'qend', 'sstart', 'send']
        
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        df_clean = df.dropna(subset=['subject_length', 'subject_sequence'])
        return df_clean
    except Exception as e:
        print(f'  ERROR: {e}')
        return None

def filter_by_length(df, query_length, threshold):
    """Κρατάει μόνο τις αλληλουχίες που είναι query_length +/- threshold."""
    if query_length is None:
        return df
    
    min_len = query_length - threshold
    max_len = query_length + threshold
    
    initial_count = len(df)
    df_filtered = df[(df['subject_length'] >= min_len) & (df['subject_length'] <= max_len)].copy()
    
    print(f"  Φιλτράρισμα μήκους (Target: {query_length}): Κρατήθηκαν {len(df_filtered)} από {initial_count} hits.")
    return df_filtered

def filter_by_organism(df):
    """Κρατάει το καλύτερο hit ανά οργανισμό."""
    df['Organism'] = df['subject_description'].str.extract(r'OS=(.*?) OX=')
    
    df_unique = df.sort_values(by='evalue', ascending=True).drop_duplicates(
        subset=['Organism'],
        keep='first'
    ).reset_index(drop=True)
    
    return df_unique

def process_single_file(input_tsv, output_tsv, col_names, threshold, protein_id):
    print(f"\n{'='*60}")
    print(f"Processing: {protein_id}")
    
    # 1. Φόρτωση δεδομένων
    df = load_blast_results(input_tsv, col_names)
    if df is None or len(df) == 0:
        return False

    # 2. Εύρεση μήκους Query
    q_len = get_query_length(protein_id)
    if q_len is None:
        print(f"  SKIPPING {protein_id} λόγω έλλειψης FASTA.")
        return False

    # 3. Φιλτράρισμα βάσει μήκους (+/- 30)
    df_len_filtered = filter_by_length(df, q_len, threshold)
    
    if len(df_len_filtered) == 0:
        print(f"Κανένα hit δεν ικανοποιεί το κριτήριο μήκους.")
        return False

    # 4. Φιλτράρισμα ανά οργανισμό
    df_final = filter_by_organism(df_len_filtered)

    # Επιλογή στηλών (Χρησιμοποιούμε την αρχική subject_sequence, όχι trimmed)
    final_columns = [
        'query_id', 'subject_id', 'Organism', 'pident', 'evalue', 
        'bitscore', 'subject_length', 'query_coverage', 'subject_coverage',
        'qstart', 'qend', 'sstart', 'send', 'has_pdb_structure', 'subject_sequence'
    ]
    
    df_final[final_columns].to_csv(output_tsv, sep='\t', index=False)
    print(f"  ✓ Αποθηκεύτηκαν {len(df_final)} αλληλουχίες.")
    return True

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    all_enriched = glob.glob(os.path.join(INPUT_DIR, '*_enriched.tsv'))
    
    if not all_enriched:
        print("No files found!")
        return

    for input_file in sorted(all_enriched):
        basename = os.path.basename(input_file)
        protein_id = basename.replace('_enriched.tsv', '')
        output_file = os.path.join(OUTPUT_DIR, f"{protein_id}_filtered.tsv")
        
        try:
            process_single_file(input_file, output_file, column_names, LEN_THRESHOLD, protein_id)
        except Exception as e:
            print(f"  ERROR processing {protein_id}: {e}")

if __name__ == '__main__':
    main()