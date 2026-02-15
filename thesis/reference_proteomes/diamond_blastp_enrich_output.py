import sys
import argparse
import pandas as pd
import requests
from pathlib import Path
from Bio import SeqIO

def parse_fasta (fasta_files):
    '''search fasta and return dict of id:{sequence,length}'''
    seq_dict = {}
    for record in SeqIO.parse(fasta_files, 'fasta'):
        seq_dict[record.id] = {
            'sequence': str(record.seq),
            'length': len(record.seq),
            'description': record.description
        }
    return seq_dict

def check_pdb_structure (protein_id):
    '''checks if protein has pdb structure using Uniprot API'''
    uniprot_id = protein_id.split('|')[1] if '|' in protein_id else protein_id.split('_')[0]
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        response = requests.get(url, timeout=5)

        if response.status_code == 200:
            data = response.json()
            # Look for PDB cross-references
            if 'uniProtKBCrossReferences' in data:
                pdb_ids = []
                for ref in data['uniProtKBCrossReferences']:  
                    if ref['database'] == 'PDB':              
                       pdb_ids.append(ref['id'])
                return ','.join(pdb_ids) if pdb_ids else 'No'
        return 'Unknown'
    except:
        return 'Unknown'


check_pdb_structure("A0A0K1YW34")

def enrich_diamond_output(diamond_tsv, database_fasta, output_file, check_pdb=False):
    columns = [
        'query_id', 'subject_id', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]

    print("Reading DIAMOND output...")
    df = pd.read_csv(diamond_tsv, sep='\t', names=columns, comment='#') 
    #Διαβάζει το TSV ως pandas DataFrame. 
    #Χρησιμοποιεί τα ονόματα στηλών που όρισε παραπάνω. 
    #Αγνοεί γραμμές που αρχίζουν με '#'
    
    print(f"Parsing database FASTA: {database_fasta}")
    seq_dict = parse_fasta(database_fasta)

    print("Adding sequence information...")
    # Add sequence information
    df['subject_sequence'] = df['subject_id'].map(lambda x: seq_dict.get(x, {}).get('sequence', 'N/A'))
    df['subject_length'] = df['subject_id'].map(lambda x: seq_dict.get(x, {}).get('length', 0))
    df['subject_description'] = df['subject_id'].map(lambda x: seq_dict.get(x, {}).get('description', 'N/A'))
    
    if check_pdb:
        print('Checking for PDB structure')
        unique_subjects = df['subject_id'].unique()
        pdb_dict = {}
        
        for i, subj_id in enumerate(unique_subjects):
            if i % 100 == 0: #Εκτυπώνει πρόοδο κάθε 100 πρωτεΐνες
                print(f"  Processed {i}/{len(unique_subjects)} proteins") 
            pdb_dict[subj_id] = check_pdb_structure(subj_id)
        df['has_pdb_structure'] = df['subject_id'].map(pdb_dict)
    else:
        df['has_pdb_structure'] = 'Not_checked'
    #calculate alignment percentage
    df['query_coverage'] = ((df['qend'] - df['qstart'] + 1) / df['length'] * 100).round(2)
    df['subject_coverage'] = ((df['send'] - df['sstart'] + 1) / df['subject_length'] * 100).round(2)

    col_order = [
        'query_id', 'subject_id', 'subject_description', 'pident', 'evalue', 'bitscore',
        'subject_length', 'query_coverage', 'subject_coverage',
        'qstart', 'qend', 'sstart', 'send', 'has_pdb_structure', 'subject_sequence'
    ]
    
    df = df[col_order]
    
    print(f"Writing enriched output to: {output_file}")
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Processed {len(df)} alignments.")
    
    return df

def main():
    parser = argparse.ArgumentParser(
        description='Enrich DIAMOND output with sequence info and PDB structures'
    )
    parser.add_argument('-i', '--input', required=True,
                        help='DIAMOND output file (TSV format, outfmt 6)')
    parser.add_argument('-d', '--database', required=True,
                        help='FASTA file used to create the DIAMOND database')
    parser.add_argument('-o', '--output', required=True,
                        help='Output enriched TSV file')
    parser.add_argument('--check-pdb', action='store_true',
                        help='Check for PDB structures (WARNING: slow for many hits)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not Path(args.input).exists():
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)
    
    if not Path(args.database).exists():
        print(f"Error: Database FASTA not found: {args.database}")
        sys.exit(1)
    
    # Run enrichment
    enrich_diamond_output(
        args.input,
        args.database,
        args.output,
        args.check_pdb
    )

if __name__ == "__main__":
    main()