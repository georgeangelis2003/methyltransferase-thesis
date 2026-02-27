import requests
import time
import csv
import os

output_file = '/home/angelis/thesis/reference_proteomes/active_sites.csv'
fasta_file = '/home/angelis/thesis/reference_proteomes/new_filtered_diamond_tsv/total_fastas_unique.fasta'

sequence_dict = {}
with open(fasta_file, 'r') as f:
    current_id = None
    current_seq = []
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if current_id:
                sequence_dict[current_id] = "".join(current_seq)
            parts = line.split('|')
            current_id = parts[1] if len(parts) > 1 else None
            current_seq = []
        else:
            current_seq.append(line)
    if current_id: 
        sequence_dict[current_id] = "".join(current_seq)

with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Accession', 'Position', 'Amino_Acid', 'Description'])
    
    accessions = list(sequence_dict.keys())

    for acc in accessions:
        print(f'Processing: {acc}')
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                features = data.get('features', [])
                found_site = False
                
                for feat in features:
                    if feat.get('type') == 'Active site':
                        pos = feat['location']['start']['value']
                        desc = feat.get('description', 'noinfo')
                        
                        protein_seq = sequence_dict.get(acc, "")
                        if protein_seq and len(protein_seq) >= pos:
                            amino_acid = protein_seq[pos - 1] 
                        else:
                            amino_acid = "Seq_Error"

                        writer.writerow([acc, pos, amino_acid, desc])
                        found_site = True
                        
                if not found_site:
                    writer.writerow([acc, 'noinfo', 'noinfo', 'noinfo'])
            else:
                writer.writerow([acc, 'API_error', 'API_error', 'API_error'])
        
        except Exception as e:
            writer.writerow([acc, 'Error', str(e), 'Error'])
            
        time.sleep(0.1)

print('Finished - check active_sites.csv')