import requests
import time
import csv
import os

PDB_FOLDER = '/home/angelis/thesis/reference_proteomes/foldseek/structures_over80pLDDT'
OUTPUT_CSV = '/home/angelis/thesis/reference_proteomes/foldseek/only_active_sites.tsv'

def get_uniprot_data(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            return response.json()
    except:
        pass
    return None

pdb_files = [f for f in os.listdir(PDB_FOLDER) if f.endswith('.pdb')]
total = len(pdb_files)

with open(OUTPUT_CSV, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    writer.writerow(['Accession', 'Position', 'Amino_Acid', 'Type', 'Description'])
    
    for i, pdb_name in enumerate(pdb_files, 1):
        acc = pdb_name.replace('.pdb', '')
        print(f"[{i}/{total}] Checking {acc}...", end='\r')
        
        data = get_uniprot_data(acc)
        found_any = False
        
        if data and 'features' in data:
            sequence = data.get('sequence', {}).get('value', '')
            for feat in data['features']:
                f_type = feat.get('type')
                
                if f_type == 'Active site':
                    pos = feat['location']['start']['value']
                    desc = feat.get('description', 'no_desc')
                    aa = sequence[pos-1] if sequence and 0 < pos <= len(sequence) else "NA"
                    
                    writer.writerow([acc, pos, aa, f_type, desc])
                    found_any = True
            
        if not found_any:
            writer.writerow([acc, 'NA', 'NA', 'noinfo', 'No active site annotation found'])
            
        csvfile.flush()
        time.sleep(0.1)

print(f"\nFinished! Results in {OUTPUT_CSV}")