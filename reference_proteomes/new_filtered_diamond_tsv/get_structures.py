import requests
import os
import time

INPUT_FILE = '/home/angelis/thesis/reference_proteomes/new_filtered_diamond_tsv/total_unique_uniprot_ids.txt'
OUTPUT_DIR = '/home/angelis/thesis/reference_proteomes/structures'

os.makedirs(OUTPUT_DIR, exist_ok=True)

def get_real_url(uniprot_id):
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        r = requests.get(api_url, timeout=10)
        if r.status_code == 200:
            data = r.json()
            if data and len(data) > 0:
                return data[0].get('pdbUrl')
    except:
        return None
    return None

def download_structure(uniprot_id):
    uniprot_id = uniprot_id.strip().upper()
    file_path = os.path.join(OUTPUT_DIR, f"{uniprot_id}.pdb")
    
    if os.path.exists(file_path):
        return "Exists"

    pdb_url = get_real_url(uniprot_id)
    
    if pdb_url:
        try:
            r = requests.get(pdb_url, timeout=15)
            if r.status_code == 200:
                with open(file_path, 'wb') as f:
                    f.write(r.content)
                return "AlphaFold OK"
        except:
            return "Download Error"
            
    return "Not in AlphaFold DB"

with open(INPUT_FILE, 'r') as f:
    ids = [line.strip() for line in f if line.strip()]

print(f"Starting API-based download...")

downloaded = 0
for i, uid in enumerate(ids, 1):
    status = download_structure(uid)
    if status == "AlphaFold OK":
        downloaded += 1
        print(f"[{i}/{len(ids)}] {uid}: {status} | Total: {downloaded}")
    elif i % 50 == 0:
        print(f"Checked {i}/{len(ids)} IDs... Found: {downloaded}")
    
    time.sleep(0.1)