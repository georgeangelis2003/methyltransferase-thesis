import requests
import time

with open("clean_ids.txt") as f:
    protein_ids = [line.strip() for line in f if line.strip()]

for pid in protein_ids:
    accession = pid.split('|')[1] if '|' in pid else pid
    data = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.json", timeout=10).json()
    refs = data.get('uniProtKBCrossReferences', [])
    go_terms = {r.get('id') for r in refs if r.get('database') == 'GO'}
    if 'GO:0016279' in go_terms or 'GO:0008757' in go_terms:
        print(f"N-MT: {pid}")
    time.sleep(0.3)

print("Done.")
