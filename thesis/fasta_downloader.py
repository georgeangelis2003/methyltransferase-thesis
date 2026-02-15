import requests
from pathlib import Path
import time
import os
def download_fasta (uniprot_id, output_dir='fasta_files', skip_existing = True):

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    #Check if file already exists
    output_path = os.path.join(output_dir, f"{uniprot_id}.fasta")
    if skip_existing and os.path.exists(output_path):
        print(f"Skipping {uniprot_id}... (already exists)")
        return True

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        print (f'Downloading {uniprot_id}', end = '')
        response = requests.get(url, timeout=10)

        if response.status_code == 200:
            output_path = os.path.join(output_dir, f'{uniprot_id}.fasta')
            with open(output_path, 'w') as f:
                f.write (response.text)
                print(f'Saved to {output_path}')
                return True
        else:
            print(f"Error {response.status_code}")
            return False
    except requests.exceptions.RequestException as e:
        print(f" Network error: {e}")
        return False
        
def main():
    uniprot_ids = [
        'O04385', 'Q6VMW2', 'Q06YR0', 'I2FFE9', 'Q9FK25', 'Q42654', 'Q6VMV9', 'Q6ZD89', 'Q6VCW3',
         'C6TAY1', 'Q6VMV8', 'A0A1S2XJF8', 'Q1EDY7', 'O24529', 'Q84KK6', 'Q84KK4', 'Q29U70', 'Q8W013',
         'Q9ZTU2', 'P28002', 'A0A125SA05', 'A0A4P8DY91', 'A0A0K1YW34', 'Q9M602', 'P28002', 'Q9ZTU2',
         'Q6ZD89', 'Q84XW5', 'A0A166U5H3'
    ]
    output_dir = 'fasta_files'
    print(f'dowloading {len(uniprot_ids)} protein sequences')

    successful = 0
    failed = 0
    for id in uniprot_ids:
        if download_fasta(id.strip(), output_dir):
            successful += 1
        else:
            failed += 1
        time.sleep(0.5)
    
    print(f"Download complete: {successful} successful, {failed} failed")
    print(f"Files saved in '{output_dir}/' directory")

if __name__ == "__main__":
    main()