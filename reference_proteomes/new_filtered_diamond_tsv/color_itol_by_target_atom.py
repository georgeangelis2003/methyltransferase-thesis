import requests
import time

with open("/home/angelis/thesis/reference_proteomes/new_filtered_diamond_tsv/clean_ids.txt", "r") as f:
    protein_ids = [line.strip() for line in f if line.strip()]

CATEGORIES = {
    "O-MT": "#e74c3c",
    "N-MT": "#3498db",
    "C-MT": "#2ecc71",
    "S-MT": "#f39c12",
}

NAME_PATTERNS = [
    ("o-methyltransferase", "O-MT"),
    ("n-methyltransferase", "N-MT"),
    ("c-methyltransferase", "C-MT"),
    ("s-methyltransferase", "S-MT"),
]

def get_enzyme_category(protein_id):
    accession = protein_id.split('|')[1] if '|' in protein_id else protein_id
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            desc = data.get('proteinDescription', {})
            rec = desc.get('recommendedName', {})
            name = rec.get('fullName', {}).get('value', '')
            if not name:
                subs = desc.get('submissionNames', [])
                if subs:
                    name = subs[0].get('fullName', {}).get('value', '')
            name = name.lower()
            for pattern, cat in NAME_PATTERNS:
                if pattern in name:
                    return cat, CATEGORIES[cat]
        return None, None
    except requests.exceptions.RequestException as e:
        print(f"  Request error: {e}")
        return None, None


header = "TREE_COLORS\nSEPARATOR TAB\nDATA\n"

with open("itol_methyltransferase_annotations.txt", "w") as out:
    out.write(header)
    total = len(protein_ids)
    found = 0

    for i, pid in enumerate(protein_ids, 1):
        category, color = get_enzyme_category(pid)

        if category:
            species = pid.split('_')[-1] if '|' not in pid else pid.split('_')[-1]
            out.write(f"{pid}\trange\t{color}\t{species}\n")
            found += 1
            print(f"[{i}/{total}] {pid[:20]} -> {category}")
        else:
            print(f"[{i}/{total}] {pid[:20]} -> uncolored")

        time.sleep(0.3)

print(f"\nDone! {found}/{total} proteins categorized.")