import os
import io
pdb_dir = '/home/angelis/thesis/reference_proteomes/foldseek/structures'
filtered_dir = '/home/angelis/thesis/reference_proteomes/foldseek/structures_over80pLDDT'
os.makedirs(filtered_dir, exist_ok=True)

for pdb_file in os.listdir(pdb_dir):
    if not pdb_file.endswith('.pdb'):
        continue
    
    plddts = []
    with open(os.path.join(pdb_dir, pdb_file), 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                try:
                    plddt = float(line[60:66].strip())
                    plddts.append(plddt)
                except:
                    continue
    if plddts:
        avg_plddt = sum(plddts) / len(plddts)
        if avg_plddt >= 80.0:
            os.system(f"cp {os.path.join(pdb_dir, pdb_file)} {os.path.join(filtered_dir, pdb_file)}")
            
print (f'filtering complete. Structures with average pLDDT >= 80 have been copied to {filtered_dir}')
