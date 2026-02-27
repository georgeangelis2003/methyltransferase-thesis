import csv

uniprot_file = '/home/angelis/thesis/reference_proteomes/foldseek/active_sites.tsv'
foldmason_file = '/home/angelis/thesis/reference_proteomes/foldseek/foldmason_active_sites.tsv'
output_file = '/home/angelis/thesis/reference_proteomes/foldseek/final_msta_table.tsv'

fold_data = {}
with open(foldmason_file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    next(reader)
    for row in reader:
        if len(row) >= 2:
            fold_data[row[0]] = row[1]

with open(uniprot_file, 'r') as f_in, open(output_file, 'w', newline='') as f_out:
    reader = csv.reader(f_in, delimiter='\t')
    writer = csv.writer(f_out, delimiter='\t')
    
    next(reader)
    writer.writerow(['Accession', 'UniProt_Pos', 'UniProt_AA', 'FoldMason_AA_4101'])
    
    for row in reader:
        acc = row[0]
        pos = row[1]
        aa = row[2]
        
        fold_aa = fold_data.get(acc, 'NOT_FOUND')

        writer.writerow([acc, pos, aa, fold_aa])

print(f"Done! Created: {output_file}")