from Bio import SeqIO
import csv
from collections import Counter

alignment_file = '/home/angelis/thesis/reference_proteomes/foldseek/foldmason_msta_aa.fa'
output_csv = '/home/angelis/thesis/reference_proteomes/foldseek/foldmason_active_sites.csv'
column_number = 4101
allign_pos_index = column_number - 1

print(f"Starting scan of alignment at column {column_number}")

stats = Counter()

with open(output_csv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    writer.writerow(['Accession', 'Amino_Acid'])
    
    for record in SeqIO.parse(alignment_file, "fasta"):
        parts = record.id.split('|')
        acc = parts[1] if len(parts) > 1 else record.id
        
        try:
            residue = record.seq[allign_pos_index].upper()
            stats[residue] += 1
        except IndexError:
            residue = 'INDEX_ERROR'
            stats['ERROR'] += 1
            
        writer.writerow([acc, residue])

print(f"\nFinished - Results saved to: {output_csv}")

print("\nResidue Distribution at Active Site:")
for res, count in stats.most_common():
    percentage = (count / sum(stats.values())) * 100
    print(f"  {res}: {count} ({percentage:.2f}%)")