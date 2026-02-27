from Bio import AlignIO
import csv
from collections import Counter

alignment_file = '/home/angelis/thesis/reference_proteomes/new_filtered_diamond_tsv/mafft_output_aligned.fasta'
output_csv = '/home/angelis/thesis/reference_proteomes/active_sites_mafft.csv'
mafft_pos_index = 14633

print(f"Starting scan of alignment at column {mafft_pos_index}...")

stats = Counter()

with open(output_csv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Accession', 'Residue_at_Site'])
    
    alignment = AlignIO.read(alignment_file, "fasta")
    
    for record in alignment:
        parts = record.id.split('|')
        acc = parts[1] if len(parts) > 1 else record.id
        
        try:
            residue = record.seq[mafft_pos_index].upper()
            stats[residue] += 1
        except IndexError:
            residue = 'INDEX_ERROR'
            stats['ERROR'] += 1
            
        writer.writerow([acc, residue])

print(f"Finished - Results saved to: {output_csv}")

print("Residue Distribution at Active Site:")
for res, count in stats.most_common():
    print(f"  {res}: {count}")