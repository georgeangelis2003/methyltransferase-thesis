import sys
from Bio import SeqIO

def find_col(fasta_path, target_id, target_pos):
    target_pos = int(target_pos)
    
    with open(fasta_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if target_id in record.id:
                seq = str(record.seq)
                current_ungapped = 0
                
                for i, char in enumerate(seq):
                    if char != "-":
                        current_ungapped += 1
                        
                    if current_ungapped == target_pos:
                        print(f"--- Result for {target_id} ---")
                        print(f"Amino Acid:       {char}")
                        print(f"Original Pos:     {target_pos}")
                        print(f"Alignment Column: {i + 1}")
                        return
        
        print(f"Error: Could not find position {target_pos} for ID {target_id}")

MY_ALIGNMENT = "/home/angelis/thesis/reference_proteomes/foldseek/foldmason_msta.fasta_aa.fa"


targets = [
    ("A0A0K1YW34", 268),
    ("A0AAD8M673", 269),
    ("A0AAD8N2I6", 271),
    ('A0AAD1Z3H2', 258),
    ('A0A2I4ECM1', 274)
]

for prot_id, pos in targets:
    find_col(MY_ALIGNMENT, prot_id, pos)