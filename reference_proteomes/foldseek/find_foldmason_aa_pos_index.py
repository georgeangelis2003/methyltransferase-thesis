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
    ("A0A251QIF7", 266),
    ("A0A5D2TVI4", 269),
    ("A0A8B8QA76", 266),
    ('A0A8I6WST8', 264),
    ('A0AAV8T7Z1', 261)
]

for prot_id, pos in targets:
    find_col(MY_ALIGNMENT, prot_id, pos)