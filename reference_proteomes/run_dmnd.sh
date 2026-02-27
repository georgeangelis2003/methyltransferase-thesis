

for f in /home/angelis/thesis/fasta_files/*.fasta; do diamond blastp -q "$f" --db uniprot_db.dmnd --out "$(basename "$f" .fasta)_dmnd_e02.tsv" --outfmt 6 --max-target-seqs 2000 --evalue 1e-02 --threads 8; done