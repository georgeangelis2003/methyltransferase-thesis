#!/bin/bash


OUTPUT_DIR="/home/angelis/thesis/reference_proteomes" 


for f in /home/angelis/thesis/reference_proteomes/dmnd_tsv_files/*.tsv; do 
    
 
    OUTPUT_FILENAME="$(basename "$f" .tsv)_enriched.tsv"
    
  
    OUTPUT_PATH="${OUTPUT_DIR}/${OUTPUT_FILENAME}"
    
    # 3. Έλεγχος αν το αρχείο εξόδου υπάρχει ήδη
    if [ -f "$OUTPUT_PATH" ]; then
        echo "Skipping $f: Output $OUTPUT_FILENAME already exists."
        continue # Παρακάμπτει την τρέχουσα επανάληψη και συνεχίζει με το επόμενο αρχείο
    fi
    
    
    echo "Processing $f..."
    python diamond_blastp_enrich_output.py \
        -i "$f" \
        -d uniprot_all_proteomes.fasta \
        -o "$OUTPUT_PATH" \
        --check-pdb
        
    echo "Completed processing for $f."
done