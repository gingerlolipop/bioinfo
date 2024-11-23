#!/bin/bash

GFF_FILE="SalmonAnnotations.gff"       
CSV_FILE="high_fst_regions.csv"        
OUTPUT_FILE="high_fst_genes.gff"      

echo "## Extracted genes overlapping high FST regions ##" > $OUTPUT_FILE

tail -n +2 $CSV_FILE | while IFS=',' read -r CHROM BIN_START BIN_END N_VARIANTS WEIGHTED_FST MEAN_FST CHR CUMPOS; do
    if [[ ! "$CHROM" =~ ^chr ]]; then
        CHROM="chr_${CHROM}"
    fi
    
    awk -v chr="$CHROM" -v start="$BIN_START" -v end="$BIN_END" \
        '$1 == chr && $3 == "gene" && $4 <= end && $5 >= start' $GFF_FILE >> $OUTPUT_FILE
done

echo "Gene extraction complete. Results saved to $OUTPUT_FILE."
