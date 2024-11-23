#!/bin/bash

GFF_FILE="SalmonAnnotations.gff"       
CLUSTERS_FILE="multi_locus_clusters.csv"  
OUTPUT_FILE="multi_locus_genes.gff"    

echo "## Genes in multi-locus FST clusters ##" > $OUTPUT_FILE

tail -n +2 $CLUSTERS_FILE | while IFS=',' read -r CHROM BIN_START BIN_END; do
    if [[ ! "$CHROM" =~ ^chr ]]; then
        CHROM="chr_${CHROM}"
    fi

    awk -v chr="$CHROM" -v start="$BIN_START" -v end="$BIN_END" \
        '$1 == chr && $3 == "gene" && $4 <= end && $5 >= start' $GFF_FILE >> $OUTPUT_FILE
done

echo "Gene extraction complete. Results saved to $OUTPUT_FILE."
