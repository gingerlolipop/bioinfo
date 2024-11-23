import pandas as pd
import numpy as np
from pathlib import Path

def find_cooccurring_genes():
    # Input file paths
    base = Path(r"C:\Users\jillb\OneDrive - UBC\CONS 503A\Assignment")
    high_fst_file = base / "vcf_analysis/fst_results/high_fst_regions.csv"
    gff_file = base / "SalmonAnnotations.gff"
    output_file = base / "vcf_analysis/fst_results/cooccurring_genes.gff"

    # Load high FST regions and add debug info
    print("Loading FST data...")
    fst_df = pd.read_csv(high_fst_file)
    print("\nFST file head:")
    print(fst_df.head())
    print("\nFST value range:", fst_df['WEIGHTED_FST'].min(), "to", fst_df['WEIGHTED_FST'].max())

    # Filter significant high FST windows
    significant_windows = fst_df[fst_df["WEIGHTED_FST"] >= 0.05]
    print("\nNumber of significant windows:", len(significant_windows))
    print("Chromosome names in FST data:", fst_df['CHROM'].unique())

    # Dictionary to store genes by chromosome
    genes_by_chr = {}

    # Debug: Check GFF chromosome names
    print("\nChecking GFF file chromosome names...")
    seen_chroms = set()
    with open(gff_file, "r", encoding='utf-8') as gff:
        for line in gff:
            if line.startswith('#'):
                continue
            chrom = line.strip().split("\t")[0]
            seen_chroms.add(chrom)
            if len(seen_chroms) >= 5:
                break
    print("First few chromosome names in GFF:", seen_chroms)

    # First pass: collect all genes in high FST regions
    print("\nProcessing GFF file for overlaps...")
    overlap_count = 0
    with open(gff_file, "r", encoding='utf-8') as gff:
        for line in gff:
            if line.startswith('#'):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
                
            chrom, start, end = cols[0], int(cols[3]), int(cols[4])
            
            # Check for overlap with significant high FST windows
            for _, window in significant_windows.iterrows():
                window_chrom = f"chr_{window['CHROM']}"
                if chrom == window_chrom and start <= window["BIN_END"] and end >= window["BIN_START"]:
                    if chrom not in genes_by_chr:
                        genes_by_chr[chrom] = []
                    genes_by_chr[chrom].append({
                        'line': line,
                        'start': start,
                        'end': end,
                        'fst': window['WEIGHTED_FST']
                    })
                    overlap_count += 1
                    break

    print(f"\nFound {overlap_count} overlapping genes")

    # Write results
    with open(output_file, "w", encoding='utf-8') as out:
        out.write("## Cooccurring genes in high FST regions ##\n")
        
        # Write genes sorted by FST value within each chromosome
        for chrom in sorted(genes_by_chr.keys()):
            genes = genes_by_chr[chrom]
            if len(genes) > 1:  # Only write if multiple genes found
                out.write(f"# Chromosome {chrom}: {len(genes)} genes #\n")
                # Sort genes by FST value
                sorted_genes = sorted(genes, key=lambda x: x['fst'], reverse=True)
                for gene in sorted_genes:
                    out.write(gene['line'])
                out.write("\n")

    print(f"Cooccurring genes in high FST regions saved to {output_file}")

    # Return summary statistics
    return {
        'total_chromosomes': len(genes_by_chr),
        'genes_per_chr': {chrom: len(genes) for chrom, genes in genes_by_chr.items()},
        'total_genes': sum(len(genes) for genes in genes_by_chr.values())
    }

if __name__ == "__main__":
    stats = find_cooccurring_genes()
    print("\nAnalysis Summary:")
    print(f"Total chromosomes with high FST genes: {stats['total_chromosomes']}")
    print("\nGenes per chromosome:")
    for chrom, count in stats['genes_per_chr'].items():
        print(f"{chrom}: {count} genes")
    print(f"\nTotal genes found: {stats['total_genes']}")