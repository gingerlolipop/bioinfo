#!/bin/bash  

# Pipeline Logic: Genome alignment -> Mark Duplicates -> Variant Calling -> Variant Filtration -> FST Analysis

### 0. Set up paths, variables, and groups
base_dir="/home/abcde/as"
fpath="${base_dir}/fastq"
ref="${base_dir}/SalmonReference.fasta"
anno="${base_dir}/SalmonAnnotations.gff"
out_dir="${base_dir}/results"
vcf_dir="${base_dir}/vcf_analysis"
dict_dir="${base_dir}"

mkdir -p $out_dir
mkdir -p "${vcf_dir}/individual_vcfs"  
mkdir -p "${vcf_dir}/merged_vcfs"      
mkdir -p "${vcf_dir}/fst_results"      

# Use readarray for better array handling
readarray -t red_r1 < <(find "$fpath" -name "*red*R1*.fastq.gz")
readarray -t red_r2 < <(find "$fpath" -name "*red*R2*.fastq.gz")
readarray -t white_r1 < <(find "$fpath" -name "*white*R1*.fastq.gz")
readarray -t white_r2 < <(find "$fpath" -name "*white*R2*.fastq.gz")

### 1. Index reference and prepare dictionary
echo "Preparing reference genome..."
bwa index $ref
samtools faidx $ref

echo "Creating sequence dictionary..."
java -jar /mnt/bin/picard.jar CreateSequenceDictionary \
    R=$ref \
    O="${dict_dir}/SalmonReference.dict"

### 2. Alignment and sorting
process_samples() {
    local r1_files=("${!1}")
    local color="$2"
    
    echo "Processing ${color} samples..."
    for r1 in "${r1_files[@]}"; do
        r2="${r1/R1/R2}"
        sample_name=$(basename "$r1" | sed 's/_R1.*//')
        echo "Aligning $sample_name..."
        bwa mem -t 4 -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA" \
            "$ref" "$r1" "$r2" | \
            samtools view -bS -q 20 - | \
            samtools sort -o "${out_dir}/${sample_name}_aligned_sorted.bam"
        samtools index "${out_dir}/${sample_name}_aligned_sorted.bam"
    done
}

process_samples "red_r1[@]" "red"
process_samples "white_r1[@]" "white"

### 3. Mark duplicates
process_duplicates() {
    local color="$1"
    echo "Marking duplicates for ${color} samples..."
    for bam_file in "${out_dir}"/*${color}*_aligned_sorted.bam; do
        sample_name=$(basename "$bam_file" "_aligned_sorted.bam")
        java -jar $picard MarkDuplicates \
            I="$bam_file" \
            O="${out_dir}/${sample_name}_duplicates.bam" \
            M="${out_dir}/${sample_name}_marked_dup_metrics.txt" \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT
        samtools index "${out_dir}/${sample_name}_duplicates.bam"
    done
}

picard=/mnt/bin/picard.jar
process_duplicates "red"
process_duplicates "white"

### 4. Variant calling and filtering
gatk="java -jar /mnt/bin/gatk-4.2.6.0/gatk-package-4.2.6.0-local.jar"

# Call variants and create initial VCF
echo "Creating merged VCF..."
$gatk HaplotypeCaller \
    -R "$ref" \
    $(printf -- '-I %s ' "${out_dir}"/*_duplicates.bam) \
    --standard-min-confidence-threshold-for-calling 30 \
    -O "${vcf_dir}/merged_vcfs/all_merged.vcf.gz"

# Basic variant filtering
$gatk VariantFiltration \
    -R "$ref" \
    -V "${vcf_dir}/merged_vcfs/all_merged.vcf.gz" \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
    --filter-name "basic_filter" \
    -O "${vcf_dir}/merged_vcfs/all_merged_filtered.vcf.gz"

# Fix chromosome names
zcat "${vcf_dir}/merged_vcfs/all_merged_filtered.vcf.gz" | \
    sed 's/^chr_//g' | \
    bgzip > "${vcf_dir}/merged_vcfs/all_merged_fixed.vcf.gz"
tabix -p vcf "${vcf_dir}/merged_vcfs/all_merged_fixed.vcf.gz"

# Convert to PLINK format
plink=/mnt/bin/plink
$plink --vcf "${vcf_dir}/merged_vcfs/all_merged_fixed.vcf.gz" \
    --make-bed \
    --out "${vcf_dir}/merged_vcfs/salmon" \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --double-id \
    --maf 0.05 \
    --geno 0.1 \
    --mind 0.1

# Create cluster file - use names to identify
echo -e "FID\tIID\tCLUSTER" > "${vcf_dir}/fst_results/clusters.txt"
for sample in $(awk '{print $2}' "${vcf_dir}/merged_vcfs/salmon.fam"); do
    if [[ $sample == *"red"* ]]; then
        echo -e "$sample\t$sample\tred" >> "${vcf_dir}/fst_results/clusters.txt"
    else
        echo -e "$sample\t$sample\twhite" >> "${vcf_dir}/fst_results/clusters.txt"
    fi
done

### 5. FST Analysis
echo "Preparing for FST analysis..."
fst_dir="${vcf_dir}/fst_results"

# Clear and recreate population files
> "${fst_dir}/red.samples"
> "${fst_dir}/white.samples"

# Create population files
for bam in "${out_dir}"/*_duplicates.bam; do
    sample=$(basename "$bam" _duplicates.bam)
    if [[ $sample == *"red"* ]]; then
        echo "$sample" >> "${fst_dir}/red.samples"
    else
        echo "$sample" >> "${fst_dir}/white.samples"
    fi
done

# Verify sample files
echo "Checking sample files..."
echo "Red samples:"
cat "${fst_dir}/red.samples"
echo -e "\nWhite samples:"
cat "${fst_dir}/white.samples"

# Ensure VCF is indexed
echo "Indexing VCF file..."
tabix -p vcf "${vcf_dir}/merged_vcfs/all_merged_fixed.vcf.gz"

# Calculate per-site FST first
echo "Calculating per-site FST..."
vcftools --gzvcf "${vcf_dir}/merged_vcfs/all_merged_fixed.vcf.gz" \
    --weir-fst-pop "${fst_dir}/red.samples" \
    --weir-fst-pop "${fst_dir}/white.samples" \
    --out "${fst_dir}/red_vs_white" \
    --min-alleles 2

# Calculate windowed FST - have half window overlap
echo "Calculating windowed FST..."
vcftools --gzvcf "${vcf_dir}/merged_vcfs/all_merged_fixed.vcf.gz" \
    --weir-fst-pop "${fst_dir}/red.samples" \
    --weir-fst-pop "${fst_dir}/white.samples" \
    --out "${fst_dir}/red_vs_white_windowed" \
    --fst-window-size 10000 \
    --fst-window-step 5000 \
    --min-alleles 2

# Check results
echo "Checking FST results..."
for file in "${fst_dir}"/red_vs_white*.fst; do
    echo -e "\nContents of $file:"
    head -n 5 "$file"
done