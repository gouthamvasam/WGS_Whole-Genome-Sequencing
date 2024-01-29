#!/bin/bash

# Define paths to input and output directories
raw_fastq_dir="/path/to/raw_fastq"
output_dir="/path/to/output"
reference_genome="/path/to/reference_genome"
known_sites="/path/to/known_sites" # e.g., dbSNP, Mills, 1000G indels
adapter_sequences="/path/to/adapter_sequences.fa"
annovar_database_dir="/path/to/annovar_database"

# Ensure output directory exists
mkdir -p "$output_dir"

# Programs needed: FastQC, Cutadapt, Trimmomatic, BWA, GATK, ANNOVAR, samtools
# Ensure these are installed and available in your PATH

# Step 1: Raw Data Quality Control with FastQC
fastqc --outdir "$output_dir/fastqc_results" "$raw_fastq_dir"/*.fastq.gz

# Step 2: Data Preprocessing with Cutadapt and Trimmomatic
# Adapter removal with Cutadapt (example command, adjust parameters as needed)
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o "$output_dir/sample_R1_trimmed.fastq.gz" -p "$output_dir/sample_R2_trimmed.fastq.gz" "$raw_fastq_dir/sample_R1.fastq.gz" "$raw_fastq_dir/sample_R2.fastq.gz"

# Quality and adapter trimming with Trimmomatic
trimmomatic PE -phred33 \
  "$output_dir/sample_R1_trimmed.fastq.gz" \
  "$output_dir/sample_R2_trimmed.fastq.gz" \
  "$output_dir/sample_R1_paired_trimmed.fastq.gz" \
  "$output_dir/sample_R1_unpaired_trimmed.fastq.gz" \
  "$output_dir/sample_R2_paired_trimmed.fastq.gz" \
  "$output_dir/sample_R2_unpaired_trimmed.fastq.gz" \
  ILLUMINACLIP:"$adapter_sequences":2:30:10:8:true \
  SLIDINGWINDOW:4:20 MINLEN:36

# Step 3: Sequence Alignment with BWA
bwa mem -t 4 -R "@RG\tID:sample\tSM:sample\tPL:ILLUMINA" \
  "$reference_genome" \
  "$output_dir/sample_R1_paired_trimmed.fastq.gz" \
  "$output_dir/sample_R2_paired_trimmed.fastq.gz" | \
  samtools sort -o "$output_dir/sample_sorted.bam" -

# Step 4: Post-alignment Processing with GATK
# Mark duplicates
gatk MarkDuplicates \
  -I "$output_dir/sample_sorted.bam" \
  -O "$output_dir/sample_marked.bam" \
  -M "$output_dir/sample_metrics.txt"

# Index the BAM file
samtools index "$output_dir/sample_marked.bam"

# Base recalibration
gatk BaseRecalibrator \
  -I "$output_dir/sample_marked.bam" \
  -R "$reference_genome" \
  --known-sites "$known_sites" \
  -O "$output_dir/recal_data.table"

gatk ApplyBQSR \
  -R "$reference_genome" \
  -I "$output_dir/sample_marked.bam" \
  --bqsr-recal-file "$output_dir/recal_data.table" \
  -O "$output_dir/sample_final.bam"

# Step 5: Variant Calling with GATK HaplotypeCaller
gatk HaplotypeCaller \
  -R "$reference_genome" \
  -I "$output_dir/sample_final.bam" \
  -O "$output_dir/sample_raw_variants.vcf"

# Step 6: Variant Annotation and Filtering with ANNOVAR
# Convert VCF to ANNOVAR input format
convert2annovar.pl -format vcf4 "$output_dir/sample_raw_variants.vcf" > "$output_dir/sample_raw_variants.avinput"

# Annotate variants using ANNOVAR
table_annovar.pl "$output_dir/sample_raw_variants.avinput" "$annovar_database_dir" \
  -buildver hg19 \
  -out "$output_dir/sample_annotated" \
  -remove \
  -protocol refGene,cytoBand,exac03,dbnsfp33a \
  -operation g,r,f,f \
  -nastring . \
  -csvout \
  -polish \
  -xref /path/to/gene_xref.txt

# The above command annotates the variants using several databases and saves the output in CSV format.
# Ensure that the gene_xref.txt file and other database files are present in the specified ANNOVAR database directory.

# Additional filtering based on ANNOVAR annotations can be done here
# For example, filtering out common variants with high allele frequency in exac03
awk -F, '$10 < 0.01 || $10 == "."' "$output_dir/sample_annotated.hg19_multianno.csv" > "$output_dir/sample_filtered_variants.csv"

# The above command filters the annotated variants to retain those with an allele frequency less than 1% in the ExAC database or not present in the database.
# Adjust the column number ('$10') based on the actual column in 