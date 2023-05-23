#!/bin/bash

# Step 1: Download fastq files from SRA
# Replace 'file.fastq' with the actual file name
fastq_file="file.fastq"
# Add commands to download the fastq file from SRA
# For example: wget https://www.ncbi.nlm.nih.gov/sra/fastq-file.fastq.gz

# Step 2: Check quality with fastqc
fastqc "$fastq_file"

# Step 3: Download and install Trimmomatic
trimmomatic_url="http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip"
trimmomatic_dir="Trimmomatic-0.39"
trimmomatic_jar="tools/trimmomatic.jar"

# Download Trimmomatic
wget "$trimmomatic_url"

# Unzip the downloaded file
unzip "$trimmomatic_dir.zip"

# Move Trimmomatic jar file to the tools folder
mv "$trimmomatic_dir/trimmomatic-0.39.jar" "$trimmomatic_jar"

# Make Trimmomatic jar file executable
chmod +x "$trimmomatic_jar"

# Test Trimmomatic installation
java -jar "$trimmomatic_jar"

# Step 4: Index the human genome
reference_genome="reference2/Homo_sapiens_assembly38.fasta"

# Step 5: Alignment
aligned_sam="aligned.sam"
bwa mem "$reference_genome" "$fastq_file" > "$aligned_sam"

# Step 6: Conversion to BAM format
aligned_bam="aligned.bam"
samtools view -b "$aligned_sam" > "$aligned_bam"

# Step 7: Sorting and indexing the BAM file
sorted_bam="sorted_reads.bam"
samtools sort "$aligned_bam" -o "$sorted_bam"
samtools index "$sorted_bam"

# Step 8: Variant calling
variants_vcf="variants.vcf"
bcftools mpileup -f "$reference_genome" "$sorted_bam" | bcftools call -mv -Ov > "$variants_vcf"

# Step 9: Annotation
output_avinput="output.avinput"
annovar_dir="tools/annovar"
convert2annovar.pl -format vcf4 "$variants_vcf" > "$output_avinput"
"$annovar_dir/annotate_variation.pl" -geneanno -dbtype refGene "$output_avinput" "$annovar_dir/humandb/"

# Step 10: Grep
lynch_positions_file="Lynch_positions_in_file.txt"
output_file="output_file.hg38_multianno.txt"
grep -E 'MLH1|MSH2|MSH6|PMS2|EPCAM' "$output_file" > "$lynch_positions_file"
