#!/bin/bash

# Step 1: Check the number of input files
if [ $# -eq 1 ]; then
  # Single-end data
  fastq_file="$1"
  trimmomatic_jar="/home/mariam/Lynch/tools/trimmomatic.jar"
  termo_out="output.fq.gz"
  adapter_file="/home/mariam/Lynch/tools/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"

  # Step 2: Check quality with FastQC
  fastqc "$fastq_file"

  # Step 3: Trimmomatic
  java -jar "$trimmomatic_jar" SE -phred33 "$fastq_file" "$termo_out" \
    ILLUMINACLIP:"$adapter_file":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

  # Step 4: Index the human genome
  reference_genome="reference2/Homo_sapiens_assembly38.fasta"

  # Step 5: Alignment
  aligned_sam="aligned.sam"
  bwa mem -t 18 "$reference_genome" "$termo_out" > "$aligned_sam"

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
  output_avinput_function="output.avinput.variant_function"
  annovar_dir="tools/annovar"
  "$annovar_dir/convert2annovar.pl" -format vcf4 "$variants_vcf" > "$output_avinput"
  "$annovar_dir/annotate_variation.pl" -geneanno -dbtype refGene "$output_avinput" "$annovar_dir/humandb/" > "$output_avinput_function"

  # Step 10: Grep
  lynch_positions_file="Lynch_positions_in_file.txt"
  grep -E 'MLH1|MSH2|MSH6|PMS2|EPCAM' "$output_avinput_function" > "$lynch_positions_file"

# Clean up by deleting unnecessary files
rm -f "$aligned_sam" "$sorted_bam" "$variants_vcf" "$output_avinput" "$output_avinput_function" \
  "$termo_out" "unpaired_forward.fq.gz" "unpaired_reverse.fq.gz" \
  "output.avinput.exonic_variant_function" "output.avinput.log" "sorted_reads.bam.bai" \
  "${fastq_file/_1.fastq/_1_fastqc.html}" "${fastq_file/_1.fastq/_1_fastqc.zip}"


elif [ $# -eq 2 ]; then
  # Paired-end data
  fastq_file1="$1"
  fastq_file2="$2"
  trimmomatic_jar="/home/mariam/Lynch/tools/trimmomatic.jar"
  termo_out1="output_1.fq.gz"
  termo_out2="output_2.fq.gz"
  unpaired_forward="unpaired_forward.fq.gz"
  unpaired_reverse="unpaired_reverse.fq.gz"
  adapter_file="/home/mariam/Lynch/tools/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"

  # Step 2: Check quality with FastQC
  fastqc "$fastq_file1"
  fastqc "$fastq_file2"

  # Step 3: Trimmomatic
  java -jar "$trimmomatic_jar" PE -phred33 "$fastq_file1" "$fastq_file2" \
    "$termo_out1" "$unpaired_forward" "$termo_out2" "$unpaired_reverse" \
    ILLUMINACLIP:"$adapter_file":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

  # Step 4: Index the human genome
  reference_genome="reference2/Homo_sapiens_assembly38.fasta"

  # Step 5: Alignment
  aligned_sam="aligned.sam"
  bwa mem -t 18 "$reference_genome" "$termo_out1" "$termo_out2" > "$aligned_sam"

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
  output_avinput_function="output.avinput.variant_function"
  annovar_dir="tools/annovar"
  "$annovar_dir/convert2annovar.pl" -format vcf4 "$variants_vcf" > "$output_avinput"
  "$annovar_dir/annotate_variation.pl" -geneanno -dbtype refGene "$output_avinput" "$annovar_dir/humandb/" > "$output_avinput_function"

  # Step 10: Grep
  lynch_positions_file="Lynch_positions_in_file.txt"
  grep -E 'MLH1|MSH2|MSH6|PMS2|EPCAM' "$output_avinput_function" > "$lynch_positions_file"

  # Clean up by deleting unnecessary files
  rm -f "$aligned_sam" "$sorted_bam" "$variants_vcf" "$output_avinput" "$output_avinput_function" \
    "$termo_out1" "$termo_out2" "$unpaired_forward" "$unpaired_reverse" \
    "output.avinput.exonic_variant_function" "output.avinput.log" "sorted_reads.bam.bai" \
    "${fastq_file1/_1.fastq/_1_fastqc.html}" "${fastq_file1/_1.fastq/_1_fastqc.zip}" \
    "${fastq_file2/_2.fastq/_2_fastqc.html}" "${fastq_file2/_2.fastq/_2_fastqc.zip}"
else
  echo "Invalid number of input files. Please provide either one file for single-end data or two files for paired-end data."
  exit 1
fi
