#!/bin/bash

#Step 1: Replace 'file.fastq' with the actual file name provided through terminal input
fastq_file="$1"


#Step 2: Check quality with fastqc
fastqc "$fastq_file"


#Step 3: Use existing Trimmomatic installation
trimmomatic_jar="/home/mariam/Lynch/tools/trimmomatic.jar"
termo_out="output.fq.gz"
adapter_file="/home/mariam/Lynch/tools/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"
java -jar "$trimmomatic_jar" SE -phred33 "$fastq_file" "$termo_out" ILLUMINACLIP:"$adapter_file":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


#Step 4: Index the human genome
reference_genome="reference2/Homo_sapiens_assembly38.fasta"


#Step 5: Alignment
aligned_sam="aligned.sam"
bwa mem "$reference_genome" "$termo_out" > "$aligned_sam"


#Step 6: Conversion to BAM format
aligned_bam="aligned.bam"
samtools view -b "$aligned_sam" > "$aligned_bam"


#Step 7: Sorting and indexing the BAM file
sorted_bam="sorted_reads.bam"
samtools sort "$aligned_bam" -o "$sorted_bam"
samtools index "$sorted_bam"


#Step 8: Variant calling
variants_vcf="variants.vcf"
bcftools mpileup -f "$reference_genome" "$sorted_bam" | bcftools call -mv -Ov > "$variants_vcf"


#Step 9: Annotation
output_avinput="output.avinput"
output_avinput_function="output.avinput.variant_function"
annovar_dir="tools/annovar"
"$annovar_dir/convert2annovar.pl" -format vcf4 "$variants_vcf" > "$output_avinput"
"$annovar_dir/annotate_variation.pl" -geneanno -dbtype refGene "$output_avinput" "$annovar_dir/humandb/" > "$output_avinput_function"


#Step 10: Grep
lynch_positions_file="Lynch_positions_in_file.txt"
grep -E 'MLH1|MSH2|MSH6|PMS2|EPCAM' "$output_avinput_function" > "$lynch_positions_file"


#Clean up by deleting unnecessary files
rm -f "$aligned_sam" "$sorted_bam" "$variants_vcf" "$output_avinput" "$output_avinput_function" "$termo_out"  "$adapter_file" "output.avinput.exonic_variant_function" "output.avinput.log" "sorted_reads.bam.bai"