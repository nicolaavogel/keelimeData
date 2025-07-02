#!/bin/bash
module load samtools
#module load bcftools


# Check if the correct number of arguments is given
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_fastq> <output_bam> <reference_fasta> <output_prefix>"
    exit 1
fi

# Assign arguments to variables
input_fastq=$1
output_bam=$2
reference_fasta=$3
output_prefix=$4

# Check if the reference file exists
if [ ! -f "$reference_fasta" ]; then
    echo "Error: Reference FASTA file does not exist."
    exit 1
fi

# Create the reference index
samtools faidx $reference_fasta

# Variables for generated file names
sorted_bam="${output_bam%.*}.sort.bam"
bam_5p_prof="${output_prefix}_5p.prof"
bam_3p_prof="${output_prefix}_3p.prof"
log_file="${output_prefix}.log"
bcf_file="${output_prefix}.calls.bcf"
filtered_vcf="${output_prefix}_filtered.vcf"
consensus_fa="${output_prefix}.cons.fa"
final_consensus_fa="${output_prefix}_final.fa"

# Command 1: Map reads to the reference (assuming mapBWA.sh does that)
./mapBWA.sh $input_fastq $sorted_bam $reference_fasta

# Command 2: Generate profiles from BAM
/projects/wintherpedersen/people/bfj994/execut/bam2prof -double -5p $bam_5p_prof -3p $bam_3p_prof $sorted_bam

# Command 3: Index the sorted BAM
samtools index $sorted_bam

# Command 4: Run endoCaller to generate BCF file
/projects/wintherpedersen/people/bfj994/execut/endoCaller -seq ${output_prefix}Cons.fa -log $log_file -name ${output_prefix}one -deam5p $bam_5p_prof -deam3p $bam_3p_prof $reference_fasta $sorted_bam
if [ $? -ne 0 ]; then
    echo "Error: endoCaller failed to run."
    exit 1
fi

##  Command 5: Annotate with depth information (for VCF creation)
# bcftools mpileup -Ou -f $reference_fasta $sorted_bam | bcftools call -mv -Ob -o $bcf_file
# if [ $? -ne 0 ]; then
#     echo "Error: bcftools mpileup failed."
#     exit 1
# fi
#
# # Index the BCF file
# bcftools index $bcf_file
# if [ $? -ne 0 ]; then
#     echo "Error: Failed to index BCF file."
#     exit 1
# fi
#
# # Annotate the BCF file
# bcftools annotate -x FORMAT,^INFO/DP -Ov -o depth_annotated.vcf $bcf_file
# if [ $? -ne 0 ]; then
#     echo "Error: bcftools annotate failed."
#     exit 1
# fi
#
# # Filter and create a mask file with positions where DP < 32
# awk '$1 !~ /^#/ && $8 ~ /DP=[0-9]+/ {split($8, arr, "DP="); if (arr[2] < 3) print $0}' depth_annotated.vcf > filtered_positions.vcf
#
# # Create a proper VCF header for the filtered positions file
# grep "^##" depth_annotated.vcf > mask_header.vcf
# echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> mask_header.vcf
# cat mask_header.vcf filtered_positions.vcf > $filtered_vcf
#
# # Command 6: Generate consensus and mask low coverage positions
# bcftools consensus -f $reference_fasta -m $filtered_vcf $bcf_file > $consensus_fa
# if [ $? -ne 0 ]; then
#     echo "Error: bcftools consensus failed."
#     exit 1
# fi
#
# # Command 7: Replace non-ACGT characters with 'N' and modify the header
# # Extract the first part of output_prefix (e.g., SC9 from SC9SqComp)
# header_prefix=$(echo $output_prefix | grep -o "^SC[0-9]\{1,2\}")
#
# # Replace non-ACGT characters with 'N' and modify the header
# sed -e "s/^>.*$/>${header_prefix}/" -e '/^>/! s/[^ACGTacgt]/N/g' $consensus_fa > $final_consensus_fa
#
#
# echo "Workflow completed."
