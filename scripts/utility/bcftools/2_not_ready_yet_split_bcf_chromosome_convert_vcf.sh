#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 17 May 2022
################################################################################
# Note

# BCFTools specializes in tasks on genotype files in the Variant Call Format
# (VCF) that describe Single Nucleotide Polymorphisms (SNPs) across samples.

# BCFTools documentation
# https://samtools.github.io/bcftools/bcftools.html

# The "concat" command in BCFTools requires that the multiple VCF files all have
# the same samples.

################################################################################

################################################################################
# Organize arguments.
path_file_list_source_vcf_bcf_files=${1} # full path to file with line-delimiter list of full paths to genotype files in VCF or BCF formats for combination
path_file_product_vcf=${2} # full path to product file in VCF format with BGZip compression
threads=${3} # count of processing threads to use
path_bcftools=${4} # full path to installation executable file of BCFTools
report=${5} # whether to print reports

################################################################################
# Organize paths.

name_file_base_product="$(basename $path_file_product_vcf .vcf.gz)"
path_product_container="$(dirname $path_vcf_product)"
name_base_product_file="$(basename $path_vcf_product .vcf.gz)"
path_vcf_source_no_compression="${path_product_container}/${name_base_source_file}_no_compression.vcf"
path_vcf_product_assembly_no_compression="${path_product_container}/${name_base_product_file}_assembly_no_compression.vcf"
path_vcf_product_sort_no_compression="${path_product_container}/${name_base_product_file}_sort_no_compression.vcf"




################################################################################

# Combine genetic records across identical samples from multiple files.
# Samples must be identical and have same sequence.
$path_bcftools \
concat \
--allow-overlaps \
--rm-dups exact \
--file-list $path_file_list_source_vcf_bcf_files \
--output $path_file_product_vcf \
--output-type u \
--threads $threads

# Sort records for SNPs or other genetic features.
$path_bcftools \
sort \
--max-mem 4G \
--output $path_file_intermediate_sort_records_chromosome \
--output-type u \
--temp-dir $path_directory_product_temporary_chromosome \
$path_file_intermediate_sort_samples_chromosome
echo "----------"
echo "$path_file_intermediate_sort_records_chromosome"
echo "----------"



# Create Tabix index for product file in VCF format.
# BCFTools sometimes requires this Tabix index to read a file.
$path_bcftools \
index \
--force \
--tbi \
--threads $threads \
$path_bcf_product



#
