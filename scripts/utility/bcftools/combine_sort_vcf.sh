#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 27 May 2022
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
path_file_list_source_vcf_files=${1} # full path to file with line-delimiter list of full paths to genotype files in VCF formats with BGZip compression and Tabix indices
path_file_product_vcf=${2} # full path to product file in VCF format with BGZip compression
threads=${3} # count of processing threads to use
path_bcftools=${4} # full path to installation executable file of BCFTools
report=${5} # whether to print reports

################################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_product_vcf .vcf.gz)"
path_directory_product="$(dirname $path_file_product_vcf)"
path_directory_product_temporary="${path_directory_product}/temporary"

path_file_temporary_combination="${path_directory_product_temporary}/${name_base_file_product}_combination.bcf"
path_file_temporary_sort="${path_directory_product_temporary}/${name_base_file_product}_sort.bcf"

# Initialize directory.
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

################################################################################

# Combine genetic records across identical samples from multiple files.
# Samples must be identical and have same sequence.
# The BCFTools "concat" command requires source files to have Tabix indices.
$path_bcftools \
concat \
--allow-overlaps \
--rm-dups exact \
--file-list $path_file_list_source_vcf_files \
--output $path_file_temporary_combination \
--output-type u \
--threads $threads

# Sort records for SNPs or other genetic features.
$path_bcftools \
sort \
--max-mem 16G \
--output $path_file_temporary_sort \
--output-type u \
--temp-dir $path_directory_product_temporary \
$path_file_temporary_combination

# Convert genotype files from VCF format with BGZip compression to working
# format.
# The BCF format allows for greater performance in BCFTools.
# Read file in VCF format with BGZip compression in BCFTools.
# Write file in BCF format without compression or other format.
$path_bcftools \
view \
--output $path_file_product_vcf \
--output-type z9 \
--threads $threads \
$path_file_temporary_sort

# Create Tabix index for file in VCF format with BGZip compression.
# BCFTools is unable to create a Tabix index for files in BCF format.
# Some commands in BCFTools require this Tabix index to read a file.
$path_bcftools \
index \
--force \
--tbi \
--threads $threads \
$path_file_product_vcf

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary

#
