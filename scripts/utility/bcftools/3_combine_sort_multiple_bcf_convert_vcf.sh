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
path_list_line_source_vcf_bcf_files=${1} # full path to file list with new line delimiter of paths to genotype files in VCF or BCF formats for combination
path_bcf_product=${2} # full path to product file in BCF format
threads=${3} # count of processing threads to use
path_bcftools=${4} # full path to installation executable file of BCFTools
report=${5} # whether to print reports

################################################################################
# Remove "chr" prefix from chromosome identifiers in VCF genotype file.
# Write to file in VCF format with BGZIP compression.

# 1. convert all genotype files from VCF to BCF format
# 2. sort samples in all BCF files
# 3.

# Convert genotype files from VCF format to BCF format.
# The BCF format allows for greater performance in BCFTools.

# Sort samples in VCF files.
# bcftools query --list-samples input.vcf | sort > samples.txt
# bcftools view --samples-file samples.txt input.vcf > output.vcf

# Combine genetic records across identical samples from multiple files.
# Samples must be identical and have same sequence.
$path_bcftools \
concat \
--allow-overlaps \
--rm-dups exact \
--file-list $path_list_line_source_vcf_files \
--output $path_bcf_product \
--output-type b9 \
--threads $threads

# Create Tabix index for product file in VCF format.
# BCFTools sometimes requires this Tabix index to read a file.
$path_bcftools \
index \
--force \
--tbi \
--threads $threads \
$path_bcf_product



#
