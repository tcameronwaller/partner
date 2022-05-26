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
path_file_source_vcf_chromosome=${1}
path_directory_product_genotype_vcf_temporary=${2} # full path to directory for temporary intermediate files
path_file_intermediate_bcf_chromosome=${3} # full path to intermediate file
path_file_intermediate_remove_duplicates_chromosome=${4} # full path to intermediate file
path_file_intermediate_list_samples_chromosome=${5} # full path to intermediate file
path_file_intermediate_sort_samples_chromosome=${6} # full path to intermediate file
path_file_intermediate_sort_records_chromosome=${7} # full path to intermediate file
threads=${8} # count of processing threads to use
path_bcftools=${9} # full path to installation executable file of BCFTools
report=${10} # whether to print reports

################################################################################
# Prepare genotype files for combination.

# Convert genotype files from VCF format to BCF format.
# The BCF format allows for greater performance in BCFTools.
# Read file in VCF format with BGZip compression in BCFTools.
# Write file in BCF format without compression.
$path_bcftools \
view \
--output $path_file_intermediate_bcf_chromosome \
--output-type u \
--threads $threads \
$path_file_source_vcf_chromosome

# Remove duplicate records for SNPs or other genetic features.
$path_bcftools \
norm \
--rm-dup exact \
--output $path_file_intermediate_remove_duplicates_chromosome \
--output-type u \
--threads $threads \
$path_file_intermediate_bcf_chromosome

# Sort samples.
# bcftools query --list-samples input.vcf | sort > samples.txt
# bcftools view --samples-file samples.txt input.vcf > output.vcf
# Extract sample identifiers from genotype file and sort these externally.
$path_bcftools \
query \
--list-samples \
--threads $threads \
$path_file_intermediate_remove_duplicates_chromosome | sort > $path_file_intermediate_list_samples_chromosome
# Sort samples within genotype file.
$path_bcftools \
view \
--samples-file $path_file_intermediate_list_samples_chromosome \
--output $path_file_intermediate_sort_samples_chromosome \
--output-type u \
--threads $threads \
$path_file_intermediate_remove_duplicates_chromosome

# Sort records for SNPs or other genetic features.
$path_bcftools \
sort \
--output $path_file_intermediate_sort_records_chromosome \
--output-type u \
--temp-dir $path_directory_product_genotype_vcf_temporary \
$path_file_intermediate_sort_samples_chromosome

#
