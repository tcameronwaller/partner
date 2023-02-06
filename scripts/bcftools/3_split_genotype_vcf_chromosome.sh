#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 31 May 2022
################################################################################
# Notes:
# This script reads information in BCFTools from a single genotype file in
# Variant Call Format (VCF).
# This script extracts records for genetic features on each chromosome and
# writes these to separate files in VCF format.

################################################################################
################################################################################
################################################################################

################################################################################
# Organize arguments.

path_file_genotype_vcf_source=${1} # full path to source genotype file in VCF format
chromosome_identifier=${2} # chromosomal identifier for query on records for genotype genetic features
path_file_vcf_product_chromosome=${3} # file name suffix for product genotype files in VCF format
threads=${4} # count of processing threads to use
path_bcftools=${5} # full path to installation executable file of BCFTools
report=${6} # whether to print reports

################################################################################
# Organize paths.

# Initialize directory.
rm $path_file_vcf_product_chromosome

################################################################################
# Extract records for genotype genetic features within specific chromosome.

# Extract genotype records for genetic features on current chromosome.
$path_bcftools \
view \
--regions $chromosome_identifier \
--output $path_file_vcf_product_chromosome \
--output-type z9 \
--threads $threads \
$path_file_genotype_vcf_source

# Create Tabix index for file in VCF format with BGZip compression.
# BCFTools is unable to create a Tabix index for files in BCF format.
# Some commands in BCFTools require this Tabix index to read a file.
$path_bcftools \
index \
--force \
--tbi \
--threads $threads \
$path_file_vcf_product_chromosome



#
