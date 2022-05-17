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

################################################################################

################################################################################
# Organize arguments.
path_chromosome_translations=${1} # full path to file for chromosome name translations
path_vcf_source=${2} # full path to source file in VCF format
path_vcf_product=${3} # full path to product file in VCF format
path_bcftools=${4} # full path to installation of BCFTools
report=${5} # whether to print reports

################################################################################
# Introduce annotation information from dbSNP reference to file in VCF format.

$path_bcftools \
annotate \
--rename_chrs $path_chromosome_translations \
--output $path_vcf_product \
--output-type b9 \
--threads 4 \
$path_vcf_source
