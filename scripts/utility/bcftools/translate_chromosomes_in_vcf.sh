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

# Examples of calls to the "annotate" function of BCFTools.
# http://samtools.github.io/bcftools/howtos/annotate.html

# BGZIP documentation
# http://www.htslib.org/doc/bgzip.html

# "Cheat sheet" on BCFTools
# https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b

# Example of simple script call format for using BCFTools Annotate to introduce rsIDs
# https://gist.github.com/obenshaindw/99440ac43de07548e453
# "/usr/bin/htslib/bcftools/bcftools annotate -a /reference/dbsnp_137.b37.vcf.gz -c ID vcf_to_add_id_to.vcf"

################################################################################

################################################################################
# Organize arguments.
path_chromosome_translations=${1} # full path to file for chromosome name translations in format for BCFTools "annotate --rename-chrs"
path_vcf_source=${2} # full path to source file in VCF format
path_vcf_product=${3} # full path to product file in VCF format
path_bcftools=${4} # full path to installation of BCFTools
report=${5} # whether to print reports

################################################################################
# Remove "chr" prefix from chromosome identifiers in VCF genotype file.
# Introduce dbSNP rsID annotations VCF genotype file.
# Write to file in VCF format with BGZIP compression.

# Only remove "chr" prefix from chromosome identifiers. Tested successfully.
$path_bcftools \
annotate \
--rename-chrs $path_chromosome_translations \
--output $path_vcf_product \
--output-type z9 \
--threads 4 \
$path_vcf_source

# Create Tabix index for product file in VCF format.
# BCFTools sometimes requires this Tabix index to read a file.
$path_bcftools \
index \
--force \
--tbi \
--threads 4 \
$path_vcf_product

# Both remove "chr" prefix from chromosome identifiers and introduce dbSNP rsID
# annotations. Only the removal of "chr" prefix works (TCW; 17 May 2022).
#$path_bcftools \
#annotate \
#--rename-chrs $path_chromosome_translations \
#--annotations $path_dbsnp_reference \
#--columns ID \
#--output $path_vcf_product \
#--output-type b9 \
#--threads 4 \
#$path_vcf_source
