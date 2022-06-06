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
path_file_vcf_source=${1} # full path to source genotype file in VCF format
path_file_vcf_product=${2} # full path to product genotype file in VCF format
path_file_chromosome_translations=${3} # full path to file for chromosome name translations in format for BCFTools "annotate --rename-chrs"
threads=${4} # count of processing threads to use
path_bcftools=${5} # full path to installation executable file of BCFTools
report=${6} # whether to print reports

################################################################################
# Remove "chr" prefix from chromosome identifiers in VCF genotype file.
# Write to file in VCF format with BGZIP compression.

$path_bcftools \
annotate \
--rename-chrs $path_file_chromosome_translations \
--output $path_file_vcf_product \
--output-type z9 \
--threads $threads \
$path_file_vcf_source

# Create Tabix index for file in VCF format with BGZip compression.
# BCFTools is unable to create a Tabix index for files in BCF format.
# Some commands in BCFTools require this Tabix index to read a file.
$path_bcftools \
index \
--force \
--tbi \
--threads $threads \
$path_file_vcf_product



#
