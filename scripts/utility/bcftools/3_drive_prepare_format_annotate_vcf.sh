#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 31 May 2022
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
path_file_vcf_source=${1} # full path to source genotype file in VCF format
path_file_vcf_product=${2} # full path to product genotype file in VCF format
path_file_genome_assembly_sequence=${3} # full path to genome assembly sequence file in FASTA format without compression that matches source and product genome assembly
path_file_chromosome_translations=${4} # full path to file for chromosome name translations in format for BCFTools "annotate --rename-chrs"
path_file_dbsnp_reference=${5} # full path to file for dbSNP reference in VCF format
threads=${6} # count of processing threads to use
path_script_decompose_align_unique_sort=${7} # full path to script for preparation of genotype VCF files
path_script_translate_chromosomes=${8} # full path to script for preparation of genotype VCF files
path_script_introduce_dbsnp_rsid=${9} # full path to script for preparation of genotype VCF files
path_bcftools=${10} # full path to installation executable file of BCFTools
report=${11} # whether to print reports

################################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_vcf_product .vcf.gz)"
path_directory_product="$(dirname $path_file_vcf_product)"
path_directory_product_temporary="${path_directory_product}/temporary_csvcf_${name_base_file_product}" # hopefully unique

path_file_temporary_combination="${path_directory_product_temporary}/${name_base_file_product}_combination.bcf"
path_file_temporary_sort="${path_directory_product_temporary}/${name_base_file_product}_sort.bcf"

# Initialize directory.
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary
rm $path_file_vcf_product

################################################################################


# ...
# ...
# ...


# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
