#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 17 May 2022
################################################################################
# Notes:
# This script extracts information about single nucleotide polymorphisms (SNPs)
# from a genotype file in Variant Call Format (VCF) and represents this
# information within a new file in PLINK2 BIM format.
# The BIM format does not represent all information from the original genotype
# file in VCF format.
################################################################################

################################################################################
# Organize arguments.
path_genotype_source_vcf=${1} # full path to source genotype file in VCF format
path_genotype_product_bim_container=${2} # full path to directory for product genotype files in BIM format
name_file_product_bim=${3} # name for product file in BIM format
threads=${4} # count of processing threads to use
path_plink2=${5} # full path to installation executable file of PLINK2
report=${6} # whether to print reports

###########################################################################
# Organize paths.

################################################################################
# Read information from VCF format file to PLINK2 and write file in BIM format.
# https://www.cog-genomics.org/plink/2.0/data#make_pgen
# Include the 'zs' modifier for '--make-just-bim' to apply Z-standard compression.

cd $path_genotype_product_bim_container
$path_plink2 \
--memory 90000 \
--threads $threads \
--vcf $path_genotype_source_vcf \
--make-just-bim \
--out "${name_file_product_bim}"

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "extract_vcf_snps_to_plink_bim.sh"
  echo "----------"
  echo "___ report:"
  #cat $path_heritability_report_suffix
  echo "----------"
fi
