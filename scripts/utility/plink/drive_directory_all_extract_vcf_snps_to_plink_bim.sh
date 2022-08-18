#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 17 May 2022
################################################################################
# Notes:
# This script finds within a parent directory all child genotype files in
# Variant Call Format (VCF).
# For each of these child genotype files in VCF format, the script calls
# another script to extract information about Single Nucleotide polymorphisms
# (SNPs) to a new file in PLINK2 BIM format.

################################################################################
################################################################################
################################################################################

################################################################################
# Organize arguments.
path_directory_genotypes_vcf=${1} # full path to parent directory with source genotype files in VCF format
pattern_genotype_source_vcf_file=${2} # string glob pattern by which to recognize source genotype files in VCF format
path_directory_genotypes_snp_bim=${3} # full path to parent directory for product genotype files in BIM format
threads=${4} # count of processing threads to use
path_plink2=${5} # full path to installation executable file of PLINK2
path_directory_promiscuity_scripts=${6} # full path to directory of general scripts
report=${7} # whether to print reports

###########################################################################
# Organize paths.

# Scripts.
path_script_extract_vcf_to_bim="${path_directory_promiscuity_scripts}/utility/plink/extract_vcf_snps_to_plink_bim.sh"

###########################################################################
# Find source genotype files in VCF format within container directory.

cd $path_directory_genotypes_vcf
for path_file in `find . -maxdepth 1 -mindepth 1 -type f -name "$pattern_genotype_source_vcf_file"`; do
  if [[ -f "$path_file" ]]; then
    # Current content item is a file.
    # Extract directory's base name.
    name_file="$(basename -- $path_file)"
    echo $name_file
    path_genotype_source_vcf="${path_directory_genotypes_vcf}/${name_file}"
    name_file_product_bim="${name_file}" # optional to add prefix or suffix here
    # Convert information from genotype files in VCF format to BIM format.
    /usr/bin/bash "${path_script_extract_vcf_to_bim}" \
    $path_genotype_source_vcf \
    $path_directory_genotypes_snp_bim \
    $name_file_product_bim \
    $threads \
    $path_plink2 \
    $report
  fi
done

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "drive_directory_all_extract_vcf_snps_to_plink_bim.sh"
  echo "----------"
  echo "___ report:"
  #cat $path_heritability_report_suffix
  echo "----------"
fi
