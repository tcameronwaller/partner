#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes...

################################################################################
################################################################################
################################################################################

################################################################################
# Organize arguments.
path_genotype_source_vcf_container=${1} # full path to directory with source genotype files in VCF format
pattern_genotype_source_vcf_file=${2} # string glob pattern by which to recognize source genotype files in VCF format
path_genotype_product_bim_container=${3} # full path to directory for product genotype files in BIM format
name_prefix_file_product_bim=${4} # name prefix for product file in BIM format
report=${5} # whether to print reports

###########################################################################
# Find source genotype files in VCF format within container directory.

cd $path_genotype_source_vcf_container
for path_file in `find . -maxdepth 1 -mindepth 1 -type f -name "$pattern_genotype_source_vcf_file"`; do
  if [[ -f "$path_file" ]]; then
    # Current content item is a file.
    # Extract directory's base name.
    name_file="$(basename -- $path_file)"
    echo $name_file

    # pass the same file name to PLINK2... will it append new suffix???
  fi
done

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "drive_convert_directory_all_vcf_to_plink_bim.sh"
  echo "----------"
  echo "___ report:"
  #cat $path_heritability_report_suffix
  echo "----------"
fi
