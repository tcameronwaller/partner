#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes:
# This script extracts information from a genotype file in VCF format and
# represents this information within a new file in BIM format.
# The BIM format does not represent all information from the original genotype
# file in VCF format.

################################################################################
################################################################################
################################################################################

################################################################################
# Organize arguments.
path_genotype_source_vcf_container=${1} # full path to directory with source genotype files in VCF format
pattern_genotype_source_vcf_file=${2} # string glob pattern by which to recognize source genotype files in VCF format
path_genotype_product_bim_container=${3} # full path to directory for product genotype files in BIM format
name_prefix_file_product_bim=${4} # name prefix for product file in BIM format
path_promiscuity_scripts=${5} # full path to directory of general scripts
report=${6} # whether to print reports

###########################################################################
# Organize paths.

# Scripts.
path_script_convert_vcf_to_bim="${path_promiscuity_scripts}/utility/convert_vcf_to_plink_bim.sh"

###########################################################################
# Find source genotype files in VCF format within container directory.

cd $path_genotype_source_vcf_container
for path_file in `find . -maxdepth 1 -mindepth 1 -type f -name "$pattern_genotype_source_vcf_file"`; do
  if [[ -f "$path_file" ]]; then
    # Current content item is a file.
    # Extract directory's base name.
    name_file="$(basename -- $path_file)"
    echo $name_file
    path_genotype_source_vcf="${path_genotype_source_vcf_container}/${name_file}"
    name_file_product_bim="${name_prefix_file_product_bim}${name_file}.bim.gz"
    # Convert information from genotype files in VCF format to BIM format.
    /usr/bin/bash "${path_script_convert_vcf_to_bim}" \
    $path_genotype_source_vcf \
    $path_genotype_product_bim_container \
    $name_file_product_bim \
    $report
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
