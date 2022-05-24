#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes...
# ...
################################################################################
################################################################################
################################################################################

################################################################################
# Organize arguments.
path_gwas_source=${1} # full path to file for source GWAS summary statistics
path_gwas_product_container=${2} # full path to directory for product GWAS summary statistics
path_script_format_gwas=${3} # full path to script for format adjustment
report=${4} # whether to print reports

################################################################################
# Paths.
path_gwas_constraint="${path_gwas_product_container}/gwas_constraint.txt"
path_gwas_format="${path_gwas_product_container}/gwas_format.txt"
path_gwas_format_compress="${path_gwas_product_container}/gwas_format.txt.gz"

################################################################################
# Format adaptation.
/usr/bin/bash "$path_script_format_gwas" \
$path_gwas_source \
$path_gwas_constraint \
$path_gwas_format \
$path_gwas_format_compress \
$report

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "drive_format_gwas_prscs.sh"
fi
