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
path_gwas_concatenation_compress=${1} # full path to file for source GWAS summary statistics
path_gwas_target_parent=${2} # full path to parent directory for target GWAS summary statistics
path_promiscuity_scripts=${3} # full path to directory of general scripts
path_script_format_gwas=${4} # full path to script for format adjustment
report=${5} # whether to print reports

################################################################################
# Paths.
path_gwas_constraint="${path_gwas_target_parent}/gwas_constraint.txt"
path_gwas_format="${path_gwas_target_parent}/gwas_format.txt"
path_gwas_format_compress="${path_gwas_target_parent}/gwas_format.txt.gz"

################################################################################
# Scripts.
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"

################################################################################
# Format adaptation.
path_gwas_source=$path_gwas_concatenation_compress
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
