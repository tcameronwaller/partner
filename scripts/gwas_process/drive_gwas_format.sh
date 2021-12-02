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
path_script_gwas_format=${4} # full path to script for format adjustment
response_standard_scale=${5} # whether to convert reponse (effect, coefficient) to z-score standard scale ("true" or "false")
report=${6} # whether to print reports

################################################################################
# Paths.
path_gwas_constraint="${path_gwas_target_parent}/gwas_constraint.txt"
path_gwas_format="${path_gwas_target_parent}/gwas_format.txt"
path_gwas_standard="${path_gwas_target_parent}/gwas_standard.txt"
path_gwas_format_compress="${path_gwas_target_parent}/gwas_format.txt.gz"

################################################################################
# Scripts.
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
path_script_calculate_z_score="${path_scripts_gwas_process}/calculate_z_score_column_5_of_6.sh"

################################################################################
# Format adaptation.
path_gwas_source=$path_gwas_concatenation_compress
/usr/bin/bash "$path_script_gwas_format" \
$path_gwas_source \
$path_gwas_constraint \
$path_gwas_format \
$path_gwas_standard \
$path_gwas_format_compress \
$path_script_calculate_z_score \
$response_standard_scale \
$report

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "drive_gwas_format.sh"
fi
