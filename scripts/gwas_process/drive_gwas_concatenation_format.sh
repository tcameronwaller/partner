#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes...
# ...
################################################################################
################################################################################
################################################################################

# TODO: TCW 26 August 2021
# TODO: I think this function ought to accept an argument for the format script to use
# TODO: need to switch between linear formatting and logistic formatting

################################################################################
# Organize arguments.
pattern_gwas_report_file=${1} # string glob pattern by which to recognize PLINK2 GWAS report files
path_gwas_source_parent=${2} # full path to parent directory for source GWAS summary statistics
path_gwas_target_parent=${3} # full path to parent directory for target GWAS summary statistics
response_standard_scale=${4} # whether to convert response (coefficient) to z-score standard scale
path_promiscuity_scripts=${5} # full path to directory of general scripts
report=${6} # whether to print reports

################################################################################
# Paths.
path_gwas_concatenation="${path_gwas_target_parent}/gwas_concatenation.txt"
path_gwas_concatenation_compress="${path_gwas_target_parent}/gwas_concatenation.txt.gz"
path_gwas_constraint="${path_gwas_target_parent}/gwas_constraint.txt"
path_gwas_format="${path_gwas_target_parent}/gwas_format.txt"
path_gwas_standard="${path_gwas_target_parent}/gwas_standard.txt"
path_gwas_format_compress="${path_gwas_target_parent}/gwas_format.txt.gz"

################################################################################
# Scripts.
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
path_script_gwas_collect_concatenate="${path_scripts_gwas_process}/concatenate_compress_gwas_chromosomes.sh"
path_script_gwas_format="${path_scripts_gwas_process}/format_gwas_ldsc/format_gwas_ldsc_plink_linear.sh"
path_script_calculate_z_score="${path_scripts_gwas_process}/calculate_z_score_column_5_of_6.sh"

################################################################################
# Concatenation across chromosomes

/usr/bin/bash "$path_script_gwas_collect_concatenate" \
$pattern_gwas_report_file \
$path_gwas_source_parent \
$path_gwas_concatenation \
$path_gwas_concatenation_compress \
$report

################################################################################
# Format adaptation
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
  echo "drive_gwas_concatenation_format.sh"
fi
