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
path_gwas_source_parent=${1} # full path to parent directory for source GWAS summary statistics
path_gwas_target_parent=${2} # full path to parent directory for target GWAS summary statistics
path_promiscuity_scripts=${3} # full path to directory of general scripts
report=${4} # whether to print reports

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
path_script_gwas_collect_concatenate="${path_scripts_gwas_process}/collect_concatenate_gwas_chromosomes.sh"
path_script_gwas_format="${path_scripts_gwas_process}/format_gwas_ldsc/format_gwas_ldsc_plink_linear.sh"
path_script_calculate_z_score="${path_scripts_gwas_process}/calculate_z_score_column_5_of_6.sh"

################################################################################
# Concatenation across chromosomes
pattern_source_file="report.*.glm.linear" # do not expand with full path yet
chromosome_start=1
chromosome_end=22
/usr/bin/bash "$path_script_gwas_collect_concatenate" \
$pattern_source_file \
$path_gwas_source_parent \
$chromosome_start \
$chromosome_end \
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
$report

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "drive_gwas_concatenation_format.sh"
fi
