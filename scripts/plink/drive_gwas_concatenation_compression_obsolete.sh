#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes...
# ...
################################################################################
################################################################################
################################################################################

# TODO: TCW; 12 July 2022
# TODO: this driver script might also need updates... for example, pass the file names instead of 'hard-coding' them here?

################################################################################
# Organize arguments.
pattern_gwas_report_file=${1} # string glob pattern by which to recognize PLINK2 GWAS report files
path_gwas_source_parent=${2} # full path to parent directory for source GWAS summary statistics
path_gwas_target_parent=${3} # full path to parent directory for target GWAS summary statistics
path_promiscuity_scripts=${4} # full path to directory of general scripts
chromosome_x=${5} # whether to collect GWAS summary statistics report for Chromosome X
report=${6} # whether to print reports

################################################################################
# Paths and names of files.
name_report_log_file="report.log"
name_report_allele_frequency_file="report.afreq"
path_gwas_concatenation="${path_gwas_target_parent}/gwas_concatenation.txt"
path_gwas_concatenation_compress="${path_gwas_target_parent}/gwas_concatenation.txt.gz"

################################################################################
# Scripts.
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
path_script_gwas_collect_concatenate="${path_scripts_gwas_process}/concatenate_compress_gwas_chromosomes.sh"

################################################################################
# Concatenation across chromosomes

/usr/bin/bash "$path_script_gwas_collect_concatenate" \
$pattern_gwas_report_file \
$name_report_log_file \
$name_report_allele_frequency_file \
$path_gwas_source_parent \
$path_gwas_target_parent \
$path_gwas_concatenation \
$path_gwas_concatenation_compress \
$chromosome_x \
$report

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "drive_gwas_concatenation.sh"
fi
