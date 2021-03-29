#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize variables.
study=${1} # unique identifier of the current GWAS study
path_source_file=${2} # full path to source file with GWAS summary statistics
path_gwas_collection=${3} # full path to temporary file for collection of GWAS summary statistics
path_gwas_format=${4} # full path to file for formatted GWAS summary statistics
path_gwas_format_compress=${5} # full path to file for formatted GWAS summary statistics after compression
path_promiscuity_scripts=${6} # complete path to directory of scripts for z-score standardization
report=${7} # whether to print reports

#path_calculate_z_score="$path_promiscuity_scripts/calculate_z_score_column_4_of_5.sh"
path_calculate_z_score="$path_promiscuity_scripts/calculate_z_score_column_5_of_6.sh"

###########################################################################
# Execute procedure.

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "study: " $study
  echo "path to original file: " $path_source_file
  echo "path to new file: " $path_gwas_format
fi

# Format of GWAS summary statistics for LDSC.
# https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics
# description: ............................ LDSC column
# variant identifier: ........................ "SNP"
# alternate allele (effect allele): .......... "A1"
# reference allele (non-effect allele): ...... "A2"
# sample size: ............................... "N"
# effect (coefficient or odds ratio): ........ "BETA" or "OR"
# probability (p-value): ..................... "P"

# Format of GWAS summary statistics for PRS-CS.
# https://github.com/getian107/PRScs
# description: ............................ PRS-CS column
# variant identifier: ........................ "SNP"
# alternate allele (effect allele): .......... "A1"
# reference allele (non-effect allele): ...... "A2"
# effect (coefficient or odds ratio): ........ "BETA" or "OR"
# probability (p-value): ..................... "P"

# Remove any previous versions of temporary files.
rm $path_gwas_collection
rm $path_gwas_format

# Organize information from linear GWAS.
echo "SNP A1 A2 N BETA P" > $path_gwas_collection
zcat $path_source_file | awk 'BEGIN { FS=","; OFS=" " } NR > 1 {print $1, toupper($4), toupper($5), $14, $10, $12}' >> $path_gwas_collection
# Calculate Z-score standardization of Beta coefficients.
/usr/bin/bash $path_calculate_z_score \
5 \
$path_gwas_collection \
$path_gwas_format \
$report

# Compress file format.
# No need in this situation, since each iteration replaces the previous file.
#gzip -cvf $path_temporary_gwas_format > $path_temporary_gwas_format_zip

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "before standardization:"
  head -10 $path_gwas_collection
  echo "after standardization:"
  head -10 $path_gwas_format
  echo "----------"
fi
