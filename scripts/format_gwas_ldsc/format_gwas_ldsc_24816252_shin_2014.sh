#!/bin/bash

###########################################################################
###########################################################################
###########################################################################

# "Organize GWAS summary statistics."
# PubMed: 24816252
# author: Shin
# date: 11 May 2014
# phenotype: metabolite abundance in blood
# Human genome version: GRCh37, hg19 <-- assume since after 2009; but article methods, data servers, and README don't specify
# variant identifier (rsID) version: ???

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
# description: ............................ LDSC column ........... source column .......... position
# variant identifier (RS ID): .............  "SNP" ................  "MarkerName" .......... 1
# alternate allele (effect allele): .......  "A1" .................  "Allele1" ............. 2
# reference allele (non-effect allele): ...  "A2" .................  "Allele2" ............. 3
# sample size: ............................  "N" ..................  "TotalSampleSize" ..... 16
# effect (coefficient or odds ratio): .....  "BETA" or "OR" .......  "Effect" .............. 8
# probability (p-value): ..................  "P" ..................  "P-value" ............. 10

# Remove any previous versions of temporary files.
rm $path_gwas_collection
rm $path_gwas_format

# Organize information from linear GWAS.
echo "SNP A1 A2 N BETA P" > $path_gwas_collection
zcat $path_source_file | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {print $1, toupper($2), toupper($3), $16, $8, $10}' >> $path_gwas_collection
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
