#!/bin/bash

###########################################################################
###########################################################################
###########################################################################

# "Organize GWAS summary statistics."
# PubMed: 29906448
# author: Ruderfer
# date: 2018
# phenotype: bipolar disorder and schizophrenia
# Human genome version: GRCh37, hg19 <-- assume since after 2009; but article methods, data servers, and README don't specify
# variant identifier (rsID) version: ???
# file: "BDSCZvsCONT.sumstats.gz", "SCZvsBD.sumstats.gz", ...
# note: study "BDSCZvsCONT.sumstats.gz" has 53,555 cases and 54,065 controls, for total 107,620 samples
# note: study "SCZvsBD.sumstats.gz" has total 38,855 samples
# note: on 6 April 2021, TCW confirmed that table columns have the same order in both files

###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize arguments.
study=${1} # unique identifier of the current GWAS study
path_source_file=${2} # full path to source file with GWAS summary statistics
path_gwas_collection=${3} # full path to temporary file for collection of GWAS summary statistics
path_gwas_format=${4} # full path to file for formatted GWAS summary statistics
path_gwas_format_compress=${5} # full path to file for formatted GWAS summary statistics after compression
path_promiscuity_scripts=${6} # complete path to directory of scripts for z-score standardization
report=${7} # whether to print reports

################################################################################
# Organize variables.

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
# variant identifier (RS ID): .............  "SNP" ................  "SNP" .................  2
# alternate allele (effect allele): .......  "A1" .................  "A1" ..................  4
# reference allele (non-effect allele): ...  "A2" .................  "A2" ..................  5
# sample size: ............................  "N" ..................  None ..................  (33,426 + 32,541)
# effect (coefficient or odds ratio): .....  "BETA" or "OR" .......  "OR" ..................  9
# probability (p-value): ..................  "P" ..................  "P" ...................  11

# Remove any previous versions of temporary files.
rm $path_gwas_collection
rm $path_gwas_format

# Organize information from linear GWAS.
echo "SNP A1 A2 N BETA P" > $path_gwas_collection
zcat $path_source_file | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {print $2, toupper($4), toupper($5), (33426 + 32541), log($9), $11}' >> $path_gwas_collection
# Calculate Z-score standardization of Beta coefficients.
/usr/bin/bash $path_calculate_z_score \
5 \
$path_gwas_collection \
$path_gwas_format \
$report

# Compress file format.
gzip -cvf $path_gwas_format > $path_gwas_format_compress

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "before standardization:"
  head -10 $path_gwas_collection
  echo "after standardization:"
  head -10 $path_gwas_format
  echo "----------"
fi
