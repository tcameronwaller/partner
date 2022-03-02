#!/bin/bash

###########################################################################
###########################################################################
###########################################################################

# "Organize GWAS summary statistics."
# PubMed: 34255042
# author: Schmitz
# date: 2021
# phenotype: Detectability of Oestradiol in UK Biobank
# Human genome version: GRCh37, hg19 <-- assume since after 2009; but article methods, data servers, and README don't specify
# variant identifier (rsID) version: ???
# file: "gwas_females.tsv.gz"

# review:
# TCW 02 March 2022

###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize variables.
path_gwas_source=${1} # full path to source file with GWAS summary statistics, with gzip compression
path_gwas_constraint=${2} # full path to temporary file for constraint of GWAS summary statistics to bounds
path_gwas_format=${3} # full path to file for formatted GWAS summary statistics
path_gwas_standard=${4} # full path to file for GWAS summary statistics with standard z-scores
path_gwas_format_compress=${5} # full path to file for formatted GWAS summary statistics after compression
path_script_calculate_z_score=${6} # full path to directory of scripts for z-score standardization
response_standard_scale=${7} # whether to convert reponse (effect, coefficient) to z-score standard scale ("true" or "false")
report=${8} # whether to print reports

###########################################################################
# Execute procedure.

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "path to source file: " $path_gwas_source
  echo "path to target file: " $path_gwas_format_compress
fi

# Here are some useful regular expressions to evaluate values in "awk".
# ( $12 ~ /^[0-9]+$/ ); ( $12 ~ /^[[:alpha:]]+$/ ); ( $12 ~ /^[[:punct:]]+$/ )

# Remove any previous versions of temporary files.
rm $path_gwas_constraint
rm $path_gwas_format
rm $path_gwas_standard
rm $path_gwas_format_compress

##########
# Constrain GWAS probability values within bounds for LDSC.
# Probability values in column "P" from PLINK2 have values of "NA" for missing
# or floating point values between zero and one.
# LDSC uses 64-bit floating point precision (double precision) to represent
# values from +/- 2.23E-308 to +/- 1.80E308 (https://github.com/bulik/ldsc/issues/144).
# Constrain probability values from 1.0E-305 to 1.0.
zcat $path_gwas_source | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_gwas_constraint
zcat $path_gwas_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 16)
    # Skip any rows with incorrect count of column fields.
    next
  else if ( ( $16 != "NA" ) && ( ($16 + 0) < 1.0E-305 ) )
    # Constrain probability value.
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, ( 1.0E-305 )
  else if ( ( $16 != "NA" ) && ( ($16 + 0) > 1.0 ) )
    # Constrain probability value.
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, ( 1.0 )
  else
    print $0
  }' >> $path_gwas_constraint

##########
# Format of GWAS summary statistics for LDSC.
# https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics

# description: ............................ LDSC column ........... source column .......... position
# variant identifier (RS ID): .............  "SNP" ................  "ID" ..................  3
# alternate allele (effect allele): .......  "A1" .................  "A1" ..................  6
# reference allele (non-effect allele): ...  "A2" .................  "AX" ..................  7
# sample size: ............................  "N" ..................  "OBS_CT" ..............  10
# effect (coefficient or odds ratio): .....  "BETA" or "OR" .......  "OR" ..................  11
# probability (p-value): ..................  "P" ..................  "PVAL" ................  16

# Organize information from linear GWAS.
if [[ "$response_standard_scale" == "true" ]]; then
  echo "SNP A1 A2 N Z P" > $path_gwas_format
else
  echo "SNP A1 A2 N OR P" > $path_gwas_format
fi
cat $path_gwas_constraint | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {print $3, toupper($6), toupper($7), $10, $11, $16}' >> $path_gwas_format

##########

# Calculate Z-score standardization of Beta coefficients.
if [[ "$response_standard_scale" == "true" ]]; then
  /usr/bin/bash $path_script_calculate_z_score \
  5 \
  $path_gwas_format \
  $path_gwas_standard \
  $report
else
  cp $path_gwas_format $path_gwas_standard
fi

# Compress file format.
gzip -cvf $path_gwas_standard > $path_gwas_format_compress

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "after bound constraints:"
  head -10 $path_gwas_constraint
  echo "before standardization:"
  head -10 $path_gwas_format
  echo "after standardization:"
  head -10 $path_gwas_standard
  echo "----------"
fi

# Remove temporary files.
rm $path_gwas_constraint
rm $path_gwas_format
rm $path_gwas_standard
