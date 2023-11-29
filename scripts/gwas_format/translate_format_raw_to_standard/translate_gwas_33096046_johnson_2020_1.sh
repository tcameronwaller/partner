#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 24 November 2023
# Date, last execution: 24 November 2023
# Review: TCW; 24 November 2023
################################################################################
# Note

# Note: TCW; 22 November 2023

# z-score = (beta - (mean of all betas)) / (standard deviation of all betas)
# beta = (z-score)*(standard deviation of all betas) + (mean of all betas)

# Here is a reference for the calculation of beta effect parameters from the
# z-score in GWAS summary statistics, but I do not trust this source.
# https://www.biostars.org/p/319584/

# If we assume that across all betas the mean is approximately zero and the
# standard deviation is approximately one, then the z-score is an adequate
# approximation of the beta effect.

# Determine maximal counts of total observations (samples), cases, and controls.
# $ zcat <file name> | awk 'BEGIN{FS=" "; OFS=" "; a=0} NR>1{if ((toupper($8) != "NA") && (($8+0)>(a+0))) a=$8} END{print a}'
# $ zcat <file name> | awk 'BEGIN{FS=" "; OFS=" "; a=0} NR>1{if ((toupper($9) != "NA") && (($9+0)>(a+0))) a=$9} END{print a}'
# $ zcat <file name> | awk 'BEGIN{FS=" "; OFS=" "; a=0} NR>1{if ((toupper($10) != "NA") && (($10+0)>(a+0))) a=$10} END{print a}'

################################################################################



################################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_product=${2} # full path to file for product GWAS summary statistics in format with GZip compression
fill_observations=${3} # logical binary indicator of whether to fill count of observations across all variants
observations=${4} # count of observations
fill_case_control=${5} # logical binary indicator of whether to fill counts of cases and controls across all variants
cases=${6} # count of cases
controls=${7} # count of controls
report=${8} # whether to print reports

################################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_product .txt.gz)"
path_directory_product="$(dirname $path_file_product)"
path_directory_product_temporary="${path_directory_product}/temporary_format_${name_base_file_product}" # hopefully unique
path_file_temporary_format_1="${path_directory_product_temporary}/${name_base_file_product}_format_1.txt"
path_file_temporary_format_2="${path_directory_product_temporary}/${name_base_file_product}_format_2.txt"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

# Remove any previous version of the product file.
rm $path_file_product

################################################################################
# Execute procedure.

# Calculate an estimate of the standard error of the effect from the Z-score of
# the effect and the probability (p-value) of the effect.
# Reference:
# PubMed:21824904
# Title: "How to obtain the confidence interval from a P value"

##########
# Translate format of GWAS summary statistics.
# Note that AWK interprets a single space delimiter (FS=" ") as any white space.
# For conciseness, only support the conditions that are relevant.

# 1. Translate column format while simplifying identifier.
echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_format_1
if [ "$fill_observations" != "1" ] && [ "$fill_case_control" != "1" ]; then
  zcat $path_file_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
    print $2, $1, $3, toupper($4), toupper($5), "NA", $6, "NA", $7, $8, "NA", (1.0), $9, $10
  }' >> $path_file_temporary_format_1
fi

# 2. Estimate standard error.
# Notice that the effect or z-score can have values less than zero.
echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_format_2
cat $path_file_temporary_format_1 | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  if ((toupper($7) != "NA") && (toupper($9) != "NA") && (($9 + 0) > 0))
    # Calculate estimate of standard error of effect from the effect and its
    # probability (p-value).
    # Use the square root of the square to take absolute value.
    print $1, $2, $3, toupper($4), toupper($5), $6, $7, sqrt((($7) / (-0.862 + sqrt(0.743 - (2.404 * log($9)))))^2), $9, $10, $11, $12, $13, $14
  else
    # Print missing value for standard error.
    print $1, $2, $3, toupper($4), toupper($5), $6, $7, "NA", $9, $10, $11, $12, $13, $14
}' >> $path_file_temporary_format_2

# Compress file format.
gzip -cvf $path_file_temporary_format_2 > $path_file_product

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "Translate format of GWAS summary statistics."
  echo "----------"
  echo "path to source GWAS file: " $path_file_source
  echo "path to product GWAS file: " $path_file_product
  echo "----------"
  echo "table before format translation:"
  zcat $path_file_source | head -5
  echo "----------"
  echo "table after format translation:"
  zcat $path_file_product | head -5
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
