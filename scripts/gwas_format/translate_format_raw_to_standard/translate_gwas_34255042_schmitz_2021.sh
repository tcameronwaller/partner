#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Notes

# This script translates the format of GWAS summary statistics from
# Schmitz et al, Journal of Clinical Endocrinology and Metabolism, 2021
# (PubMed:34255042).
# Host: https://zenodo.org/record/4926701

# Source Format
# Human Genome Assembly: GRCh37 (UK Biobank)
# Effect allele: "A1"
# Delimiter: white space
# Columns: CHR POS ID REF ALT A1 AX A1_FREQ TEST OBS_CT OR LOG(OR)_SE L95 U95 Z_STAT P
#          1   2   3  4   5   6  7  8       9    10     11 12         13  14  15     16 (TCW; 15 February 2023)
# Note: It seems that the authors already converted the standard error of the odds ratio to the standard error for a beta coefficient.

# Format Translation
# columns: $3, $1, $2, toupper($6), toupper($7), $8, log($11), $12, $16, $10, $15, (1), "NA", "NA"

# Product Format (Team Standard)
# effect allele: "A1"
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

# Review: TCW; 17 February 2023

###########################################################################
###########################################################################
###########################################################################



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
path_file_temporary_format="${path_directory_product_temporary}/${name_base_file_product}_format.txt"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

# Remove any previous version of the product file.
rm $path_file_product

###########################################################################
# Execute procedure.

# Note: TCW; 17 February 2023
# The logarithm of a negative number or zero is undefined.

##########
# Translate format of GWAS summary statistics.
# Note that AWK interprets a single space delimiter (FS=" ") as any white space.
echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_format
# For conciseness, only support the conditions that are relevant.
if [ "$fill_observations" != "1" ] && [ "$fill_case_control" == "1" ]; then
  zcat $path_file_source | awk -v cases=$cases -v controls=$controls 'BEGIN {FS = " "; OFS = " "} NR > 1 {
    if ((toupper($11) != "NA") && (($11 + 0) > 0))
      print $3, $1, $2, toupper($6), toupper($7), $8, log($11), $12, $16, $10, $15, (1.0), (cases), (controls)
    else
      print $3, $1, $2, toupper($6), toupper($7), $8, "NA", $12, $16, $10, $15, (1.0), (cases), (controls)
  }' >> $path_file_temporary_format
fi

# Compress file format.
gzip -cvf $path_file_temporary_format > $path_file_product

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Translate format of GWAS summary statistics."
  echo "path to source GWAS file: " $path_file_source
  echo "path to product GWAS file: " $path_file_product
  echo "table after format:"
  head -10 $path_file_temporary_format
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
