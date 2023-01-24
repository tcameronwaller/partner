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
# effect allele: _____
# delimiter: white space
# columns: CHR POS ID REF ALT A1 AX A1_FREQ TEST OBS_CT OR LOG(OR)_SE L95 U95 Z_STAT P

# Format Translation
# columns: $3, $1, $2, toupper($6), toupper($7), $8, log($11), $12, $16, $10, $15, (1), "NA", "NA"

# Product Format (Team Standard)
# effect allele: "A1"
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

# review: TCW; 23 December 2022
# check: Standard Format Columns [TCW; 22 December 2022]
# check: Study citation, PubMed, and Host website [TCW; 22 December 2022]
# check: Study field delimiters [TCW; 23 December 2022]
# check: Study source columns [TCW; 23 December 2022]
# check: Translation column order [TCW; 23 December 2022]


###########################################################################
###########################################################################
###########################################################################



################################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_product=${2} # full path to file for product GWAS summary statistics in format with GZip compression
report=${3} # whether to print reports

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

##########
# Translate format of GWAS summary statistics.
# Note that AWK interprets a single space delimiter (FS=" ") as any white space.

echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_format
zcat $path_file_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print $3, $1, $2, toupper($6), toupper($7), $8, log($11), $12, $16, $10, $15, (1), "NA", "NA"
}' >> $path_file_temporary_format

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
