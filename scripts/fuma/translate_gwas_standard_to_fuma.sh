#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Notes

# This script translates the format of linear GWAS summary statistics from
# PLINK2.

# Source Format (Team Standard)
# effect allele: "A1"
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

# Product Format (Concise for FUMA)
# effect allele: "A1"
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N

# Review: TCW; __ April 2023

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
echo "SNP CHR BP A1 A2 A1AF BETA SE P N" > $path_file_temporary_format
# For conciseness, only support the conditions that are relevant.
zcat $path_file_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print $1, $2, $3, toupper($4), toupper($5), $6, $7, $8, $9, $10
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
