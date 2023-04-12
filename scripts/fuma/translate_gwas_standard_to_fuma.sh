#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 12 April 2023
# Date, last execution: 12 April 2023
# Review: TCW; 12 April 2023
################################################################################
# Note

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

################################################################################



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
path_file_temporary_format_1="${path_directory_product_temporary}/${name_base_file_product}_format_1.txt"
path_file_temporary_format_2="${path_directory_product_temporary}/${name_base_file_product}_format_2.txt"

# Initialize directories.
mkdir -p $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

# Initialize files.
rm $path_file_product

###########################################################################
# Execute procedure.

##########
# Translate format of GWAS summary statistics.
# Note that AWK interprets a single space delimiter (FS=" ") as any white space.

# 1. Remove columns irrelevant to FUMA.
echo "SNP CHR BP A1 A2 A1AF BETA SE P N" > $path_file_temporary_format_1
# For conciseness, only support the conditions that are relevant.
zcat $path_file_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print $1, $2, $3, toupper($4), toupper($5), $6, $7, $8, $9, $10
}' >> $path_file_temporary_format_1

# 2. Remove rows with missing information about effects.
echo "SNP CHR BP A1 A2 A1AF BETA SE P N" > $path_file_temporary_format_2
cat $path_file_temporary_format_1 | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  if ((toupper($7) != "NA") && (($7 + 0) > 0) && (toupper($8) != "NA") && (($8 + 0) > 0) && (toupper($9) != "NA") && (($9 + 0) > 0))
    # Columns "BETA", "SE", and "P" have non-missing values.
    print $1, $2, $3, toupper($4), toupper($5), $6, $7, $8, $9, $10
  else
    # Column "BETA", "SE", or "P" has missing value.
    next
}' >> $path_file_temporary_format_2

# Compress file format.
gzip -cvf $path_file_temporary_format_2 > $path_file_product

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Translate format of GWAS summary statistics."
  echo "path to source GWAS file: " $path_file_source
  echo "path to product GWAS file: " $path_file_product
  echo "----------"
  echo "table after format:"
  #head -10 $path_file_temporary_format_2
  zcat $path_file_product | head -10
  echo "----------"
  echo "rows in table before format: "
  zcat $path_file_source | wc -l
  echo "rows in table after format: "
  zcat $path_file_product | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
