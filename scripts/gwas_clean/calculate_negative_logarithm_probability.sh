#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 29 March 2023
# Date, last execution: 29 March 2023
# Review: TCW; 29 March 2023
################################################################################
# Note

# This script calculates the negative, base-ten logarithm (-log10) of the
# probabilities in a table of GWAS summary statistics.

# Source Format
# Effect allele: "A1"
# Delimiter: white space
# Columns: SNP CHR BP A1 A2 A1AF BETA SE P N  Z  INFO NCASE NCONT
#          1   2   3  4  5  6    7    8  9 10 11 12   13    14

# Product Format: Team Standard
# Effect allele: "A1"
# Delimiter: white space
# Columns: SNP CHR BP A1 A2 A1AF BETA SE P_NEG_LOG_10 N  Z  INFO NCASE NCONT
#          1   2   3  4  5  6    7    8  9            10 11 12   13    14

################################################################################



################################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_product=${2} # full path to file for product GWAS summary statistics with GZip compression
report=${3} # whether to print reports

################################################################################
# Organize paths.

# Directories.
name_base_file_product="$(basename $path_file_product .txt.gz)"
path_directory_product="$(dirname $path_file_product)"
path_directory_product_temporary="${path_directory_product}/temporary_${name_base_file_product}" # hopefully unique

# Files.
path_file_temporary_calculation="${path_directory_product_temporary}/${name_base_file_product}_calculation.txt"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

# Initialize files.
rm $path_file_product



###########################################################################
# Execute procedure.



##########
# Note that AWK interprets a single space delimiter (FS=" ") as any white space.
# Calculate negative, base-ten logarithm (-log10) of probabilities.
echo "SNP CHR BP A1 A2 A1AF BETA SE P_NEG_LOG_10 N Z INFO NCASE NCONT" > $path_file_temporary_calculation
zcat $path_file_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  if ((toupper($9) != "NA") && (($9 + 0) > 0))
    # Probability has a non-missing, numeric value greater than zero.
    # Calculate negative, base-ten logarithm (-log10) of probability.
    print $1, $2, $3, $4, $5, $6, $7, $8, (-1*(log($9)/log(10))), $10, $11, $12, $13, $14
  else
    # Probability has a missing value.
    # Skip the row entirely.
    next
}' >> $path_file_temporary_calculation



##########

# Compress file format.
gzip -cvf $path_file_temporary_calculation > $path_file_product

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Calculate negative, base-ten logarithm (-log10)"
  echo "of probabilities in GWAS summary statistics."
  echo "path to source GWAS file: " $path_file_source
  echo "path to product GWAS file: " $path_file_product
  echo "table after calculation:"
  head -10 $path_file_temporary_calculation
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
