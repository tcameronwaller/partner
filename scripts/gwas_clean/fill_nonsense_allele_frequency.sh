#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 15 November 2023
# Date, last execution: 27 November 2023
# Date, review: 27 November 2023
################################################################################
# Note

# Some analyses on GWAS summary statistics, such as LDSC, do not require or use
# the frequencies of effect alleles.
# The purpose of this script is to prevent SNPs from failing filters on the
# basis of allele frequency.

################################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_product=${2} # full path to file for product GWAS summary statistics with GZip compression
report=${3} # whether to print reports

################################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_product .txt.gz)"
path_directory_product="$(dirname $path_file_product)"
path_directory_product_temporary="${path_directory_product}/temporary_${name_base_file_product}" # hopefully unique
path_file_temporary="${path_directory_product_temporary}/${name_base_file_product}_temporary.txt"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

# Remove any previous version of the product file.
rm $path_file_product

###########################################################################
# Execute procedure.

#else if ( ( ($6 + 0) > 0 ) && ( ($6 + 0) < 1.0E-307 ) )
#  # Constrain allele frequency.
#  print $1, $2, $3, $4, $5, (1.0E-307), $7, $8, $9, $10, $11, $12, $13, $14

#echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_check
zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_file_temporary
zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( ($6 == "") || (toupper($6) == "NA") || (toupper($6) == "NAN") || ( ($6 + 0) < 0 ) || ( ($6 + 0) > 1.0 ) )
    # Current row has missing or nonsense allele frequency.
    # Fill with nonsense allele frequency that will not trip filters.
    print $1, $2, $3, $4, $5, (0.5), $7, $8, $9, $10, $11, $12, $13, $14
  else
    # Row passes all checks, filters, and constraints.
    # Print the row entirely as it is.
    print $0
  }' >> $path_file_temporary

# Compress file format.
gzip -cvf $path_file_temporary > $path_file_product

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "Filter and constrain values in GWAS summary statistics."
  echo "----------"
  echo "path to source GWAS file: " $path_file_source
  echo "path to product GWAS file: " $path_file_product
  echo "----------"
  echo "table before transformation:"
  zcat $path_file_source | head -5
  echo "- - Count of lines in source GWAS summary statistics:"
  zcat $path_file_source | wc -l
  echo "----------"
  echo "table after transformation:"
  zcat $path_file_product | head -5
  echo "- - Count of lines in product GWAS summary statistics:"
  zcat $path_file_product | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
