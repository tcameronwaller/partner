#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 21 September 2023
# Date, last execution: 15 November 2023
# Date, review: 15 November 2023
################################################################################
# Note

# This script accompanies script file "check_gwas_summary_values.sh".

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

# Here are some useful regular expressions to evaluate values in "awk".
# ( $12 ~ /^[0-9]+$/ ); ( $12 ~ /^[[:alpha:]]+$/ ); ( $12 ~ /^[[:punct:]]+$/ )

##########
# The 64-bit floating point precision (double precision) can represent values
# from +/- 2.23E-308 to +/- 1.80E308 (https://github.com/bulik/ldsc/issues/144).

# ( (toupper($4) != "T") && (toupper($4) != "C") && (toupper($4) != "G") && (toupper($4) != "A") )

#else if ( (toupper($4) !~ /^[TCGA]+$/) )
#  # Check effect allele.
#  print $0

#else if ( (toupper($5) !~ /^[TCGA]+$/) )
#  # Check other allele.
#  print $0


#echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_check
zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_file_temporary
zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 14)
    # Skip any of table rows with incorrect count of column fields, indicating empty cells.
    next
  else if ( ($3 == "") || (toupper($3) == "NA") || (toupper($3) == "NAN") || ( ($3 + 0) < 1.0 ) )
    # Skip any rows with missing base position coordinate.
    next
  else if ( ($4 == "") || (toupper($4) == "NA") || (toupper($4) == "NAN") )
    # Skip any rows with missing effect allele.
    next
  else if ( (toupper($4) !~ /[TCGAtcga]*/) )
    # Skip any rows with nonsense effect allele.
    next
  else if ( ($5 == "") || (toupper($5) == "NA") || (toupper($5) == "NAN") )
    # Skip any rows with missing other allele.
    next
  else if ( (toupper($5) !~ /[TCGAtcga]*/) )
    # Skip any rows with nonsense other allele.
    next
  else if ( ($6 == "") || (toupper($6) == "NA") || (toupper($6) == "NAN") || ( ($6 + 0) < 0 ) || ( ($6 + 0) > 1.0 ) )
    # Skip any rows with missing or nonsense allele frequency.
    next
  else if ( ( ($6 + 0) > 0 ) && ( ($6 + 0) < 1.0E-307 ) )
    # Constrain allele frequency.
    print $1, $2, $3, $4, $5, (1.0E-307), $7, $8, $9, $10, $11, $12, $13, $14
  else if ( ($7 == "") || (toupper($7) == "NA") || (toupper($7) == "NAN") || ( ($7 + 0) < 0 ) )
    # Skip any rows with missing effect parameter.
    next
  else if ( ($8 == "") || (toupper($8) == "NA") || (toupper($8) == "NAN") || ( ($8 + 0) < 0 ) )
    # Skip any rows with missing standard error.
    next
  else if ( ($9 == "") || (toupper($9) == "NA") || (toupper($9) == "NAN") || ( ($9 + 0) < 0 ) )
    # Skip any rows with missing or nonsense probability.
    next
  else if ( ( ($9 + 0) > 0 ) && ( ($9 + 0) < 1.0E-307 ) )
    # Constrain probability.
    print $1, $2, $3, $4, $5, $6, $7, $8, (1.0E-307), $10, $11, $12, $13, $14
  else if ( ( ($9 + 0) > 1.0 ) )
    # Constrain probability.
    print $1, $2, $3, $4, $5, $6, $7, $8, (1.0), $10, $11, $12, $13, $14
  else if ( (toupper($10) == "NA") || (toupper($10) == "NAN") || ( ($10 + 0) < 1.0 ) )
    # Skip any rows with missing count of observations (sample size).
    next
  else
    # Row passes all checks, filters, and constraints.
    # Print the row entirely.
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
  echo "table before format translation:"
  zcat $path_file_source | head -5
  echo "- - Count of lines in source GWAS summary statistics:"
  zcat $path_file_source | wc -l
  echo "----------"
  echo "table after format translation:"
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
