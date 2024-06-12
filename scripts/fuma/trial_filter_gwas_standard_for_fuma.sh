#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 12 June 2024
# Date, last execution: 12 June 2024
# Review: TCW; 12 June 2024
################################################################################
# Note

# Trial script to prepare summary statistics that TCW ran on data from UK
# Biobank for upload and analyses in FUMA.

# Source Format (standard for Psychiatric Genetics team)
# effect allele: "A1"
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

# TODO: TCW; 12 June 2024
# TODO: Still need to remove excess columns "Z", "INFO", "NCAS", and "NCONT".

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

################################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_product .txt.gz)"
path_directory_product="$(dirname $path_file_product)"
path_directory_product_temporary="${path_directory_product}/temporary_${name_base_file_product}" # hopefully unique
path_file_temporary_1="${path_directory_product_temporary}/${name_base_file_product}_temporary_1.txt"
path_file_temporary_2="${path_directory_product_temporary}/${name_base_file_product}_temporary_2.txt"
path_file_temporary_3="${path_directory_product_temporary}/${name_base_file_product}_temporary_3.txt"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

# Remove any previous version of the product file.
rm $path_file_product

################################################################################
# Organize parameters.

set +x
set +v

################################################################################
# Execute procedure.

# Source format of GWAS summary statistics was the standard of the Psychiatric
# Genetics team in 2022.
# columns: "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT"

##########
# 1. Filter records in GWAS summary statistics.
# Note that AWK interprets a single space delimiter (FS=" ") as any white space.
#echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_1
zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_file_temporary_1
zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 14)
    # Skip any of table rows with incorrect count of column fields, indicating empty cells.
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
  else if ( ($6 == "") || (toupper($6) == "NA") || (toupper($6) == "NAN") || ( ($6 + 0) < 0 ) )
    # Skip any rows with missing or nonsense allele frequency.
    # Some subsequent analyses do not require this information, but in that case
    # it is possible to fill missing allele frequency with value 0.5.
    next
  else if ( ($7 == "") || (toupper($7) == "NA") || (toupper($7) == "NAN") )
    # Skip any rows with missing effect parameter.
    # Remember that the effect parameter (beta) can be less than zero.
    next
  else if ( ($8 == "") || (toupper($8) == "NA") || (toupper($8) == "NAN") || ( ($8 + 0) < 0 ) )
    # Skip any rows with missing or nonsense less than zero standard error.
    next
  else if ( ($9 == "") || (toupper($9) == "NA") || (toupper($9) == "NAN") || ( ($9 + 0) < 0 ) )
    # Skip any rows with missing or nonsense probability.
    next
  else if ( ($10 == "") || (toupper($10) == "NA") || (toupper($10) == "NAN") || ( ($10 + 0) < 1.0 ) )
    # Skip any rows with missing count of observations (sample size).
    next
  else if ( (toupper(substr($1, 1, 2)) !~ /[RS]/) )
    # Skip any rows with invalid dbSNP reference SNP cluster identifiers (rsIDs).
    next
  else
    # Row passes all checks, filters, and constraints.
    # Print the row entirely.
    print $0
  }' >> $path_file_temporary_1




##########
# Compression

# Compress file format.
gzip -cvf $path_file_temporary_1 > $path_file_product



##########
# Report.

count_source=$(zcat $path_file_source | wc -l)
count_product=$(zcat $path_file_product | wc -l)
count_difference=$(($count_source - $count_product))
proportion_difference=$(echo "scale=10; $count_difference / $count_source" | bc)

if [ "$report" == "true" ]; then
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
fi

if [ "$report" == "true" ] && (( $(echo "$count_difference > 0.0" | bc -l) )); then
  echo "----------"
  echo "There are differences in source and product after filters."
  echo "----------"
  echo "table before transformation:"
  zcat $path_file_source | head -5
  echo "- - Count of lines in source GWAS summary statistics:"
  echo $count_source
  echo "----------"
  echo "table after transformation:"
  zcat $path_file_product | head -5
  echo "- - Count of lines in product GWAS summary statistics:"
  echo $count_product
  echo "----------"
  echo "loss count (difference): " $count_difference
  echo "loss proportion: " $proportion_difference
  echo "----------"
  echo "----------"
  echo "----------"
fi

if [ "$report" == "true" ] && (( $(echo "$proportion_difference >= 0.0001" | bc -l) )); then
  echo "**************************************************"
  echo "**************************************************"
  echo "**************************************************"
  echo "WARNING:"
  echo "Current GWAS summary statistics lost 0.01% or more of SNPs to filters!"
  echo "path to source GWAS file: " $path_file_source
  echo "path to product GWAS file: " $path_file_product
  echo "**************************************************"
  echo "**************************************************"
  echo "**************************************************"
fi


##########
# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
