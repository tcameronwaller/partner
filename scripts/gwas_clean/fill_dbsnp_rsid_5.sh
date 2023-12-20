#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 11 December 2023
# Date, last execution: 19 December 2023
# Review: 19 December 2023
################################################################################
# Notes:



################################################################################
# Organize arguments.

path_file_gwas_source=${1} #
path_file_gwas_product=${2} #
path_directory_temporary=${3} #
strict=${4} #
report=${5} #

################################################################################
# Organize paths.

# Temporary files.
path_ftemp_merge_priority_clean_strict="${path_directory_temporary}/merge_priority_clean_strict.txt"
path_ftemp_merge_priority_check="${path_directory_temporary}/merge_priority_check.txt"

###########################################################################
# Execute procedure.



# $1               $2  $3  $4 $5 $6 $7   $8   $9 $10 $11 $12 $13  $14   $15   $16 $17   $18 $19 $20
# identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P   N   Z   INFO NCASE NCONT ID  CHROM POS ALT REF



if true; then
  # Test the procedure by determining the proportion of SNP rsIDs from dbSNP that match those in original GWAS summary statistics.
  echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID CHROM POS ALT REF" > $path_ftemp_merge_priority_check
  cat $path_ftemp_merge_priority_clean_strict | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
    if ( $2 != $16)
      # Keep rows for which the source SNP rsID does not match the dbSNP rsID.
      print $0
    else
      # Skip rows for which source SNP rsID matches the dbSNP rsID.
      next
  }' >> $path_ftemp_merge_priority_check
  # Report.
  if [[ "$report" == "true" ]]; then
    # Calculations.
    # https://www.gnu.org/software/bc/manual/html_mono/bc.html
    count_check=$(cat $path_ftemp_merge_priority_check | wc -l)
    count_total=$(cat $path_ftemp_merge_priority_clean_strict | wc -l)
    proportion_check=$(echo "scale=5; ($count_check / $count_total)" | bc -l) # order of operations and scale rounding matters
    percentage_check=$(echo "scale=5; ($proportion_check * 100)" | bc -l) # order of operations and scale rounding matters
    echo "----------"
    echo "----------"
    echo "----------"
    echo "For special situations where the source GWAS summary statistics"
    echo "already included rsIDs for SNPs, it is possible to test the"
    echo "proportion of the rsIDs extracted from dbSNP that match the"
    echo "originals."
    echo "----------"
    echo "This check uses the strict version that only considers the SNPs that"
    echo "matched and merged with dbSNP rsIDs successfully."
    echo "----------"
    echo "Table that only includes SNPs for which the dbSNP rsID does not match"
    echo "the source rsID:"
    echo "----------"
    head -5 $path_ftemp_merge_priority_check
    echo "----------"
    echo "Lines that do not match between original and dbSNP: " $count_check
    echo "Lines total (strict): " $count_total
    echo "Proportion that do not match: " $proportion_check
    echo "Percentage that do not match: " $percentage_check "%"
    echo "----------"
    echo "----------"
    echo "----------"
  fi
fi



if [[ "$report" == "true" ]]; then
  # Calculations.
  # https://www.gnu.org/software/bc/manual/html_mono/bc.html
  count_source=$(zcat $path_file_gwas_source | wc -l)
  count_strict=$(cat $path_ftemp_merge_priority_clean_strict | wc -l)
  count_product=$(zcat $path_file_gwas_product | wc -l)
  #percentage_poduct=$(echo "100 * ($count_product / $count_source)" | bc -l) # order of operations and scale rounding matters
  proportion_strict=$(echo "scale=5; ($count_strict / $count_source)" | bc -l) # order of operations and scale rounding matters
  percentage_strict=$(echo "scale=5; ($proportion_strict * 100)" | bc -l) # order of operations and scale rounding matters
  count_difference=$(($count_source - $count_strict))
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "Fill dbSNP rsIDs for SNPs in GWAS summary statistics."
  echo "----------"
  echo "Count of original lines in source: " $count_source
  echo "Count of lines that matched and merged with dbSNP (strict): " $count_strict
  echo "Count of lines difference between source and strict: " $count_difference
  echo "Proportion of lines that matched and merged: " $proportion_strict
  echo "Percentage of lines that matched and merged: " $percentage_strict "%"
  echo "----------"
  echo "Use strict filter to SNPs that matched and merged with dbSNP: " $strict
  echo "If active this filter determines count of lines in product."
  echo "----------"
  echo "path to source GWAS file: " $path_file_gwas_source
  echo "path to product GWAS file: " $path_file_gwas_product
  echo "----------"
  echo "source table before transformation:"
  zcat $path_file_gwas_source | head -5
  echo "- - Count of lines in source table: " $count_source
  echo "----------"
  echo "product table after transformation:"
  zcat $path_file_gwas_product | head -5
  echo "- - Count of lines in product table: " $count_product
  echo "----------"
  echo "----------"
  echo "----------"
fi



#
