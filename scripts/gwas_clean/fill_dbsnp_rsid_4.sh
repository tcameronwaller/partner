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

path_file_gwas_product=${1} #
path_directory_temporary=${2} #
strict=${3} #
report=${4} #


################################################################################
# Organize paths.

# Temporary files.
path_ftemp_merge_ref_alt_clean="${path_directory_temporary}/merge_ref_alt_clean.txt"
path_ftemp_merge_priority="${path_directory_temporary}/merge_priority.txt"
path_ftemp_merge_priority_clean="${path_directory_temporary}/merge_priority_clean.txt"
path_ftemp_merge_priority_clean_strict="${path_directory_temporary}/merge_priority_clean_strict.txt"
path_ftemp_product_format="${path_directory_temporary}/gwas_product_format.txt"


###########################################################################
# Execute procedure.



# $1               $2  $3  $4 $5 $6 $7   $8   $9 $10 $11 $12 $13  $14   $15   $16   $17      $18    $19    $20    $21   $22      $23    $24    $25
# identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P   N   Z   INFO NCASE NCONT ID_ar CHROM_ar POS_ar ALT_ar REF_ar ID_ra CHROM_ra POS_ra ALT_ra REF_ra



# Note: TCW; 12 December 2023
# Where there are multiple rsIDs for a SNP, dbSNP records use semicolon ";" delimited lists.
# Split by semicolon and take the first instance in the array to accommodate.

# Determine whether to keep information from Merge 1 or Merge 2.
# First print conditional block: A1 is ALT, A2 is REF, use Merge 1
# Second print conditional block: A1 is REF, A2 is ALT, use Merge 2
# https://www.gnu.org/software/gawk/manual/html_node/String-Functions.html
echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID CHROM POS ALT REF" > $path_ftemp_merge_priority
cat $path_ftemp_merge_ref_alt_clean | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 25)
    # Skip any rows with incorrect count of column fields.
    next
  else if ( ($3 == $17) && ($4 == $18) && ($5 == $19) && ($6 == $20) )
    # A1 is ALT
    # A2 is REF
    # Use Merge 1: ALT_REF
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20
  else if ( ($3 == $22) && ($4 == $23) && ($5 == $25) && ($6 == $24) )
    # A1 is REF
    # A2 is ALT
    # Use Merge 2: REF_ALT
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $21, $22, $23, $24, $25
  else
    # The record did not match or merge with information from dbSNP.
    # Fill with missing information.
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, "NA", "NA", "NA", "NA", "NA"
}' >> $path_ftemp_merge_priority
# GNU Awk can only handle a few operations at a time.
# Split delimited lists of identifiers and keep only the first.
echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID CHROM POS ALT REF" > $path_ftemp_merge_priority_clean
cat $path_ftemp_merge_priority | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  (a = $16); split(a, b, ";"); print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, b[1], $17, $18, $19, $20
}' >> $path_ftemp_merge_priority_clean

# Strict.
echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID CHROM POS ALT REF" > $path_ftemp_merge_priority_clean_strict
cat $path_ftemp_merge_priority_clean | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( ($16 == "NA") && ($17 == "NA") && ($18 == "NA") && ($19 == "NA") && ($20 == "NA") )
    # Skip any row records with missing information from merge.
    next
  else
    # Keep the row record.
    print $0
}' >> $path_ftemp_merge_priority_clean_strict


# $1               $2  $3  $4 $5 $6 $7   $8   $9 $10 $11 $12 $13  $14   $15   $16 $17   $18 $19 $20
# identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P   N   Z   INFO NCASE NCONT ID  CHROM POS ALT REF



##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Table after Merge 1, Merge 2, and prioritization:"
  head -10 $path_ftemp_merge_priority_clean
  cat $path_ftemp_merge_priority_clean | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi



##########
# 5. Adjust format of product GWAS summary statistics.

if [ "$strict" == "true" ]; then

  # The strict version has already been filtered to non-missing matches with
  # dbSNP.

  echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_product_format
  cat $path_ftemp_merge_priority_clean_strict | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
    print $16, $3, $4, toupper($5), toupper($6), $7, $8, $9, $10, $11, $12, $13, $14, $15
  }' >> $path_ftemp_product_format

elif [ "$strict" == "false" ]; then

  # The non-strict version has not yet been filtered to non-missing matches with
  # dbSNP.
  # Only replace the original SNP identifier if the match from dbSNP is not
  # missing.

    echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_product_format
    cat $path_ftemp_merge_priority_clean | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
      if ( ($16 != "NA") && ($17 != "NA") && ($18 != "NA") && ($19 != "NA") && ($20 != "NA") )
        print $16, $3, $4, toupper($5), toupper($6), $7, $8, $9, $10, $11, $12, $13, $14, $15
      else
        print $2, $3, $4, toupper($5), toupper($6), $7, $8, $9, $10, $11, $12, $13, $14, $15
    }' >> $path_ftemp_product_format

fi

# Compress file format.
gzip -cvf $path_ftemp_product_format > $path_file_gwas_product



#
