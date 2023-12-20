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

path_directory_temporary=${1} #
report=${2} #


################################################################################
# Organize paths.


# Files.
if true; then
  path_ftemp_dbsnp_extraction_ref_alt="${path_directory_reference_dbsnp}/dbsnp_extraction_ref_alt.txt"
fi

# Temporary files.
if false; then
  path_ftemp_dbsnp_extraction_ref_alt="${path_directory_temporary}/dbsnp_extraction_ref_alt.txt"
fi
path_ftemp_merge_alt_ref_clean="${path_directory_temporary}/merge_alt_ref_clean.txt"
path_ftemp_merge_ref_alt="${path_directory_temporary}/merge_ref_alt.txt"
path_ftemp_merge_ref_alt_clean="${path_directory_temporary}/merge_ref_alt_clean.txt"


###########################################################################
# Execute procedure.



##########
# 4.2. Merge 2.
# Table 1: dbSNP extraction with identifier format <CHROM>_<POS>_<REF>_<ALT>.
# Table 1: 7 total columns with merge identifier in column 1.
# Table 2: Table from Merge 1.
# Table 2: 20 total columns with merge identifier in column 1.
# Merged table: 26 total columns with merge identifier in column 1.
# Delimiter: Space
# It is not necessary to print the header row separately.
#echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID_ar CHROM_ar POS_ar ALT_ar REF_ar ID CHROM POS ALT REF RS_ID" > $path_ftemp_merge_ref_alt
awk 'FNR==NR{a[$1]=$2FS$3FS$4FS$5FS$6FS$7; next} {
  if(a[$1]==""){a[$1]="NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"}; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, a[$1]
} END {
  delete a
}' $path_ftemp_dbsnp_extraction_ref_alt $path_ftemp_merge_alt_ref_clean > $path_ftemp_merge_ref_alt



##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "This is your chance to check merge accuracy!"
  echo "Table immediately after Merge 2:"
  head -10 $path_ftemp_merge_ref_alt
  cat $path_ftemp_merge_ref_alt | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi

# $1               $2  $3  $4 $5 $6 $7   $8   $9 $10 $11 $12 $13  $14   $15   $16   $17      $18    $19    $20    $21 $22   $23 $24 $25 $26
# identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P   N   Z   INFO NCASE NCONT ID_ar CHROM_ar POS_ar ALT_ar REF_ar ID  CHROM POS ALT REF RS_ID

# Check and filter the information from the merge.
echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID_ar CHROM_ar POS_ar ALT_ar REF_ar ID_ra CHROM_ra POS_ra ALT_ra REF_ra" > $path_ftemp_merge_ref_alt_clean
#cat $path_ftemp_merge_ref_alt | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_ftemp_merge_ref_alt_clean
cat $path_ftemp_merge_ref_alt | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 26)
    # Skip any rows with incorrect count of column fields.
    next
  else if ( ($1 == "NA") )
    # Missing identifier indicates that records did not match or other problem.
    # Skip the record.
    next
  else if ( ($3 == "NA") && ($4 == "NA") && ($5 == "NA") && ($6 == "NA") )
    # Missingness of match criteria from Table 2 should not occur, but this
    # would indicate that the merge procedure had kept records from Table 1
    # that did not match records in Table 2.
    # Skip the record.
    next
  else
    # Keep record from successful merge.
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25
}' >> $path_ftemp_merge_ref_alt_clean



##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Table after Merge 2 and subsequent filters:"
  head -10 $path_ftemp_merge_ref_alt_clean
  cat $path_ftemp_merge_ref_alt_clean | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi



#
