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
path_directory_temporary=${2} #
report=${3} #

################################################################################
# Organize paths.


# Files.
if true; then
  path_ftemp_dbsnp_extraction_alt_ref="${path_directory_reference_dbsnp}/dbsnp_extraction_alt_ref.txt"
  path_ftemp_dbsnp_extraction_ref_alt="${path_directory_reference_dbsnp}/dbsnp_extraction_ref_alt.txt"
fi

# Temporary files.
if false; then
  path_ftemp_dbsnp_extraction_alt_ref="${path_directory_temporary}/dbsnp_extraction_alt_ref.txt"
  path_ftemp_dbsnp_extraction_ref_alt="${path_directory_temporary}/dbsnp_extraction_ref_alt.txt"
fi
path_ftemp_gwas_identifier_a1_a2="${path_directory_temporary}/gwas_source_identifier_a1_a2.txt"
path_ftemp_merge_alt_ref="${path_directory_temporary}/merge_alt_ref.txt"
path_ftemp_merge_alt_ref_clean="${path_directory_temporary}/merge_alt_ref_clean.txt"


###########################################################################
# Execute procedure.



##########
# 3. For GWAS summary statistics, assemble information with unique identifiers
# specific to site (chromosome, position) and both effect and other alleles.
# Identifier format: <CHROM>_<POS>_<ALT>_<REF>
# GWAS summary statistics.
echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_gwas_identifier_a1_a2
zcat $path_file_gwas_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print ($2"_"$3"_"$4"_"$5), $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14
}' >> $path_ftemp_gwas_identifier_a1_a2



##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "dbSNP extraction identifier: <CHROM>_<POS>_<ALT>_<REF>"
  echo "- - table before merge:"
  head -5 $path_ftemp_dbsnp_extraction_alt_ref
  echo "----------"
  echo "dbSNP extraction identifier: <CHROM>_<POS>_<REF>_<ALT>"
  echo "- - table before merge:"
  head -5 $path_ftemp_dbsnp_extraction_ref_alt
  echo "----------"
  echo "GWAS summary statistics identifier: <CHROM>_<POS>_<A1>_<A2>"
  echo "- - table before merge:"
  head -5 $path_ftemp_gwas_identifier_a1_a2
  echo "----------"
  echo "----------"
  echo "----------"
fi



##########
# 4. Merge dbSNP information to GWAS summary statistics.
# Merge text tables by a common identifier in "awk".
# Bash command "join" might also be capable of this type of merge or join, but
# it might require the same identifiers in sort order.
# These merges treat Table 2 as the priority table.
# The merge of Table 1 with Table 2 will include all records from Table 2 but
# only those records from Table 1 that have matching identifiers.
# These merges print each row of Table 2 regardless of whether there is a record
# with matching identifier in Table 1.

# Process explanation:
# 1. GNU awk reads lines from the first file (FNR==NR) into an array in memory
# (array "a") that uses values from the first column of the first file (column
# "identifier_merge") as a hashable index (a[$1]=...).
# 2. GNU awk then proceeds to the second file and prints line by line while
# including any of the hashable lines from the first file that have a matching
# index, this time from the first column of the second file (column
# "identifier_merge").

# Reference:
# https://stackoverflow.com/questions/32481877/what-are-nr-and-fnr-and-what-does-nr-fnr-imply


##########
# 4.1. Merge 1.
# Table 1: dbSNP extraction with identifier format <CHROM>_<POS>_<ALT>_<REF>.
# Table 1: 7 total columns with merge identifier in column 1.
# Table 2: GWAS summary statistics
# Table 2: 15 total columns with merge identifier in column 1.
# Merged table: 21 total columns with merge identifier in column 1.
# Delimiter: Space
# It is not necessary to print the header row separately.
#echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID CHROM POS ALT REF RS_ID" > $path_ftemp_merge_alt_ref
awk 'FNR==NR{a[$1]=$2FS$3FS$4FS$5FS$6FS$7; next} {
  if(a[$1]==""){a[$1]="NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"}; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, a[$1]
} END {
  delete a
}' $path_ftemp_dbsnp_extraction_alt_ref $path_ftemp_gwas_identifier_a1_a2 > $path_ftemp_merge_alt_ref

##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "This is your chance to check merge accuracy!"
  echo "Table immediately after Merge 1:"
  head -10 $path_ftemp_merge_alt_ref
  cat $path_ftemp_merge_alt_ref | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi

# $1               $2  $3  $4 $5 $6 $7   $8   $9 $10 $11 $12 $13  $14   $15   $16 $17   $18 $19 $20  $21
# identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P   N   Z   INFO NCASE NCONT ID  CHROM POS ALT REF RS_ID

# Check and filter the information from the merge.
echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID_ar CHROM_ar POS_ar ALT_ar REF_ar" > $path_ftemp_merge_alt_ref_clean
#cat $path_ftemp_merge_alt_ref | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_ftemp_merge_alt_ref_clean
cat $path_ftemp_merge_alt_ref | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 21)
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
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20
}' >> $path_ftemp_merge_alt_ref_clean

##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Table after Merge 1 and subsequent filters:"
  head -10 $path_ftemp_merge_alt_ref_clean
  cat $path_ftemp_merge_alt_ref_clean | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi



#
