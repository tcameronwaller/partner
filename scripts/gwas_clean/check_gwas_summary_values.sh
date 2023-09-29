#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 21 September 2023
# Date, last execution: 29 September 2023
# Date, review: 29 September 2023
################################################################################
# Note

# This script accompanies script file "filter_constrain_gwas_summary_values.sh".

################################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_directory_parent_temporary=${2} # full path to parent directory for temporary child directory and files
report=${3} # whether to print reports

################################################################################
# Organize paths.

name_base_file_temporary="$(basename $path_file_source .txt.gz)"
path_directory_temporary="${path_directory_parent_temporary}/temporary_${name_base_file_temporary}" # hopefully unique

path_file_temporary_check="${path_directory_temporary}/${name_base_file_temporary}_check.txt"

# Initialize directory.
rm -r $path_directory_temporary
mkdir -p $path_directory_temporary

# Remove any previous version of the product file.
rm $path_file_temporary_check

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

##########
# TCW; 29 September 2023
# Test regular expressions with an online tool (https://chortle.ccsu.edu/finiteautomata/Section07/sect07_12.html).
#else if ( (toupper($4) ~ /[^T]/) && (toupper($4) ~ /[^C]/) && (toupper($4) ~ /[^G]/) && (toupper($4) ~ /[^A]/) )
#  # Check effect allele for any characters other than "T", "C", "G", or "A".
#  # This condition is too inclusive because it checks for each character individually.
#  # Example match: TCATC (effect allele)
#  print $0
#else if ( (toupper($5) ~ /[^T]/) && (toupper($5) ~ /[^C]/) && (toupper($5) ~ /[^G]/) && (toupper($5) ~ /[^A]/) )
#  # Check other allele for any characters other than "T", "C", "G", or "A".
#  # This condition is too inclusive because it checks for each character individually.
#  # Example match: TATATGTGT (other allele)
#  print $0
#else if ( (toupper($4) ~ /[^T][^C][^G][^A]/) )
#  # Check effect allele for any characters other than "T", "C", "G", or "A".
#  # I do not know for sure what this is doing.
#  print $0
#else if ( (toupper($4) !~ /T*/) || (toupper($4) !~ /C*/) || (toupper($4) !~ /G*/) || (toupper($4) !~ /A*/) )
#  # This condition would only match strings that consisted entirely of one or more "T", "C", "G", or "A".
#  print $0
#else if ( (toupper($4) ~ /[^TCGAtcga]*/) )
#  # This condition will only match strings that do not have any of the characters "T", "C", "G", "A", etc.
#  # Without the asterisk, this condition would only perform properly on strings with a single character.
#  # Example matches: "D", "E", "DEX"
#  # Example not matches: "TAT", "TAD", "DAT", "DEXT"
#  print $0

##########

#echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_check
zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_file_temporary_check
zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 14)
    # Check for any rows with incorrect count of column fields, indicating empty cells.
    print $0
  else if ( (toupper($3) == "NA") || (toupper($3) == "NAN") || ( ($3 + 0) < 1.0 ) )
    # Check base position coordinate for missingness.
    print $0
  else if ( ($4 == "") || (toupper($4) == "NA") || (toupper($4) == "NAN") )
    # Check effect allele for missing values.
    print $0
  else if ( (toupper($4) !~ /[TCGAtcga]*/) )
    # Check effect allele for any characters other than "T", "C", "G", or "A".
    # This condition matches a string that consists of any characters other than "T", "C", "G", or "A".
    # The "toupper" transformation is unnecessary.
    # Without the asterisk, this condition would only perform properly on strings with a single character.
    # Example matches: "D", "DEX", "TAD", "TATDCA"
    # Example not matches: "T", "TAT", "TATGCA"
    print $0
  else if ( ($5 == "") || (toupper($5) == "NA") || (toupper($5) == "NAN") )
    # Check other allele for missing values.
    print $0
  else if ( (toupper($5) !~ /[TCGAtcga]*/) )
    # Check other allele for any characters other than "T", "C", "G", or "A".
    print $0
  else if ( (toupper($6) == "NA") || (toupper($6) == "NAN") || ( ($6 + 0) < 0 ) || ( ($6 + 0) > 1.0 ) )
    # Check allele frequency value for missingness or out of range.
    print $0
  else if ( ( ($6 + 0) > 0 ) && ( ($6 + 0) < 1.0E-307 ) )
    # Check allele frequency value for missingness or out of range.
    print $0
  else if ( (toupper($7) == "NA") || (toupper($7) == "NAN") )
    # Check effect parameter value for missingness.
    print $0
  else if ( (toupper($8) == "NA") || (toupper($8) == "NAN") )
    # Check standard error value for missingness.
    print $0
  else if ( (toupper($9) == "NA") || (toupper($9) == "NAN") || ( ($9 + 0) < 1.0E-307 ) || ( ($9 + 0) > 1.0 ) )
    # Check probability value for missingness or out of range.
    print $0
  else if ( (toupper($10) == "NA") || (toupper($10) == "NAN") || ( ($10 + 0) < 1.0 ) )
    # Check count of observations (sample size) for missingness.
    print $0
  else
    # Skip any rows that pass all checks.
    next
  }' >> $path_file_temporary_check



# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "Check values in GWAS summary statistics."
  echo "path to source GWAS file: " $path_file_source
  echo "path to product GWAS file: " $path_file_temporary_check
  echo "- - Count of lines in source GWAS summary statistics:"
  zcat $path_file_source | wc -l
  echo "- - Count of lines in product GWAS summary statistics:"
  cat $path_file_temporary_check | wc -l
  echo "----------"
  echo "first 25 lines of product GWAS summary statistics:"
  cat $path_file_temporary_check | head -25
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Remove temporary, intermediate files.
#rm -r $path_directory_temporary



#
