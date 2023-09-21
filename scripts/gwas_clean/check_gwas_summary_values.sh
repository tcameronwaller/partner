#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 21 September 2023
# Date, last execution: ___ 2023
# Date, review: ___ 2023
################################################################################
# Note



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

#echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_check
zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_file_temporary_check
zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 14)
    # Check for any rows with incorrect count of column fields, indicating empty cells.
    print $0
  else if ( (toupper($4) != "T") && (toupper($4) != "C") && (toupper($4) != "G") && (toupper($4) != "A") )
    # Check effect allele.
    print $0
  else if ( (toupper($4) !~ /^[TCGA]+$/) )
    # Check effect allele.
    print $0
  else if ( (toupper($5) != "T") && (toupper($5) != "C") && (toupper($5) != "G") && (toupper($5) != "A") )
    # Check other allele.
    print $0
  else if ( (toupper($5) !~ /^[TCGA]+$/) )
    # Check other allele.
    print $0
  else if ( (toupper($7) == "NA") )
    # Check effect parameter value for missingness.
    print $0
  else if ( (toupper($8) == "NA") )
    # Check standard error value for missingness.
    print $0
  else if ( (toupper($9) == "NA") || ( ($9 + 0) < 1.0E-307 ) || ( ($9 + 0) > 1.0 ) )
    # Check probability value for missingness or out of range.
    print $0
  else if ( (toupper($10) == "NA") || ( ($10 + 0) < 1.0 ) )
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
  echo "Check values in GWAS summary statistics."
  echo "path to source GWAS file: " $path_file_source
  echo "path to product GWAS file: " $path_file_temporary_check
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Remove temporary, intermediate files.
#rm -r $path_directory_temporary



#
