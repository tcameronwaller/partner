#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Notes

# This script combines reports of SNP effect weights from GCTB SBayesR across
# chromosomes.

# Format of SNP effect weights from SBayesR.
# file suffix: ".snpRes"
# delimiter: white space
# columns: Id Name Chrom Position A1 A2 A1Frq A1Effect SE PIP (TCW; 2023-03-10)
#          1  2    3     4        5  6  7     8        9  10

# review: TCW; 10 March 2023

###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize arguments.

path_directory_source=${1} # full path to parent directory for source files from GCTB SBayesR
name_file_source_prefix=${2} # prefix of name of file for source files from GCTB SBayesR
name_file_source_suffix=${3} # suffix of name of file for source files from GCTB SBayesR
path_file_product=${4} # full path to file for product combination of information from GCTB SBayesR
chromosome_x=${5} # whether to include Chromosome X
report=${6} # whether to print reports

################################################################################
# Organize paths.

path_directory_product="$(dirname $path_file_product)"

# Initialize directory.
mkdir -p $path_directory_product

# Initialize files.
rm $path_file_product

###########################################################################
# Execute procedure.

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Collection and combination of rows from multiple files for chromosomes."
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Initialize cumulative variable.
count_lines_cumulative=0

# Collect information.
echo "Id Name Chrom Position A1 A2 A1Frq A1Effect SE PIP" > $path_file_product
# Iterate on relevant chromosomes.
if [[ "$chromosome_x" == "true" ]]; then
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
else
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
fi
for chromosome in "${chromosomes[@]}"; do
  # Define path to source file for current chromosome.
  name_file_source="${name_file_source_prefix}${chromosome}${name_file_source_suffix}"
  path_file_source="${path_directory_source}/${name_file_source}"
  # Collect and combine all lines except for first (header).
  cat $path_file_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
    # Print the row entirely.
    print $0
  } ' >> $path_file_product
  # Report.
  if [[ "$report" == "true" ]]; then
    count_lines=($(cat $path_file_source | wc -l))
    ((count_lines_cumulative += $count_lines))
    #count_lines_cumulative=$((count_lines_cumulative + $count_lines))
    echo "Chromosome " $chromosome " original lines: " $count_lines
  fi
done

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Collection complete."
  echo "Cumulative count of lines: " $count_lines_cumulative
  count_lines=($(cat $path_file_product | wc -l))
  echo "Total lines in collection product: " $count_lines
  echo "----------"
  echo "----------"
  echo "----------"
fi



#
