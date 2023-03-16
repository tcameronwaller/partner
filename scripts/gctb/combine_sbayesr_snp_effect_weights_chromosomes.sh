#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 10 March 2023
################################################################################
# Notes

# This script combines reports of SNP effect weights from GCTB SBayesR across
# chromosomes.

# Source Format and Product Format
# Description: Format of SNP effect weights from SBayesR.
# Documentation site: https://cnsgenomics.com/software/gctb/#Tutorial
# File suffix: ".snpRes"
# File type: text
# File compression: none
# Delimiter: Tab
# Chromosome base position coordinate system: base 1
#   Site: https://www.biostars.org/p/84686/
#   Note: Coordinates designate 1-based integer position of each base
# Columns: Id Name Chrom Position A1 A2 A1Frq A1Effect SE PIP (TCW; 2023-03-10)
#          1  2    3     4        5  6  7     8        9  10
# Column Descriptions
# 1. "Id": sequential integer designator of SNP that SBayesR creates without general significance
#   - Not unique after combining SNPs across chromosomes
# 2. "Name": reference SNP cluster identifier (rsID) of SNP
# 3. "Chrom": chromosome
# 4. "Position": base position in human genome assembly (GRCh37 as of 10 March 2023)
# 5. "A1": effect allele
# 6. "A2": other allele
# 7. "A1Frq": frequency of effect allele
# 8. "A1Effect": LD-adjusted effect
# 9. "SE": standard error of LD-adjusted effect
# 10. "PIP": posterior inclusion probability
#   - "frequency of the SNP being fitted as a non-zero effect in the model across MCMC samples"

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
echo -e "Id\tName\tChrom\tPosition\tA1\tA2\tA1Frq\tA1Effect\tSE\tPIP" > $path_file_product
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
  cat $path_file_source | awk 'BEGIN {FS = " "; OFS ="\t"} NR > 1 {
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
  echo "(Expect that cumulative count to be about 21 greater due to header lines.)"
  count_lines=($(cat $path_file_product | wc -l))
  echo "Total lines in collection product: " $count_lines
  echo "----------"
  echo "----------"
  echo "----------"
fi



#
