#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 17 June 2024
# Date, last execution: 17 June 2024
# Review: TCW; 17 June 2024
################################################################################
# Note

# The Next Generation Sequencing Core Facility at Mayo Clinic delivers raw
# sequence data along with a brief report of notes and sample descriptions in a
# text file. This script parses the text file and extracts the sample
# descriptions, which it then writes to a table in text format with delimiters
# of tabls and new lines (tsv).

# These sample descriptions will be necessary to map RNA sequence data to the
# identifiers of the original samples from the experiment.

# Product Format (table in text format with delimiters of tabs and new lines)
# columns: identifier_sequence molecule type identifier_experiment

################################################################################



################################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source report and sample descriptions in text format
path_file_product=${2} # full path to file for product sample descriptions in text format
report=${3} # whether to print reports

################################################################################
# Organize paths.

################################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_product .txt.gz)"
path_directory_product="$(dirname $path_file_product)"

# Initialize directory.
mkdir -p $path_directory_product

# Remove any previous version of the product file.
rm $path_file_product

################################################################################
# Organize parameters.

set +x
set +v

################################################################################
# Execute procedure.

##########
# 1. Find the line with the header immediately before the start of the table of
# sample descriptions.
# awk "/Samples:/{print NR}" $path_file_source
# sed -n "/Samples:/=" $path_file_source
line_header="$(grep -n 'Samples:' $path_file_source | cut -d: -f1)"
line_start=$(($line_header + 1))

##########
# 2. Find the line with the footer after the end of the table of sample
# descriptions.
line_footer="$(grep -n '\----------' $path_file_source | cut -d: -f1)"
line_end=$(($line_footer - 3))

##########
# 1. Parse text file and extract table of sample descriptions.
# Note that AWK interprets a single space delimiter (FS=" ") as any white space.
echo "identifier_sequence molecule type identifier_experiment" | awk 'BEGIN { FS=" "; OFS="\t" } {print $1, $2, $3, $4}' > $path_file_product
#echo "identifier_sequence\tmolecule\ttype\tidentifier_experiment" > $path_file_product
#cat $path_file_source | awk "BEGIN { FS=" "; OFS=" " } NR == $line_start, NR == $line_end {}"
cat $path_file_source | awk -v line_start=$line_start -v line_end=$line_end 'BEGIN { FS=" "; OFS="\t" } NR==line_start, NR==line_end {
    # Print information from the row.
    print $1, $2, $3, $4
  }' >> $path_file_product


################################################################################
# Report.

count_source=$(cat $path_file_source | wc -l)
count_product=$(cat $path_file_product | wc -l)

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "----------"
  echo "Lines of start and end for table of sample descriptions:"
  echo "line start: " $line_start
  echo "line end: " $line_end
  #echo "line footer: " $line_footer
  echo "----------"
  echo "Count of lines in source file:"
  echo $count_source
  echo "Count of lines in product file:"
  echo $count_product
  echo "----------"
  echo "Extract table of sample descriptions:"
  echo "head"
  cat $path_file_product | head -10
  echo "tail"
  cat $path_file_product | tail -10
  echo "----------"
fi



#
