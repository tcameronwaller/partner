#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 16 July 2024
# Date, last execution or modification: 16 July 2024
# Review: TCW; 16 July 2024
###############################################################################
# Note

# The goal of this script is to extract the original identifiers or names
# across the headers of columns in a table report from quantification of
# RNA sequence reads in HTSeq.

###############################################################################
# Organize arguments.

path_file_table_source=${1} # full path to source file for original table
path_file_table_product=${2} # full path to product file for novel table
report=${3} # whether to print reports to terminal

###############################################################################
# Organize paths.

# Directories.
path_directory_product="$(dirname $path_file_table_product)"

# Files.

# Initialize directory.
mkdir -p $path_directory_product

# Initialize file.
rm $path_file_table_product

###############################################################################
# Organize parameters.

###############################################################################
# Execute procedure.

# Collect lines for product table.
lines_product=()
line="name_full;name_base"
lines_product+=($line)
# Read first line of header names for columns in table.
IFS=$'\t' read -r -a array_line < $path_file_table_source
# Iterate across header names of columns in table.
for name_full in "${array_line[@]}"; do
  # Extract base name of file.
  name_base="$(basename $name_full .tsv)"
  # Define line.
  line="${name_full};${name_base}"
  lines_product+=($line)
done

# Write to file lines for table.
for line in "${lines_product[@]}"; do
  echo $line >> $path_file_table_product
done

##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Source original table:"
  head -10 $path_file_table_source
  echo "----------"
  echo "Product novel table:"
  head -10 $path_file_table_product
  echo "----------"
  echo "----------"
  echo "----------"
fi


###############################################################################
# End.
