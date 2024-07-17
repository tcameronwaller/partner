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
stamp_date=$(date +%Y-%m-%d)
name_base_file_product="$(basename $path_file_table_product .tsv)"
path_directory_product="$(dirname $path_file_table_product)"
path_directory_temporary="${path_directory_product}/temporary_${name_base_file_product}_${stamp_date}" # hopefully unique

# Files.
path_file_temporary_1="${path_directory_temporary}/${name_base_file_product}_temporary_1.txt"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_temporary
mkdir -p $path_directory_temporary

# Initialize file.
rm $path_file_table_product

###############################################################################
# Organize parameters.

###############################################################################
# Execute procedure.

# Copy first line of source file.
head -1 $path_file_table_source > $path_file_temporary_1

# Replace problematic delimiters.
#sed -i 's%\\t%;%g' $path_file_temporary_1

# Collect lines for product table.
lines_product=()
line="identifier_original;identifier_novel;identifier_original_simple"
lines_product+=($line)
# Collect unique novel identifiers.
declare -A identifiers_novel_unique # initialize an associative array
# Read first line of header names for columns in table.
IFS=$'\t' read -r -a array_line < $path_file_temporary_1
# Iterate across header names of columns in table.
for name_full in "${array_line[@]}"; do
  # Initialize and reset the novel identifier.
  identifier_novel=""
  # Extract keys of associative array the unique original header names.
  string_identifiers_check=${!identifiers_novel_unique[@]} # transfers a space-delimited list of keys
  # Extract base name of file.
  name_base="$(basename $name_full .bam)"
  # Generate novel identifier.
  identifier_novel_prefix=`cat /dev/urandom | tr -dc [:upper:] | head -c3`
  identifier_novel_suffix=`cat /dev/urandom | tr -dc [:digit:] | head -c4`
  identifier_novel="${identifier_novel_prefix}-${identifier_novel_suffix}"
  # Make sure that the novel identifier does not already exist in the collection.
  if [[ "${string_identifiers_check}" =~ "${identifier_novel}" ]]; then
    # Generate another novel identifier.
    # The probability of generating an existing identifier twice must be small.
    identifier_novel_prefix=`cat /dev/urandom | tr -dc [:upper:] | head -c3`
    identifier_novel_suffix=`cat /dev/urandom | tr -dc [:digit:] | head -c4`
    identifier_novel="${identifier_novel_prefix}-${identifier_novel_suffix}"
  fi
  if [[ ! "${string_identifiers_check}" =~ "${identifier_novel}" ]]; then
    identifiers_novel_unique[$identifier_novel]=0 # assign a key-value pair
    # Define line.
    line="${name_full};${identifier_novel};${name_base}"
    lines_product+=($line)
  fi
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
  head -1 $path_file_table_source
  echo "----------"
  echo "Temporary table:"
  head -1 $path_file_temporary_1
  echo "----------"
  echo "Product novel table:"
  head -10 $path_file_table_product
  echo "----------"
  echo "----------"
  echo "----------"
fi

##########
# Remove directory of temporary, intermediate files.
rm -r $path_directory_temporary

###############################################################################
# End.
