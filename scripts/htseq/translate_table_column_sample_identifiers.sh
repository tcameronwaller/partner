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



###############################################################################
# Organize arguments.

path_file_column_translations=${1} # full path to parameter file of translations for names of columns
path_file_table_source=${2} # full path to source file for original table
path_file_table_product=${3} # full path to product file for novel table
report=${4} # whether to print reports to terminal

###############################################################################
# Organize paths.

# Directories.
stamp_date=$(date +%Y-%m-%d)
name_base_file_product="$(basename $path_file_table_product .tsv)"
path_directory_product="$(dirname $path_file_table_product)"
path_directory_temporary="${path_directory_product}/temporary_${name_base_file_product}_${stamp_date}" # hopefully unique

# Files.
path_file_temporary_1="${path_directory_temporary}/${name_base_file_product}_temporary_1.tsv"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_temporary
mkdir -p $path_directory_temporary

# Initialize file.
rm $path_file_temporary_1
rm $path_file_table_product

###############################################################################
# Organize parameters.

###############################################################################
# Execute procedure.

# TODO: IF the column name appears in the reference translations

# Extract translations of header names for columns in table.
declare -A translations # initialize an associative array
input=$path_file_column_translations
while IFS=$'\n' read -r -a array_lines; do
  for line in "${array_lines[@]}"; do
    # Separate segments within current line.
    IFS=$'\t' read -r -a array_segments <<< "${line}"
    # Select original and novel identifiers for samples from segments in current line.
    name_original="${array_segments[0]}"
    name_novel="${array_segments[1]}"
    # Skip the header line.
    if [[ $name_original == "identifier_original" ]]; then
      continue
    fi
    # Organize identifiers within associative array.
    translations[$name_original]=$name_novel # assign a key-value pair
  done
done < "${input}"
# Extract keys of associative array the unique original header names.
#string_names_original=${!translations[@]} # transfers a space-delimited list of keys
IFS=$' ' read -r -a array_names_original <<< "${!translations[@]}"
count_names_original=${#array_names_original[@]}
IFS=$' '; string_names_original="${array_names_original[*]}"

# Copy table from source file.
cp $path_file_table_source $path_file_temporary_1

# The command call to "sed" cannot use a delimiter that also occurs in the
# original or novel strings for replacement.
# The first option is to use a delimiter in the "sed" command that does not
# appear in the strings for replacement.
# The second option is to replace the problematic delimiter temporarily and
# then revert back after the main replacement.
# For example: sed -i -e 's/\/foo-bar-blargh/g' $path_file_temporary

# Read first line of header names for columns in table.
IFS=$'\t' read -r -a array_line < $path_file_table_source
# Iterate across header names of columns in table.
# By iterating across the header names of columns in the source table rather
# than those from the parameter table, there is some control over which
# translations to use.
for name_original in "${array_line[@]}"; do
  if [[ "${string_names_original}" =~ "${name_original}" ]]; then
    name_novel="${translations[$name_original]}"
    #echo "original: $name_original"
    #echo "novel: $name_novel"
    sed -i -e "s%${name_original}%${name_novel}%g" $path_file_temporary_1 # replace name throughout file
  fi
done

# Copy table to product file.
cp $path_file_temporary_1 $path_file_table_product

##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Source original table:"
  head -10 $path_file_table_source
  echo "----------"
  echo "Count of names for translation: $count_names_original"
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
