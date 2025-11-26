#!/bin/bash

################################################################################
# Author: T. Cameron Waller, Ph.D.
# Date, initialization: 17 November 2025
# Date, review or revision: 17 November 2025
################################################################################
# Note

# The purpose of this procedural script is to concatenate together sets of
# items from separate files into new group sets. The script reads and
# interprets a table of parameters that specifies which individual sets,
# corresponding to individual files, belong in each group set. # The individual
# files encode items in sets using text format with new line delimiters between
# items.

# Recent example of usage:

###############################################################################
# Organize arguments.

path_directory_source=${1} # full path to source file in text format from which to read identifiers or names of genes
path_directory_product=${2} # full path to product file to which to write in text format the identifiers or names of genes
path_file_source_table=${3} # full path to source file in text format with tab and newline delimiters
suffix_file_source=${4} # text suffix in names of source files
suffix_file_product=${5} # text suffix in names of product files
delimiter_source=${6} # text delimiter between items in source file, not white space
delimiter_product=${7} # text delimiter between items in product file, not white space
report=${8} # whether to print reports to terminal

################################################################################
# Organize paths.

# Directories.

# Files.

# Scripts.

# Initialize directories.
#rm -r $path_directory_temporary # caution
#rm -r $path_directory_product # caution
mkdir -p $path_directory_product
#mkdir -p $path_directory_temporary
cd $path_directory_product

# Initialize files.

###############################################################################
# Organize parameters.

# Parameters.
#delimiter_source="newline" # "newline", "tab", "\n", "\t", ";", ":", ",", not " "
#delimiter_product="newline" # "newline", "tab", "\n", "\t", ";", ":", ",", not " "
report="true"
#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error

################################################################################
# Execute procedure.

##########
# Read and collect names of files from source directory.
#cd $path_directory_source
# Bash version 4.4 introduced the "-d" option for "readarray".
#readarray -d "" -t paths_files_source < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.txt.gz" -print0)
paths_files_source=()
while IFS= read -r -d $'\0'; do
  paths_files_source+=("$REPLY")
done < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*${suffix_file_source}" -print0)
count_paths_files_source=${#paths_files_source[@]}

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Source directory:"
  echo $path_directory_source
  echo "count of files: " $count_paths_files_source
  echo "first file: " ${paths_files_source[0]} # notice base-zero indexing
  echo "last file: " ${paths_files_source[$count_paths_files_source - 1]}
  echo "Product directory:"
  echo $path_directory_product
  echo "----------"
fi

##########
# Read table of parameters and collect unique names for groups of sets.

# Initialize array for collection.
names_groups=()

# Read lines from file and split fields within each line by space, tab, or new-line delimiters.
input=$path_file_table_parameter
while IFS=$' \t\n' read -r -a array
do

  # Extract values from individual columns within table's row.
  # These parameters match the format of the table as of 13 December 2023.
  raw_execution="${array[0]}"
  raw_sequence="${array[1]}"
  raw_category="${array[2]}"
  raw_name_group="${array[3]}"
  raw_name_set="${array[4]}"
  raw_size="${array[5]}"
  raw_date_accession="${array[6]}"
  raw_note="${array[7]}"

  # Report.
  if [ $raw_execution == "1" ] && [ "$report" == "true" ]; then
    echo "----------"
    echo "field 0, execution: ${raw_execution}"
    echo "field 1, sequence: ${raw_sequence}"
    echo "field 2, category: ${raw_category}"
    echo "field 3, name_group: ${raw_name_group}"
    echo "field 4, name_set: ${raw_name_set}"
    echo "field 5, size: ${raw_size}"
    echo "field 6, date_accession: ${raw_date_accession}"
    echo "field 7, note: ${raw_note}"
    echo "----------"
  fi
  # Execute procedure for current record's parameters.
  if [ $raw_execution == "1" ]; then
    # Collect name of group.
    names_groups+=("${raw_name_group}")
  fi
done < "${input}"

# Collect unique names of groups.

IFS=" " read -r -a names_groups_unique <<< "$(echo "${names_groups[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')"



######################3

# first iteration across rows of table
# collect names of groups (check if "execution == 1")
# collect unique names of groups

# iteration across names of groups
# iterate on each name of groups
# iterate on each row of the table
# for each row of the table...
# check if "execution == 1"
# check if name of group matches current iteration group
# check if matching file exists for the gene set's name of the row
# # c



###########################
# SCRAP
################################


##########
# Iterate on files.
for path_file_source in "${paths_files_source[@]}"; do

  ##########
  # Organize file names and paths.
  # Extract base name of file.
  name_base_file_product="$(basename $path_file_source $suffix_file_source)"
  name_file_product="${name_base_file_product}${suffix_file_product}"
  path_file_product="${path_directory_product}/${name_file_product}" # hopefully unique
  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo "Source file path:"
    echo $path_file_source
    echo "Product file path:"
    echo $path_file_product
    echo "----------"
  fi

  ##########
  # Read text items from file to array.
  # Initialize array.
  items_source=() # lines
  # Read text items from file using delimiters such as new line.
  input=$path_file_source
  if [[ "$delimiter_source" == "newline" ]]; then
    while IFS=$'\n' read -r -a item
    do
    # Report.
    #if [ "$report" == "true" ]; then
    #  echo "----------"
    #  echo "item: ${item}"
    #  echo "----------"
    #fi
    # Collect.
    items_source+=("${item}")
  done < <(tail -n +0 "${input}"; echo) # append new line to tail to ensure read of last line
  fi

  ##########
  # Write text items to file.
  # Skip first two items in array.
  count_items_source=${#items_source[@]}
  for item in "${items_source[@]:2}"; do
    # Determine whether item is non-empty.
    if [[ -n "$item" ]]; then
      # Write text item to file.
      echo $item >> $path_file_product
    fi
  done
  # Read items from product file.
  readarray -t items_product < $path_file_product
  count_items_product=${#items_product[@]}
  #index_array_maximum=$(($count_items_product))
  index_array_maximum=$(($count_items_product - 1))
  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo "source file"
    echo "count of items in source file: " $count_items_source
    echo "first item: " ${items_source[0]} # notice base-zero indexing
    echo "last item: " ${items_source[$count_items_source - 1]}
    #echo "last item: " ${items_source[$count_items_source]}
    echo "----------"
    echo "count of items in product file: " $count_items_product
    echo "maximum array index: " $index_array_maximum
    echo "first item: " ${items_product[0]} # notice base-zero indexing
    echo "last item: " ${items_product[$index_array_maximum]}
    echo "----------"
  fi


done



################################################################################
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "script: clean_simplify_gene_sets_msigdb_directory.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
fi


#
