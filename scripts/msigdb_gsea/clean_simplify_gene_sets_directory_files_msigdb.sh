#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 12 February 2025
# Date, last execution or modification: 24 September 2025
# Review: TCW; 24 September 2025
################################################################################
# Note

# Recent example of usage:
# /.../pails_process/omega3/2025-09-22_heterogeneity_candidate_adipose_fibrosis

###############################################################################
# Organize arguments.

path_directory_source=${1} # full path to source file in text format from which to read identifiers or names of genes
path_directory_product=${2} # full path to product file to which to write in text format the identifiers or names of genes
suffix_file_source=${3} # text suffix in names of source files
suffix_file_product=${4} # text suffix in names of product files
delimiter_source=${5} # text delimiter between items in source file, not white space
delimiter_product=${6} # text delimiter between items in product file, not white space
report=${7} # whether to print reports to terminal

################################################################################
# Organize paths.

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_tools=$(<"$path_directory_paths/path_directory_tools.txt")
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_scripts="$path_directory_process/scripts"
path_directory_package="$path_directory_process/package"

# Files.

# Scripts.

# Initialize directories.
#rm -r $path_directory_temporary # caution
#rm -r $path_directory_product # caution
#mkdir -p $path_directory_product
#mkdir -p $path_directory_temporary
cd $path_directory_product

# Initialize files.
#rm $path_file_product

###############################################################################
# Organize parameters.

# Parameters.
delimiter_source="newline" # "newline", "tab", "\n", "\t", ";", ":", ",", not " "
delimiter_product="newline" # "newline", "tab", "\n", "\t", ";", ":", ",", not " "
report="true"
#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error

################################################################################
# Execute procedure.

##########
# Collect files from source directory.
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
