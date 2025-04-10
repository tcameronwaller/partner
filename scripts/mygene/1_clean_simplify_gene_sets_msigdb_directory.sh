#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 12 February 2025
# Date, last execution or modification: 12 February 2025
# Review: TCW; 12 February 2025
################################################################################
# Note


################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_tools=$(<"$path_directory_paths/path_directory_tools.txt")
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_scripts="$path_directory_process/scripts"

path_directory_dock="$path_directory_process/dock"
path_directory_data="$path_directory_dock/in_data" # restore script does not modify "in_data" for efficiency
path_directory_parameters="$path_directory_dock/in_parameters"
path_directory_parameters_private="$path_directory_dock/in_parameters_private"

path_directory_source="${path_directory_dock}/in_data/mygene_2025-02-14/sets_gene_msigdb_2025-02-14_raw"
path_directory_product="${path_directory_dock}/in_data/mygene_2025-02-14/sets_gene_msigdb_2025-02-14_clean_symbol"
path_directory_temporary="${path_directory_product}/temporary_process" # hopefully unique

# File suffix.
suffix_file_source=".v2024.1.Hs.grp"
suffix_file_product=".txt"

# Scripts.

# Initialize directories.
rm -r $path_directory_temporary
rm -r $path_directory_product
mkdir -p $path_directory_product
mkdir -p $path_directory_temporary
cd $path_directory_product

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

##########
# Remove temporary, intermediate files.
rm -r $path_directory_temporary

#
