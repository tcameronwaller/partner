#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 11 September 2025
# Date, last execution or modification: 11 September 2025
# Review: TCW; 11 September 2025
################################################################################
# Note

# MSigDB format
# delimiter between sets in a collection: newline
# delimiter between genes in a set: tab
# file name suffix: ".gmt"

# Recent example of usage:
# /.../pails_process/omega3/2025-09-08_genes_sets_hypotheses_custom

################################################################################
# Organize parameters.

parent_pail="omega3"
pail_process="2025-09-08_genes_sets_hypotheses_custom"
dock_pail_process="pail_${pail_process}"
report="true"

name_file_product="collection_fibrosis"
delimiter_source="newline" # "newline", "tab", "\n", "\t", ";", ":", ",", not " "
delimiter_product_set="tab" # "newline", "tab", "\n", "\t", ";", ":", ",", not " "
delimiter_product_collection="newline" # "newline", "tab", "\n", "\t", ";", ":", ",", not " "
suffix_file_source=".txt"
suffix_file_product=".gmt"
report="true"
#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error

################################################################################
# Organize paths.

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_tools=$(<"$path_directory_paths/path_directory_tools.txt")
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_scripts="$path_directory_process/scripts"
path_directory_package="$path_directory_process/package"
path_directory_package_partner="$path_directory_package/partner"

path_directory_dock="$path_directory_process/dock"
path_directory_data="$path_directory_dock/in_data" # restore script does not modify "in_data" for efficiency
path_directory_demonstration="$path_directory_dock/in_demonstration"
path_directory_parameters="$path_directory_dock/in_parameters"
path_directory_parameters_private="$path_directory_dock/in_parameters_private"

path_directory_endocrinology=$(<"$path_directory_paths/path_directory_mayo_clinic_endocrinology.txt")
path_directory_pails_process="$path_directory_endocrinology/privacy_security/pails_process"
path_directory_pail="$path_directory_pails_process/${parent_pail}/$pail_process"
path_directory_pail_scripts="$path_directory_pail/scripts"
path_directory_pail_parameters="$path_directory_pail/parameters"

#path_directory_source="${path_directory_dock}/${dock_pail_process}/genes_sets_msigdb_raw"
#path_directory_product="${path_directory_dock}/${dock_pail_process}/genes_sets_msigdb_clean"

path_directory_source="${path_directory_dock}/${dock_pail_process}/genes_sets_for_collection"
path_directory_product="${path_directory_dock}/${dock_pail_process}/collection_genes_sets"

stamp_date=$(date +%Y-%m-%d)
path_directory_temporary="${path_directory_product}/temporary_${stamp_date}" # hopefully unique

# Scripts.

# Initialize directories.
#rm -r $path_directory_temporary # caution
#rm -r $path_directory_product # caution
mkdir -p $path_directory_product
mkdir -p $path_directory_temporary
cd $path_directory_product


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


path_file_product="${path_directory_product}/${name_file_product}${suffix_file_product}" # hopefully unique


##########
# Iterate on files.

# Initialize counter.
counter_collection=0

for path_file_source in "${paths_files_source[@]}"; do

  ##########
  # Organize file names and paths.
  # Extract base name of file.
  name_base_file_source="$(basename $path_file_source $suffix_file_source)"
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

  # Read items in array.
  readarray -t items_source < $path_file_source
  count_items=${#items_source[@]}
  index_array_maximum=$(($count_items - 1))
  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo "path: " $path_file_source
    echo "name: " $name_base_file_source
    echo "count of items: " $count_items
    echo "first item: " ${items_source[0]} # notice base-zero indexing
    echo "last item: " ${items_source[$index_array_maximum]}
    echo "----------"
    for item in "${items_source[@]}"; do
      echo $item
    done
  fi

  ##########
  # Write text items to file.
  count_items_source=${#items_source[@]}
  # Initialize counter.
  counter_set=0
  for item in "${items_source[@]}"; do
    # Determine whether item is non-empty.
    if [[ -n "$item" ]]; then
      if [[ $counter_collection -eq "0" ]] && [[ $counter_set -eq "0" ]]; then
        # Initialize entry for the set.
        printf "${name_base_file_source}\t" > $path_file_product
        printf "${name_base_file_source}\t" >> $path_file_product
        # Write item to text in file with appropriate delimiter.
        #echo -n "${item}" >> $path_file_product
        printf "${item}" >> $path_file_product
      elif [[ $counter_collection -eq "0" ]] && [[ $counter_set -gt "0" ]]; then
        # Write item to text in file with appropriate delimiter.
        #echo -n "\t${item}" >> $path_file_product
        printf "\t${item}" >> $path_file_product
      elif [[ $counter_collection -gt "0" ]] && [[ $counter_set -eq "0" ]]; then
        # Initialize entry for the set.
        printf "${name_base_file_source}\t" >> $path_file_product
        printf "${name_base_file_source}\t" >> $path_file_product
        # Write item to text in file with appropriate delimiter.
        #echo -n "${item}" >> $path_file_product
        printf "${item}" >> $path_file_product
      elif [[ $counter_collection -gt "0" ]] && [[ $counter_set -gt "0" ]]; then
        # Write item to text in file with appropriate delimiter.
        #echo -n "\t${item}" >> $path_file_product
        printf "\t${item}" >> $path_file_product
      fi
      # Increment counter.
      ((counter_set++))
    fi
  done

  # Create new line.
  printf "\n" >> $path_file_product
  # Increment counter.
  ((counter_collection++))

done



################################################################################
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "script: assemble_collection_genes_sets_directory.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
fi

##########
# Remove temporary, intermediate files.
rm -r $path_directory_temporary

#
