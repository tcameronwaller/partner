#!/bin/bash

################################################################################
# Author: T. Cameron Waller, Ph.D.
# Date, initialization: 17 November 2025
# Date, review or revision: 20 November 2025
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
#cd ~
#path_directory_paths="./Downloads/paths_process_local"
#path_directory_tools=$(<"$path_directory_paths/path_directory_tools.txt")
#path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
#path_directory_scripts="$path_directory_process/scripts"
#path_directory_package="$path_directory_process/package"

# Files.

# Scripts.

# Initialize directories.
#rm -r $path_directory_temporary # caution
#rm -r $path_directory_product # caution
mkdir -p $path_directory_product
#mkdir -p $path_directory_temporary
cd $path_directory_product

# Initialize files.
#rm $path_file_product

###############################################################################
# Organize parameters.

# Parameters.

# Delimiters.
# "space", "newline", "tab", "\n", "\t", ";", ":", ","
if [ $delimiter_source == "space" ]; then
  #IFS=$' \t\n'
  delimiter_source=" \t\n"
fi
if [ $delimiter_product == "newline" ]; then
  delimiter_product="\n"
fi

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
input=$path_file_source_table
while IFS=$' \t\n' read -r -a array_line
do
  # Extract values from individual columns within table's row.
  raw_execution="${array_line[0]}"
  raw_sequence="${array_line[1]}"
  raw_category="${array_line[2]}"
  raw_name_group="${array_line[3]}"
  raw_name_set="${array_line[4]}"
  raw_size="${array_line[5]}"
  raw_date_accession="${array_line[6]}"
  raw_note="${array_line[7]}"
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

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "unique names of groups for sets"
  echo "----------"
fi
for name_group in "${names_groups_unique[@]}"; do
  # Report.
  if [[ "$report" == "true" ]]; then
    echo $name_group
  fi
done
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Report.
if [[ "$report" == "true" ]]; then
  echo -e "\n"
  echo -e "\n"
  echo -e "\n"
  echo "----------"
  echo "----------"
  echo "----------"
  echo "collection of items from sets belonging to each group"
  echo "----------"
  echo "----------"
  echo "----------"
  echo -e "\n"
  echo -e "\n"
  echo -e "\n"
fi

# Iterate on names for new groups of sets.
for name_group in "${names_groups_unique[@]}"; do

  # Initialize array for collection.
  names_sets_group=()
  # Read table of parameters and collect names and paths to files for sets in current group.
  input=$path_file_source_table
  while IFS=$' \t\n' read -r -a array_line
  do
    # Extract values from individual columns within table's row.
    raw_execution="${array_line[0]}"
    raw_sequence="${array_line[1]}"
    raw_category="${array_line[2]}"
    raw_name_group="${array_line[3]}"
    raw_name_set="${array_line[4]}"
    raw_size="${array_line[5]}"
    raw_date_accession="${array_line[6]}"
    raw_note="${array_line[7]}"
    # Execute procedure for current record's parameters.
    if [ $raw_execution == "1" ] && [ "$raw_name_group" == "$name_group" ]; then
      # Collect name of set in current group.
      names_sets_group+=("${raw_name_set}")
    fi
  done < "${input}"

  # Collect unique names of sets in current group.
  IFS=" " read -r -a names_sets_group_unique <<< "$(echo "${names_sets_group[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')"

  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo "current group: $name_group"
    echo "----------"
    echo "unique names of sets in current group"
    echo "----------"
  fi
  for name_set_group in "${names_sets_group_unique[@]}"; do
    # Report.
    if [[ "$report" == "true" ]]; then
      echo $name_set_group
    fi
  done
  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo -e "\n"
    echo -e "\n"
    echo -e "\n"
    echo "----------"
    echo "counts of items in each set belonging to current group"
    echo "set name: set size"
    echo "----------"
    echo -e "\n"
  fi

  # Initialize array for collection.
  items_sets_group=()
  # Iterate on names of sets in current group.
  # Read files for each individual set and collect items.
  for name_set_group in "${names_sets_group_unique[@]}"; do
    # Initialize array for filter and collection of items in current set.
    items_set=()
    # Determine name and path of file.
    name_file_set="${name_set_group}${suffix_file_source}"
    path_file_set="${path_directory_source}/${name_file_set}"
    # Determine whether file exists.
    if [[ -f "$path_file_set" ]]; then
      # Read items from file.
      input=$path_file_set
      while IFS=$delimiter_source read -r -a line
      do
        # Determine whether the line is empty or only contains white space.
        #[ -z "$item" ] && continue
        if [[ ! -z "$line" ]] && \
           [[ -n "$line" ]] && \
           [[ $line =~ [^[:space:]] ]]; then
          # Current line is not empty and contains characters other than white
          # space.
          # Collect item.
          items_set+=("${line}")
        else
          # Skip any empty lines with only white space.
          continue
        fi
      done < <(tail -n +0 "${input}"; echo) # append new line to tail to ensure read of last line
    fi
    # Collect items from all sets.
    items_sets_group+=( "${items_sets_group[@]}" "${items_set[@]}" )
    # Report.
    if [[ "$report" == "true" ]]; then
      # Organize information.
      count_items=${#items_set[@]}
      # Print information.
      echo "${name_set_group}: ${count_items}"
    fi
  done

  # Collect unique items for current group of sets.
  IFS=" " read -r -a items_sets_group_unique <<< "$(echo "${items_sets_group[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')"

  # Report.
  if [[ "$report" == "true" ]]; then
    # Organize information.
    IFS=" "
    count_items=${#items_sets_group_unique[@]}
    # Print information.
    echo "----------"
    echo -e "\n"
    echo "----------"
    echo "count of unique items from all sets in group:"
    echo $count_items
    echo "----------"
    echo "first item: " ${items_sets_group_unique[0]} # notice base-zero indexing
    echo "last item: " ${items_sets_group_unique[$count_items - 1]}
    #echo "${items_sets_group_unique[*]}"
    echo "----------"
    echo -e "\n"
    echo -e "\n"
    echo -e "\n"
  fi

  # Determine name and path of file for current group of sets.
  name_file_group="${name_group}${suffix_file_product}"
  path_file_group="${path_directory_product}/${name_file_group}"
  # Initialize file.
  rm $path_file_group # caution
  # Write to file unique items from sets in current group.
  for item in "${items_sets_group_unique[@]}"; do
    # Write text item to file.
    echo $item >> $path_file_group
  done
  #IFS=$delimiter_product
  #echo "${items_sets_group_unique[*]}" > $path_file_group
  #printf '%s' "${items_sets_group_unique[*]}" > $path_file_group
done

################################################################################
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "script: concatenate_union_sets_from_table.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
fi


#
