#!/bin/bash

################################################################################
# Author: T. Cameron Waller, Ph.D.
# Date, initialization: 19 November 2025
# Date, review or revision: 19 November 2025
################################################################################
# Note

# The purpose of this procedural script is to access from the Molecular
# Signatures Database (MSigDB) a collection of gene sets. The script reads and
# interprets a table of parameters that specifies which individual sets to
# access.

# Recent example of usage:

###############################################################################
# Organize arguments.

path_directory_source=${1} # full path to source file in text format from which to read identifiers or names of genes
path_directory_product=${2} # full path to product file to which to write in text format the identifiers or names of genes
path_file_source_table=${3} # full path to source file in text format with tab and newline delimiters
prefix_accession=${4} # prefix of path and query for accession
suffix_accession=${5} # suffix of path and query for accession
prefix_file_product=${6} # prefix for name of file product
suffix_file_product=${7} # suffix for name of file product
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
if [ $prefix_file_product == "none" ]; then
  prefix_file_product=""
fi
if [ $suffix_file_product == "none" ]; then
  suffix_file_product=""
fi

#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error

################################################################################
# Execute procedure.

##########
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
# Read table of parameters and collect unique names for sets of genes.

# Initialize array for collection.
names_sets=()

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
    names_sets+=("${raw_name_set}")
  fi
done < "${input}"

# Collect unique names of groups.
IFS=" " read -r -a names_sets_unique <<< "$(echo "${names_sets[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')"

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "unique names of sets for accession"
  echo "----------"
fi
for name_set in "${names_sets_unique[@]}"; do
  # Report.
  if [[ "$report" == "true" ]]; then
    echo $name_set
  fi
done
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo -e "\n"
  echo -e "\n"
  echo -e "\n"
fi



# Iterate on names of sets and manage their accession and write to file.
for name_set in "${names_sets_unique[@]}"; do

  # Determine full path and query for accession.
  path_query_set="${prefix_accession}${name_set}${suffix_accession}"

  # Determine name and path of file to which to write information.
  name_file_product="${prefix_file_product}${name_set}${suffix_file_product}"
  path_file_product="${path_directory_product}/${name_file_product}"

  # Initialize file.
  rm $path_file_product # caution

  # Access information and write to file.
  # wget path_query --output-document path_and_or_name_for_file # specify full path and name for file
  # wget path_query --directory-prefix # specify only directory and preserve original name of file
  # Accession from the Molecular Signatures Database (MSigDB) operates as a
  # query, and it does not return the name of the file as works normally with
  # wget.
  wget "${path_query_set}" --output-document "${path_file_product}"
done



##########
# Report.
if [[ "$report" == "true" ]]; then
  # Organize information.
  # Read and collect names of files from source directory.
  #cd $path_directory_source
  # Bash version 4.4 introduced the "-d" option for "readarray".
  #readarray -d "" -t paths_files_source < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.txt.gz" -print0)
  paths_files_product=()
  while IFS= read -r -d $'\0'; do
    paths_files_product+=("$REPLY")
  done < <(find $path_directory_product -maxdepth 1 -mindepth 1 -type f -name "*${suffix_file_product}" -print0)
  count_paths_files_product=${#paths_files_product[@]}
  # Print information.
  echo "----------"
  echo "accession complete"
  echo "----------"
  echo -e "\n"
  echo "source directory:"
  echo $path_directory_source
  echo "product directory:"
  echo $path_directory_product
  echo "table of parameters for accession:"
  echo $path_file_source_table
  echo "----------"
  echo "count of files: " $count_paths_files_product
  echo "first file: " ${paths_files_product[0]} # notice base-zero indexing
  echo "last file: " ${paths_files_product[$count_paths_files_product - 1]}
  echo "----------"
fi




################################################################################
# Report.
if [ "$report" == "true" ]; then
  echo -e "\n"
  echo -e "\n"
  echo -e "\n"
  echo "----------"
  echo "script: access_msigdb_gene_sets_from_table.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
fi


#
