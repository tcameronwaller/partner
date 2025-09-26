#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 12 February 2025
# Date, last execution or modification: 24 September 2025
# Review: TCW; 24 September 2025
# Review: TCW; 11 September 2025
# Review: TCW; 12 February 2025
################################################################################
# Note

# Recent example of usage:
# /.../pails_process/omega3/2025-09-22_heterogeneity_candidate_adipose_fibrosis

###############################################################################
# Organize arguments.

path_directory_source=${1} # full path to source file in text format from which to read identifiers or names of genes
path_directory_product=${2} # full path to product file to which to write in text format the identifiers or names of genes
suffix_file_source=${3} # text suffix in names of source files
suffix_file_product_list=${4} # text suffix in names of product files
suffix_file_product_table=${5} # text suffix in names of product files
delimiter_source=${6} # text delimiter between items in source file, not white space
delimiter_product=${7} # text delimiter between items in product file, not white space
type_source=${8} # type of name or identifier for genes
type_product=${9} # type of name or identifier for genes
species=${10} # species reference for genes
report=${11} # whether to print reports to terminal

################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_tools=$(<"$path_directory_paths/path_directory_tools.txt")
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_scripts="$path_directory_process/scripts"
path_directory_package="$path_directory_process/package"

path_directory_product_lists="${path_directory_product}/lists"
path_directory_product_tables="${path_directory_product}/tables"

stamp_date=$(date +%Y-%m-%d)
path_directory_temporary="${path_directory_product}/temporary_${stamp_date}" # hopefully unique


# Files.

# Scripts.
path_file_script="${path_directory_scripts}/partner/mygene/convert_gene_identifiers_names.sh"

# Initialize directories.
#rm -r $path_directory_temporary # caution
#rm -r $path_directory_product # caution
mkdir -p $path_directory_product
mkdir -p $path_directory_product_lists
mkdir -p $path_directory_product_tables
mkdir -p $path_directory_temporary
cd $path_directory_product

###############################################################################
# Organize parameters.

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
  name_file_product_list="${name_base_file_product}${suffix_file_product_list}"
  name_file_product_table="${name_base_file_product}${suffix_file_product_table}"
  path_file_product_list="${path_directory_product_lists}/${name_file_product_list}" # hopefully unique
  path_file_product_table="${path_directory_product_tables}/${name_file_product_table}" # hopefully unique
  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo "Path to source file:"
    echo $path_file_source
    echo "Path to product files:"
    echo $path_file_product_list
    echo $path_file_product_table
    echo "----------"
  fi
  ##########
  # Translate identifiers or names.
  /usr/bin/bash $path_file_script \
  $path_file_source \
  $path_file_product_list \
  $path_file_product_table \
  $delimiter_source \
  $delimiter_product \
  $type_source \
  $type_product \
  $species \
  $report

done

################################################################################
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "script: convert_gene_identifiers_names_directory_files.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
fi

##########
# Remove temporary, intermediate files.
rm -r $path_directory_temporary


#
