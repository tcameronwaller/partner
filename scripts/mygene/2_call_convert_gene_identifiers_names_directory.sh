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

# dock/out_age_exercise/transcriptomics/operate_sets/lists
#path_directory_source="${path_directory_dock}/in_data/mygene_2025-02-14/sets_gene_msigdb_2025-02-14_clean_symbol"
#path_directory_product="${path_directory_dock}/in_data/mygene_2025-02-14/sets_gene_msigdb_2025-02-14_clean_ensembl"

path_directory_source="${path_directory_dock}/temp_mygene_2025-07-17/input"
path_directory_product="${path_directory_dock}/temp_mygene_2025-07-17/output"
path_directory_temporary="${path_directory_product}/temporary_process" # hopefully unique

# File suffix.
suffix_file_source=".txt"
suffix_file_product=".txt"

# Scripts.
path_file_script="${path_directory_scripts}/partner/mygene/convert_gene_identifiers_names.sh"

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
type_source="ensembl.gene" # "entrezgene", "ensembl.gene", "symbol",
type_product="symbol" # "entrezgene", "ensembl.gene", "symbol",
species="human" # "human"
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
  # Translate identifiers or names.
  /usr/bin/bash $path_file_script \
  $path_file_source \
  $path_file_product \
  $delimiter_source \
  $delimiter_product \
  $type_source \
  $type_product \
  $species

done


################################################################################
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "script: call_convert_gene_identifiers_names_directory.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
fi

##########
# Remove temporary, intermediate files.
rm -r $path_directory_temporary


#
