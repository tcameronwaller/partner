#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 31 October 2024
# Date, last execution or modification: 31 October 2024
# Review: TCW; 31 Octoboer 2024
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

# Files.
path_file_source="${path_directory_dock}/in_data/source_msigdb_hp_insulin_resistance.txt"
path_file_product="${path_directory_dock}/in_data/msigdb_hp_insulin_resistance.txt"

# Scripts.
path_file_script="${path_directory_scripts}/partner/mygene/convert_gene_identifiers_names.sh"

# Initialize directories.

###############################################################################
# Organize parameters.

# Parameters.
delimiter_source="," # "newline", "tab", "\n", "\t", ";", ":", ",", not " "
delimiter_product="\n" # "newline", "tab", "\n", "\t", ";", ":", ",", not " "
type_source="entrezgene" # "entrezgene", "ensembl.gene", "symbol",
type_product="ensembl.gene" # "entrezgene", "ensembl.gene", "symbol",
species="human" # "human"
report="true"
#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error

################################################################################
# Execute procedure.

/usr/bin/bash $path_file_script \
$path_file_source \
$path_file_product \
$delimiter_source \
$delimiter_product \
$type_source \
$type_product \
$species

################################################################################
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "script: example_call_convert_gene_identifiers_names.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
fi


#
