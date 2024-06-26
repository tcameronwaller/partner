#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 27 March 2023
# Date, last execution or modification: 26 June 2024
# Review: TCW; 26 June 2024
################################################################################
# Note


################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_aai.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_source="${path_directory_dock}/in_data/gwas_priority_collection_format_team_2022-08-10"
path_directory_product="${path_directory_dock}/plots_qq_gwas"

# Files.
name_file_source_prefix="tcw_ukb_" # must not be empty string
name_file_source_suffix=".txt.gz" # must not be empty string
name_file_source_not=".place_holder_blarg" # exclude any files that include this character string in file name

# Scripts.
path_file_script_create_plots="${path_directory_process}/partner/scripts/gwas_clean/create_gwas_qq_plots.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product

################################################################################
# Organize parameters.

# Report.
report="true"

################################################################################
# Execute procedure.

# Collect and standardize polygenic scores.

/usr/bin/bash $path_file_script_create_plots \
$path_directory_source \
$name_file_source_prefix \
$name_file_source_suffix \
$name_file_source_not \
$path_directory_product


################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "create_qq_plots.sh"
  echo "----------"
fi



#
