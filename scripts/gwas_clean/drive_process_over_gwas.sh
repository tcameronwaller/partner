#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 15 November 2023
# Date, last execution: 13 December 2023
# Review: TCW; 13 December 2023
################################################################################
# Note

################################################################################


################################################################################
# Organize arguments.

path_file_table_parameter=${1} # full path to file of parameter table in text format with information about translations
path_directory_source=${2} # full path to source directory within which to find files of source GWAS summary statistics
path_directory_product=${3} # full path to product directory within which to save files of product GWAS summary statistics
path_script_process=${4} # full path to file of script to apply to each file of GWAS summary statistics in the parameter table
report=${5} # whether to print reports

################################################################################
# Procedure.


##########
# Read lines from file and split fields within each line by space, tab, or new-line delimiters.
input=$path_file_table_parameter
while IFS=$' \t\n' read -r -a array
do

  # Extract values from individual columns within table's row.
  # These parameters match the format of the table as of 13 December 2023.
  raw_availability="${array[0]}"
  raw_inclusion="${array[1]}"
  raw_directory="${array[2]}"
  raw_name_study="${array[3]}"
  raw_phenotype="${array[4]}"
  raw_sex="${array[5]}"
  raw_name_file_source="${array[6]}"
  raw_suffix_file_source="${array[7]}"
  raw_bgzip="${array[8]}"
  raw_gzip="${array[9]}"
  raw_type="${array[10]}"
  raw_fill_observations="${array[11]}"
  raw_observations_total="${array[12]}"
  raw_fill_case_control="${array[13]}"
  raw_cases="${array[14]}"
  raw_controls="${array[15]}"
  raw_observations_effective="${array[16]}"
  raw_prevalence_sample="${array[17]}"
  raw_prevalence_population="${array[18]}"
  raw_script="${array[19]}"
  raw_note="${array[20]}"

  # Report.
  if [ $raw_inclusion == "1" ] && [ "$report" == "true" ]; then
    echo "----------"
    echo "field 0, availability: ${raw_availability}"
    echo "field 1, inclusion: ${raw_inclusion}"
    echo "field 2, directory: ${raw_directory}"
    echo "field 3, name_study: ${raw_name_study}"
    echo "field 4, phenotype: ${raw_phenotype}"
    echo "field 5, sex: ${raw_sex}"
    echo "field 6, file: ${raw_name_file_source}"
    echo "field 7, suffix: ${raw_suffix_file_source}"
    echo "field 8, bgzip: ${raw_bgzip}"
    echo "field 9, gzip: ${raw_gzip}"
    echo "field 10, type: ${raw_type}"
    echo "field 11, fill_observations: ${raw_fill_observations}"
    echo "field 12, observations_total: ${raw_observations_total}"
    echo "field 13, fill_case_control: ${raw_fill_case_control}"
    echo "field 14, cases: ${raw_cases}"
    echo "field 15, controls: ${raw_controls}"
    echo "field 16, observations_effective: ${raw_observations_effective}"
    echo "field 17, prevalence_sample: ${raw_prevalence_sample}"
    echo "field 18, prevalence_population: ${raw_prevalence_population}"
    echo "field 19, script: ${raw_script}"
    echo "field 20, note: ${raw_note}"
    echo "----------"
  fi
  # Execute procedure for current record's parameters.
  if [ $raw_inclusion == "1" ]; then
    # Organize paths.
    #path_directory_temporary="${path_directory_product}/temporary_${name_base_file_source}_${raw_name_study}" # hopefully unique
    path_file_gwas_source="${path_directory_source}/${raw_name_study}.txt.gz"
    path_file_gwas_product="${path_directory_product}/${raw_name_study}.txt.gz"
    # Call script for translation.
    /usr/bin/bash $path_script_process \
    $path_file_gwas_source \
    $path_file_gwas_product \
    $report
    # Remove temporary directory and files.
    #rm -r $path_directory_temporary
  fi
done < "${input}"

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "drive_process_over_gwas.sh"
  echo "----------"
fi



#
