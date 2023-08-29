#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 22 December 2022
# Date, last execution: 2 August 2023
# Review: TCW; 2 August 2023
################################################################################
# Note

################################################################################



################################################################################
# Organize arguments.

path_file_translation=${1} # full path to file of parameter table in text format with information about translations
path_directory_parent_source=${2} # full path to parent directory within which to find child directories and files of source GWAS summary statistics
path_directory_script=${3} # full path to directory within which to find scripts for format translations
path_directory_product=${4} # full path to directory within which to save files of product GWAS summary statistics
path_bgzip=${5} # full path to installation executable file of BGZip
report=${6} # whether to print reports

################################################################################
# Procedure.


# Read lines from file and split fields within each line by space, tab, or new-line delimiters.
input=$path_file_translation
while IFS=$' \t\n' read -r -a array
do

  # Extract values from individual columns within table's row.
  raw_inclusion="${array[0]}"
  raw_directory="${array[1]}"
  raw_name_study="${array[2]}"
  raw_phenotype="${array[3]}"
  raw_sex="${array[4]}"
  raw_name_file_source="${array[5]}"
  raw_suffix_file_source="${array[6]}"
  raw_bgzip="${array[7]}"
  raw_gzip="${array[8]}"
  raw_type="${array[9]}"
  raw_fill_observations="${array[10]}"
  raw_observations="${array[11]}"
  raw_fill_case_control="${array[12]}"
  raw_cases="${array[13]}"
  raw_controls="${array[14]}"
  raw_prevalence_sample="${array[15]}"
  raw_prevalence_population="${array[16]}"
  raw_script="${array[17]}"
  raw_note="${array[18]}"

  # Report.
  if [ $raw_inclusion == "1" ] && [ "$report" == "true" ]; then
    echo "----------"
    echo "field 0, inclusion: ${array[0]}"
    echo "field 1, directory: ${array[1]}"
    echo "field 2, name: ${array[2]}"
    echo "field 3, phenotype: ${array[3]}"
    echo "field 4, sex: ${array[4]}"
    echo "field 5, file: ${array[5]}"
    echo "field 6, suffix: ${array[6]}"
    echo "field 7, bgzip: ${array[7]}"
    echo "field 8, gzip: ${array[8]}"
    echo "field 9, type: ${array[9]}"
    echo "field 10, fill_observations: ${array[10]}"
    echo "field 11, observations: ${array[11]}"
    echo "field 12, fill_case_control: ${array[12]}"
    echo "field 13, cases: ${array[13]}"
    echo "field 14, controls: ${array[14]}"
    echo "field 15, prevalence_sample: ${array[15]}"
    echo "field 16, prevalence_population: ${array[16]}"
    echo "field 17, script: ${array[17]}"
    echo "field 18, note: ${array[18]}"
    echo "----------"
  fi
  # Execute procedure for current record's parameters.
  if [[ $raw_inclusion == "1" ]]; then
    # Organize paths.
    path_directory_child_source="${path_directory_parent_source}/$raw_directory"
    name_base_file_source=$(echo $raw_name_file_source | sed "s/$raw_suffix_file_source//")
    path_file_source="${path_directory_child_source}/${raw_name_file_source}"
    path_directory_temporary="${path_directory_product}/temporary_${name_base_file_source}_${raw_name_study}" # hopefully unique
    path_file_source_standard="${path_directory_temporary}/${name_base_file_source}.txt.gz"
    path_file_script="${path_directory_script}/$raw_script"
    path_file_product="${path_directory_product}/$raw_name_study.txt.gz"
    mkdir -p $path_directory_temporary
    # Manage compression formats and file suffixes.
    if [ $raw_bgzip == "1" ] && [ $raw_gzip == "1" ]; then
      # 1. Decompress from BGZip format (http://www.htslib.org/doc/bgzip.html).
      $path_bgzip --decompress "$path_file_source" --stdout > "${path_directory_temporary}/${name_base_file_source}.txt"
      # 2. Compress to GZip format.
      gzip -cvf "${path_directory_temporary}/${name_base_file_source}.txt" > $path_file_source_standard
    fi
    if [ $raw_bgzip != "1" ] && [ $raw_gzip == "1" ]; then
      # 1. Compress to Gzip format
      gzip -cvf "$path_file_source" > $path_file_source_standard
    fi
    if [ $raw_bgzip != "1" ] && [ $raw_gzip != "1" ]; then
      # Manage suffix.
      cp "$path_file_source" $path_file_source_standard
    fi
    # Call script for translation.
    /usr/bin/bash $path_file_script \
    $path_file_source_standard \
    $path_file_product \
    $raw_fill_observations \
    $raw_observations \
    $raw_fill_case_control \
    $raw_cases \
    $raw_controls \
    $report
    # Remove temporary directory and files.
    rm -r $path_directory_temporary
  fi
done < "${input}"

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "drive_translations_gwas_to_standard_format.sh"
  echo "----------"
fi



#
