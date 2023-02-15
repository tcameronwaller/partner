#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 22 December 2022
################################################################################
# Note

# Format of parameter table:

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
  # Report.
  if [[ "$report" == "true" ]]; then
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
    echo "field 15, script: ${array[15]}"
    echo "field 16, note: ${array[16]}"
    echo "----------"
  fi
  # Execute procedure for current record's parameters.
  if [[ "${array[0]}" == "1" ]]; then
    # Organize paths.
    path_directory_child_source="${path_directory_parent_source}/${array[1]}"
    name_file_source="${array[5]}"
    suffix_file_source="${array[6]}"
    name_base_file_source=$(echo $name_file_source | sed "s/$suffix_file_source//")
    path_file_source="${path_directory_child_source}/${name_file_source}"
    path_directory_temporary="${path_directory_product}/temporary_tcw_73216"
    path_file_source_standard="${path_directory_temporary}/${name_base_file_source}.txt.gz"
    path_file_script="${path_directory_script}/${array[15]}"
    path_file_product="${path_directory_product}/${array[2]}.txt.gz"
    mkdir -p $path_directory_temporary
    # Manage compression formats and file suffices.
    if [ "${array[7]}" == "1" ] && [ "${array[8]}" == "1" ]; then
      # 1. Decompress from BGZip format (http://www.htslib.org/doc/bgzip.html).
      $path_bgzip --decompress "$path_file_source" --stdout > "${path_directory_temporary}/${name_base_file_source}.txt"
      # 2. Compress to GZip format.
      gzip -cvf "${path_directory_temporary}/${name_base_file_source}.txt" > $path_file_source_standard
    fi
    if [ "${array[7]}" != "1" ] && [ "${array[8]}" == "1" ]; then
      # 1. Compress to Gzip format
      gzip -cvf "$path_file_source" > $path_file_source_standard
    fi
    if [ "${array[7]}" != "1" ] && [ "${array[8]}" != "1" ]; then
      # Manage suffix.
      cp "$path_file_source" $path_file_source_standard
    fi
    # Call script for translation.
    /usr/bin/bash $path_file_script \
    $path_file_source_standard \
    $path_file_product \
    ${array[10]} \
    ${array[11]} \
    ${array[12]} \
    ${array[13]} \
    ${array[14]} \
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
  echo "drive_translations_gwas_to_standard_format.sh"
  echo "----------"
fi



#
