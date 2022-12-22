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
    echo "field 9, count: ${array[9]}"
    echo "field 10, script: ${array[10]}"
    echo "----------"
  fi
  # Execute procedure for current record's parameters.
  if [[ "${array[0]}" == "1" ]]; then
    # Organize paths.
    path_directory_child_source="${path_directory_parent_source}/${array[1]}"
    path_file_script="${path_directory_script}/${array[10]}"
    path_file_product="${path_directory_product}/${array[2]}.txt.gz"
    # Organize name and suffix of source file.
    name_file_source="${array[5]}"
    suffix_file_source="${array[6]}"
    name_base_file_source=$(echo $name_file_source | sed "s/$suffix_file_source//")
    path_file_source_standard="${path_directory_child_source}/${name_base_file_source}_temporary_tcw_73216.txt.gz"
    #name_suffix_file_progress="${name_file_source}"
    #path_file_source="${path_directory_child_source}/${array[5]}"
    #name_base_file_source="$(basename $path_file_source "${array[6]}")"

    # Manage compression formats and file suffices.
    if [ "${array[7]}" == "1" ] && [ "${array[8]}" == "1" ]; then
      # 1. Decompress from BGZip format (http://www.htslib.org/doc/bgzip.html).
      $path_bgzip --decompress "${path_directory_child_source}/${name_file_source}" --stdout > "${path_directory_child_source}/${name_base_file_source}_temporary_tcw_73216.txt"
      # 2. Compress to GZip format.
      gzip -cvf "${path_directory_child_source}/${name_base_file_source}_temporary_tcw_73216.txt" > $path_file_source_standard
      # Collect paths to temporary files for removal after completion of the remainder of procedure.
      removals=()
      removals+=("${path_directory_child_source}/${name_base_file_source}_temporary_tcw_73216.txt")
      removals+=($path_file_source_standard)
    fi
    if [ "${array[7]}" != "1" ] && [ "${array[8]}" == "1" ]; then
      # 1. Compress to Gzip format
      gzip -cvf "${path_directory_child_source}/${name_file_source}" > $path_file_source_standard
      # Collect paths to temporary files for removal after completion of the remainder of procedure.
      removals=()
      removals+=($path_file_source_standard)
    fi
    if [ "${array[7]}" != "1" ] && [ "${array[8]}" != "1" ]; then
      # Manage suffix.
      cp "${path_directory_child_source}/${name_file_source}" $path_file_source_standard
      # Collect paths to temporary files for removal after completion of the remainder of procedure.
      removals=()
      removals+=($path_file_source_standard)
    fi
    # Call script for translation.
    if [[ "${array[9]}" == "NA" ]]; then
      /usr/bin/bash "${path_directory_script}/${array[10]}" \
      $path_file_source_standard \
      $path_file_product \
      $report
    elif [[ "${array[9]}" != "NA" ]]; then
      /usr/bin/bash "${path_directory_script}/${array[10]}" \
      $path_file_source_standard \
      $path_file_product \
      ${array[9]} \
      $report
    fi

    # Remove temporary files.
    for path_file_temporary in "${removals[@]}"; do
      #echo "remove this file: "
      #echo $path_file_temporary
      rm $path_file_temporary
    done
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
