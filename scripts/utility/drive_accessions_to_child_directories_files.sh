#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 20 December 2022
# Date, last execution: 18 September 2023
# Review: TCW; 18 September 2023
################################################################################
# Note

# Format of parameter table:
# - File in text format
# - New-line characters as delimiters between lines for each row
# - Tab characters as delimiters between fields for columns within each row
# - column 0: "inclusion"; Logical binary indicator of whether to execute procedure for current row's record
# - - The text character string must be "1" to execute procedure for current row's record.
# - - This conditional statement also excludes the first header row.
# - column 1: "directory"; Name of child directory to create within parent directory and within which to store the file or files
# - column 2: "phenotype"; Dependent variable, only for readability and reference
# - column 3: "path_accession"; Full path on internet (Uniform Resource Locator) to the file for download
# - column 4: "zip"; Logical binary indicator of whether after accession to decompress the file from the Zip format, using command "unzip".
# - - The text character string must be "1" to enable this decompression.
# - There can be additional columns for reference, but the program will not read these.

################################################################################
# Organize arguments.

path_file_accession=${1} # full path to file of parameter table in text format with information about accessions
path_directory_parent=${2} # full path to parent directory within which to create child directories and save files
report=${3} # whether to print reports

################################################################################
# Procedure.


# Read lines from file and split fields within each line by space, tab, or new-line delimiters.
input=$path_file_accession
while IFS=$' \t\n' read -r -a array
do
  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo "field 0, inclusion: ${array[0]}"
    echo "field 1, directory: ${array[1]}"
    echo "field 2, description: ${array[2]}"
    echo "field 3, path_accession: ${array[3]}"
    echo "field 4, zip: ${array[4]}"
    echo "field 5, Tarball Gzip archive: ${array[5]}"
    echo "----------"
  fi
  # Execute procedure for current record's parameters.
  if [[ "${array[0]}" == "1" ]]; then
    # Create child directory if it does not already exist.
    path_directory_child="${path_directory_parent}/${array[1]}"
    mkdir -p $path_directory_child
    #cd $path_directory_child
    # Access the specific file and save within the child directory.
    wget "${array[3]}" --directory-prefix "${path_directory_child}" --content-disposition --no-check-certificate --show-progress
    # Decompress content from Zip archive format.
    if [[ "${array[4]}" == "1" ]]; then
      # The argument flag "--extract-dir" (synonym for "-d") does not seem to work properly.
      unzip "${path_directory_child}/*.zip" -d "${path_directory_child}"
    fi
    # Decompress content from Tarball GZip archive format.
    if [[ "${array[5]}" == "1" ]]; then
      tar -xzvf "${path_directory_child}/*.tar.gz" --directory $path_directory_child
      tar -xzvf "${path_directory_child/*.tgz}" --directory $path_directory_child
    fi
  fi
done < "${input}"

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "drive_accessions_to_child_directories_files.sh"
  echo "----------"
fi



#
