#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 20 December 2022
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

path_directory_parent=${1} # full path to parent directory within which to create child directories and save files
path_file_accession=${2} # full path to file of parameter table in text format with information about accessions
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
    echo "field 2, phenotype: ${array[2]}"
    echo "field 3, path_accession: ${array[3]}"
    echo "field 4, zip: ${array[4]}"
    echo "----------"
  fi
  # Execute procedure for current record's parameters.
  if [[ "${array[0]}" == "1" ]]; then
    # Create child directory if it does not already exist.
    path_directory_child="${path_directory_parent}/${array[1]}"
    mkdir -p $path_directory_child
    #cd $path_directory_child
    # Access the specific file and save within the child directory.
    wget "${array[3]}" --directory-prefix "${path_directory_child}" --content-disposition --show-progress
    # Decompress file from Zip format.
    if [[ "${array[4]}" == "1" ]]; then
      # The argument flag "--extract-dir" (synonym for "-d") does not seem to work properly.
      unzip "${path_directory_child}/*.zip" -d "${path_directory_child}"
    fi
  fi
done < "${input}"

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "drive_access_child_directories_files.sh"
  echo "----------"
fi



#
