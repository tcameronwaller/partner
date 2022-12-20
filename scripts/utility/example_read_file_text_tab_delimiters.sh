#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: 19 December 2022
################################################################################
# Note

# Delimiter of each line must be either Tab or a single Space.

################################################################################
# Organize arguments.

path_file_source=${1} # full path to source file in text format
report=${2} # whether to print reports

# Read lines from file split fields by space, tab, or new-line delimiters.
if true; then
  input=$path_file_source
  while IFS=$' \t\n' read -r -a array
  do
    # Report.
    if [[ "$report" == "true" ]]; then
      echo "----------"
      echo "field 0: ${array[0]}"
      echo "field 1: ${array[1]}"
      echo "----------"
    fi
  done < "${input}"
fi

# First read lines from file and then second split fields from each line separately.
if false; then
  input=$path_file_source
  while IFS="" read -r line
  do
    IFS=$' \t\n' read -r -a array <<< "${line}"
    #field_0="${array[0]}"
    #field_1="${array[1]}"
    # Report.
    if [[ "$report" == "true" ]]; then
      echo "----------"
      echo "line:"
      echo "${line}"
      echo "field 0: ${array[0]}"
      echo "field 1: ${array[1]}"
      echo "----------"
    fi
  done < "${input}"
fi



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "test.sh"
  echo "----------"
fi



#
