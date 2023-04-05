#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 4 April 2023
# Date, last execution: 4 April 2023
# Review: TCW; 4 April 2023
################################################################################
# Note

# This script executes a very simple procedure to demonstrate management of
# batch jobs on the Oracle Sun Grid Engine scheduler.


################################################################################



################################################################################
# Organize arguments.

chromosome=${1} # identifier of a single chromosome
name_file_prefix=${2} # prefix of name of file
name_file_suffix=${3} # suffix of name of file
message_common=${4} # message parameter common to all jobs
path_directory_product=${5} # full path to parent directory for product files
threads=${6} # count of concurrent or parallel process threads on node cores
report=${7} # whether to print reports

################################################################################
# Execute procedure.

path_file_product="${path_directory_product}/${name_file_prefix}${chromosome}${name_file_suffix}"
echo "$message_common" > $path_file_product

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "3_execute_example.sh"
  echo "----------"
fi



#
