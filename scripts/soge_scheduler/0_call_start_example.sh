#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 4 April 2023
# Date, last execution: 4 April 2023
# Review: TCW; 4 April 2023
################################################################################
# Note


################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_promiscuity="${path_directory_process}/promiscuity"
path_directory_product="${path_directory_dock}/test_soge_scheduler"

# Scripts.
path_script_define_submit="${path_directory_promiscuity}/scripts/soge_scheduler/1_define_submit_batch_example.sh"

# Initialize directories.
mkdir -p $path_directory_product

################################################################################
# Organize parameters.

name_file_prefix="test_test_test_chr_"
name_file_suffix=".txt"
message_common="Hello_World!" # Any white space disrupts the handling of arguments.
threads=1
report="true"



################################################################################
# Call.

/usr/bin/bash $path_script_define_submit \
$path_directory_product \
$name_file_prefix \
$name_file_suffix \
$message_common \
$threads \
$report



#
