#!/bin/bash

################################################################################
# Notes:

# Plan: TCW; 7 February 2023
# Test the pipeline on the NCSA computational cluster.
# Use 16 threads for job on GWAS summary statistics from each study.
# Use incremental memory allocations to find the optimum.



################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_tools=$(<"./waller_tools.txt")

path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_product="${path_directory_dock}/test_pipe_gwas_clean"

# Files.

#identifier_gwas="30367059_teumer_2018_hyperthyroidism"
#identifier_gwas="30367059_teumer_2018_hypothyroidism"
#identifier_gwas="30367059_teumer_2018_tsh_female"
#identifier_gwas="30367059_teumer_2018_tsh_male"
#identifier_gwas="32769997_zhou_2020_tsh"
#identifier_gwas="32042192_ruth_2020_testosterone_female"
identifier_gwas="34255042_schmitz_2021_estradiol_female"
path_file_gwas_source="${path_directory_dock}/hormone_genetics/gwas_format_standard/${identifier_gwas}.txt.gz"
path_file_gwas_product="${path_directory_product}/${identifier_gwas}.txt.gz"


# Scripts.
path_file_script_pipe_gwas_clean="${path_directory_process}/promiscuity/scripts/gwas_clean/pipe_gwas_clean.sh"

################################################################################
# Organize parameters.

report="true"
threads=16

################################################################################
# Execute procedure.


# Call pipe script clean procedure on GWAS summary statistics.
/usr/bin/bash $path_file_script_pipe_gwas_clean \
$path_file_gwas_source \
$path_file_gwas_product \
$threads \
$report



#
