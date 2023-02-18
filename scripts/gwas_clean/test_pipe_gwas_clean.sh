#!/bin/bash

################################################################################
# Notes:

# Plan: TCW; 7 February 2023
# Test the pipeline on the NCSA computational cluster.
# Use 16 threads for job on GWAS summary statistics from each study.
# Use incremental memory allocations to find the optimum.



################################################################################

# Strategy below is for the driving script that submits the batch.

# Find element in array...
# Look-up parameters for each GWAS within the parameter table
# Find the correct record in the array by matching the study identifier

#array=("a b" "c d")

#for ((i=0; i<${#array[@]}; i++)); do
#  if [[ ${array[$i]} == "a b" ]]; then
#    echo "Element $i matched"
#  fi
#done

#Output:

#Element 0 matched

#${#array[@]} contains number of last element in array.




################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_tools=$(<"./waller_tools.txt")

path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_product="${path_directory_dock}/test_pipe_gwas_clean"
path_directory_temporary="${path_directory_product}/temporary_tcw_test_20230208" # hopefully unique

# Files.

#identifier_gwas="BMI_GIANTUKB_EUR"
#identifier_gwas="BMI_GIANT_EUR"
#identifier_gwas="BD_PGC3_EUR"
#path_file_gwas_source="${path_directory_dock}/hormone_genetics/gwas_from_team_collection/${identifier_gwas}.txt.gz"

identifier_gwas="32581359_saevarsdottir_2020_aitd" # logistic

#identifier_gwas="30367059_teumer_2018_hyperthyroidism"
#identifier_gwas="30367059_teumer_2018_hypothyroidism" # <-- test logistic
#identifier_gwas="30367059_teumer_2018_tsh_female" # <-- test linear
#identifier_gwas="30367059_teumer_2018_tsh_male"
#identifier_gwas="32769997_zhou_2020_tsh"
#identifier_gwas="32042192_ruth_2020_testosterone_female"
#identifier_gwas="32042192_ruth_2020_shbg_female"
#identifier_gwas="34255042_schmitz_2021_estradiol_female"
path_file_gwas_source="${path_directory_dock}/hormone_genetics_tcw_2023-02-17/gwas_format_standard/${identifier_gwas}.txt.gz"
path_file_gwas_product="${path_directory_product}/${identifier_gwas}.txt.gz"

# Temporary files.
path_ftemp_gwas_source_decomp="${path_directory_temporary}/${identifier_gwas}_source.txt"
path_ftemp_gwas_source="${path_directory_temporary}/${identifier_gwas}_source.txt.gz"

# Scripts.
path_file_script_pipe_gwas_clean="${path_directory_process}/promiscuity/scripts/gwas_clean/pipe_gwas_clean.sh"

# Initialize directories.
rm -r $path_directory_product
rm -r $path_directory_temporary
mkdir -p $path_directory_product
mkdir -p $path_directory_temporary
cd $path_directory_product


################################################################################
# Organize parameters.

#type="linear"
type="logistic"
count_cases=532
threads=1
report="true"

################################################################################
# Execute procedure.

# Keep only the first 10,000 records.
zcat $path_file_gwas_source | awk 'NR < 10002 {
  print $0
}' >> $path_ftemp_gwas_source_decomp

# Compress file.
gzip -cvf $path_ftemp_gwas_source_decomp > $path_ftemp_gwas_source

# Call pipe script clean procedure on GWAS summary statistics.
/usr/bin/bash $path_file_script_pipe_gwas_clean \
$path_ftemp_gwas_source \
$path_file_gwas_product \
$type \
$count_cases \
$threads \
$report



##########
# Remove temporary directories and files.
# Suppress this block for debugging.
if true; then
  rm -r $path_directory_temporary
fi



#
