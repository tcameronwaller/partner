#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date: ___ 2023
################################################################################
# Note

# TODO: TCW; 8 February 2023
# TODO: This is an EARLY draft of this script... not yet functional.

################################################################################
# Organize arguments.

path_file_translation=${1} # full path to file of parameter table in text format with information about translations
path_directory_source=${2} # full path to directory within which to find files of source GWAS summary statistics
path_directory_product=${3} # full path to directory within which to save files of product GWAS summary statistics
report=${4} # whether to print reports

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes

# Files.
path_file_batch_instances="${path_directory_product}/batch_instances.txt"

# Scripts.
path_file_script_run_batch="${path_directory_process}/promiscuity/scripts/gwas_clean/2_run_batch_pipe_gwas_clean.sh"
path_file_script_pipe_gwas_clean="${path_directory_process}/promiscuity/scripts/gwas_clean/pipe_gwas_clean.sh"

# Initialize files.
rm $path_file_batch_instances

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product

################################################################################
# Organize parameters.

threads=1
#report="true"

################################################################################
# Execute procedure.

# 5. pass the "type" variable on to the "run_batch" script.


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
    # Define variables.
    name="${array[2]}"
    type="${array[9]}"
    count_cases="${array[13]}"
    # Organize paths.
    path_file_gwas_source="${path_directory_source}/${name}.txt.gz"
    path_file_gwas_product="${path_directory_product}/${name}.txt.gz"
    # Define and append a new batch instance.
    instance="${path_file_gwas_source};${path_file_gwas_product};${type};${count_cases}"
    echo $instance >> $path_file_batch_instances
  fi
done < "${input}"



################################################################################
# Submit batch instances to cluster scheduler.

# Read batch instances.
readarray -t batch_instances < $path_file_batch_instances
batch_instances_count=${#batch_instances[@]}
echo "----------"
echo "count of batch instances: " $batch_instances_count
echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
echo "last batch instance: " ${batch_instances[$batch_instances_count - 1]}

# Execute batch with grid scheduler.
if false; then
  # Submit array batch to Sun Grid Engine.
  # Array batch indices must start at one (not zero).
  qsub -t 1-${batch_instances_count}:1 \
  -o "${path_directory_product}/batch_out.txt" \
  -e "${path_directory_product}/batch_error.txt" \
  $path_file_script_run_batch \
  $path_file_batch_instances \
  $batch_instances_count \
  $path_file_script_pipe_gwas_clean \
  $threads \
  $report
fi



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "1_submit_batch_pipe_gwas_clean.sh"
  echo "----------"
fi



#
