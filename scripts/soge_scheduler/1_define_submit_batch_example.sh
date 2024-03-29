#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 4 April 2023
# Date, last execution: 4 April 2023
# Review: TCW; 4 April 2023
################################################################################
# Note

# This script initiates a sequence of scripts that execute a batch array of
# example jobs on the Oracle Sun Grid Engine grid scheduler.
# documentation: "https://docs.oracle.com/cd/E19279-01/820-3257-12/n1ge.html"
# documentation: "http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html"

################################################################################



################################################################################
# Organize arguments.

path_directory_product=${1} # full path to parent directory for product files
name_file_prefix=${2} # prefix of name of file
name_file_suffix=${3} # suffix of name of file
message_common=${4} # message parameter common to all jobs
threads=${5} # count of concurrent or parallel process threads on node cores
report=${6} # whether to print reports

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_promiscuity="${path_directory_process}/promiscuity"

# Files.
path_file_batch_instances="${path_directory_product}/batch_instances.txt"
path_file_batch_out="${path_directory_product}/batch_out.txt"
path_file_batch_error="${path_directory_product}/batch_error.txt"

# Scripts.
path_script_run_batch_job="${path_directory_promiscuity}/scripts/soge_scheduler/2_run_batch_job_example.sh"
path_script_execute_procedure="${path_directory_promiscuity}/scripts/soge_scheduler/3_execute_example.sh"

# Initialize directories.
mkdir -p $path_directory_product

# Initialize files.
rm $path_file_batch_instances
rm $path_file_batch_out
rm $path_file_batch_error

################################################################################
# Organize batch job instances.

# Iterate on autosomal chromosomes.
chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
for chromosome in "${chromosomes[@]}"; do
  # Define parameters in array instance for batch job.
  instance="${chromosome};${name_file_prefix};${name_file_suffix}"
  echo $instance >> $path_file_batch_instances
done

# Read batch instances.
readarray -t batch_instances < $path_file_batch_instances
batch_instances_count=${#batch_instances[@]}



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo "1_define_submit_batch_example.sh"
  echo "----------"
  echo "count of batch instances: " $batch_instances_count
  echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
  echo "last batch instance: " ${batch_instances[$batch_instances_count - 1]}
  echo "----------"
fi



################################################################################
# Submit batch of jobs to grid cluster scheduler for processing.
# Submit to Oracle Sun Grid Engine.
# Indices in array of batch jobs start at one, not zero.

if true; then
  qsub -t 1-${batch_instances_count}:1 \
  -o $path_file_batch_out \
  -e $path_file_batch_error \
  $path_script_run_batch_job \
  $path_file_batch_instances \
  $batch_instances_count \
  $message_common \
  $path_directory_product \
  $path_script_execute_procedure \
  $threads \
  $report
fi



#
