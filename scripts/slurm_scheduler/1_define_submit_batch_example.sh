#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 4 April 2023
# Date, last execution: 4 April 2023
# Review: TCW; 4 April 2023
################################################################################
# Note

# This script initiates a sequence of scripts that execute a batch array of
# example jobs on the SLURM grid scheduler.
# documentation: "https://slurm.schedmd.com/quickstart.html"
# documentation: "https://wiki.hpc.rug.nl/peregrine/advanced_job_management/passing_parameters_to_a_job_script"

################################################################################


# TODO: TCW; 4 April 2023
# TODO: Confirm that the indices are correct and do not miss any instances from the array.
# TODO: I'm having difficulty pass the arguments to the run script...
# See reference "https://wiki.hpc.rug.nl/peregrine/advanced_job_management/passing_parameters_to_a_job_script"


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
path_script_run_batch_job="${path_directory_promiscuity}/scripts/slurm_scheduler/2_run_batch_job_example.sh"
path_script_execute_procedure="${path_directory_promiscuity}/scripts/slurm_scheduler/3_execute_example.sh"

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
# Indices in array of batch jobs start at zero.

if true; then
  sbatch --array 0-${batch_instances_count}:1 --chdir $path_directory_product \
  $path_script_run_batch_job \
  $path_file_batch_instances \
  $batch_instances_count \
  $message_common \
  $path_directory_product \
  $path_script_execute_procedure \
  $threads \
  $report
fi

if false; then
  sbatch --array 0-${batch_instances_count}:1 --chdir $path_directory_product \
  --export=path_file_batch_instances=$path_file_batch_instances \
  --export=batch_instances_count=$batch_instances_count \
  --export=message_common=$message_common \
  --export=path_directory_product=$path_directory_product \
  --export=path_script_execute_procedure=$path_script_execute_procedure \
  --export=threads=$threads \
  --export=report=$report \
  $path_script_run_batch_job
fi


#
