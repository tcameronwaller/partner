#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 6 August 2022
# Date, last execution: 30 November 2023
# Review: TCW; 30 November 2023
################################################################################
# Note


################################################################################



################################################################################
# Organize arguments.

path_file_batch_instances=${1} # full directory path and file name for batch instances in a text file
path_directory_batch=${2} # full path to directory for batch files
path_directory_product=${3} # full path to parent directory within which to write product directories and files
path_directory_disequilibrium=${4} # full directory path to directory for LDSC reference on linkage disequilibrium
path_directory_process=${5} # full path to directory for all processes relevant to current project
threads=${6} # count of processing threads to use
report=${7} # whether to print reports

################################################################################
# Organize paths.

# Directories.
cd ~/paths
#path_directory_process=$(<"./process_psychiatric_metabolism.txt")
#path_directory_dock="${path_directory_process}/dock"
path_directory_partner="${path_directory_process}/partner"
path_directory_batch_1="${path_directory_batch}/batch_1"
path_directory_batch_2="${path_directory_batch}/batch_2"
path_directory_batch_3="${path_directory_batch}/batch_3"
path_directory_batch_4="${path_directory_batch}/batch_4"

# Files.
#path_file_batch_instances="${path_directory_product}/batch_instances.txt"
#path_file_batch_out="${path_directory_product}/batch_out.txt"
#path_file_batch_error="${path_directory_product}/batch_error.txt"

# Scripts.
path_script_batch_2="${path_directory_partner}/scripts/ldsc/ldsc_correlation_batch_2.sh"
path_script_batch_3="${path_directory_partner}/scripts/ldsc/estimate_gwas_genetic_correlation_ldsc.sh"

# Initialize directories.
cd $path_directory_batch # execute batch from within this directory
mkdir -p $path_directory_batch_1
mkdir -p $path_directory_batch_2
mkdir -p $path_directory_batch_3
mkdir -p $path_directory_batch_4

################################################################################
# Organize batch job instances.

# Read batch instances.
readarray -t batch_instances < $path_file_batch_instances
batch_instances_count=${#batch_instances[@]}
index_array_maximum=$(($batch_instances_count - 1))

################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "ldsc_correlation_batch_1.sh"
  echo "----------"
  echo "count of batch instances: " $batch_instances_count
  echo "maximum array index: " $index_array_maximum
  echo "----------"
  echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
  echo "last batch instance: " ${batch_instances[$index_array_maximum]}
  echo "----------"
fi



################################################################################
# Submit batch of jobs to grid cluster scheduler for processing.
# Submit to Slurm Scheduler.
# Indices in array of batch jobs start at zero.

if false; then
  sbatch --array 0-${index_array_maximum}:1 --chdir $path_directory_batch_1 \
  $path_script_batch_2 \
  $path_file_batch_instances \
  $batch_instances_count \
  0 \
  $path_directory_product \
  $path_directory_disequilibrium \
  $threads \
  $report \
  $path_script_batch_3
fi


# NCSA implementation of SLURM has value of "MaxArraySize" of 4,000.
# Find the value of "MaxArraySize" within the text file at path "/etc/slurm/slurm.conf".
# SLURM will not allow an array index greater than "MaxArraySize".
# It is necessary to split larger jobs.

if true; then
  sbatch --array 0-3000:1 --chdir $path_directory_batch_1 \
  $path_script_batch_2 \
  $path_file_batch_instances \
  $batch_instances_count \
  0 \
  $path_directory_product \
  $path_directory_disequilibrium \
  $threads \
  $report \
  $path_script_batch_3

  sbatch --array 0-$((index_array_maximum - 3000)):1 --chdir $path_directory_batch_2 \
  $path_script_batch_2 \
  $path_file_batch_instances \
  $batch_instances_count \
  3000 \
  $path_directory_product \
  $path_directory_disequilibrium \
  $threads \
  $report \
  $path_script_batch_3

  sbatch --array 0-$((index_array_maximum - 6000)):1 --chdir $path_directory_batch_3 \
  $path_script_batch_2 \
  $path_file_batch_instances \
  $batch_instances_count \
  6000 \
  $path_directory_product \
  $path_directory_disequilibrium \
  $threads \
  $report \
  $path_script_batch_3

fi




#
