#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 27 September 2023
# Date, last execution: 27 September 2023
# Review: TCW; 27 September 2023
################################################################################
# Note



################################################################################
# Organize arguments.

path_file_batch_instances=${1} # full directory path and file name for batch instances in a text file
path_directory_batch=${2} # full path to directory for batch files
path_directory_process=${3} # full path to process directory
report=${4} # whether to print reports

################################################################################
# Organize paths.

# Scripts.
path_file_script_batch_2="${path_directory_process}/partner/scripts/gwas_clean/gwas2vcf_gwas_clean_batch_2.sh"
path_file_script_batch_3="${path_directory_process}/partner/scripts/gwas_clean/process_gwas_gwas2vcf.sh"

# Initialize directories.
cd $path_directory_batch # it is important to execute batch from within directory


################################################################################
# Organize parameters.

threads=1
#report="true"

################################################################################
# Execute procedure.


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
  echo "gwas2vcf_gwas_clean_batch_1.sh"
  echo "----------"
  echo "count of batch instances: " $batch_instances_count
  echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
  echo "last batch instance: " ${batch_instances[$index_array_maximum]}
  echo "----------"
fi

################################################################################
# Submit batch of jobs to grid cluster scheduler for processing.
# Submit to Slurm Scheduler.
# Indices in array of batch jobs start at zero.

if true; then
  sbatch --array 0-${index_array_maximum}:1 --chdir $path_directory_batch \
  $path_file_script_batch_2 \
  $path_file_batch_instances \
  $batch_instances_count \
  $path_file_script_batch_3 \
  $threads \
  $report
fi



#
