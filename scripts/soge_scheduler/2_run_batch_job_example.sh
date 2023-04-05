#!/bin/bash

################################################################################
# Specify arguments for qsub command.
# Note that bash does not interpret qsub parameters, which are bash comments.
# Bash will not expand variables in qsub parameters.
# Shell.
#$ -S /bin/bash
# Name of job.
#$ -N tcw_example
# Contact.
# "b": beginning, "e": end, "a": abortion, "s": suspension, "n": never
#$ -M waller.tcameron@mayo.edu
#$ -m as
# Standard output and error.
# Specify as arguments when calling qsub.
### -o "./out"
### -e "./error"
# Queue.
# "1-hour", "1-day", "4-day", "7-day", "30-day", "lg-mem"
#$ -q 1-hour
# Priority 0-15.
### -p -10
# Memory per iteration.
# Segmentation errors commonly indicate a memory error.
#$ -l h_vmem=1G
# Concurrent threads; assigns value to variable NSLOTS.
#$ -pe threaded 1
# Range of indices.
# Specify as argument when calling qsub.
# Array batch indices cannot start at zero.
### -t 1-100:1
# Limit on concurrent processes.
#$ -tc 25

################################################################################
# Note.

# TODO: TCW; 4 April 2023
# TODO: it seems like it ought to be possible to pass the "threads" parameter
# to SGE when submitting the batch.

################################################################################
# Organize arguments.

path_file_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
message_common=${3} # message parameter common to all jobs
path_directory_product=${4} # full path to parent directory for product files
path_script_execute_procedure=${5} # full path fo file of script for execution of job procedure
threads=${6} # count of concurrent or parallel process threads on node cores
report=${7} # whether to print reports



###########################################################################
# Organize parameters.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1)) # Indices in array of batch jobs start at one, not zero.
readarray -t batch_instances < $path_file_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
chromosome="${array[0]}"
name_file_prefix="${array[1]}"
name_file_suffix="${array[2]}"


################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo "2_run_batch_job_example.sh"
  echo "----------"
  echo "Common message: " $message_common
  echo "chromosome: " $chromosome
  echo "prefix: " $name_file_prefix
  echo "suffix: " $name_file_suffix
  echo "----------"
fi



###########################################################################
# Execute procedure.

if true; then
  /usr/bin/bash $path_script_execute_procedure \
  $chromosome \
  $name_file_prefix \
  $name_file_suffix \
  $message_common \
  $path_directory_product \
  $threads \
  $report
fi



#
