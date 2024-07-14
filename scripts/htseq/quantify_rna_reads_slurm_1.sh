#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 13 July 2024
# Date, last execution or modification: 13 July 2024
# Review: TCW; 13 July 2024
###############################################################################
# Note


###############################################################################
# Organize arguments.

path_file_parallel_instances=${1} # full path to file in text format with list of information for each instance in parallel batch
path_directory_parallel=${2} # full path to directory for files relating to management of parallel batch
path_file_annotation_gtf_gzip=${3} # full path to file for annotation of reference genome with GZip compression
threads=${4} # count of concurrent or parallel process threads on node cores
report=${5} # whether to print reports to terminal
path_script_parallel_2=${6} # full path to file of script for execution of job procedure
path_script_parallel_3=${7} # full path to file of script for execution of job procedure
path_environment_htseq=${8} # full path to Python virtual environment

###############################################################################
# Organize paths.

# Initialize directories.
cd $path_directory_parallel # execute parallel batch procedures from within this directory

###############################################################################
# Organize job instances in parallel batch.

# Read job instances in parallel batch.
readarray -t instances_parallel < $path_file_parallel_instances
count_parallel_instances=${#instances_parallel[@]}
index_array_maximum=$(($count_parallel_instances - 1)) # indices start at zero in array of parallel batch job instances

###############################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "script:"
  echo $0 # Print full file path to script.
  echo "quantify_rna_reads_slurm_1.sh"
  echo "----------"
  echo "count of instances in parallel batch: " $count_parallel_instances
  echo "maximum index in array of instances: " $index_array_maximum
  echo "----------"
  echo "first instance: " ${instances_parallel[0]} # notice base-zero indexing
  echo "last instance: " ${instances_parallel[$index_array_maximum]}
  echo "----------"
fi

###############################################################################
# Submit batch of parallel job instances to Slurm Scheduler to manage
# processes on computational grid cluster.
# Indices start at zero in array of parallel batch job instances.

# Determine shift in indices if necessary.
# NCSA mForge implementation of SLURM has value of "MaxArraySize" of 10,000
# as confirmed by TCW on 14 May 2024.
# Find the value of "MaxArraySize" within the text file at path
# "/etc/slurm/slurm.conf".
# SLURM will not allow an array index greater than "MaxArraySize".
# It is necessary to split larger jobs.
index_shift=0

if true; then
  sbatch \
  --exclude=mforge160 \
  --array 0-${index_array_maximum}:1 \
  --chdir $path_directory_parallel \
  $path_script_parallel_2 \
  $path_file_parallel_instances \
  $count_parallel_instances \
  $path_directory_parallel \
  $path_file_annotation_gtf_gzip \
  $index_shift \
  $threads \
  $report \
  $path_script_parallel_3 \
  $path_environment_htseq
fi


###############################################################################
# End.
