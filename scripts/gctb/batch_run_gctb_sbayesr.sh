#!/bin/bash

###########################################################################
# Specify arguments for qsub command.
# Note that bash does not interpret qsub parameters, which are bash comments.
# Bash will not expand variables in qsub parameters.
# Shell.
#$ -S /bin/bash
# Name of job.
#$ -N tcw_gctb_sbayesr
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
#$ -l h_vmem=4G
# Concurrent threads; assigns value to variable NSLOTS.
#$ -pe threaded 1
# Range of indices.
# Specify as argument when calling qsub.
# Array batch indices cannot start at zero.
### -t 1-100:1
# Limit on concurrent processes.
#$ -tc 25

# http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

################################################################################
# Note.

# TCW; 17 January 2023
# Each batch instance job ran the GCTB SBayesR procedure on a single autosomal
# chromosome.
# This batch completed successfully in the queue with duration of 1 hour
# ("-q 1-hour"), with 1 node core thread ("-pe threaded 1"), and 4 Gigabytes of
# memory ("-l h_vmem=4G").



################################################################################
# Organize arguments.

path_file_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
path_file_gwas=${3} # full path and name to file for source GWAS summary statistics in GCTB and GCTA-COJO ".ma" format without compression
path_script_run_sbayesr=${4} # full path to directory and file of script for direct run of GCTB SBayesR
path_gctb=${5} # full path to directory and file for local executable installation of GCTB SBayesR
threads=${6} # count of concurrent or parallel process threads on node cores
report=${7} # whether to print reports



###########################################################################
# Organize parameters.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_file_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
chromosome="${array[0]}"
path_file_base_ld_matrix="${array[1]}"
path_file_base_product="${array[2]}"


################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo "batch_run_gctb_sbayesr.sh"
  echo "----------"
  echo "GWAS summary statistics: " $path_file_gwas
  echo "chromosome: " $chromosome
  echo "path_file_base_ld_matrix: " $path_file_base_ld_matrix
  echo "path_file_base_product: " $path_file_base_product
  echo "----------"
fi



###########################################################################
# Execute procedure.

if true; then
  /usr/bin/bash $path_script_run_sbayesr \
  $path_file_gwas \
  $path_file_base_ld_matrix \
  $path_file_base_product \
  $path_gctb \
  $threads \
  $report
fi





#
