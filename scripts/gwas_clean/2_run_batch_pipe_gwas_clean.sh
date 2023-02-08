#!/bin/bash

###########################################################################
# Specify arguments for qsub command.
# Note that bash does not interpret qsub parameters, which are bash comments.
# Bash will not expand variables in qsub parameters.
# Shell.
#$ -S /bin/bash
# Name of job.
#$ -N tcw_gwas_clean
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
#$ -q 1-day
# Priority 0-15.
### -p -10
# Memory per iteration.
# Segmentation errors commonly indicate a memory error.
#$ -l h_vmem=4G
# Concurrent threads; assigns value to variable NSLOTS.
#$ -pe threaded 16
# Range of indices.
# Specify as argument when calling qsub.
# Array batch indices cannot start at zero.
### -t 1-100:1
# Limit on concurrent processes.
#$ -tc 10

# http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

################################################################################
# Note.

################################################################################
# Organize argument variables.

path_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
path_script_format_gwas_linear=${3} # full path to script for format on linear GWAS
path_script_format_gwas_logistic=${4} # full path to script for format on logistic GWAS
report=${5} # whether to print reports

###########################################################################
# Organize variables.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
type_regression="${array[0]}"
path_file_gwas_source="${array[1]}"
path_file_gwas_product="${array[2]}"

###########################################################################
# Execute procedure.

# Determine appropriate script.
if [[ "$type_regression" == "linear" ]]; then
  path_script_format_gwas="${path_script_format_gwas_linear}"
elif [[ "$type_regression" == "logistic" ]]; then
  path_script_format_gwas="${path_script_format_gwas_logistic}"
else
  echo "invalid specification of regression type"
fi

# Adjust format of GWAS summary statistics.
if true; then
  /usr/bin/bash "${path_script_format_gwas}" \
  $path_file_gwas_source \
  $path_file_gwas_product \
  $report
fi

#
