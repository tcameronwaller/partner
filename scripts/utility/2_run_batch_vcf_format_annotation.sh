#!/bin/bash

###########################################################################
# Specify arguments for qsub command.
# Note that bash does not interpret qsub parameters, which are bash comments.
# Bash will not expand variables in qsub parameters.
# Shell.
#$ -S /bin/bash
# Name of job.
#$ -N waller_gwas
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
# Important to specify 32 threads to avoid inconsistency with interactive
# calculations.
#$ -pe threaded 8
# Range of indices.
# Specify as argument when calling qsub.
# Array batch indices cannot start at zero.
### -t 1-100:1
# Limit on concurrent processes.
#$ -tc 50

# http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

################################################################################
# Note.

# LDSC Munge of GWAS summary statistics throws a memory error with 4 thread
# slots ("-pe threaded 4") and 1 Gigabyte of memory ("-l h_vmem=1G").

# LDSC Munge of GWAS summary statistics runs with 8 thread slots
# ("-pe threaded 8") and 1 Gigabyte of memory ("-l h_vmem=1G").

################################################################################
# Organize argument variables.

path_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
regression_type=${3} # type of GWAS regression, either "linear" or "logistic"
response=${4} # whether GWAS response is beta coefficient ("coefficient"), odds ratio ("odds_ratio"), or z-scores ("z_score")
response_standard_scale=${5} # whether to convert response (coefficient) to z-score standard scale
path_gwas_format_container=${6} # full path to parent directories of GWAS summary statistics for each study
path_gwas_munge_container=${7} # full path to parent directories of GWAS summary statistics for each study
path_heritability_container=${8} # full path to parent directory for heritability reports
path_scripts_record=${9} # full path to directory of scripts for a specific analysis report date
path_process=${10} # full path to directory for all processes relevant to current project
restore_target_study_directories=${11} # whether to delete any previous directories for each study's format and munge GWAS ("true" or "false")

###########################################################################
# Organize variables.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
study="${array[0]}"
path_gwas_concatenation_compress="${array[1]}"

###########################################################################
# Execute procedure.

# Concatenate GWAS across chromosomes.
/usr/bin/bash "${path_scripts_record}/4_drive_gwas_format_munge_heritability.sh" \
$study \
$path_gwas_concatenation_compress \
$regression_type \
$response \
$response_standard_scale \
$path_gwas_format_container \
$path_gwas_munge_container \
$path_heritability_container \
$path_process \
$restore_target_study_directories
