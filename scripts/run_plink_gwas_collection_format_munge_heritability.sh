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
#$ -M tcameronwaller@gmail.com
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

# Collection and concatenation of GWAS summary statistics runs with 4 thread
# slots ("-pe threaded 4") and 1 Gigabyte of memory ("-l h_vmem=1G"), even for
# GWAS on > 300,000 persons.

# LDSC Munge of GWAS summary statistics throws a memory error with 4 thread
# slots ("-pe threaded 4") and 1 Gigabyte of memory ("-l h_vmem=1G").

# LDSC Munge of GWAS summary statistics runs with 8 thread slots
# ("-pe threaded 8") and 1 Gigabyte of memory ("-l h_vmem=1G").


################################################################################
# Organize argument variables.

path_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
path_gwas_parent=${3}
path_heritability_parent=${4}
path_genetic_reference=${5}
path_promiscuity_scripts=${6}
path_ldsc=${7} # full path to executable LDSC
report=${8}

###########################################################################
# Organize variables.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
directory=${batch_instances[$batch_index]}

path_script_gwas_collect_concatenate="${path_promiscuity_scripts}/collect_concatenate_gwas_chromosomes.sh"
path_promiscuity_scripts_ldsc_heritability="${path_promiscuity_scripts}/ldsc_genetic_heritability_correlation"
path_script_format_munge_heritability="${path_promiscuity_scripts_ldsc_heritability}/format_munge_gwas_heritability_ldsc.sh"
path_scripts_format="${path_promiscuity_scripts}/format_gwas_ldsc"
path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_plink_linear.sh"

###########################################################################
# Execute procedure.

# Concatenate GWAS across chromosomes.
if true; then
  # Organize variables.
  pattern_source_file="report.*.glm.linear" # do not expand with full path yet
  path_source_directory="${path_gwas_parent}/${directory}"
  chromosome_start=1
  chromosome_end=22
  path_gwas_concatenation="${path_source_directory}/gwas_concatenation.txt"
  path_gwas_concatenation_compress="${path_source_directory}/gwas_concatenation.txt.gz"
  /usr/bin/bash "$path_script_gwas_collect_concatenate" \
  $pattern_source_file \
  $path_source_directory \
  $chromosome_start \
  $chromosome_end \
  $path_gwas_concatenation \
  $path_gwas_concatenation_compress \
  $report
fi

# Format and munge GWAS summary statistics.
# Estimate genotype heritability.
if true; then
  # Organize paths.
  path_gwas_study="${path_gwas_parent}/${directory}"
  path_heritability_study="${path_heritability_parent}/${directory}"
  # Initialize directories.
  #mkdir -p $path_gwas_study
  mkdir -p $path_heritability_study
  # Organize variables.
  study="${directory}"
  name_prefix="null"
  path_source_file="${path_gwas_study}/gwas_concatenation.txt.gz"
  report="true" # "true" or "false"
  /usr/bin/bash "$path_script_format_munge_heritability" \
  $study \
  $name_prefix \
  $path_source_file \
  $path_genetic_reference \
  $path_gwas_study \
  $path_heritability_study \
  $path_script_gwas_format \
  $path_promiscuity_scripts \
  $path_ldsc \
  $report
fi
