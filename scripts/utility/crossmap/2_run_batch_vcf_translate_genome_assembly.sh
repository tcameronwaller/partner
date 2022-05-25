#!/bin/bash

###########################################################################
# Specify arguments for qsub command.
# Note that bash does not interpret qsub parameters, which are bash comments.
# Bash will not expand variables in qsub parameters.
# Shell.
#$ -S /bin/bash
# Name of job.
#$ -N tcw_crossmap
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
#$ -pe threaded 16
# Range of indices.
# Specify as argument when calling qsub.
# Array batch indices cannot start at zero.
### -t 1-100:1
# Limit on concurrent processes.
#$ -tc 30

# http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

################################################################################
# Note.

################################################################################
# Organize argument variables.

path_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
path_assembly_translation_chain=${3} # full path to chain file for assembly translation
path_product_genome_assembly_sequence=${4} # full path to product genome assembly sequence file in FASTA format without compression
path_script_translate_genome_assembly_vcf=${5} # full path to script for translation of genome assembly in VCF format
path_environment_crossmap=${6} # full path to Python 3 environment with installation of CrossMap
report=${7} # whether to print reports

###########################################################################
# Organize variables.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
path_vcf_source="${array[0]}"
path_vcf_product="${array[1]}"

###########################################################################
# Execute procedure.

if true; then
  # Genome assembly.
  /usr/bin/bash "${path_script_translate_genome_assembly_vcf}" \
  $path_vcf_source \
  $path_vcf_product \
  $path_assembly_translation_chain \
  $path_product_genome_assembly_sequence \
  $path_environment_crossmap \
  $report
fi

#