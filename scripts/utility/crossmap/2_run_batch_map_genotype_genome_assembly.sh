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
#$ -q 4-day
# Priority 0-15.
### -p -10
# Memory per iteration.
# Segmentation errors commonly indicate a memory error.
#$ -l h_vmem=10G
# Concurrent threads; assigns value to variable NSLOTS.
# Important to specify 32 threads to avoid inconsistency with interactive
# calculations.
#$ -pe threaded 32
# Range of indices.
# Specify as argument when calling qsub.
# Array batch indices cannot start at zero.
### -t 1-100:1
# Limit on concurrent processes.
#$ -tc 5

# http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

################################################################################
# Note.

################################################################################
# Organize argument variables.

path_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
path_assembly_translation_chain=${3} # full path to chain file for assembly translation
path_genome_assembly_sequence=${4} # full path to genome assembly sequence file in FASTA format without compression that matches product genome assembly
threads=${5} # count of processing threads to use
path_script_translate_assembly=${6} # full path to script for combination and sort of genotype files in VCF format
path_environment_crossmap=${7} # full path to Python 3 environment with installation of CrossMap
path_bcftools=${8} # full path to installation executable file of BCFTools
report=${9} # whether to print reports

###########################################################################
# Organize variables.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
path_file_source_vcf="${array[0]}"
path_file_product_vcf="${array[1]}"

###########################################################################
# Execute procedure.

if true; then
  # Genome assembly.
  /usr/bin/bash "${path_script_translate_assembly}" \
  $path_file_source_vcf \
  $path_file_product_vcf \
  $path_assembly_translation_chain \
  $path_genome_assembly_sequence \
  $threads \
  $path_environment_crossmap \
  $path_bcftools \
  $report
fi



#
