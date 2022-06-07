#!/bin/bash

###########################################################################
# Specify arguments for qsub command.
# Note that bash does not interpret qsub parameters, which are bash comments.
# Bash will not expand variables in qsub parameters.
# Shell.
#$ -S /bin/bash
# Name of job.
#$ -N tcw_bcftools
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
#$ -l h_vmem=5G
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
path_file_genome_assembly_sequence=${3} # full path to genome assembly sequence file in FASTA format without compression that matches source and product genome assembly
path_file_chromosome_translations=${4} # full path to file for chromosome name translations in format for BCFTools "annotate --rename-chrs"
path_file_dbsnp_reference=${5} # full path to file for dbSNP reference in VCF format
threads=${6} # count of processing threads to use
path_script_drive_preparation_format_annotation=${7} # full path to script for preparation of genotype VCF files
path_script_decompose_align_unique_sort=${8} # full path to script for preparation of genotype VCF files
path_script_translate_chromosomes=${9} # full path to script for preparation of genotype VCF files
path_script_introduce_dbsnp_rsid=${10} # full path to script for preparation of genotype VCF files
path_bcftools=${11} # full path to installation executable file of BCFTools
report=${12} # whether to print reports

###########################################################################
# Organize variables.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
path_file_vcf_source_chromosome="${array[0]}"
path_file_vcf_product_chromosome="${array[1]}"

###########################################################################
# Execute procedure.

if true; then
  # Prepare genotype files.
  /usr/bin/bash "${path_script_drive_preparation_format_annotation}" \
  $path_file_vcf_source_chromosome \
  $path_file_vcf_product_chromosome \
  $path_file_genome_assembly_sequence \
  $path_file_chromosome_translations \
  $path_file_dbsnp_reference \
  $threads \
  $path_script_decompose_align_unique_sort \
  $path_script_translate_chromosomes \
  $path_script_introduce_dbsnp_rsid \
  $path_bcftools \
  $report
fi



#
