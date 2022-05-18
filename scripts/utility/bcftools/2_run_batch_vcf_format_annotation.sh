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
#$ -q 1-day
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
path_chromosome_translations=${3} # full path to file for chromosome name translations in format for BCFTools "annotate --rename-chrs"
path_dbsnp_reference=${4} # full path to file for dbSNP reference in VCF format
path_script_chromosome_in_vcf=${5} # full path to script for chromosome identifier
path_script_dbsnp_rsid_to_vcf=${6} # full path to script for SNP identifier
threads=${7} # count of processing threads to use
path_bcftools=${8} # full path to installation of BCFTools
report=${9} # whether to print reports

###########################################################################
# Organize variables.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
path_vcf_source_chromosome="${array[0]}"
path_vcf_product_chromosome="${array[1]}"
path_vcf_source_snp="${array[2]}"
path_vcf_product_snp="${array[3]}"

###########################################################################
# Execute procedure.

if true; then
  # Chromosome identifiers.
  /usr/bin/bash "${path_script_chromosome_in_vcf}" \
  $path_chromosome_translations \
  $path_vcf_source_chromosome \
  $path_vcf_product_chromosome \
  $threads \
  $path_bcftools \
  $report
  # SNP identifiers.
  /usr/bin/bash "${path_script_dbsnp_rsid_to_vcf}" \
  $path_dbsnp_reference \
  $path_vcf_source_snp \
  $path_vcf_product_snp \
  $threads \
  $path_bcftools \
  $report
fi



#
