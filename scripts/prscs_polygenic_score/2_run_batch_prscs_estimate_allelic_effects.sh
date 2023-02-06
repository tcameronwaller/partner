#!/bin/bash

###########################################################################
# Specify arguments for qsub command.
# Note that bash does not interpret qsub parameters, which are bash comments.
# Bash will not expand variables in qsub parameters.
# Shell.
#$ -S /bin/bash
# Name of job.
#$ -N tcw_prscs
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
#$ -l h_vmem=16G
# Concurrent threads; assigns value to variable NSLOTS.
# Important to specify 32 threads to avoid inconsistency with interactive
# calculations.
#$ -pe threaded 1
# Range of indices.
# Specify as argument when calling qsub.
# Array batch indices cannot start at zero.
### -t 1-100:1
# Limit on concurrent processes.
#$ -tc 30

# http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

################################################################################
# Note.

# TCW; 18 August 2022
# On genotypes from the Mayo Clinic Bipolar Biobank, PRS-CSX ran successfully
# with allocations of 1 thread ("-pe threaded 1"), 16 Gigabytes of memory
# ("-l h_vmem=16G"), and 1 hour of process time ("-q 1-hour").
# With these allocations, PRS-CSX ran to completion on genotypes for most
# chromosomes.
# The SGE manager killed instance jobs for chromosomes 1, 2, 3, 6, and 8
# presumably because they required more than 1 hour to run.

# Recommendation of allocations for approximately 3,000 genotypes with
# stratification by chromosome:
# 24 hours of process time ("-q 1-day")
# 1 thread ("-pe threaded 1")
# 16 Gigabytes of memory ("-l h_vmem=16G")

################################################################################
# Organize argument variables.

path_file_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
path_file_gwas_summary=${3} # full path to source GWAS summary statistics in format for PRS-CS
count_gwas_samples=${4} # integer count of samples for the GWAS study
path_directory_genetic_reference_prscs=${5} # full path to directory for genetic references
population_ancestry=${6} # character code of ancestral population of GWAS study: 'AFR', 'AMR', 'EAS', 'EUR', or 'SAS'
path_directory_allele_effect=${7} # full path to directory for product reports on posterior allele effect size estimates
threads=${8} # count of processing threads to use
path_script_prscs_estimate_allelic_effects=${9} # full path to script for estimation of allelic posterior effects in PRS-CSX
path_environment_prscs=${10} # full path to Python 3 environment with installation of CrossMap
path_prscsx=${11} # full path to installation executable file of PRS-CSX
report=${12} # whether to print reports

###########################################################################
# Organize variables.

# Determine batch instance.
batch_index=$((SGE_TASK_ID-1))
readarray -t batch_instances < $path_file_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
chromosome="${array[0]}"
path_file_genotype_snp_bim_prefix="${array[1]}"
name_file_allele_effect="${array[2]}"

###########################################################################
# Execute procedure.

if true; then
  # Estimate allelic effects in PRS-CSX.
  /usr/bin/bash "${path_script_prscs_estimate_allelic_effects}" \
  $path_file_gwas_summary \
  $count_gwas_samples \
  $path_file_genotype_snp_bim_prefix \
  $path_directory_genetic_reference_prscs \
  $population_ancestry \
  $path_directory_allele_effect \
  $name_file_allele_effect \
  $chromosome \
  $threads \
  $path_environment_prscs \
  $path_prscsx \
  $report
fi



#
