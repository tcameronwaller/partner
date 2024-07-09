#!/bin/bash

#SBATCH --job-name=samtools                  # name of job
#SBATCH --mail-user=waller.tcameron@mayo.edu # email address
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50   # situations in which to send email
#SBATCH --partition=cpu-short                # queue: cpu-short, cpu-med, cpu-long
#SBATCH --nodes=1                            # count of cluster nodes (CPUs)
#SBATCH --ntasks-per-node=2                  # count of cores or threads on node
#SBATCH --mem=4G                             # memory per node (per CPU)
#SBATCH --time=0-01:00:00                    # time allocation request (days-hours:minutes:seconds)
#SBATCH --output ./%x.%A.%N.%j.%a.stdout
#SBATCH --error ./%x.%A.%N.%j.%a.stderr
#SBATCH --signal=USR1@60

################################################################################
# Note.

# This script must execute within a working directory set up as the batch
# directory.

# NCSA mForge implementation of SLURM has value of "MaxArraySize" of 10,000
# as confirmed by TCW on 14 May 2024.
# Find the value of "MaxArraySize" within the text file at path
# "/etc/slurm/slurm.conf".
# SLURM will not allow an array index greater than "MaxArraySize".
# It is necessary to split larger jobs.

################################################################################
# Organize arguments.

path_file_parallel_instances=${1} # full path to file in text format with list of information for each instance in parallel batch
count_parallel_instances=${2} # count of instances in parallel batch
shift_index=${3} # shift index necessary when count of instances exceeds variable "MaxArraySize"

path_file_reference_genome=${3} # full path to file for reference genome sequence
path_execution_samtools=${4} # full path to executable file for SamTools
report=${5} # whether to print reports



path_directory_product=${4} # full path to parent directory within which to write product directories and files
path_directory_disequilibrium=${5} # full directory path to directory for LDSC reference on linkage disequilibrium
threads=${6} # count of concurrent or parallel process threads on node cores
report=${7} # whether to print reports
path_script_parallel_3=${8} # full path to file of script for execution of job procedure
