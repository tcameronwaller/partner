#!/bin/bash

#SBATCH --job-name=samtools                  # name of job
#SBATCH --mail-user=waller.tcameron@mayo.edu # email address
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50   # situations in which to send email
#SBATCH --partition=cpu-short                # queue: cpu-short, cpu-med, cpu-long
#SBATCH --nodes=1                            # count of cluster nodes (CPUs)
#SBATCH --ntasks-per-node=2                  # count of cores or threads on node
#SBATCH --mem=0.5G                           # memory per node (per CPU)
#SBATCH --time=0-01:00:00                    # time allocation request (days-hours:minutes:seconds)
#SBATCH --output ./%x.%A.%N.%j.%a.stdout
#SBATCH --error ./%x.%A.%N.%j.%a.stderr
#SBATCH --signal=USR1@60

###############################################################################
# Note.

# This script must execute within a working directory set up as the parent
# directory of child directories and files for management of parallel batch.

# NCSA mForge implementation of SLURM has value of "MaxArraySize" of 10,000
# as confirmed by TCW on 14 May 2024.
# Find the value of "MaxArraySize" within the text file at path
# "/etc/slurm/slurm.conf".
# SLURM will not allow an array index greater than "MaxArraySize".
# It is necessary to split larger jobs.

###############################################################################
# Organize arguments.

path_file_parallel_instances=${1} # full path to file in text format with list of information for each instance in parallel batch
count_parallel_instances=${2} # count of instances in parallel batch
path_directory_parallel=${3} # full path to directory for files relating to management of parallel batch
index_shift=${4} # shift index necessary when count of instances exceeds system variable "MaxArraySize"
path_file_reference_genome=${5} # full path to file for reference genome sequence
path_file_reference_genome_index=${6} # full path to file for reference genome sequence index
threads=${7} # count of concurrent or parallel process threads on node cores
report=${8} # whether to print reports to terminal
path_script_parallel_3=${9} # full path to file of script for execution of job procedure
path_execution_samtools=${10} # full path to executable file for SamTools

###############################################################################
# Organize parameters.

# Determine current job instance of parallel batch.
index_parallel=$((SLURM_ARRAY_TASK_ID + index_shift)) # indices start at zero in array of parallel batch job instances
readarray -t instances_parallel < $path_file_parallel_instances
instance=${instances_parallel[$index_parallel]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
path_file_source="${array[0]}"
path_file_product="${array[1]}"

###############################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "script:"
  echo $0 # Print full file path to script.
  echo "convert_cram_to_bam_slurm_2.sh"
  echo "----------"
  echo "instance in parallel batch"
  echo "Slurm job id: " $SLURM_JOB_ID
  echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
  echo "----------"
  echo "path to source file:"
  echo $path_file_source
  echo "path to product file:"
  echo $path_file_product
  echo "----------"
fi


###############################################################################
# Execute procedure.

# Execute main procedure.
if true; then
  /usr/bin/bash $path_script_parallel_3 \
  $path_file_source \
  $path_file_product \
  $path_file_reference_genome \
  $path_file_reference_genome_index \
  $threads \
  $report \
  $path_execution_samtools
fi

###############################################################################
# End.
