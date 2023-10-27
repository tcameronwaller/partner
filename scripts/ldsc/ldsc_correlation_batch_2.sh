#!/bin/bash

#SBATCH --job-name=ldsc_rg                   # name of job
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



################################################################################
# Organize arguments.

path_file_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
shift_index=${3} # shift index
path_directory_product=${4} # full path to parent directory within which to write product directories and files
path_directory_disequilibrium=${5} # full directory path to directory for LDSC reference on linkage disequilibrium
threads=${6} # count of concurrent or parallel process threads on node cores
report=${7} # whether to print reports
path_script_batch_3=${8} # full path fo file of script for execution of job procedure



################################################################################
# Organize parameters.

# Determine batch instance.
batch_index=$((SLURM_ARRAY_TASK_ID + shift_index)) # Indices in array of batch jobs start at zero.
readarray -t batch_instances < $path_file_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
path_file_base_product="${array[0]}"
path_file_source_primary="${array[1]}"
path_file_source_secondary="${array[2]}"

################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "ldsc_correlation_batch_2.sh"
  echo "----------"
  echo "Product file path:"
  echo $path_file_base_product
  echo "Primary source file path:"
  echo $path_file_source_primary
  echo "Secondary source file path:"
  echo $path_file_source_secondary
  echo "----------"
  echo "Slurm job id: " $SLURM_JOB_ID
  echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
fi



###########################################################################
# Execute procedure.

if true; then
  /usr/bin/bash $path_script_batch_3 \
  $path_file_source_primary \
  $path_file_source_secondary \
  $path_file_base_product \
  $path_directory_disequilibrium \
  $threads \
  $report
fi



#
