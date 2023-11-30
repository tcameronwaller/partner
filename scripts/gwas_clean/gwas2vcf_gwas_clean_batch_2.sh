#!/bin/bash

#SBATCH --job-name=gwas_clean                # name of job
#SBATCH --mail-user=waller.tcameron@mayo.edu # email address
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50   # situations in which to send email
#SBATCH --partition=cpu-med                  # queue: cpu-short, cpu-med, cpu-long
#SBATCH --nodes=1                            # count of cluster nodes (CPUs)
#SBATCH --ntasks-per-node=1                  # count of cores or threads on node
#SBATCH --mem=16G                            # memory per node (per CPU)
#SBATCH --time=3-00:00:00                    # time allocation request (days-hours:minutes:seconds)
#SBATCH --output ./%x.%A.%N.%j.%a.stdout
#SBATCH --error ./%x.%A.%N.%j.%a.stderr
#SBATCH --signal=USR1@60

################################################################################
# Note.



################################################################################
# Organize argument variables.

path_file_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
path_file_script_batch_3=${3} # full path to script for format on logistic GWAS
threads=${4} # count of processing threads to use
report=${5} # whether to print reports

################################################################################
# Organize parameters.

# Determine batch instance.
batch_index=$SLURM_ARRAY_TASK_ID # Indices in array of batch jobs start at zero.
readarray -t batch_instances < $path_file_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
path_file_gwas_source="${array[0]}"
path_file_gwas_product="${array[1]}"
type="${array[2]}"

################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "gwas2vcf_gwas_clean_batch_2.sh"
  echo "----------"
  echo "Slurm job id: " $SLURM_JOB_ID
  echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
  echo "----------"
  echo "source file GWAS:"
  echo $path_file_gwas_source
  echo "product file GWAS:"
  echo $path_file_gwas_product
  echo "----------"
fi

###########################################################################
# Execute procedure.

# Adjust format of GWAS summary statistics.
if true; then
  /usr/bin/bash $path_file_script_batch_3 \
  $path_file_gwas_source \
  $path_file_gwas_product \
  $type \
  $threads \
  $report
fi



#
