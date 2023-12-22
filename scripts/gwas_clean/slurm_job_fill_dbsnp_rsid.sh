#!/bin/bash

#SBATCH --job-name=fill_rsid                 # name of job
#SBATCH --mail-user=waller.tcameron@mayo.edu # email address
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50   # situations in which to send email
#SBATCH --partition=lg-mem                   # queue: cpu-short, cpu-med, cpu-long, lg-mem
#SBATCH --nodes=1                            # count of cluster nodes (CPUs)
#SBATCH --ntasks-per-node=1                  # count of cores or threads on node
#SBATCH --mem=320G                           # memory per node (per CPU)
#SBATCH --time=0-05:00:00                    # time allocation request (days-hours:minutes:seconds)
#SBATCH --output ./%x.%A.%N.%j.%a.stdout
#SBATCH --error ./%x.%A.%N.%j.%a.stderr
#SBATCH --signal=USR1@60

################################################################################
# Note.

# Review: TCW; 19 December 2023

# Note: TCW; 19 December 2023
# The process of the full "fill_dbsnp_rsid.sh" script, without any division into
# smaller sections (subscripts), completed for the study
# "37872160_williams_2023", using 278.20 Gigabytes of memory with management by
# the SLURM scheduler.
# Subsequently, the process of the "fill_dbsnp_rsid.sh" script, after division
# into five smaller sections (subscripts) in an attempt to save memory,
# completed for the study "37872160_williams_2023", using 274.89 Gigabytes of
# memory with management by the SLURM scheduler.
# It appears that division of the process into smaller subscripts did not offer
# any advantage of clearing system memory between subscripts.

################################################################################
# Organize arguments.

path_file_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
strict=${3} # whether to return GWAS summary statistics filtered to SNPs with successful match to dbSNP rsID
report=${4} # whether to print reports
path_file_script=${5} # full path to script file

################################################################################
# Organize parameters.

# Determine batch instance.
batch_index=$((SLURM_ARRAY_TASK_ID)) # Indices in array of batch jobs start at zero.
readarray -t batch_instances < $path_file_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
path_file_source="${array[0]}"
path_file_product="${array[1]}"

################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "slurm_job_fill_dbsnp_rsid.sh"
  echo "----------"
  echo "Source file path:"
  echo $path_file_source
  echo "Product file path:"
  echo $path_file_product
  echo "----------"
  echo "Slurm job id: " $SLURM_JOB_ID
  echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
  echo "----------"
fi

###########################################################################
# Execute procedure.


##########
# Call script.
/usr/bin/bash $path_file_script \
$path_file_source \
$path_file_product \
$strict \
$report



#
