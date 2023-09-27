#!/bin/bash

#SBATCH --job-name=tcw_ldsc                  # name of job
#SBATCH --mail-user=waller.tcameron@mayo.edu # email address
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50   # situations in which to send email
#SBATCH --partition=cpu-short                # queue: cpu-short, cpu-med, cpu-long
#SBATCH --nodes=1                            # count of cluster nodes (CPUs)
#SBATCH --tasks-per-node=4                   # count of cores or threads on node
#SBATCH --time=0-00:30:00                    # time allocation request (days-hours:minutes:seconds)
#SBATCH --mem=4G                             # memory per node (per CPU)
#SBATCH --output logs/%x.%A.%N.%j.%a.stdout
#SBATCH --output logs/%x.%A.%N.%j.%a.stderr
#SBATCH --signal=USR1@60

# --mem=1G                             # total memory per node (per CPU)
# --mem-per-cpu=1G                     # memory per task (per core or thread on CPU node)

# Slurm shortcut variables.
# x: Job name
# A: Slurm array Job identifier
# j: Slurm job number
# a: Slurm array task ID
# N: Node name

# Slurm is not able to create the "logs" child directory within the current
# working parent directory. It is necessary to create this "logs" child
# directory before executing the job script.

# It is also possible to pass any of the parameters above when calling the
# SLURM run script.
# Especially consider calling "--array" and "--chdir" when calling SLURM script.

# squeue -u [User LANID]
# scontrol show job [job identifier] # Only lasts for 5 minutes.
# sacct -j [job identifier] --format=[names of variables "sacct -e"] # "%25" allows to expand a specific column
# sacct -j [job identifier] --format=User,JobID%15,Jobname%25,partition,nodelist,SystemComment,reason$20,exitCode
# seff -f [job identifier] # Efficiency of job utilization of allocated resources.

# salloc # Request an allocation of a single node for interactive use.
# srun # In some ways similar to "salloc".

# https://slurm.schedmd.com/sbatch.html
# https://researchcomputing.princeton.edu/support/knowledge-base/memory

################################################################################
# Note.



################################################################################
# Organize arguments.

path_file_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
path_directory_product=${3} # full path to parent directory within which to write product directories and files
path_directory_disequilibrium=${4} # full directory path to directory for LDSC reference on linkage disequilibrium
threads=${5} # count of concurrent or parallel process threads on node cores
report=${6} # whether to print reports
path_script_batch_3=${7} # full path fo file of script for execution of job procedure



################################################################################
# Organize parameters.

# Determine batch instance.
batch_index=$SLURM_ARRAY_TASK_ID # Indices in array of batch jobs start at zero.
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
