#!/bin/bash

#SBATCH --job-name=tcw_score                 # name of job
#SBATCH --mail-user=waller.tcameron@mayo.edu # email address
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50   # situations in which to send email
#SBATCH --partition=cpu-med                  # queue: cpu-short, cpu-med, cpu-long
#SBATCH --nodes=1                            # count of cluster nodes (CPUs)
#SBATCH --tasks=4                            # count of cores or threads on node
#SBATCH --time=0-07:00:00                    # time allocation request (days-hours:minutes:seconds)
#SBATCH --mem=8G                             # memory per node (per CPU)
#SBATCH --output logs/%x.%A.%N.%j.%a.stdout
#SBATCH --output logs/%x.%A.%N.%j.%a.stderr
#SBATCH --signal=USR1@60

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

################################################################################
# Note.



################################################################################
# Organize arguments.

path_file_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
path_directory_product=${3} # full path to parent directory for product files
threads=${4} # count of concurrent or parallel process threads on node cores
report=${5} # whether to print reports
path_script_calculate_score=${6} # full path fo file of script for execution of job procedure



###########################################################################
# Organize parameters.

# Determine batch instance.
batch_index=$SLURM_ARRAY_TASK_ID # Indices in array of batch jobs start at zero.
readarray -t batch_instances < $path_file_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
path_file_source_effects="${array[0]}"
path_file_source_genotypes="${array[1]}"
name_base_file_product="${array[2]}"



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "2_run_batch_polygenic_score.sh"
  echo "----------"
  echo "Slurm job id: " $SLURM_JOB_ID
  echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
  echo "path file effects: " $path_file_source_effects
  echo "path file genotypes: " $path_file_source_genotypes
  echo "base name scores: " $name_base_file_product
  echo "----------"
fi



###########################################################################
# Execute procedure.

if true; then
  /usr/bin/bash $path_script_calculate_score \
  $path_file_source_effects \
  $path_file_source_genotypes \
  $path_directory_product \
  $name_base_file_product \
  $threads \
  $report
fi



#
