#!/bin/bash

#SBATCH --job-name=tcw_example               # name of job
#SBATCH --mail-user=waller.tcameron@mayo.edu # email address
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50   # situations in which to send email
#SBATCH --partition=cpu-short                # queue: cpu-short, cpu-med, cpu-long
#SBATCH --nodes=1                            # count of cluster nodes (CPUs)
#SBATCH --ntasks-per-node=1                  # count of cores or threads on node
#SBATCH --mem=1G                             # memory per node (per CPU)
#SBATCH --time=0-00:01:00                    # time allocation request (days-hours:minutes:seconds)
#SBATCH --output logs/%x.%A.%N.%j.%a.stdout
#SBATCH --error logs/%x.%A.%N.%j.%a.stderr
#SBATCH --signal=USR1@60

################################################################################
# Note.

# Use "sbatch --help" to see descriptions of parameters.
# Use syntax "--<option>=<value>" for full name of parameter.
# Use syntax "-<option> <value>" for abbreviations of parameter.

# --mem=1G                             # total memory per node (per CPU)
# --mem-per-cpu=1G                     # memory per task (per core or thread on CPU node)

# cpu-short default time 2 minutes and default max is 4 days (must set time parameter)
# cpu-med default time 2 minutes and default max is 4 days (must set time parameter)
# cpu-long default time 2 minutes and default max is 4 days (must set time parameter)

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
# scancel {job identifier}
# scontrol show job [job identifier] # Only lasts for 5 minutes.
# sacct -j [job identifier] --format=[names of variables "sacct -e"] # "%25" allows to expand a specific column
# sacct -j [job identifier] --format=User,JobID%15,Jobname%25,partition,nodelist,SystemComment,reason$20,exitCode
# seff -f [job identifier] # Efficiency of job utilization of allocated resources.

# salloc # Request an allocation of a single node for interactive use.
# srun # In some ways similar to "salloc".

# https://slurm.schedmd.com/sbatch.html
# https://researchcomputing.princeton.edu/support/knowledge-base/memory


################################################################################
# Organize arguments.

path_file_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
message_common=${3} # message parameter common to all jobs
path_directory_product=${4} # full path to parent directory for product files
path_script_execute_procedure=${5} # full path fo file of script for execution of job procedure
threads=${6} # count of concurrent or parallel process threads on node cores
report=${7} # whether to print reports



################################################################################
# Organize parameters.

# Determine batch instance.
batch_index=$SLURM_ARRAY_TASK_ID # Indices in array of batch jobs start at zero.
readarray -t batch_instances < $path_file_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
chromosome="${array[0]}"
name_file_prefix="${array[1]}"
name_file_suffix="${array[2]}"


################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "2_run_batch_job_example.sh"
  echo "----------"
  echo "Common message: " $message_common
  echo "chromosome: " $chromosome
  echo "prefix: " $name_file_prefix
  echo "suffix: " $name_file_suffix
  echo "----------"
  echo "Slurm job id: " $SLURM_JOB_ID
  echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
fi



###########################################################################
# Execute procedure.

if true; then
  /usr/bin/bash $path_script_execute_procedure \
  $chromosome \
  $name_file_prefix \
  $name_file_suffix \
  $message_common \
  $path_directory_product \
  $threads \
  $report
fi



#
