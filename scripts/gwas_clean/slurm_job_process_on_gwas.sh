#!/bin/bash

#SBATCH --job-name=clean_gwas                # name of job
#SBATCH --mail-user=waller.tcameron@mayo.edu # email address
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50   # situations in which to send email
#SBATCH --partition=cpu-short                # queue: cpu-short, cpu-med, cpu-long
#SBATCH --nodes=1                            # count of cluster nodes (CPUs)
#SBATCH --ntasks-per-node=1                  # count of cores or threads on node
#SBATCH --mem=8G                             # memory per node (per CPU)
#SBATCH --time=0-05:00:00                    # time allocation request (days-hours:minutes:seconds)
#SBATCH --output ./%x.%A.%N.%j.%a.stdout
#SBATCH --error ./%x.%A.%N.%j.%a.stderr
#SBATCH --signal=USR1@60

################################################################################
# Note.

# Review: TCW; 19 December 2023



################################################################################
# Organize arguments.

path_file_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
path_directory_source=${3} # full path to source directory within which to find files of source GWAS summary statistics
path_directory_product=${4} # full path to product directory within which to save files of product GWAS summary statistics
report=${5} # whether to print reports
path_file_script_process=${6} # full path to file of script to apply to each file of GWAS summary statistics in the parameter table

################################################################################
# Organize parameters.

# Determine batch instance.
batch_index=$((SLURM_ARRAY_TASK_ID)) # Indices in array of batch jobs start at zero.
readarray -t batch_instances < $path_file_batch_instances
instance=${batch_instances[$batch_index]}

# Separate fields from instance.
IFS=";" read -r -a array <<< "${instance}"
raw_directory="${array[0]}"
raw_name_file_source="${array[1]}"
raw_suffix_file_source="${array[2]}"
raw_name_study="${array[3]}"
raw_type="${array[4]}"
raw_fill_observations="${array[5]}"
raw_fill_case_control="${array[6]}"
raw_observations_total="${array[7]}"
raw_cases="${array[8]}"
raw_controls="${array[9]}"
raw_observations_effective="${array[10]}"
raw_prevalence_sample="${array[11]}"
raw_prevalence_population="${array[12]}"
raw_bgzip="${array[13]}"
raw_gzip="${array[14]}"
raw_script="${array[15]}"

################################################################################
# Organize paths.

path_file_gwas_source="${path_directory_source}/${raw_name_study}.txt.gz"
path_file_gwas_product="${path_directory_product}/${raw_name_study}.txt.gz"

################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "slurm_job_process_on_gwas.sh"
  echo "----------"
  echo "Process script:"
  echo $path_file_script_process
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
# Call script for process.
/usr/bin/bash $path_file_script_process \
$path_file_gwas_source \
$path_file_gwas_product \
$report



#
