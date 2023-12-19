#!/bin/bash

#SBATCH --job-name=format_gwas               # name of job
#SBATCH --mail-user=waller.tcameron@mayo.edu # email address
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50   # situations in which to send email
#SBATCH --partition=cpu-short                # queue: cpu-short, cpu-med, cpu-long
#SBATCH --nodes=1                            # count of cluster nodes (CPUs)
#SBATCH --ntasks-per-node=1                  # count of cores or threads on node
#SBATCH --mem=8G                             # memory per node (per CPU)
#SBATCH --time=0-03:00:00                    # time allocation request (days-hours:minutes:seconds)
#SBATCH --output ./%x.%A.%N.%j.%a.stdout
#SBATCH --error ./%x.%A.%N.%j.%a.stderr
#SBATCH --signal=USR1@60

################################################################################
# Note.

# Review: TCW; 26 November 2023

################################################################################
# Organize arguments.

path_file_batch_instances=${1} # text list of information for each instance in batch
batch_instances_count=${2} # count of instances in batch
path_directory_parent_source=${3} # full path to parent directory within which to find child directories and files of source GWAS summary statistics
path_directory_product=${4} # full path to directory within which to save files of product GWAS summary statistics
path_directory_parent_scripts_format=${5} # full path to directory within which to find scripts
path_bgzip=${6} # full path to installation executable file of BGZip
report=${7} # whether to print reports

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
raw_fill_observations="${array[4]}"
raw_fill_case_control="${array[5]}"
raw_observations="${array[6]}"
raw_cases="${array[7]}"
raw_controls="${array[8]}"
raw_bgzip="${array[9]}"
raw_gzip="${array[10]}"
raw_script="${array[11]}"

################################################################################
# Organize paths.

path_directory_child_source="${path_directory_parent_source}/$raw_directory"
name_base_file_source=$(echo $raw_name_file_source | sed "s/$raw_suffix_file_source//")
path_file_source="${path_directory_child_source}/${raw_name_file_source}"
path_directory_temporary="${path_directory_product}/temporary_${name_base_file_source}_${raw_name_study}" # hopefully unique
path_file_source_standard="${path_directory_temporary}/${name_base_file_source}.txt.gz"
path_file_script="${path_directory_parent_scripts_format}/$raw_script"
path_file_product="${path_directory_product}/$raw_name_study.txt.gz"

# Initialize directories.
mkdir -p $path_directory_temporary

################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "translate_format_batch_instance.sh"
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
# Manage compression formats and file suffixes.
if [ $raw_bgzip == "1" ] && [ $raw_gzip == "1" ]; then
  # 1. Decompress from BGZip format (http://www.htslib.org/doc/bgzip.html).
  $path_bgzip --decompress "$path_file_source" --stdout > "${path_directory_temporary}/${name_base_file_source}.txt"
  # 2. Compress to GZip format.
  gzip -cvf "${path_directory_temporary}/${name_base_file_source}.txt" > $path_file_source_standard
fi
if [ $raw_bgzip != "1" ] && [ $raw_gzip == "1" ]; then
  # 1. Compress to Gzip format
  gzip -cvf "$path_file_source" > $path_file_source_standard
fi
if [ $raw_bgzip != "1" ] && [ $raw_gzip != "1" ]; then
  # Manage suffix.
  cp "$path_file_source" $path_file_source_standard
fi

##########
# Call script for translation.
/usr/bin/bash $path_file_script \
$path_file_source_standard \
$path_file_product \
$raw_fill_observations \
$raw_observations \
$raw_fill_case_control \
$raw_cases \
$raw_controls \
$report

##########
# Remove temporary directory and files.
rm -r $path_directory_temporary



#
