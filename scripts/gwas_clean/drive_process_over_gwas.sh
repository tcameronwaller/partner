#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 22 December 2022
# Date, last execution: 19 December 2023
# Review: TCW; 19 December 2023
################################################################################
# Note



################################################################################
# Organize arguments.

path_file_table_parameter=${1} # full path to file of parameter table in text format with information about translations
path_directory_source=${2} # full path to source directory within which to find files of source GWAS summary statistics
path_directory_product=${3} # full path to product directory within which to save files of product GWAS summary statistics
path_file_script_process=${4} # full path to file of script to apply to each file of GWAS summary statistics in the parameter table
path_directory_partner_scripts=${5} # full path to directory within which to find scripts
report=${6} # whether to print reports

################################################################################
# Organize paths.

# Directories.
path_directory_batch="${path_directory_product}/batch"

# Files.
path_file_batch_instances="${path_directory_batch}/batch_instances.txt"

# Scripts.
path_file_script_slurm_job="${path_directory_partner_scripts}/gwas_clean/slurm_job_process_on_gwas.sh"

# Initialize directories.
rm -r $path_directory_batch # caution
mkdir -p $path_directory_batch
cd $path_directory_batch # execute batch from within this directory

# Initialize files.
rm $path_file_batch_instances

################################################################################
# Procedure.

##########
# Initialize array of instance parameters.
instances=()

##########
# Read lines from file and split fields within each line by space, tab, or new-line delimiters.
input=$path_file_table_parameter
while IFS=$' \t\n' read -r -a array
do

  # Extract values from individual columns within table's row.
  # These parameters match the format of the table as of 13 December 2023.
  raw_availability="${array[0]}"
  raw_inclusion="${array[1]}"
  raw_directory="${array[2]}"
  raw_name_study="${array[3]}"
  raw_phenotype="${array[4]}"
  raw_sex="${array[5]}"
  raw_name_file_source="${array[6]}"
  raw_suffix_file_source="${array[7]}"
  raw_bgzip="${array[8]}"
  raw_gzip="${array[9]}"
  raw_type="${array[10]}"
  raw_fill_observations="${array[11]}"
  raw_observations_total="${array[12]}"
  raw_fill_case_control="${array[13]}"
  raw_cases="${array[14]}"
  raw_controls="${array[15]}"
  raw_observations_effective="${array[16]}"
  raw_prevalence_sample="${array[17]}"
  raw_prevalence_population="${array[18]}"
  raw_script="${array[19]}"
  raw_note="${array[20]}"

  # Report.
  if [ $raw_inclusion == "1" ] && [ "$report" == "true" ]; then
    echo "----------"
    echo "field 0, availability: ${raw_availability}"
    echo "field 1, inclusion: ${raw_inclusion}"
    echo "field 2, directory: ${raw_directory}"
    echo "field 3, name_study: ${raw_name_study}"
    echo "field 4, phenotype: ${raw_phenotype}"
    echo "field 5, sex: ${raw_sex}"
    echo "field 6, file: ${raw_name_file_source}"
    echo "field 7, suffix: ${raw_suffix_file_source}"
    echo "field 8, bgzip: ${raw_bgzip}"
    echo "field 9, gzip: ${raw_gzip}"
    echo "field 10, type: ${raw_type}"
    echo "field 11, fill_observations: ${raw_fill_observations}"
    echo "field 12, observations_total: ${raw_observations_total}"
    echo "field 13, fill_case_control: ${raw_fill_case_control}"
    echo "field 14, cases: ${raw_cases}"
    echo "field 15, controls: ${raw_controls}"
    echo "field 16, observations_effective: ${raw_observations_effective}"
    echo "field 17, prevalence_sample: ${raw_prevalence_sample}"
    echo "field 18, prevalence_population: ${raw_prevalence_population}"
    echo "field 19, script: ${raw_script}"
    echo "field 20, note: ${raw_note}"
    echo "----------"
  fi
  # Execute procedure for current record's parameters.
  if [ $raw_inclusion == "1" ]; then
    # Assemble parameters for instance.
    part_one="${raw_directory};${raw_name_file_source};${raw_suffix_file_source};${raw_name_study}"
    part_two="${raw_type};${raw_fill_observations};${raw_fill_case_control}"
    part_three="${raw_observations_total};${raw_cases};${raw_controls};${raw_observations_effective}"
    part_four="${raw_prevalence_sample};${raw_prevalence_population}"
    part_five="${raw_bgzip};${raw_gzip};${raw_script}"
    instances+=("${part_one};${part_two};${part_three};${part_four};${part_five}")
  fi
done < "${input}"



##########
# Organize batch instances.

count_instances=${#instances[@]}

for instance in "${instances[@]}"; do
  # Define parameters in array instance for batch job.
  echo $instance >> $path_file_batch_instances
done

# Read batch instances.
readarray -t batch_instances < $path_file_batch_instances
batch_instances_count=${#batch_instances[@]}
index_array_maximum=$(($batch_instances_count - 1))

################################################################################
# Report.

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Source directory:"
  echo $path_directory_parent_source
  echo "Product directory:"
  echo $path_directory_product
  echo "count of instances: " $count_instances
  echo "first instance: " ${instances[0]} # notice base-zero indexing
  echo "last instance: " ${instances[$count_instances - 1]}
  echo "----------"
  echo "----------"
  echo "count of batch instances: " $batch_instances_count
  echo "maximum array index: " $index_array_maximum
  echo "----------"
  echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
  echo "last batch instance: " ${batch_instances[$index_array_maximum]}
  echo "----------"

fi

sleep 5s

################################################################################
# Submit batch of jobs to grid cluster scheduler for processing.
# Submit to Slurm Scheduler.
# Indices in array of batch jobs start at zero.

if true; then
  sbatch --array 0-${index_array_maximum}:1 --chdir $path_directory_batch \
  $path_file_script_slurm_job \
  $path_file_batch_instances \
  $batch_instances_count \
  $path_directory_source \
  $path_directory_product \
  $report \
  $path_file_script_process
fi

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "drive_process_over_gwas.sh"
  echo "----------"
fi



#
