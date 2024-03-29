#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 20 March 2023
# Date, last execution: 21 March 2023
# Review: TCW; 21 March 2023
################################################################################
# Note

# This Bash script calls a Python script to read all matching files from a
# source directory for polygenic scores. These files are probably for polygenic
# scores of the same trait on separate chromosomes. The Python script reads
# these files as text tables, assigns unique names to the columns for the
# polygenic scores, merges the tables together by matching identifiers for
# genotypes, calculates the sum of all polygenic scores for each identifier,
# calculates the mean of polygenic scores by division by count of variants, and
# writes the table as a text table to the product file path.

# It is necessary to call the Python script within a package directory where it
# has access to appropriate subpackages and modules.

################################################################################
# Organize arguments.

path_directory_source=${1} # full path to source directory from which to read files of polygenic scores across chromosomes
name_file_source_prefix=${2} # prefix of name of files of polygenic scores from each chromosome
name_file_source_suffix=${3} # suffix of name of files of polygenic scores from each chromosome
name_file_source_not=${4} # character string within file names for exclusion
path_file_product=${5} # full path to product file in UCSC Browser Extensible Data (BED) format

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_waller_tools=$(<"./waller_tools.txt")
path_environment_main="${path_waller_tools}/python/environments/main"
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_package="${path_directory_process}/psychiatric_metabolism/psychiatric_metabolism"
path_directory_product="$(dirname $path_file_product)"

# Scripts.
path_file_script_source="${path_directory_process}/promiscuity/scripts/python/script_combine_sum_polygenic_scores.py"
path_file_script_product="${path_directory_package}/script_combine_sum_polygenic_scores.py"

# Initialize directories.
#rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product

# Initialize files.
rm $path_file_product

################################################################################
# Organize parameters.

#name_file_source_suffix=".sscore"
#name_file_source_not=".vars"
name_column_identifier="#IID"
name_column_allele_total="ALLELE_CT"
name_column_allele_dosage="NAMED_ALLELE_DOSAGE_SUM"
name_column_score_sum="SCORE1_SUM"
name_column_score_mean="SCORE1_AVG"
threads=2
report="true"

################################################################################
# Execute procedure.

# Copy Python script to package directory.
cp $path_file_script_source $path_file_script_product

# Activate Virtual Environment.
source "${path_environment_main}/bin/activate"
echo "confirm Python Virtual Environment path..."
which python3
sleep 5s
# Regulate concurrent or parallel process threads on node cores.
# Force Python program (especially SciPy) not to use all available cores on a
# cluster computation node.
export MKL_NUM_THREADS=$threads
export NUMEXPR_NUM_THREADS=$threads
export OMP_NUM_THREADS=$threads

python3 $path_file_script_product $path_directory_source $name_file_source_prefix $name_file_source_suffix $name_file_source_not $name_column_identifier $name_column_allele_total $name_column_allele_dosage $name_column_score_sum $name_column_score_mean $path_file_product

# Deactivate Virtual Environment.
deactivate
which python3

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Collection and combination of polygenic scores complete."
  echo "path to source directory: " $path_directory_source
  echo "file name prefix: " $name_file_source_prefix
  echo "file name suffix: " $name_file_source_suffix
  echo "file name exclusion: " $name_file_source_not
  echo "path to product file: " $path_file_product
  echo "----------"
  echo "----------"
  echo "----------"
fi



#
