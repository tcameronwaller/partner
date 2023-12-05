#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 5 December 2023
# Date, last execution: 5 December 2023
# Review: 5 December 2023
################################################################################
# Note

# This Bash script calls a Python script to convert LDSC estimates of SNP
# heritability from the observed scale to the liability scale.

# The format of the source file is a text table with tab delimiters
# (tab-separated values; tsv). The table has a first column of labels for the
# names of the studies, a second column of LDSC estimates of SNP heritability on
# the observed scale, a third column of the standard errors, a fourth column of
# the sample prevalence, and a fifth column of the population prevalence.

# The format of the product file is a text table with tab delimiters
# (tab-separated values; tsv) that includes all columns from the source file
# along with sixth and seventh columns of the SNP heritability and standard
# error on the liability scale.

# It is necessary to call the Python script within a package directory where it
# has access to appropriate subpackages and modules.


################################################################################
# Organize arguments.

path_file_source=${1} # full path to source file
path_file_product=${2} # full path to product file

################################################################################
# Organize paths.

# Directories.
path_tools="/media/tcameronwaller/primary/data/local/work/affiliation/general/tools"
path_environment_main="${path_tools}/python/environments/main"
path_directory_process="/home/tcameronwaller/Downloads/process_scratch_local"

path_directory_package="${path_directory_process}/scratch/scratch"

# Scripts.
path_file_script_source="${path_directory_process}/partner/scripts/python/script_convert_heritability_observed_liability.py"
path_file_script_product="${path_directory_package}/script_convert_heritability_observed_liability.py"

# Initialize directories.
rm $path_file_product # caution
cd $path_directory_process

################################################################################
# Organize parameters.

threads=1
report="true"

################################################################################
# Execute procedure.

##########
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

python3 $path_file_script_product $path_file_source $path_file_product

# Deactivate Virtual Environment.
deactivate
which python3

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Conversion of SNP heritabilities from observed to liability scales."
  echo "path to source file: " $path_file_source
  echo "path to product file: " $path_file_product
  echo "----------"
  echo "----------"
  echo "----------"
fi



#
