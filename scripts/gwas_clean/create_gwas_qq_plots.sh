#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 27 March 2023
# Date, last execution: 27 March 2023
# Review: TCW; __ March 2023
################################################################################
# Note

# This Bash script calls a Python script to create QQ Plots for sets of GWAS
# summary statistics in files within a source directory.

# It is necessary to call the Python script within a package directory where it
# has access to appropriate subpackages and modules.

# Source Format
# Description: Team Standard
# File suffix: ".txt.gz"
# File type: text
# File compression: GZip
# Delimiter: white space
# Effect allele: A1
# Probability (p-value): P
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

# References on QQ Plots for GWAS summary statistics:
# https://www.broadinstitute.org/diabetes-genetics-initiative/plotting-genome-wide-association-results
# https://www.broadinstitute.org/files/shared/diabetes/scandinavs/qqplot.R


################################################################################
# Organize arguments.

path_directory_source=${1} # full path to source directory from which to read files of GWAS summary statistics
name_file_source_prefix=${2} # prefix of name of files of polygenic scores from each chromosome
name_file_source_suffix=${3} # suffix of name of files of polygenic scores from each chromosome
name_file_source_not=${4} # character string within file names for exclusion
path_directory_product=${5} # full path to product directory in which to write files for plots

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_waller_tools=$(<"./waller_tools.txt")
path_environment_main="${path_waller_tools}/python/environments/main"
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_package="${path_directory_process}/psychiatric_metabolism/psychiatric_metabolism"

# Scripts.
path_file_script_source="${path_directory_process}/promiscuity/scripts/python/script_create_gwas_qq_plots.py"
path_file_script_product="${path_directory_package}/script_create_gwas_qq_plots.py"

# Initialize directories.
#rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product

################################################################################
# Organize parameters.

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

python3 $path_file_script_product $path_directory_source $name_file_source_prefix $name_file_source_suffix $name_file_source_not $path_directory_product

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
  echo "path to product directory: " $path_directory_product
  echo "----------"
  echo "----------"
  echo "----------"
fi



#
