#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 17 May 2024
# Date, last execution: 21 May 2024
# Review: TCW; __ May 2024
################################################################################
# Note

# This Bash script calls a Python script to extract information from the raw
# text report logs that the LDSC tool creates for estimates of SNP heritability
# (h2) and genetic correlation (rg). This procedure accommodates a large count
# of individual files and assembles and organizes these within a text table.

# It is necessary to call the Python script within a package directory where it
# has access to appropriate subpackages and modules. Call this script within a
# context in which it has access to the rest of the partner package.

################################################################################
# Organize arguments.

type_analysis=${1} # type of analysis, heritability or correlation, for which to extract information
path_directory_source=${2} # full path to source directory from which to read files
traversal=${3} # whether to extract from all files in child source directories, preserving names of child directories
name_file_source_prefix=${4} # prefix in name of files to read
name_file_source_suffix=${5} # suffix in name of files to read
name_file_source_not=${6} # character string within name of files to exclude
name_file_product=${7} # either none for traversal of child directories or the name to use for product file
path_directory_product=${8} # full path to product directory within which to write files
path_directory_process=${9} # full path to process directory
path_directory_temporary=${10} # full path to temporary directory for working space
path_directory_environment=${11} # full path to directory of Python virtual environment
report=${12} # whether to print reports

################################################################################
# Organize paths.

# Directories.
path_directory_package_source="${path_directory_process}/partner/package"
path_directory_package_source_copy="${path_directory_temporary}/package/package"
path_directory_package_product_parent="${path_directory_temporary}/package"
path_directory_package_product="${path_directory_temporary}/package/partner"
#path_directory_product="$(dirname $path_file_product)"

# Scripts.
path_file_script_source="${path_directory_process}/partner/scripts/python/extract_ldsc_heritability_correlation.py"
path_file_script_product="${path_directory_package_product_parent}/extract_ldsc_heritability_correlation.py"

# Initialize directories.
mkdir -p $path_directory_package_product_parent
cd $path_directory_product

################################################################################
# Organize parameters.

threads=2

################################################################################
# Execute procedure.

# Organize package directory.
cp -r $path_directory_package_source $path_directory_package_product_parent
mv $path_directory_package_source_copy $path_directory_package_product
cp "${path_directory_package_product}/__init__.py" "${path_directory_package_product_parent}/__init__.py"

# Copy Python script to package directory.
cp $path_file_script_source $path_file_script_product

# Activate Python Virtual Environment.
source "${path_directory_environment}/bin/activate"
echo "confirm Python Virtual Environment path..."
which python3
sleep 5s
# Regulate concurrent or parallel process threads on node cores.
# Force Python program (especially SciPy) not to use all available cores on a
# cluster computation node.
export MKL_NUM_THREADS=$threads
export NUMEXPR_NUM_THREADS=$threads
export OMP_NUM_THREADS=$threads

python3 $path_file_script_product \
$type_analysis \
$path_directory_source \
$traversal \
$name_file_source_prefix \
$name_file_source_suffix \
$name_file_source_not \
$name_file_product \
$path_directory_product \
$path_directory_temporary \
$report

# Deactivate Python Virtual Environment.
deactivate
which python3

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Extraction complete."
  echo "path to source directory: " $path_directory_source
  echo "path to product directory: " $path_directory_product
  echo "----------"
  echo "----------"
  echo "----------"
fi



#
