#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes...

# Review: TCW; 22 November 2022

# Review: TCW; 8 September 2023
# This script does not yet call the liability-scale functionality of LDSC!!!

################################################################################
################################################################################
################################################################################

################################################################################
# Organize arguments.

path_file_source=${1} # full directory path and file name for source GWAS summary statistics in format with compression
path_file_base_product=${2} # full directory path and base file name for LDSC heritability log report
path_directory_disequilibrium=${3} # full directory path to directory for LDSC reference on linkage disequilibrium
threads=${4} # count of processing threads to use
report=${5} # whether to print reports

################################################################################
# Organize paths.

# Read private, local file paths.
cd ~/paths
path_tools=$(<"./waller_tools.txt")
#path_ldsc=$(<"./tools_ldsc.txt")
path_ldsc=$(<"./tools_ldsc_biotools.txt")

path_environment_ldsc="${path_tools}/python/environments/ldsc"
path_file_product_suffix="${path_file_base_product}.log"

# Remove any previous version of the product files.
rm $path_file_product_suffix

################################################################################
# Activate Virtual Environment.

source "${path_environment_ldsc}/bin/activate"
echo "confirm Python Virtual Environment path..."
which python2
sleep 5s

# Force SciPy not to use all available cores on a cluster computation node.
export MKL_NUM_THREADS=$threads
export NUMEXPR_NUM_THREADS=$threads
export OMP_NUM_THREADS=$threads

################################################################################
# Estimate SNP heritability in LDSC from GWAS summary statistics.

$path_ldsc/ldsc.py \
--h2 $path_file_source \
--ref-ld-chr ${path_directory_disequilibrium}/ \
--w-ld-chr ${path_directory_disequilibrium}/ \
--out $path_file_base_product

################################################################################
# Deactivate Virtual Environment.
deactivate
which python2

################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "LDSC estimate of SNP heritability on GWAS summary statistics."
  echo "path to source GWAS file: " $path_file_source
  echo "path to product report file: " $path_file_product_suffix
  #head -10 $path_file_product_suffix
  #echo "LDSC Munge Log"
  #cat $path_file_product_log
  echo "----------"
  echo "----------"
  echo "----------"
fi



#
