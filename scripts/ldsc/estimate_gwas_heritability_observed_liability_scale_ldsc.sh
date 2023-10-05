#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 22 November 2022
# Date, last execution: 5 October 2023
# Date, review: 8 September 2023
################################################################################
# Note

# TODO: TCW; 19 July 2023
# TODO: attempt to decrease standard error by constraining the LDSC intercept to
# 1.
# https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation


################################################################################
# Organize arguments.

path_file_source=${1} # full directory path and file name for source GWAS summary statistics in format with compression
path_file_base_product=${2} # full directory path and base file name for LDSC heritability log report
path_directory_disequilibrium=${3} # full directory path to directory for LDSC reference on linkage disequilibrium
scale=${4} # whether to estimate SNP heritability on the observed or liability scales
prevalence_sample=${5} # estimate of phenotype prevalence in the sample
prevalence_population=${6} # estimate of phenotype prevalence in the population
threads=${7} # count of processing threads to use
report=${8} # whether to print reports

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

#$path_ldsc/ldsc.py \
#--h2 $path_file_source \
#--ref-ld-chr ${path_directory_disequilibrium}/ \
#--w-ld-chr ${path_directory_disequilibrium}/ \
#--out $path_file_base_product

if [[ "$scale" == "observed" ]]; then
  $path_ldsc/ldsc.py \
  --h2 $path_file_source \
  --ref-ld-chr ${path_directory_disequilibrium}/ \
  --w-ld-chr ${path_directory_disequilibrium}/ \
  --out $path_file_base_product
elif [[ "$scale" == "liability" ]]; then
  $path_ldsc/ldsc.py \
  --h2 $path_file_source \
  --ref-ld-chr ${path_directory_disequilibrium}/ \
  --w-ld-chr ${path_directory_disequilibrium}/ \
  --samp-prev $prevalence_sample \
  --pop-prev $prevalence_population \
  --out $path_file_base_product
fi




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
