#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 1 August 2022
# Date, last execution: 4 August 2023
# Review: TCW; 4 August 2023
################################################################################
# Note


################################################################################



################################################################################
# Organize arguments.

path_file_source=${1} # full directory path and file name for source GWAS summary statistics in format with compression
path_file_base_product=${2} # full directory path and base file name GWAS summary statistics after munge
path_file_alleles=${3} # full directory path to LDSC reference file for HapMap3 SNPs
response=${4} # whether GWAS response is beta coefficient ("coefficient"), odds ratio ("odds_ratio"), or z-scores ("z_score")
threads=${5} # count of processing threads to use
report=${6} # whether to print reports

################################################################################
# Organize paths.

# Read private, local file paths.
cd ~/paths
path_tools=$(<"./waller_tools.txt")
#path_ldsc=$(<"./tools_ldsc.txt")
path_ldsc=$(<"./tools_ldsc_biotools.txt")

path_environment_ldsc="${path_tools}/python/environments/ldsc"
path_file_product_suffix="${path_file_base_product}.sumstats.gz"
path_file_product_log="${path_file_base_product}.log"

# Remove any previous version of the product files.
rm $path_file_product_suffix
rm $path_file_product_log

################################################################################
# Activate Virtual Environment.

source "${path_environment_ldsc}/bin/activate"
echo "confirm Python Virtual Environment path..."
which python2
sleep 1s

# Force SciPy not to use all available cores on a cluster computation node.
export MKL_NUM_THREADS=$threads
export NUMEXPR_NUM_THREADS=$threads
export OMP_NUM_THREADS=$threads

################################################################################
# Munge GWAS summary statistics for use in LDSC.
# Argument "--a1-inc" tells LDSC that GWAS sumstats do not have signed statistics but are coded so that A1 allele always increases effect.
# Argument ""--signed-sumstats" designates the column name for signed statistics.
# Examples:
# - "--signed-sumstats BETA,0"
# - "--signed-sumstats OR,1"
# - "--signed-sumstats Z,0"

if [[ "$response" == "coefficient" ]]; then
  $path_ldsc/munge_sumstats.py \
  --sumstats $path_file_source \
  --signed-sumstats BETA,0 \
  --merge-alleles $path_file_alleles \
  --out $path_file_base_product
elif [[ "$response" == "z_score" ]]; then
  $path_ldsc/munge_sumstats.py \
  --sumstats $path_file_source \
  --signed-sumstats Z,0 \
  --merge-alleles $path_file_alleles \
  --out $path_file_base_product
elif [[ "$response" == "odds_ratio" ]]; then
  $path_ldsc/munge_sumstats.py \
  --sumstats $path_file_source \
  --signed-sumstats OR,1 \
  --merge-alleles $path_file_alleles \
  --out $path_file_base_product
else
  echo "invalid specification of GWAS effect"
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
  echo "LDSC Munge on GWAS summary statistics."
  echo "path to source GWAS file: " $path_file_source
  echo "path to product GWAS file: " $path_file_product_suffix
  echo "table after munge:"
  #head -10 $path_file_product_suffix
  echo "LDSC Munge Log"
  #cat $path_file_product_log
  echo "----------"
  echo "----------"
  echo "----------"
fi



#
