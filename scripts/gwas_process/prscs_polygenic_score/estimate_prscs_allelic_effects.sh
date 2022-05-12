#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes...

# This script calls PRS-CS to generate posterior effects across SNP alleles.

################################################################################
################################################################################
################################################################################

# After PRS-CS, there is some post-processing...
# 1. concatenate PRS-CS reports for each chromosome
# 2. use PLINK "--score" command

################################################################################
# Organize arguments.

path_source_gwas_summary=${1} # full path to source GWAS summary statistics in format for PRS-CS
count_gwas_samples=${2} # integer count of samples for the GWAS study
path_snp_relevance_bim_prefix=${3} # full path to file name prefix for '.bim' format list of relevant target SNPs
path_target_directory_score=${4} # full path to directory for product reports on posterior effect size estimates
name_file_product=${5} # name of product report file
path_genetic_reference_prscs=${6} # full path to directory for genetic references
population_ancestry=${7} # character code of ancestral population of GWAS study: 'AFR', 'AMR', 'EAS', 'EUR', or 'SAS'
chromosome=${8} # chromosome
report=${9} # whether to print reports

###########################################################################
# Organize paths.

################################################################################
# Activate Virtual Environment.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_tools=$(<"./waller_tools.txt")
path_prs_cs="${path_tools}/prs_cs/"
path_environment_prs_cs="${path_tools}/python/environments/prs_cs"
source "${path_environment_prs_cs}/bin/activate"
echo "confirm Python Virtual Environment path..."
which python3
sleep 5s

################################################################################
# Force SciPy not to use all available cores on a cluster computation node.

count_threads=4
export MKL_NUM_THREADS=$count_threads
export NUMEXPR_NUM_THREADS=$count_threads
export OMP_NUM_THREADS=$count_threads

################################################################################
# Test PRS-CSX.
# https://github.com/getian107/PRScsx

python3 "${path_prs_cs}/PRScsx/PRScsx.py" \
--ref_dir=$path_genetic_reference_prscs \
--bim_prefix=$path_snp_relevance_bim_prefix \
--sst_file=$path_source_gwas_summary \
--n_gwas=$count_gwas_samples \
--pop=$population_ancestry \
--out_dir=$path_target_directory_score \
--out_name=$name_file_product \
--a=1.0 \
--b=0.5 \
--phi=1e-3 \
--n_iter=1000 \
--n_burnin=500 \
--thin=5 \
--chrom=$chromosome \
--meta=False \
--seed=777

################################################################################
# Deactivate Virtual Environment.
deactivate
which python3

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "estimate_prscs_allelic_effects.sh"
  echo "----------"
  echo "GWAS: ${path_source_gwas_summary}"
  echo "chromosome: ${chromosome}"
  #cat $path_heritability_report_suffix
  echo "----------"
fi
