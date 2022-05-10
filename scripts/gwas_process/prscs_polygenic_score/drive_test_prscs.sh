#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes...

################################################################################
################################################################################
################################################################################

# After PRS-CS, there is some post-processing...
# 1. concatenate PRS-CS reports for each chromosome
# 2. use PLINK "--score" command

################################################################################
# Organize arguments.
path_genetic_reference_prscs=${1} # full path to directory for genetic references
path_source_gwas_summary=${2} # full path to source GWAS summary statistics in format for PRS-CS
count_gwas_samples=${3} # integer count of samples for the GWAS study
population_ancestry=${4} # character code of ancestral population of GWAS study: 'AFR', 'AMR', 'EAS', 'EUR', or 'SAS'
path_snp_relevance_bim_prefix=${5} # full path to file name prefix for '.bim' format list of relevant target SNPs
path_target_directory_score=${6} # full path to directory for product reports on posterior effect size estimates
file_name_prefix=${7} # prefix for names of product report files
chromosome=${8} # chromosome
report=${9} # whether to print reports

###########################################################################
# Organize paths.

path_1000_genomes="$path_genetic_reference_prscs/1000_genomes"
path_uk_biobank="$path_genetic_reference_prscs/uk_biobank"

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
--ref_dir=$path_1000_genomes \
--bim_prefix=$path_snp_relevance_bim_prefix \
--sst_file=$path_source_gwas_summary \
--n_gwas=$count_gwas_samples \
--pop=$population_ancestry \
--out_dir=$path_target_directory_score \
--out_name=$file_name_prefix \
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
  echo "drive_prscs_test.sh"
  echo "----------"
  echo "PRS-CS report:"
  #cat $path_heritability_report_suffix
  echo "----------"
fi
