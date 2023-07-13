#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 17 August 2022
# Date, last execution:
# Review: TCW; 13 July 2023
################################################################################
# Note

# This script calls PRS-CS to generate posterior effects across SNP alleles.

# PRS-CSX documentation
# https://github.com/getian107/PRScsx

# Install PRS-CS and PRS-CSX.
#cd ./prs_cs/ # Navigate to the directory in which to install program.
#git clone https://github.com/getian107/PRScs.git
#git clone https://github.com/getian107/PRScsx.git

# After PRS-CS, there is some post-processing similar to SBayesR.
# 1. concatenate PRS-CS reports for each chromosome
# 2. use PLINK "--score" command
# 3. refer to pipeline for SBayesR

################################################################################



################################################################################
# Organize arguments.

path_file_gwas_summary=${1} # full path to source GWAS summary statistics in format for PRS-CS
count_gwas_samples=${2} # integer count of samples for the GWAS study
path_file_genotype_snp_bim_prefix=${3} # full path to file name prefix for list of relevant target SNPs in BIM format without '.bim' suffix
path_genetic_reference_prscs=${4} # full path to directory for genetic references
population_ancestry=${5} # character code of ancestral population of GWAS study: 'AFR', 'AMR', 'EAS', 'EUR', or 'SAS'
path_directory_allele_effect=${6} # full path to directory for product reports on posterior allele effect size estimates
name_file_allele_effect=${7} # name of product report file without file suffix
chromosome=${8} # chromosome
threads=${9} # count of processing threads to use
path_environment_prscs=${10} # full path to Python 3 environment with installation of dependencies
path_prscsx=${11} # full path to installation executable file of PRS-CSX
report=${12} # whether to print reports

###########################################################################
# Organize paths.

# Activate Virtual Environment.
source "${path_environment_prscs}/bin/activate"
echo "confirm Python Virtual Environment path..."
which python3
sleep 5s

# Force SciPy not to use all available cores on a cluster computation node.
export MKL_NUM_THREADS=$threads
export NUMEXPR_NUM_THREADS=$threads
export OMP_NUM_THREADS=$threads

# Calculate posterior effects.
python3 $path_prscsx \
--ref_dir=$path_genetic_reference_prscs \
--bim_prefix=$path_file_genotype_snp_bim_prefix \
--sst_file=$path_file_gwas_summary \
--n_gwas=$count_gwas_samples \
--pop=$population_ancestry \
--out_dir=$path_directory_allele_effect \
--out_name=$name_file_allele_effect \
--a=1.0 \
--b=0.5 \
--phi=1e-3 \
--n_iter=1000 \
--n_burnin=500 \
--thin=5 \
--chrom=$chromosome \
--meta=False \
--seed=777

# Deactivate Virtual Environment.
deactivate
which python3

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "estimate_prscs_allelic_effects.sh"
  echo "----------"
  echo "GWAS: ${path_file_gwas_summary}"
  echo "chromosome: ${chromosome}"
  #cat $path_heritability_report_suffix
  echo "----------"
fi
