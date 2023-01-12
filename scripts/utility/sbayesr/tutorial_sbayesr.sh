#!/bin/bash

################################################################################
# Notes:
# This script includes notes and code to work through an introductory tutorial
# on using the GCTB tool package.
# GCTB is a package of tools for Genome-wide Complex Trait Bayesian (GCTB)
# analysis that includes the tool SBayesR for calculation of Polygenic Scores
# (PGS).
# Documentation on GCTB.
# site: https://cnsgenomics.com/software/gctb/#Overview
# tutorial: https://cnsgenomics.com/software/gctb/#Tutorial
# educational presentation slides and video: https://www.cnsgenomics.com/data/teaching/PCTG/SBayes/

# Author: T. Cameron Waller
# Date, first execution: 10 January 2023
# Date, review: ____

################################################################################

##########
# Organize paths.

# Directories.
cd ~/paths
path_directory_tools=$(<"./waller_tools.txt")
path_plink2=$(<"./tools_plink2.txt")
path_sbayesr=$(<"./tools_waller_sbayesr.txt")
path_directory_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parent="${path_directory_dock}/test_sbayesr"

# Files.
path_file_gwas_source="${path_directory_dock}/hormone_genetics/gwas_format_standard/32042192_ruth_2020_testosterone_bioavailable_female.txt.gz"
path_file_gwas_product="${path_directory_parent}/32042192_ruth_2020_testosterone_bioavailable_female.ma"
path_file_gwas_tutorial="${path_directory_parent}/gctb_2.0_tutorial/ma/sim_1.ma"
path_file_ld_matrix_tutorial="${path_directory_parent}/gctb_2.0_tutorial/ldm/sparse/chr22/1000G_eur_chr22.ldm.sparse"
#path_file_ld_matrix_sparse_europe="${path_directory_parent}/..."
path_file_base_product="${path_directory_parent}/test_sbayesr_out_tcw_723"

# Scripts.
path_script_gwas_format="${path_directory_process}/promiscuity/scripts/utility/sbayesr/constrain_translate_gwas_standard_to_sbayesr.sh"
path_script_run_sbayesr="${path_directory_process}/promiscuity/scripts/utility/sbayesr/run_sbayesr.sh"

# Uniform Resource Locators (URLs).
url_gctb="https://cnsgenomics.com/software/gctb/download/gctb_2.04.3_Linux.zip"
url_gctb_tutorial="https://cnsgenomics.com/software/gctb/download/gctb_2.0_tutorial.zip"
url_ld_matrix_sparse_shrinkage_hapmap3="https://zenodo.org/record/3350914/files/ukbEURu_hm3_sparse.zip"
url_ld_matrix_sparse_chisquare="https://cnsgenomics.com/data/GCTB/ukbEURu_imp_v3_HM3_n50k.chisq10.zip"
url_ld_matrix_banded="https://cnsgenomics.com/data/GCTB/band_ukb_10k_hm3.zip"

# Initialize directories.
#rm -r $path_directory_parent
mkdir -p $path_directory_parent
cd $path_directory_parent



##########
# Installation.
# Authors provided a "statically linked 64-bit Linux executable, gctb".
# It is practical to run GCTB by a simple execution call to this handle.
# The authors also gave more thorough instructions for installation from source
# with local custom, local compilation.
# site: https://cnsgenomics.com/software/gctb/download/README.html
#cd ./ # navigate to directory in which to install GCTB.
if false; then
  # Access the specific file and save within the child directory.
  wget \
  "${url_gctb}" \
  --directory-prefix "${path_directory_tools}" \
  --content-disposition \
  --show-progress
  # Decompress content from Zip archive format.
  unzip "${path_directory_tools}/gctb_2.04.3_Linux.zip" -d "${path_directory_tools}"
  "${path_directory_tools}/gctb_2.04.3_Linux/gctb"
  $path_sbayesr
fi



##########
# Accession of data.

# Data for tutorial.
if false; then
  # Access the specific file and save within the child directory.
  wget \
  "${url_gctb_tutorial}" \
  --directory-prefix "${path_directory_parent}" \
  --content-disposition \
  --show-progress
  # Decompress content from Zip archive format.
  unzip "${path_directory_parent}/gctb_2.0_tutorial.zip" -d "${path_directory_parent}"
fi

# Linkage Disequilibrium (LD) reference matrix.
# Sparse (shrinkage) LD reference matrix on 50,000 persons of European ancestral
# background from the UK Biobank.
if false; then
  # Access the specific file and save within the child directory.
  wget \
  "${url_ld_matrix_sparse_shrinkage_hapmap3}" \
  --directory-prefix "${path_directory_parent}" \
  --content-disposition \
  --show-progress
  # Decompress content from Zip archive format.
  unzip "${path_directory_parent}/ukbEURu_hm3_sparse.zip" -d "${path_directory_parent}"
fi



##########
# Translate GWAS summary statistics to format for GCTB.
# columns: SNP   A1   A2   freq   b   se   p   N
if false; then
  /usr/bin/bash "${path_script_gwas_format}" \
  $path_file_gwas_source \
  $path_file_gwas_product \
  $report
fi

##########
# Apply SBayesR to adjust weights of effect sizes across SNPs.
# The path to the LD matrix actually points to three separate files with
# different suffixes: ".bin", ".info", ".log".
# Extra commands for SBayesR:
# --unscale-genotype
# --exclude-mhc
# --exclude-region
# --impute-n

if true; then
  $path_sbayesr \
  --sbayes R \
  --exclude-mhc \
  --ldm $path_file_ld_matrix_tutorial \
  --pi 0.95,0.02,0.02,0.01 \
  --gamma 0.0,0.01,0.1,1 \
  --gwas-summary $path_file_gwas_tutorial \
  --chain-length 10000 \
  --burn-in 2000 \
  --out-freq 10 \
  --out $path_file_base_product
fi

if false; then
  /usr/bin/bash $path_script_run_sbayesr \
  $path_file_gwas_tutorial \
  $path_file_ld_matrix_tutorial \
  $path_file_base_product \
  $path_sbayesr \
  $report
fi



##########
# Apply PLINK to calculate scores across SNPs in genotypes.






#
