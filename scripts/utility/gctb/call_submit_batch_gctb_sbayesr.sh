#!/bin/bash

################################################################################
# Notes:

# This script offers an example of a full execution of the SBayesR procedure
# from the GCTB tool package across all autosomal chromosomes (1-22).
# GCTB is a package of tools for Genome-wide Complex Trait Bayesian (GCTB)
# analysis that includes the tool SBayesR for calculation of Polygenic Scores
# (PGS).
# Documentation on GCTB.
# site: https://cnsgenomics.com/software/gctb/#Overview
# tutorial: https://cnsgenomics.com/software/gctb/#Tutorial
# educational presentation slides and video: https://www.cnsgenomics.com/data/teaching/PCTG/SBayes/

# Author: T. Cameron Waller
# Date, first execution: 17 January 2023
# Date, review: ____

# It is important that Genome-wide Association Study (GWAS) summary statistics
# and Linkage Disequilibrium (LD) reference matrices to use names (rs
# identifiers) and coordinates for Single Nucleotide Polymorphisms (SNPs) that
# are in the same assembly of the human genome, such as GRCh37 or GRCh38.
# As of year 2023, LD reference matrices for GCTB SBayesR used GRCh37.

################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_tools=$(<"./waller_tools.txt")
#path_plink2=$(<"./tools_plink2.txt")
path_gctb=$(<"./tools_waller_gctb.txt")
path_directory_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_product="${path_directory_dock}/test_sbayesr_full"
path_directory_ld_matrix="${path_directory_product}/ukbEURu_hm3_shrunk_sparse"

# Files.
path_file_gwas_source="${path_directory_dock}/hormone_genetics/gwas_format_standard/32042192_ruth_2020_testosterone_female.txt.gz"
path_file_gwas_product="${path_directory_product}/32042192_ruth_2020_testosterone_female.ma"
#path_file_ld_matrix_source="${path_directory_reference}/gctb/ukbEURu_hm3_sparse.zip"
path_file_ld_matrix_source="${path_directory_dock}/test_sbayesr/ukbEURu_hm3_sparse.zip"
path_file_ld_matrix_product="${path_directory_product}/ukbEURu_hm3_sparse.zip"
name_file_ld_matrix_prefix="ukbEURu_hm3_chr"
name_file_ld_matrix_suffix="_v3_50k.ldm.sparse"
name_file_product_prefix="sbayesr_female_testosterone_chromosome_"
name_file_product_suffix=""

# Scripts.
path_script_gwas_format="${path_directory_process}/promiscuity/scripts/utility/gctb/constrain_translate_gwas_standard_to_gctb.sh"
path_script_submit_batch="${path_directory_process}/promiscuity/scripts/utility/gctb/submit_batch_gctb_sbayesr_chromosomes.sh"
path_script_batch_run_sbayesr="${path_directory_process}/promiscuity/scripts/utility/gctb/batch_run_gctb_sbayesr.sh"
path_script_run_sbayesr="${path_directory_process}/promiscuity/scripts/utility/gctb/run_gctb_sbayesr.sh"

# Initialize directories.
#rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product



###########################################################################
# Organize parameters.

chromosome_x="false"
report="true"



###########################################################################
# Execute procedure.

##########
# 1. Prepare GWAS summary statistics.

if false; then
  /usr/bin/bash $path_script_gwas_format \
  $path_file_gwas_source \
  $path_file_gwas_product \
  $report
fi



##########
# 2. Prepare LD matrices.

if false; then
  cp $path_file_ld_matrix_source $path_file_ld_matrix_product
  unzip $path_file_ld_matrix_product -d $path_directory_product
fi



##########
# 3. Prepare and submit batch of jobs for processing on each chromosome.

if true; then
  /usr/bin/bash $path_script_submit_batch \
  $path_file_gwas_product \
  $path_directory_ld_matrix \
  $name_file_ld_matrix_prefix \
  $name_file_ld_matrix_suffix \
  $path_directory_product \
  $name_file_product_prefix \
  $name_file_product_suffix \
  $chromosome_x \
  $path_script_batch_run_sbayesr \
  $path_script_run_sbayesr \
  $path_gctb \
  $report
fi



##########
# 4. Remove temporary files.

# rm -r $path_directory_ld_matrix



#
