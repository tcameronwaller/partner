#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes:
# This script functions as a customizable call to the SBayesR procedure of the
# GCTB tool package.
# GCTB is a package of tools for Genome-wide Complex Trait Bayesian (GCTB)
# analysis that includes the tool SBayesR for calculation of Polygenic Scores
# (PGS).
# Documentation on GCTB.
# site: https://cnsgenomics.com/software/gctb/#Overview
# tutorial: https://cnsgenomics.com/software/gctb/#Tutorial
# educational presentation slides and video: https://www.cnsgenomics.com/data/teaching/PCTG/SBayes/

# It is advisable to use this script to run the SBayesR procedure on a single
# chromosome at a time via an upstream driver script or a submission to a
# computational cluster (such as the Sun Grid Engine) for parallel processing.

# Author: T. Cameron Waller
# Date, first execution: 12 January 2023
# Date, review: 16 January 2023

################################################################################
################################################################################
################################################################################



################################################################################
# Organize arguments.

path_file_gwas=${1} # full directory path and file name for source GWAS summary statistics in GCTB and GCTA-COJO ".ma" format without compression
path_file_base_ld_matrix=${2} # full directory path and base file name for Linkage Disequilibrium (LD) reference matrix in GCTB format
path_file_base_product=${3} # full directory path and base file name for product files from GCTB SBayesR
path_gctb=${4} # full directory path and file name for local executable installation of GCTB SBayesR
report=${5} # whether to print reports

################################################################################
# Run SBayesR.

##########
# Apply SBayesR to adjust weights of effect sizes across SNPs.
# The path to the LD matrix actually points to three separate files with
# different suffixes: ".bin", ".info", ".log".
# Extra commands for SBayesR:
# --unscale-genotype
# --exclude-mhc
# --exclude-region
# --impute-n

$path_gctb \
--sbayes R \
--exclude-mhc \
--ldm $path_file_base_ld_matrix \
--pi 0.95,0.02,0.02,0.01 \
--gamma 0.0,0.01,0.1,1 \
--gwas-summary $path_file_gwas \
--chain-length 10000 \
--burn-in 2000 \
--out-freq 10 \
--out $path_file_base_product 2>&1 | tee "${path_file_base_product}.log"



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "run_gctb_sbayesr.sh"
  echo "----------"
fi



#
