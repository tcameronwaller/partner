#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 12 January 2023
################################################################################
# Notes:
# This script functions as a customizable call to the SBayesR procedure of the
# GCTB tool package.
# GCTB is a package of tools for Genome-wide Complex Trait Bayesian (GCTB)
# analysis that includes the tool SBayesR for calculation of Polygenic Scores
# (PGS).
# Publication: Lloyd-Jones, Nature Communications, 2019 (PubMed:31704910)
# Documentation on GCTB.
# Site: http://cnsgenomics.com/software/gctb/
# Site: https://cnsgenomics.com/software/gctb/#Overview
# Site: https://cnsgenomics.com/software/gctb/#Tutorial
# Tutorial: https://cnsgenomics.com/software/gctb/#Tutorial
# educational presentation slides and video: https://www.cnsgenomics.com/data/teaching/PCTG/SBayes/

# It is advisable to use this script to run the SBayesR procedure on a single
# chromosome at a time via an upstream driver script or a submission to a
# computational cluster (such as the Sun Grid Engine) for parallel processing.

# Date, review: 16 January 2023

################################################################################
################################################################################
################################################################################

# TODO: TCW; 20 Februayr 2023
# New argument:
# observations_variant <-- whether observations are variant-specific

################################################################################
# Organize arguments.

path_file_gwas=${1} # full directory path and file name for source GWAS summary statistics in GCTB and GCTA-COJO ".ma" format without compression
path_file_base_ld_matrix=${2} # full directory path and base file name for Linkage Disequilibrium (LD) reference matrix in GCTB format
path_file_base_product=${3} # full directory path and base file name for product files from GCTB SBayesR
observations_variant=${4} # logical binary indicator of whether counts of observations are reliable and specific to each variant (SNP)
path_gctb=${5} # full directory path and file name for local executable installation of GCTB SBayesR
threads=${6} # count of concurrent or parallel process threads on node cores
report=${7} # whether to print reports

################################################################################
# Execute procedure.

# Regulate concurrent or parallel process threads on node cores.
export MKL_NUM_THREADS=$threads
export NUMEXPR_NUM_THREADS=$threads
export OMP_NUM_THREADS=$threads

##########
# Apply SBayesR to adjust weights of effect sizes across SNPs.
# The file path to the LD matrix actually points to three separate files with
# different suffixes: ".bin", ".info", ".log".
# Extra commands for SBayesR:
# --unscale-genotype # remove default assumption that genotypes on unit variance scale (not sure what this means; TCW; 20 February 2023)
# --exclude-mhc # exclude Human Leukocyte Antigen (HLA) Mayor Histocompatibility Complex (MHC) region (extremely high heterogeneity)
# --exclude-region
# --impute-n # estimate (impute) variant-specific counts of observations
# --robust # apply an alternative parameterisation for SNP effect variance

# Extra commands for SBayesR:
# --unscale-genotype # <-- Recommended!
# --exclude-mhc # <-- Recommended!
# --exclude-region
# --impute-n # <-- use this if the count of observations is unreliable
# --robust # apply a more robust algorithm for convergence
# File with suffix ".snpRes" gives the new effect sizes across SNPs after
# adjustment of weights for LD (I think; TCW; 12 January 2023).


if [[ "$observations_variant" == "1" ]]; then
  $path_gctb \
  --sbayes R \
  --gwas-summary $path_file_gwas \
  --ldm $path_file_base_ld_matrix \
  --pi 0.95,0.02,0.02,0.01 \
  --gamma 0.0,0.01,0.1,1.0 \
  --chain-length 10000 \
  --burn-in 2000 \
  --out-freq 10 \
  --unscale-genotype \
  --exclude-mhc \
  --robust \
  --out $path_file_base_product 2>&1 | tee "${path_file_base_product}.log"
elif [[ "$observations_variant" == "0" ]]; then
  $path_gctb \
  --sbayes R \
  --gwas-summary $path_file_gwas \
  --ldm $path_file_base_ld_matrix \
  --pi 0.95,0.02,0.02,0.01 \
  --gamma 0.0,0.01,0.1,1.0 \
  --chain-length 10000 \
  --burn-in 2000 \
  --out-freq 10 \
  --unscale-genotype \
  --exclude-mhc \
  --robust \
  --impute-n \
  --out $path_file_base_product 2>&1 | tee "${path_file_base_product}.log"
fi



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "run_gctb_sbayesr.sh"
  echo "----------"
fi



#
